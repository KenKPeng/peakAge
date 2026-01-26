#' Fit latent peak-age model
#'
#' Main entry point: prepares Stan data, compiles the Stan model, and runs MCMC.
#'
#' Prior design (v1):
#' - We encourage users to standardize y to roughly 0-1. If not, adjust p prior accordingly.
#' - Users typically only need to set p prior (and optionally k prior).
#' - Other priors default to weak / conservative values and are not the primary user surface.
#'
#' @param df Long-format data frame. Must include athlete id, age, and y columns.
#' @param curve Curve family: "symmetric" or "asymmetric".
#' @param peak_formula Athlete-level covariate formula for log peak age (AFT),
#'   e.g. ~ x1 + x2. Use ~ 1 for intercept-only. Categorical covariates are automatically handled using R’s default
#'   treatment contrasts. The first level of each factor is treated as the
#'   reference group, and regression coefficients represent differences
#'   relative to that reference.
#' @param p_formula Athlete-level covariate formula for p (peak height) submodel.
#'   Default ~ 1 (intercept handled by p_pop; if ~1, treated as no covariates).
#' @param k_formula Athlete-level covariate formula for curvature submodel.
#'   For symmetric: affects logk via X_k * beta_k.
#'   For asymmetric: the same X_k affects both logkL and logkR via beta_kL / beta_kR.
#'   Default ~ 1 (treated as no covariates).
#' @param prior List specifying priors. Minimal recommended:
#'   prior = list(
#'     p = list(pop = c(mean, sd)),
#'     k = list(pop = c(mean, sd))   # symmetric
#'     # or asymmetric:
#'     # k = list(pop = list(L = c(meanL, sdL), R = c(meanR, sdR)))
#'   )
#' All other priors have weak defaults (see code).
#' @param formula Deprecated/compatibility alias for `peak_formula`.
#' @param id,age,y Column names in df.
#' @param grainsize grainsize for Stan reduce_sum.
#' @param iter_warmup,iter_sampling MCMC iterations.
#' @param chains,parallel_chains,threads_per_chain Sampling parallelism.
#' @param adapt_delta,max_treedepth Stan NUTS controls.
#' @param seed Random seed passed to Stan.
#' @param refresh Stan sampling refresh.
#' @param stan_file Optional path to Stan file. If NULL, uses inst/stan based on `curve`.
#' @param ... Passed to cmdstanr::CmdStanModel$sample().
#'
#' @return An object of class `peakage_fit`.
#' @export
fit_peakage <- function(df,
                        curve = c("symmetric", "asymmetric"),
                        peak_formula = ~ 1,
                        p_formula    = ~ 1,
                        k_formula    = ~ 1,
                        prior = NULL,
                        formula = NULL,
                        id = "Athlete_id",
                        age = "age",
                        y = "y",
                        grainsize = 100,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        chains = 4,
                        parallel_chains = chains,
                        threads_per_chain = 1,
                        adapt_delta = 0.95,
                        max_treedepth = 12,
                        seed = NULL,
                        refresh = 0,
                        stan_file = NULL,
                        ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package `cmdstanr` is required. Please install it first.")
  }

  curve <- match.arg(curve)

  # ---- backward compatibility: formula -> peak_formula ----
  if (!is.null(formula)) peak_formula <- formula

  # --- locate Stan file inside installed package (preferred) ---
  if (is.null(stan_file)) {
    pkg <- utils::packageName()

    stan_name <- switch(
      curve,
      symmetric  = "peakage_sym.stan",
      asymmetric = "peakage_asym.stan"
    )

    stan_file <- system.file(file.path("stan", stan_name), package = pkg)
    if (!nzchar(stan_file)) stan_file <- file.path("inst", "stan", stan_name)
  }
  if (!file.exists(stan_file)) stop("Stan file not found at: ", stan_file)

  # --- build shared Stan data ---
  stan_data <- build_stan_data_sym(
    df = df,
    peak_formula = peak_formula,
    p_formula = p_formula,
    k_formula = k_formula,
    id = id,
    age = age,
    y = y,
    grainsize = grainsize
  )

  # ---- store athlete id mapping (aligned with Stan index i=1..I) ----
  if (!id %in% names(df)) stop("`df` must contain id column: ", id)
  if (anyNA(df[[id]])) stop("Missing values in athlete id column.")
  athlete_ids <- sort(unique(df[[id]]))

  # ---- build prior defaults + apply user overrides ----
  prior_used <- .peakage_build_priors(prior = prior, curve = curve)

  # ---- inject prior hyperparameters into Stan data (names must match Stan `data {}`) ----
  # common
  stan_data$p_pop_loc   <- prior_used$p$pop[1]
  stan_data$p_pop_scale <- prior_used$p$pop[2]

  stan_data$tau_u_rate    <- prior_used$hidden$tau_u_rate
  stan_data$sigma_y_scale <- prior_used$hidden$sigma_y_sd

  stan_data$beta_T_scale  <- prior_used$hidden$beta_T_sd
  stan_data$beta_p_scale  <- prior_used$hidden$beta_p_sd
  stan_data$sigma_T_scale <- prior_used$hidden$sigma_T_sd

  if (curve == "symmetric") {
    stan_data$logk_pop_loc   <- prior_used$k$pop[1]
    stan_data$logk_pop_scale <- prior_used$k$pop[2]
    stan_data$beta_k_scale   <- prior_used$hidden$beta_k_sd

    # Stan expects vector[2]
    if (length(stan_data$tau_u_rate) != 2L) stan_data$tau_u_rate <- rep(stan_data$tau_u_rate[1], 2)
  } else {
    # asymmetric
    kp <- prior_used$k$pop
    if (is.numeric(kp) && length(kp) == 2L) {
      stan_data$logkL_pop_loc   <- kp[1]
      stan_data$logkL_pop_scale <- kp[2]
      stan_data$logkR_pop_loc   <- kp[1]
      stan_data$logkR_pop_scale <- kp[2]
    } else if (is.list(kp) && !is.null(kp$L) && !is.null(kp$R)) {
      stan_data$logkL_pop_loc   <- kp$L[1]
      stan_data$logkL_pop_scale <- kp$L[2]
      stan_data$logkR_pop_loc   <- kp$R[1]
      stan_data$logkR_pop_scale <- kp$R[2]
    } else {
      stop("For curve='asymmetric', prior$k$pop must be c(mean, sd) or list(L=c(mean,sd), R=c(mean,sd)).")
    }

    stan_data$beta_kL_scale <- prior_used$hidden$beta_k_sd
    stan_data$beta_kR_scale <- prior_used$hidden$beta_k_sd

    # Stan expects vector[3]
    if (length(stan_data$tau_u_rate) != 3L) stan_data$tau_u_rate <- rep(stan_data$tau_u_rate[1], 3)
  }

  # --- compile model (cache) ---
  mod <- .peakage_get_model(stan_file)

  # --- init function (simple stable defaults; not data-dependent priors) ---
  init_fun <- function() {
    # crude center for logT_i init: log(mean(age)) but don't use it as prior
    logT0 <- log(mean(stan_data$age))

    base <- list(
      p_pop = prior_used$p$pop[1],
      sigma_y = max(1e-3, prior_used$hidden$sigma_y_sd / 2),

      beta_T  = rep(0, stan_data$K_T),
      sigma_T = max(1e-3, prior_used$hidden$sigma_T_sd / 2),
      logT_i  = rep(logT0, stan_data$I),

      beta_p = rep(0, stan_data$K_p)
    )

    if (curve == "symmetric") {
      c(base, list(
        logk_pop = prior_used$k$pop[1],
        tau_u = rep(0.3, 2),
        z_u = array(0, dim = c(stan_data$I, 2)),
        beta_k = rep(0, stan_data$K_k)
      ))
    } else {
      # asymmetric
      # if user gave L/R distinct, init at their means; else shared
      kp <- prior_used$k$pop
      logkL0 <- if (is.list(kp)) kp$L[1] else kp[1]
      logkR0 <- if (is.list(kp)) kp$R[1] else kp[1]

      c(base, list(
        logkL_pop = logkL0,
        logkR_pop = logkR0,
        tau_u = rep(0.3, 3),
        z_u = array(0, dim = c(stan_data$I, 3)),
        beta_kL = rep(0, stan_data$K_k),
        beta_kR = rep(0, stan_data$K_k)
      ))
    }
  }

  fit <- mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    refresh = refresh,
    init = init_fun,
    ...
  )

  out <- list(
    fit = fit,
    stan_data = stan_data,
    athlete_ids = athlete_ids,
    curve = curve,
    prior = prior_used,
    peak_formula = peak_formula,
    p_formula = p_formula,
    k_formula = k_formula,
    id_col = id,
    age_col = age,
    y_col = y,
    stan_file = stan_file,
    call = match.call()
  )
  class(out) <- "peakage_fit"
  out
}

# ---- internal: prior builder ----
.peakage_build_priors <- function(prior = NULL, curve = c("symmetric", "asymmetric")) {
  curve <- match.arg(curve)

  # v1 defaults assuming y roughly in [0, 1]
  # user should adjust p prior if not standardized
  defaults <- list(
    p = list(pop = c(0.5, 1.0)),        # p_pop ~ N(0.5, 1)
    k = list(pop = c(0.0, 2.0)),        # logk_pop ~ N(0, 2) or used for both sides
    hidden = list(
      # tau_u ~ Exp(rate): smaller rate => heavier tail
      tau_u_rate = if (curve == "symmetric") c(1.0, 1.0) else c(1.0, 1.0, 1.0),

      # sigma_y ~ N(0, scale) with <lower=1e-6>, effectively half-normal
      sigma_y_sd = 0.5,

      # regression coefficients (shared scale)
      beta_T_sd = 2.0,
      beta_p_sd = 2.0,
      beta_k_sd = 2.0,

      sigma_T_sd = 1.0
    )
  )

  if (is.null(prior)) return(defaults)
  if (!is.list(prior)) stop("`prior` must be a list.")

  out <- defaults

  # helpers
  is_num2 <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x)) && x[2] > 0
  is_pos <- function(x) is.numeric(x) && length(x) >= 1L && all(is.finite(x)) && all(x > 0)

  # p prior
  if (!is.null(prior$p)) {
    if (!is.list(prior$p) || is.null(prior$p$pop) || !is_num2(prior$p$pop)) {
      stop("prior$p$pop must be numeric c(mean, sd) with sd > 0.")
    }
    out$p$pop <- as.numeric(prior$p$pop)
  }

  # k prior
  if (!is.null(prior$k)) {
    if (!is.list(prior$k) || is.null(prior$k$pop)) {
      stop("prior$k must be a list with element `pop`.")
    }

    if (curve == "symmetric") {
      if (!is_num2(prior$k$pop)) stop("For symmetric, prior$k$pop must be numeric c(mean, sd).")
      out$k$pop <- as.numeric(prior$k$pop)
    } else {
      kp <- prior$k$pop
      if (is_num2(kp)) {
        out$k$pop <- as.numeric(kp)
      } else if (is.list(kp) && is_num2(kp$L) && is_num2(kp$R)) {
        out$k$pop <- list(L = as.numeric(kp$L), R = as.numeric(kp$R))
      } else {
        stop("For asymmetric, prior$k$pop must be c(mean, sd) or list(L=c(mean,sd), R=c(mean,sd)).")
      }
    }
  }

  # allow power users to override hidden pieces (optional; not documented heavily)
  if (!is.null(prior$hidden)) {
    if (!is.list(prior$hidden)) stop("prior$hidden must be a list.")
    h <- prior$hidden

    if (!is.null(h$tau_u_rate)) {
      if (!is_pos(h$tau_u_rate)) stop("prior$hidden$tau_u_rate must be positive numeric vector.")
      out$hidden$tau_u_rate <- as.numeric(h$tau_u_rate)
    }
    if (!is.null(h$sigma_y_sd)) {
      if (!is.numeric(h$sigma_y_sd) || length(h$sigma_y_sd) != 1L || !is.finite(h$sigma_y_sd) || h$sigma_y_sd <= 0) {
        stop("prior$hidden$sigma_y_sd must be a single positive number.")
      }
      out$hidden$sigma_y_sd <- as.numeric(h$sigma_y_sd)
    }
    if (!is.null(h$beta_T_sd)) {
      if (!is.numeric(h$beta_T_sd) || length(h$beta_T_sd) != 1L || !is.finite(h$beta_T_sd) || h$beta_T_sd <= 0) {
        stop("prior$hidden$beta_T_sd must be a single positive number.")
      }
      out$hidden$beta_T_sd <- as.numeric(h$beta_T_sd)
    }
    if (!is.null(h$beta_p_sd)) {
      if (!is.numeric(h$beta_p_sd) || length(h$beta_p_sd) != 1L || !is.finite(h$beta_p_sd) || h$beta_p_sd <= 0) {
        stop("prior$hidden$beta_p_sd must be a single positive number.")
      }
      out$hidden$beta_p_sd <- as.numeric(h$beta_p_sd)
    }
    if (!is.null(h$beta_k_sd)) {
      if (!is.numeric(h$beta_k_sd) || length(h$beta_k_sd) != 1L || !is.finite(h$beta_k_sd) || h$beta_k_sd <= 0) {
        stop("prior$hidden$beta_k_sd must be a single positive number.")
      }
      out$hidden$beta_k_sd <- as.numeric(h$beta_k_sd)
    }
    if (!is.null(h$sigma_T_sd)) {
      if (!is.numeric(h$sigma_T_sd) || length(h$sigma_T_sd) != 1L || !is.finite(h$sigma_T_sd) || h$sigma_T_sd <= 0) {
        stop("prior$hidden$sigma_T_sd must be a single positive number.")
      }
      out$hidden$sigma_T_sd <- as.numeric(h$sigma_T_sd)
    }
  }

  # normalize tau length
  if (curve == "symmetric") {
    if (length(out$hidden$tau_u_rate) != 2L) out$hidden$tau_u_rate <- rep(out$hidden$tau_u_rate[1], 2)
  } else {
    if (length(out$hidden$tau_u_rate) != 3L) out$hidden$tau_u_rate <- rep(out$hidden$tau_u_rate[1], 3)
  }

  out
}

# ---- internal: model cache ----
.peakage_model_cache <- new.env(parent = emptyenv())

.peakage_get_model <- function(stan_file) {
  key <- normalizePath(stan_file, winslash = "/", mustWork = TRUE)
  if (exists(key, envir = .peakage_model_cache, inherits = FALSE)) {
    return(get(key, envir = .peakage_model_cache, inherits = FALSE))
  }
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = TRUE),
    force_recompile = FALSE
  )
  assign(key, mod, envir = .peakage_model_cache)
  mod
}
