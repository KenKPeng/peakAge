# Build separate athlete and competition design matrices for the extended model.
.peakage_extended_data <- function(df, formulas, observation_formula, peak_sd_group,
                                   id, age, y, grainsize) {
  if (!is.null(peak_sd_group) &&
      (!is.character(peak_sd_group) || length(peak_sd_group) != 1L ||
        is.na(peak_sd_group) || !nzchar(peak_sd_group))) {
    stop("peak_sd_group must be one column name.")
  }
  columns <- c(id, age, y, peak_sd_group)
  if (!is.data.frame(df) || !nrow(df) || !all(columns %in% names(df))) {
    stop("Provide a nonempty data frame with columns: ", paste(columns, collapse = ", "), ".")
  }
  if (anyNA(df[[id]]) || !is.atomic(df[[id]]) || !is.null(dim(df[[id]]))) {
    stop("Athlete IDs must be an atomic vector without missing values.")
  }
  if (!is.numeric(df[[age]]) || !is.numeric(df[[y]]) ||
      !is.null(dim(df[[age]])) || !is.null(dim(df[[y]])) ||
      any(!is.finite(df[[age]]) | df[[age]] <= 0) || any(!is.finite(df[[y]]))) {
    stop("Ages must be positive and finite. Scores must be finite numeric values.")
  }
  if (length(grainsize) != 1L || !is.numeric(grainsize) || !is.finite(grainsize) ||
      grainsize < 1 || grainsize > .Machine$integer.max || grainsize != floor(grainsize)) {
    stop("grainsize must be a positive integer.")
  }

  athlete_ids <- sort(unique(df[[id]]))
  record_ids <- match(df[[id]], athlete_ids)
  first_row <- match(seq_along(athlete_ids), record_ids)
  athlete_df <- df[first_row, , drop = FALSE]

  one_sided <- function(formula) {
    formula <- stats::as.formula(formula)
    if (length(formula) != 2L || "." %in% all.names(formula)) {
      stop("Use one-sided formulas with explicitly named covariates.")
    }
    vars <- all.vars(formula)
    if (!all(vars %in% names(df))) {
      stop("Covariates missing from df: ", paste(setdiff(vars, names(df)), collapse = ", "), ".")
    }
    formula
  }
  formulas <- lapply(formulas, one_sided)
  observation_formula <- one_sided(observation_formula)
  athlete_vars <- unique(c(unlist(lapply(formulas, all.vars), use.names = FALSE), peak_sd_group))
  for (name in athlete_vars) {
    value <- df[[name]]
    if (!is.atomic(value) || !is.null(dim(value)) || anyNA(value) ||
        any(value != value[first_row[record_ids]])) {
      stop("Athlete-level covariate '", name, "' must be complete and constant within athlete.")
    }
  }
  frame <- function(formula, data) {
    stats::model.frame(formula, data, na.action = stats::na.fail, drop.unused.levels = TRUE)
  }
  validate_design <- function(x, rows, label) {
    if (nrow(x) != rows || any(!is.finite(x)) || anyDuplicated(colnames(x))) {
      stop("Invalid design matrix for ", label, ".")
    }
    if (ncol(x) && qr(x)$rank != ncol(x)) {
      stop("Design matrix for ", label, " is not full rank.")
    }
    x
  }
  X <- lapply(names(formulas), function(component) {
    formula <- formulas[[component]]
    x <- stats::model.matrix(formula, frame(formula, athlete_df))
    if (!ncol(x)) stop("The ", component, " formula must produce at least one column.")
    validate_design(x, nrow(athlete_df), component)
  })
  names(X) <- names(formulas)

  observation_frame <- frame(observation_formula, df)
  observation_terms <- attr(observation_frame, "terms")
  observation_matrix <- stats::model.matrix(observation_terms, observation_frame)
  observation_contrasts <- attr(observation_matrix, "contrasts")
  # Record character levels as well as factor levels to reconstruct this matrix.
  categorical <- vapply(observation_frame, function(x) is.factor(x) || is.character(x),
                        logical(1))
  observation_xlevels <- lapply(observation_frame[categorical], function(x) levels(as.factor(x)))
  Z <- observation_matrix[, colnames(observation_matrix) != "(Intercept)", drop = FALSE]
  Z <- validate_design(Z, nrow(df), "observation")
  if (qr(cbind(1, Z))$rank != ncol(Z) + 1L) {
    stop("The observation design must not include a constant. Use reference coding for factors.")
  }

  group <- if (is.null(peak_sd_group)) factor(rep("All", length(athlete_ids))) else
    factor(athlete_df[[peak_sd_group]])
  stan_data <- list(
    I = length(athlete_ids), N = nrow(df), K = as.integer(vapply(X, ncol, integer(1))),
    X_alpha = unname(X$alpha), X_L = unname(X$L), X_R = unname(X$R), X_tau = unname(X$tau),
    P = ncol(Z), Z = unname(Z), y = as.numeric(df[[y]]), age = as.numeric(df[[age]]),
    id = record_ids, G = nlevels(group), group = as.integer(group),
    grainsize = as.integer(grainsize)
  )
  list(
    stan_data = stan_data, athlete_ids = athlete_ids, group_levels = levels(group),
    design = c(lapply(X, colnames), list(observation = as.character(colnames(Z)))),
    observation_design = Z, observation_formula = observation_formula,
    observation_terms = observation_terms, observation_xlevels = observation_xlevels,
    observation_contrasts = observation_contrasts
  )
}

# Fixed prior locations and scales are passed as Stan data.
.peakage_extended_priors <- function(prior, design) {
  means <- c(alpha = .5, L = log(.25), R = log(.25), tau = log(28))
  intercept_sds <- c(alpha = 1, L = 1.5, R = 1.5, tau = .5)
  slope_sds <- c(alpha = 1, L = .5, R = .5, tau = .1)
  components <- names(means)
  defaults <- list(
    beta_mean = stats::setNames(lapply(components, function(k) {
      ifelse(design[[k]] == "(Intercept)", means[[k]], 0)
    }), components),
    beta_sd = stats::setNames(lapply(components, function(k) {
      ifelse(design[[k]] == "(Intercept)", intercept_sds[[k]], slope_sds[[k]])
    }), components),
    delta_sd = rep(1, length(design$observation)),
    scale_sd = c(alpha = 1, L = 1, R = 1, tau = .5),
    residual_scale_sd = .5, residual_heterogeneity_sd = .5, lkj_eta = 2,
    error_shape_sd = 3, error_nu_shape = 2, error_nu_rate = .1
  )
  check_names <- function(x, allowed, label) {
    if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
        anyDuplicated(names(x)) || !all(names(x) %in% allowed)))) {
      stop("Invalid ", label, ". Allowed names: ", paste(allowed, collapse = ", "), ".")
    }
  }
  if (is.null(prior)) prior <- list()
  check_names(prior, names(defaults), "extended prior")
  values <- function(x, keys, label, positive = TRUE) {
    if (!is.numeric(x) || any(!is.finite(x)) || (positive && any(x <= 0))) {
      stop(label, " must contain finite", if (positive) ", positive" else "", " values.")
    }
    if (!is.null(names(x))) {
      if (length(x) != length(keys) || anyDuplicated(names(x)) ||
          !setequal(names(x), keys)) stop(label, " names must match the design columns.")
      x <- x[keys]
    } else if (length(x) == 1L) {
      x <- rep(x, length(keys))
    } else if (length(x) != length(keys)) {
      stop(label, " must be scalar or have ", length(keys), " entries.")
    }
    as.numeric(x)
  }
  out <- defaults
  for (type in c("beta_mean", "beta_sd")) {
    overrides <- prior[[type]]
    if (is.null(overrides)) next
    check_names(overrides, components, type)
    for (k in names(overrides)) {
      out[[type]][[k]] <- values(overrides[[k]], design[[k]], paste(type, k),
                                positive = type == "beta_sd")
    }
  }
  keys <- list(delta_sd = design$observation, scale_sd = components)
  for (name in setdiff(names(prior), c("beta_mean", "beta_sd"))) {
    out[[name]] <- values(prior[[name]], if (name %in% names(keys)) keys[[name]] else name, name)
  }
  out
}

.fit_peakage_extended <- function(df, peak_formula, p_formula, left_formula, right_formula,
                                  observation_formula, peak_sd_group, prior, id, age, y,
                                  grainsize, stan_file, ...) {
  formulas <- list(alpha = p_formula, L = left_formula, R = right_formula, tau = peak_formula)
  prepared <- .peakage_extended_data(df, formulas, observation_formula, peak_sd_group,
                                     id, age, y, grainsize)
  prior <- .peakage_extended_priors(prior, prepared$design)
  stan_data <- prepared$stan_data
  for (component in names(formulas)) {
    stan_data[[paste0("beta_prior_mean_", component)]] <- prior$beta_mean[[component]]
    stan_data[[paste0("beta_prior_sd_", component)]] <- prior$beta_sd[[component]]
  }
  stan_data$delta_prior_sd <- prior$delta_sd
  stan_data$scale_prior_sd <- unname(prior$scale_sd)
  stan_data$residual_scale_prior_sd <- prior$residual_scale_sd
  stan_data$residual_heterogeneity_prior_sd <- prior$residual_heterogeneity_sd
  stan_data$lkj_eta <- prior$lkj_eta
  stan_data$error_shape_prior_sd <- prior$error_shape_sd
  stan_data$error_nu_prior_shape <- prior$error_nu_shape
  stan_data$error_nu_prior_rate <- prior$error_nu_rate
  if (is.null(stan_file)) {
    stan_file <- system.file("stan", "peakage_extended.stan", package = "peakAge")
    if (!nzchar(stan_file)) stan_file <- file.path("inst", "stan", "peakage_extended.stan")
  }
  if (!file.exists(stan_file)) stop("Stan file not found at: ", stan_file)
  model <- .peakage_get_model(stan_file)
  stan_data_cmdstan <- stan_data
  vector_fields <- c("id", "group", "age", "y", "delta_prior_sd",
    paste0("beta_prior_mean_", names(formulas)), paste0("beta_prior_sd_", names(formulas)))
  for (name in vector_fields) {
    if (length(stan_data_cmdstan[[name]]) == 1L) {
      stan_data_cmdstan[[name]] <- array(stan_data_cmdstan[[name]], dim = 1L)
    }
  }
  args <- list(...)
  if (!"init" %in% names(args)) {
    args$init <- function() {
      list(
        b_alpha = array(prior$beta_mean$alpha), b_L = array(prior$beta_mean$L),
        b_R = array(prior$beta_mean$R), b_tau = array(prior$beta_mean$tau),
        delta = array(rep(0, stan_data$P)),
        sigma_alpha_group = array(rep(prior$scale_sd[1] / 2, stan_data$G)),
        sigma_L = prior$scale_sd[2] / 2, sigma_R = prior$scale_sd[3] / 2,
        sigma_tau = prior$scale_sd[4] / 2, L_Omega = diag(4),
        z_b = matrix(0, nrow = 4, ncol = stan_data$I),
        sigma_y = max(2e-6, prior$residual_scale_sd / 2),
        omega_y = prior$residual_heterogeneity_sd / 2, z_sigma = array(rep(0, stan_data$I)),
        error_shape = 0, error_nu = max(4, prior$error_nu_shape / prior$error_nu_rate)
      )
    }
  }
  fit <- do.call(model$sample, c(list(data = stan_data_cmdstan), args))
  prepared$stan_data <- stan_data
  structure(c(prepared, list(
    fit = fit, model = "extended", curve = "asymmetric", prior = prior,
    peak_formula = peak_formula, p_formula = p_formula, left_formula = left_formula,
    right_formula = right_formula, peak_sd_group = peak_sd_group,
    id_col = id, age_col = age, y_col = y, stan_file = stan_file,
    call = match.call()
  )), class = "peakage_fit")
}
