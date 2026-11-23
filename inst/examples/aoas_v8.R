# Source this file, then call prepare_aoas_v8() with your own copy of the
# v8 application directory. It requires config.R, prepare_data.R,
# prior_profiles.R and the original data/data_clean.csv. The canonical data
# constructor verifies the cohort checksum before preparing the full dataset.
#
# Example (with peakAge installed):
#   source(system.file("examples", "aoas_v8.R", package = "peakAge"))
#   prepared <- prepare_aoas_v8("path/to/v8", profile = "primary")
#   fit <- fit_aoas_v8("path/to/v8", profile = "primary", parallel_chains = 2)
# Use profile = "weak_zero" for the wider-prior analysis.
#
# The package fit uses the same scientific likelihood and priors. The package
# samples fixed-effect coefficients directly. The original workflow samples
# standardized coordinates for the interaction and within-sport event effects.
# The coordinates induce the same fixed Normal priors on the effects. Starts
# and adaptation can differ, so draws need not match the original chains.

prepare_aoas_v8 <- function(root, profile = c("primary", "weak_zero")) {
  profile <- match.arg(profile)
  root <- normalizePath(root, mustWork = TRUE)
  workflow <- new.env(parent = baseenv())
  workflow$read.csv <- utils::read.csv
  workflow$setNames <- stats::setNames
  workflow$ave <- stats::ave
  for (file in c("config.R", "prior_profiles.R", "prepare_data.R")) {
    sys.source(file.path(root, file), envir = workflow)
  }
  cfg <- workflow$cfg
  if (!identical(cfg$protocol, "peakage_Mgroup_two_priors_v8_split_t_vectorized_ncp") ||
      !identical(cfg$scientific_model, "Mgroup") || isTRUE(cfg$fixture)) {
    stop("Use the original full-data v8 Mgroup workflow directory.")
  }
  prepared <- workflow$prepare_two_prior_data(root)
  priors <- workflow$peakage_prior_profile(profile)
  d <- prepared$stan_data

  # Match the original sport reference, sex and birth-year columns, sex-by-
  # sport interactions, and full-cohort centered within-sport event contrasts.
  sport_index <- seq.int(2L, d$S)
  event_index <- seq_len(d$K_event)
  indoor_index <- seq_len(d$K_indoor)
  sport <- 1 * outer(d$sport_id, sport_index, "==")
  colnames(sport) <- paste0("sport_", sport_index)
  interaction <- sport * d$female
  colnames(interaction) <- paste0("female_sport_", sport_index)
  event <- d$X_event
  colnames(event) <- paste0("event_", event_index)
  X <- cbind("(Intercept)" = 1, sport, female = d$female,
             birth_year = d$birth_year, interaction, event)
  Z <- d$Z_indoor
  colnames(Z) <- paste0("indoor_", indoor_index)
  athlete_formula <- stats::reformulate(colnames(X)[-1L])
  observation_formula <- stats::reformulate(colnames(Z), intercept = FALSE)

  records <- data.frame(athlete_index = d$id, age = d$age, y = d$y,
                        peak_group = d$alpha_group[d$id])
  records <- cbind(records, X[d$id, -1L, drop = FALSE], Z)
  athlete_rows <- match(seq_len(d$I), records$athlete_index)
  same_design <- function(actual, expected) {
    identical(dim(actual), dim(expected)) && identical(colnames(actual), colnames(expected)) &&
      identical(as.vector(actual), as.vector(expected))
  }
  stopifnot(
    same_design(stats::model.matrix(athlete_formula, records[athlete_rows, ]), X),
    same_design(stats::model.matrix(observation_formula, records), Z),
    identical(records$peak_group[athlete_rows], d$alpha_group),
    qr(X)$rank == ncol(X), qr(cbind(1, Z))$rank == ncol(Z) + 1L
  )

  # Map each design column to its original scientific coefficient. Locations
  # and SDs come from the canonical table, including both prior profiles.
  coefficient_keys <- function(intercept, sport, sex, birth, interaction, event) {
    stats::setNames(c(intercept, paste0(sport, "[", sport_index, "]"), sex, birth,
      paste0(interaction, "[", sport_index, "]"),
      paste0("delta_event_", event, "[", event_index, "]")), colnames(X))
  }
  keys <- list(
    alpha = coefficient_keys("alpha0", "eta_alpha", "beta_S_alpha", "beta_Y_alpha",
                             "eta_alphaS", "alpha"),
    L = coefficient_keys("logrho_L0", "eta_L", "beta_S_L", "beta_Y_L", "eta_LS", "L"),
    R = coefficient_keys("logrho_R0", "eta_R", "beta_S_R", "beta_Y_R", "eta_RS", "R"),
    tau = coefficient_keys("beta0", "eta_tau", "beta_S", "beta_Y", "eta_tauS", "tau")
  )
  table <- priors$table
  lookup <- function(key, field, family) {
    index <- match(unname(key), table$variable)
    stopifnot(!anyNA(index), all(table$family[index] == family))
    stats::setNames(table[[field]][index], names(key))
  }
  indoor_keys <- stats::setNames(paste0("delta_indoor[", indoor_index, "]"), colnames(Z))
  height_scales <- lookup(paste0("sigma_alpha_group[", seq_len(d$G_alpha), "]"),
                          "scale", "half_normal")
  stopifnot(all(height_scales == height_scales[1L]),
            all(lookup(indoor_keys, "location", "normal") == 0))
  scale_keys <- c(L = "sigma_L", R = "sigma_R", tau = "sigma_tau")
  prior <- list(
    beta_mean = lapply(keys, lookup, field = "location", family = "normal"),
    beta_sd = lapply(keys, lookup, field = "scale", family = "normal"),
    delta_sd = lookup(indoor_keys, "scale", "normal"),
    scale_sd = c(alpha = height_scales[1L], lookup(scale_keys, "scale", "half_normal")),
    residual_scale_sd = lookup("sigma_y", "scale", "half_normal"),
    residual_heterogeneity_sd = lookup("omega_y", "scale", "half_normal"),
    lkj_eta = lookup("Omega", "scale", "LKJ"),
    error_shape_sd = lookup("error_shape", "scale", "normal"),
    error_nu_shape = lookup("error_nu1", "location", "gamma_truncated"),
    error_nu_rate = lookup("error_nu1", "scale", "gamma_truncated")
  )
  sampling <- cfg$sampling[c("chains", "iter_warmup", "iter_sampling", "adapt_delta",
                            "max_treedepth", "step_size", "metric",
                            "threads_per_chain", "refresh")]
  sampling$seed <- unname(cfg$sampling$seeds[[profile]])
  list(
    fit_args = list(df = records, id = "athlete_index", curve = "asymmetric",
                    model = "extended", peak_formula = athlete_formula,
                    p_formula = athlete_formula, left_formula = athlete_formula,
                    right_formula = athlete_formula, observation_formula = observation_formula,
                    peak_sd_group = "peak_group", grainsize = cfg$grainsize, prior = prior),
    sampling = sampling, profile = profile, input_md5 = prepared$input_md5,
    athlete_ids = prepared$athletes[, c("id", "Athlete_id")],
    peak_group_map = prepared$height_group_map,
    coefficients = c(keys, list(observation = indoor_keys))
  )
}

fit_aoas_v8 <- function(root, profile = c("primary", "weak_zero"), ...) {
  prepared <- prepare_aoas_v8(root, match.arg(profile))
  controls <- list(...)
  if (length(controls) && (is.null(names(controls)) || any(!nzchar(names(controls))) ||
      anyDuplicated(names(controls)) ||
      any(names(controls) %in% c(names(prepared$fit_args), "stan_file", "formula", "data")))) {
    stop("Supply uniquely named sampling controls in ... .")
  }
  sampling <- utils::modifyList(prepared$sampling, controls, keep.null = TRUE)
  do.call(peakAge::fit_peakage, c(prepared$fit_args, sampling))
}
