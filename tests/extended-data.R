library(peakAge)

prepare <- getFromNamespace(".peakage_extended_data", "peakAge")
make_priors <- getFromNamespace(".peakage_extended_priors", "peakAge")
expect_error <- function(expr, text) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(text, message, fixed = TRUE))
}

athletes <- data.frame(
  Athlete_id = c("A", "B", "C", "D"),
  event = factor(c("track", "track", "field", "field")),
  cohort = c(-0.5, 0.5, -0.5, 0.5)
)
records <- athletes[rep(seq_len(4L), each = 3L), ]
records$age <- rep(c(23, 27, 31), times = 4L)
records$indoor <- rep(c(0, 1, 0), times = 4L)
records$y <- 0.8 - 0.01 * (records$age - 27)^2
records <- records[c(12, 1, 5, 8, 2, 11, 3, 10, 9, 7, 6, 4), ]
formulas <- list(alpha = ~ event, L = ~ cohort, R = ~ 1, tau = ~ event + cohort)
build <- function(df = records, observation_formula = ~ indoor) {
  prepare(df, formulas, observation_formula, "event", "Athlete_id", "age", "y", 2L)
}
prepared <- build()
d <- prepared$stan_data
stopifnot(
  identical(prepared$athlete_ids, athletes$Athlete_id),
  identical(d$id, match(records$Athlete_id, athletes$Athlete_id)),
  identical(unname(d$K), c(2L, 2L, 1L, 3L)),
  identical(d$group, as.integer(athletes$event)),
  identical(d$Z[, 1L], records$indoor),
  identical(d$X_alpha, unname(model.matrix(~ event, athletes))),
  identical(d$X_L, unname(model.matrix(~ cohort, athletes))),
  identical(d$X_R, unname(model.matrix(~ 1, athletes))),
  identical(d$X_tau, unname(model.matrix(~ event + cohort, athletes)))
)
no_adjustment <- build(observation_formula = ~ 0)
stopifnot(identical(dim(no_adjustment$stan_data$Z), c(12L, 0L)))

alpha_columns <- prepared$design$alpha
alpha_mean <- setNames(c(0.8, 0.1), alpha_columns)
prior <- make_priors(list(
  beta_mean = list(alpha = rev(alpha_mean)),
  beta_sd = list(tau = 0.2),
  scale_sd = c(tau = 0.1, R = 0.4, L = 0.3, alpha = 0.2)
), prepared$design)
stopifnot(
  identical(prior$beta_mean$alpha, unname(alpha_mean)),
  identical(prior$beta_sd$tau, rep(0.2, 3L)),
  identical(prior$scale_sd, c(0.2, 0.3, 0.4, 0.1))
)
expect_error(make_priors(list(beta_sd = list(tau = 0)), prepared$design), "positive")
expect_error(make_priors(list(unknown = 1), prepared$design), "Invalid extended prior")

bad <- records
bad$cohort[1L] <- 2
expect_error(build(bad), "constant within athlete")
bad <- records
bad$indoor[1L] <- NA_real_
expect_error(build(bad), "missing values")
bad <- records
bad$constant <- 1
expect_error(build(bad, ~ constant), "must not include a constant")
