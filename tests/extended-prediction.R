library(peakAge)

if (requireNamespace("posterior", quietly = TRUE)) {
  # CmdStanR returns a draws_matrix. Its single-entry slices can retain
  # dimensions, which must not be carried into the residual generator.
  draws <- cbind(
    "p_hat[1]" = c(0.8, 0.9, 0.85),
    "logkL_hat[1]" = log(c(0.01, 0.02, 0.015)),
    "logkR_hat[1]" = log(c(0.02, 0.03, 0.025)),
    "T_hat[1]" = c(26, 27, 28),
    "sigma_y_i[1]" = c(0.03, 0.04, 0.05),
    "error_shape" = c(0.7, 0.3, 0.5),
    "error_nu" = c(8, 12, 10)
  )
  fit <- list(
    fit = structure(list(draws = function(variables, format) {
      posterior::as_draws_matrix(draws[, variables, drop = FALSE])
    }), class = "CmdStanMCMC"),
    athlete_ids = "A1", model = "extended", curve = "asymmetric",
    design = list(observation = character())
  )
  records <- data.frame(Athlete_id = "A1", age = c(24, 27, 30), y = c(0.7, 0.9, 0.65))
  set.seed(124)
  plot <- plot_athlete_curve(fit, records, "A1", band_type = "pred")
  stopifnot(inherits(plot, "ggplot"))
  ribbon <- ggplot2::ggplot_build(plot)$data[[1L]]
  stopifnot(nrow(ribbon) == 200L,
            all(is.finite(ribbon$ymin)), all(is.finite(ribbon$ymax)),
            all(ribbon$ymin < ribbon$ymax))
}
