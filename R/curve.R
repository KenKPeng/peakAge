#' peak curve mean function
#'
#' Computes the mean performance curve under a symmetric quadratic peak model:
#' mu(age) = p - exp(logk) * (age - T)^2.
#'
#' @param age Numeric vector of ages.
#' @param p Numeric vector (or scalar): performance level at the peak.
#' @param logk Numeric vector (or scalar): log curvature parameter (k = exp(logk) > 0).
#' @param T Numeric vector (or scalar): peak age.
#' @return Numeric vector of mu with length equal to the recycled common length of inputs.
#' @export
curve_mu_sym <- function(age, p, logk, T) {
  .check_recycle_ok(age, p, logk, T)

  age  <- as.numeric(age)
  p    <- as.numeric(p)
  logk <- as.numeric(logk)
  T    <- as.numeric(T)

  if (any(!is.finite(age)))  stop("`age` must be finite.")
  if (any(!is.finite(p)))    stop("`p` must be finite.")
  if (any(!is.finite(logk))) stop("`logk` must be finite.")
  if (any(!is.finite(T)))    stop("`T` must be finite.")

  if (any(logk > 50)) {
    warning("Some `logk` values are very large; exp(logk) may overflow.")
  }

  k <- exp(logk)
  p - k * (age - T)^2
}

# ---- internal helpers ----

.check_recycle_ok <- function(...) {
  args <- list(...)
  lens <- vapply(args, length, integer(1))
  mx <- max(lens)
  if (!all(lens %in% c(1L, mx))) {
    stop(
      "Lengths must be 1 or equal after recycling. Got lengths: ",
      paste(lens, collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' Peak summary for symmetric quadratic curve
#'
#' @param p,logk,T As in [curve_mu_sym()].
#' @return A data.frame with columns: T, peak_mu, k, second_deriv.
#' @export
curve_peak_sym <- function(p, logk, T) {
  .check_recycle_ok(p, logk, T)

  p    <- as.numeric(p)
  logk <- as.numeric(logk)
  T    <- as.numeric(T)

  if (any(!is.finite(p)) || any(!is.finite(logk)) || any(!is.finite(T))) {
    stop("`p`, `logk`, `T` must be finite.")
  }

  if (any(logk > 50)) {
    warning("Some `logk` values are very large; exp(logk) may overflow.")
  }

  k <- exp(logk)
  data.frame(
    T = T,
    peak_mu = p,          # mu(T) = p
    k = k,
    second_deriv = -2 * k # d2/dage2 mu = -2k
  )
}

#' Predictive distribution under symmetric curve
#'
#' @param age,p,logk,T As in [curve_mu_sym()].
#' @param sigma_y Observation noise SD (>=0).
#' @param n_draw Optional integer. If provided, simulate y_rep with `n_draw` draws per age.
#' @return A list with mu, sigma_y, and optionally y_rep (matrix of dim length(mu) x n_draw).
#' @export
curve_predict_sym <- function(age, p, logk, T, sigma_y, n_draw = NULL) {
  sigma_y <- as.numeric(sigma_y)
  if (length(sigma_y) != 1L || !is.finite(sigma_y) || sigma_y < 0) {
    stop("`sigma_y` must be a single non-negative finite number.")
  }

  mu <- curve_mu_sym(age = age, p = p, logk = logk, T = T)
  out <- list(mu = mu, sigma_y = sigma_y)

  if (!is.null(n_draw)) {
    n_draw <- as.integer(n_draw)
    if (length(n_draw) != 1L || is.na(n_draw) || n_draw <= 0) {
      stop("`n_draw` must be a positive integer.")
    }

    # y_rep[i, ] are draws for age[i]
    out$y_rep <- matrix(
      stats::rnorm(length(mu) * n_draw, mean = rep(mu, each = n_draw), sd = sigma_y),
      nrow = n_draw,
      ncol = length(mu),
      byrow = TRUE
    )
    out$y_rep <- t(out$y_rep)  # length(mu) x n_draw
  }

  out
}

#' Evaluate symmetric curves on an age grid for multiple athletes
#'
#' @param age_grid Numeric vector of ages to evaluate.
#' @param params A data.frame with columns: id (optional), p, logk, T.
#' @return A data.frame with columns: id (if present), age, mu.
#' @export
curve_grid_sym <- function(age_grid, params) {
  age_grid <- as.numeric(age_grid)
  if (any(!is.finite(age_grid))) stop("`age_grid` must be finite.")

  if (!is.data.frame(params)) stop("`params` must be a data.frame.")
  need <- c("p", "logk", "T")
  if (!all(need %in% names(params))) {
    stop("`params` must contain columns: ", paste(need, collapse = ", "))
  }

  p_col    <- as.numeric(params$p)
  logk_col <- as.numeric(params$logk)
  T_col    <- as.numeric(params$T)

  if (any(!is.finite(p_col)) || any(!is.finite(logk_col)) || any(!is.finite(T_col))) {
    stop("`params$p`, `params$logk`, `params$T` must be finite.")
  }

  has_id <- "id" %in% names(params)
  ids <- if (has_id) params$id else seq_len(nrow(params))

  G <- length(age_grid)
  J <- nrow(params)

  age_out  <- rep(age_grid, times = J)
  p_out    <- rep(p_col,    each = G)
  logk_out <- rep(logk_col, each = G)
  T_out    <- rep(T_col,    each = G)

  mu_out <- curve_mu_sym(age_out, p_out, logk_out, T_out)

  out <- data.frame(age = age_out, mu = mu_out)
  if (has_id) {
    out$id <- rep(ids, each = G)
    out <- out[, c("id", "age", "mu")]
  }
  out
}




#' Asymmetric peak curve mean function (piecewise quadratic)
#'
#' Computes the mean performance curve under an asymmetric quadratic peak model:
#' mu(age) = p - exp(logkL) * (T - age)^2  for age <= T
#'         = p - exp(logkR) * (age - T)^2  for age >  T
#'
#' @param age Numeric vector of ages.
#' @param p Numeric vector (or scalar): performance level at the peak.
#' @param logkL Numeric vector (or scalar): log left curvature parameter.
#' @param logkR Numeric vector (or scalar): log right curvature parameter.
#' @param T Numeric vector (or scalar): peak age.
#' @return Numeric vector of mu with length equal to max(length(age), length(p), length(logkL), length(logkR), length(T)).
#' @export
curve_mu_asym <- function(age, p, logkL, logkR, T) {
  age   <- as.numeric(age)
  p     <- as.numeric(p)
  logkL <- as.numeric(logkL)
  logkR <- as.numeric(logkR)
  T     <- as.numeric(T)

  if (any(!is.finite(age)))   stop("`age` must be finite.")
  if (any(!is.finite(p)))     stop("`p` must be finite.")
  if (any(!is.finite(logkL))) stop("`logkL` must be finite.")
  if (any(!is.finite(logkR))) stop("`logkR` must be finite.")
  if (any(!is.finite(T)))     stop("`T` must be finite.")

  kL <- exp(logkL)
  kR <- exp(logkR)

  # vectorized, recycling rules apply
  d <- age - T
  # left side: age <= T => d <= 0, use kL
  # right side: age > T => d > 0, use kR
  p - ifelse(d <= 0, kL * d^2, kR * d^2)
}

#' Peak summary for asymmetric quadratic curve
#'
#' @param p,logkL,logkR,T As in [curve_mu_asym()].
#' @return A data.frame with columns: T, peak_mu, kL, kR, second_deriv_left, second_deriv_right.
#' @export
curve_peak_asym <- function(p, logkL, logkR, T) {
  p     <- as.numeric(p)
  logkL <- as.numeric(logkL)
  logkR <- as.numeric(logkR)
  T     <- as.numeric(T)

  if (any(!is.finite(p)) || any(!is.finite(logkL)) || any(!is.finite(logkR)) || any(!is.finite(T))) {
    stop("`p`, `logkL`, `logkR`, `T` must be finite.")
  }

  kL <- exp(logkL)
  kR <- exp(logkR)

  data.frame(
    T = T,
    peak_mu = p,                 # mu(T) = p on both sides
    kL = kL,
    kR = kR,
    second_deriv_left  = -2 * kL,
    second_deriv_right = -2 * kR
  )
}

#' Predictive distribution under asymmetric curve
#'
#' @param age,p,logkL,logkR,T As in [curve_mu_asym()].
#' @param sigma_y Observation noise SD (>=0).
#' @param n_draw Optional integer. If provided, simulate y_rep with `n_draw` draws per age.
#' @return A list with mu, sigma_y, and optionally y_rep (matrix).
#' @export
curve_predict_asym <- function(age, p, logkL, logkR, T, sigma_y, n_draw = NULL) {
  sigma_y <- as.numeric(sigma_y)
  if (length(sigma_y) != 1L || !is.finite(sigma_y) || sigma_y < 0) {
    stop("`sigma_y` must be a single non-negative finite number.")
  }

  mu <- curve_mu_asym(age = age, p = p, logkL = logkL, logkR = logkR, T = T)
  out <- list(mu = mu, sigma_y = sigma_y)

  if (!is.null(n_draw)) {
    n_draw <- as.integer(n_draw)
    if (n_draw <= 0) stop("`n_draw` must be a positive integer.")

    out$y_rep <- matrix(
      stats::rnorm(length(mu) * n_draw, mean = rep(mu, times = n_draw), sd = sigma_y),
      nrow = length(mu),
      ncol = n_draw
    )
  }

  out
}

#' Evaluate asymmetric curves on an age grid for multiple athletes
#'
#' @param age_grid Numeric vector of ages to evaluate.
#' @param params A data.frame with columns: id (optional), p, logkL, logkR, T.
#' @return A data.frame with columns: id (if present), age, mu.
#' @export
curve_grid_asym <- function(age_grid, params) {
  age_grid <- as.numeric(age_grid)
  if (any(!is.finite(age_grid))) stop("`age_grid` must be finite.")

  if (!is.data.frame(params)) stop("`params` must be a data.frame.")
  need <- c("p", "logkL", "logkR", "T")
  if (!all(need %in% names(params))) {
    stop("`params` must contain columns: ", paste(need, collapse = ", "))
  }
  if (any(!is.finite(params$p)) ||
      any(!is.finite(params$logkL)) ||
      any(!is.finite(params$logkR)) ||
      any(!is.finite(params$T))) {
    stop("`params$p`, `params$logkL`, `params$logkR`, `params$T` must be finite.")
  }

  has_id <- "id" %in% names(params)
  ids <- if (has_id) params$id else seq_len(nrow(params))

  res_list <- lapply(seq_len(nrow(params)), function(j) {
    mu <- curve_mu_asym(
      age = age_grid,
      p = params$p[j],
      logkL = params$logkL[j],
      logkR = params$logkR[j],
      T = params$T[j]
    )
    data.frame(
      id = ids[j],
      age = age_grid,
      mu = mu
    )
  })

  out <- do.call(rbind, res_list)
  if (!has_id) out$id <- NULL
  out
}

