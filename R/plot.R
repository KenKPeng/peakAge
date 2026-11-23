utils::globalVariables(c("lo", "hi", "mu", "T_hat"))

#' Plot fitted performance curve vs observed data for one or more athletes
#'
#' Draws observed performance points and the fitted posterior curve for selected athletes.
#' Adds a vertical dashed line at the posterior point estimate of peak age.
#' Optionally draws a credible/predictive band.
#'
#' If `athlete_id` is a vector, returns:
#' - a single faceted ggplot if `facet = TRUE`
#' - a named list of ggplot objects if `facet = FALSE`
#'
#' Band types:
#' - band_type = "mu": posterior credible band for latent mean curve mu(age)
#' - band_type = "pred": posterior predictive band for new observations y(age)
#' For extended fits, predictive bands use the athlete's residual SD and the
#' fitted two-piece Student-t error distribution. If the model has competition
#' covariates, curves and bands use a zero competition adjustment. Observed
#' scores are shifted to this reference using the posterior mean adjustment.
#'
#' @param fit A `peakage_fit` returned by [fit_peakage()].
#' @param df Long-format data used for fitting (or compatible), containing id/age/y.
#'   When plotting observed scores from an extended fit, include its competition
#'   covariates as well.
#' @param athlete_id One or more athlete IDs (e.g. "A2" or c("A2","A7")).
#' @param id,age,y Column names in df. If NULL, uses fit$id_col / fit$age_col / fit$y_col,
#'   falling back to "Athlete_id"/"age"/"y".
#' @param age_grid Optional numeric vector of ages to evaluate fitted curve.
#'   If NULL, uses an athlete-specific grid spanning observed ages (with small padding).
#' @param point "mean" or "median" for the fitted curve and peak-age vertical line.
#' @param show_band Logical; if TRUE, draw a band.
#' @param band_level Credible level for band, e.g. 0.95. Default 0.95.
#' @param band_type "mu" (latent mean credible band) or "pred" (posterior predictive band).
#' @param ndraws Optional integer; subsample posterior draws for speed. Default 1000.
#' @param show_points Logical; draw observed points.
#' @param title Optional plot title. If multiple athletes, this is used as the global title
#'   when `facet=TRUE`, otherwise each plot gets an automatic title if NULL.
#' @param facet Logical; if TRUE and multiple athletes, facet into one plot.
#' @param facet_ncol Optional number of columns when faceting.
#' @param gg_opts A named list of convenient ggplot options. Supported keys:
#'   - point_size, point_alpha, point_color
#'   - line_size,  curve_color
#'   - band_alpha, band_fill
#'   - vline_color, vline_linetype
#'   - theme (a ggplot2 theme object, e.g. ggplot2::theme_bw())
#'
#' @return If single athlete: a ggplot. If multiple athletes: ggplot (facet=TRUE) or a named list of ggplots (facet=FALSE).
#' @export

plot_athlete_curve <- function(fit,
                               df,
                               athlete_id,
                               id = NULL,
                               age = NULL,
                               y = NULL,
                               age_grid = NULL,
                               point = c("mean", "median"),
                               show_band = TRUE,
                               band_level = 0.95,
                               band_type = c("mu", "pred"),
                               ndraws = 1000,
                               show_points = TRUE,
                               title = NULL,
                               facet = TRUE,
                               facet_ncol = NULL,
                               gg_opts = list()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.")
  }

  point <- match.arg(point)
  band_type <- match.arg(band_type)

  if (!is.list(fit) || is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    stop("`fit` must be a peakage_fit with CmdStanMCMC.")
  }
  if (is.null(fit$athlete_ids)) {
    stop("`fit$athlete_ids` is missing. Please refit using updated fit_peakage() that stores mapping.")
  }
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")

  # ---- columns ----
  id  <- id  %||% (fit$id_col  %||% "Athlete_id")
  age <- age %||% (fit$age_col %||% "age")
  y   <- y   %||% (fit$y_col   %||% "y")
  if (!all(c(id, age, y) %in% names(df))) {
    stop("`df` must contain columns: ", paste(c(id, age, y), collapse = ", "))
  }

  # ---- gg options ----
  default_opts <- list(
    point_size = 2,
    point_alpha = 0.85,
    point_color = "black",
    line_size = 1,
    curve_color = "black",
    band_alpha = 0.20,
    band_fill = "grey70",
    vline_color = "red",
    vline_linetype = "dashed",
    theme = ggplot2::theme_bw(),
    xlim = NULL,
    ylim = NULL,
    free_x = TRUE,
    free_y = FALSE
  )
  gg_opts <- utils::modifyList(default_opts, gg_opts)

  # ---- athletes to plot ----
  athlete_id <- as.character(athlete_id)
  athlete_id <- unique(athlete_id)

  i_vec <- match(athlete_id, as.character(fit$athlete_ids))
  bad <- is.na(i_vec)
  if (any(bad)) {
    stop("These athlete_id are not found in fit$athlete_ids: ",
         paste(athlete_id[bad], collapse = ", "))
  }

  curve <- fit$curve %||% "symmetric"
  extended <- identical(fit$model, "extended")
  predictive <- isTRUE(show_band) && band_type == "pred"
  adjusted <- extended && length(fit$design$observation) > 0L
  adjusted_points <- adjusted && isTRUE(show_points)
  y_label <- if (adjusted_points) paste0(y, " (adjusted)") else y
  caption <- if (adjusted) {
    if (adjusted_points) {
      "Curves use zero competition adjustment. Scores are adjusted to the same reference."
    } else {
      "Curves use zero competition adjustment."
    }
  }

  # ---- variables from draws ----
  vars <- if (curve == "symmetric") {
    c("T_hat", "p_hat", "logk_hat")
  } else {
    c("T_hat", "p_hat", "logkL_hat", "logkR_hat")
  }
  if (extended) {
    if (predictive) vars <- c(vars, "sigma_y_i")
    vars <- as.vector(outer(vars, i_vec, function(var, i) sprintf("%s[%d]", var, i)))
    if (predictive) vars <- c(vars, "error_shape", "error_nu")
    if (adjusted_points) vars <- c(vars, "delta")
  } else if (predictive) {
    vars <- c(vars, "sigma_y")
  }

  dr <- fit$fit$draws(variables = vars, format = "draws_matrix")
  if (ncol(dr) == 0) stop("No draws returned; check variables exist in fit.")
  if (adjusted_points) {
    delta_names <- sprintf("delta[%d]", seq_along(fit$design$observation))
    delta_mean <- colMeans(dr[, delta_names, drop = FALSE])
  }

  # ---- subsample draws for speed ----
  ndraws <- as.integer(ndraws)
  if (length(ndraws) != 1L || is.na(ndraws) || ndraws <= 0) ndraws <- 1000L
  S <- nrow(dr)
  if (ndraws < S) {
    set.seed(1)
    dr <- dr[sample.int(S, ndraws), , drop = FALSE]
  }
  S <- nrow(dr)

  # helper: grab vector[I] element at index i from draws_matrix
  grab_vec_i <- function(prefix, i) {
    nm <- paste0(prefix, "[", i, "]")
    if (!nm %in% colnames(dr)) stop("Missing draw column: ", nm, " (check Stan generated quantities).")
    dr[, nm]
  }

  # For pred band
  if (predictive && !extended) {
    if (!"sigma_y" %in% colnames(dr)) stop("sigma_y draws missing; cannot draw predictive band.")
    sig_draw <- dr[, "sigma_y"]
  }

  # ---- per-athlete compute: return list(plot_df, band_df, obs_df, T_hat) ----
  build_one <- function(aid, i) {
    df_i <- df[df[[id]] == aid, , drop = FALSE]
    if (nrow(df_i) == 0) stop("No rows found in df for athlete_id = ", aid)

    a_obs <- as.numeric(df_i[[age]])
    y_obs <- as.numeric(df_i[[y]])
    if (any(!is.finite(a_obs)) || any(!is.finite(y_obs))) stop("Observed age/y must be finite for athlete ", aid)
    if (adjusted_points) {
      y_obs <- y_obs - .peakage_observation_shift(fit, df_i, delta_mean)
    }

    # athlete-specific grid unless user supplies a common one
    if (is.null(age_grid)) {
      lo <- min(a_obs); hi <- max(a_obs)
      pad <- 0.5 * (hi - lo + 1e-6)
      lo <- lo - 0.10 * pad
      hi <- hi + 0.10 * pad
      ag <- seq(lo, hi, length.out = 200)
    } else {
      ag <- as.numeric(age_grid)
      if (any(!is.finite(ag))) stop("`age_grid` must be finite.")
    }

    # draws for this athlete
    T_draw <- grab_vec_i("T_hat", i)
    p_draw <- grab_vec_i("p_hat", i)

    G <- length(ag)
    mu_mat <- matrix(NA_real_, nrow = G, ncol = S)

    if (curve == "symmetric") {
      logk_draw <- grab_vec_i("logk_hat", i)
      # Use curve.R implementation (single source of truth)
      for (s in seq_len(S)) {
        mu_mat[, s] <- curve_mu_sym(age = ag, p = p_draw[s], logk = logk_draw[s], T = T_draw[s])
      }
    } else {
      logkL_draw <- grab_vec_i("logkL_hat", i)
      logkR_draw <- grab_vec_i("logkR_hat", i)
      for (s in seq_len(S)) {
        mu_mat[, s] <- curve_mu_asym(age = ag, p = p_draw[s], logkL = logkL_draw[s], logkR = logkR_draw[s], T = T_draw[s])
      }
    }

    # point curve + peak age
    T_hat <- if (point == "mean") mean(T_draw) else stats::median(T_draw)
    mu_point <- if (point == "mean") rowMeans(mu_mat) else apply(mu_mat, 1, stats::median)

    plot_df <- data.frame(athlete_id = aid, age = ag, mu = mu_point)
    obs_df  <- data.frame(athlete_id = aid, age = a_obs, y = y_obs)

    # band (optional)
    band_df <- NULL
    if (isTRUE(show_band)) {
      band_level2 <- as.numeric(band_level)
      if (length(band_level2) != 1L || !is.finite(band_level2) || band_level2 <= 0 || band_level2 >= 1) {
        stop("`band_level` must be in (0,1).")
      }
      alpha <- (1 - band_level2) / 2
      lo_p <- alpha
      hi_p <- 1 - alpha

      if (band_type == "mu") {
        lo <- apply(mu_mat, 1, stats::quantile, probs = lo_p, names = FALSE)
        hi <- apply(mu_mat, 1, stats::quantile, probs = hi_p, names = FALSE)
      } else {
        # posterior predictive
        yrep_mat <- mu_mat
        sigma <- if (extended) grab_vec_i("sigma_y_i", i) else sig_draw
        if (extended) {
          shape <- as.numeric(dr[, "error_shape"])
          nu <- as.numeric(dr[, "error_nu"])
        }
        for (s in seq_len(S)) {
          error <- if (extended) {
            sigma[s] * .peakage_split_t_rng(G, shape[s], nu[s])
          } else {
            stats::rnorm(G, mean = 0, sd = sigma[s])
          }
          yrep_mat[, s] <- mu_mat[, s] + error
        }
        lo <- apply(yrep_mat, 1, stats::quantile, probs = lo_p, names = FALSE)
        hi <- apply(yrep_mat, 1, stats::quantile, probs = hi_p, names = FALSE)
      }
      band_df <- data.frame(athlete_id = aid, age = ag, lo = lo, hi = hi)
    }

    list(plot_df = plot_df, band_df = band_df, obs_df = obs_df, T_hat = T_hat)
  }

  pieces <- Map(build_one, athlete_id, i_vec)

  # ---- single athlete -> single ggplot ----
  if (length(athlete_id) == 1L) {
    one <- pieces[[1]]
    p <- ggplot2::ggplot()

    if (!is.null(one$band_df)) {
      p <- p + ggplot2::geom_ribbon(
        data = one$band_df,
        ggplot2::aes(x = age, ymin = lo, ymax = hi),
        alpha = gg_opts$band_alpha,
        fill  = gg_opts$band_fill
      )
    }

    p <- p + ggplot2::geom_line(
      data = one$plot_df,
      ggplot2::aes(x = age, y = mu),
      linewidth = gg_opts$line_size,
      color = gg_opts$curve_color
    )

    if (isTRUE(show_points)) {
      p <- p + ggplot2::geom_point(
        data = one$obs_df,
        ggplot2::aes(x = age, y = y),
        size = gg_opts$point_size,
        alpha = gg_opts$point_alpha,
        color = gg_opts$point_color
      )
    }

    p <- p + ggplot2::geom_vline(
      xintercept = one$T_hat,
      color = gg_opts$vline_color,
      linetype = gg_opts$vline_linetype
    )

    if (is.null(title)) title <- paste0("Athlete ", athlete_id, " (", curve, ")")
    p <- p + ggplot2::labs(title = title, x = age, y = y_label, caption = caption) + gg_opts$theme
    if (!is.null(gg_opts$xlim) || !is.null(gg_opts$ylim)) {
      p <- p + ggplot2::coord_cartesian(
        xlim = gg_opts$xlim,
        ylim = gg_opts$ylim
      )
    }
    return(p)
  }

  # ---- multiple athletes ----
  if (isTRUE(facet)) {
    plot_df_all <- do.call(rbind, lapply(pieces, `[[`, "plot_df"))
    obs_df_all  <- do.call(rbind, lapply(pieces, `[[`, "obs_df"))
    band_df_all <- if (isTRUE(show_band)) do.call(rbind, lapply(pieces, `[[`, "band_df")) else NULL

    # vline data per athlete (different T_hat)
    vline_df <- data.frame(
      athlete_id = athlete_id,
      T_hat = vapply(pieces, `[[`, numeric(1), "T_hat")
    )

    p <- ggplot2::ggplot()

    if (!is.null(band_df_all)) {
      p <- p + ggplot2::geom_ribbon(
        data = band_df_all,
        ggplot2::aes(x = age, ymin = lo, ymax = hi),
        alpha = gg_opts$band_alpha,
        fill  = gg_opts$band_fill
      )
    }

    p <- p + ggplot2::geom_line(
      data = plot_df_all,
      ggplot2::aes(x = age, y = mu),
      linewidth = gg_opts$line_size,
      color = gg_opts$curve_color
    )

    if (isTRUE(show_points)) {
      p <- p + ggplot2::geom_point(
        data = obs_df_all,
        ggplot2::aes(x = age, y = y),
        size = gg_opts$point_size,
        alpha = gg_opts$point_alpha,
        color = gg_opts$point_color
      )
    }

    p <- p + ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = T_hat),
      color = gg_opts$vline_color,
      linetype = gg_opts$vline_linetype
    )

    free_x <- isTRUE(gg_opts$free_x)
    free_y <- isTRUE(gg_opts$free_y)
    if (!is.null(gg_opts$xlim)) free_x <- FALSE
    if (!is.null(gg_opts$ylim)) free_y <- FALSE

    scales_opt <- if (free_x && free_y) "free" else if (free_x) "free_x" else if (free_y) "free_y" else "fixed"

    p <- p + ggplot2::facet_wrap(~ athlete_id, ncol = facet_ncol, scales = scales_opt)

    if (is.null(title)) title <- paste0("Fitted curves (", curve, ")")
    p <- p + ggplot2::labs(title = title, x = age, y = y_label, caption = caption) + gg_opts$theme

    if (!is.null(gg_opts$xlim) || !is.null(gg_opts$ylim)) {
      p <- p + ggplot2::coord_cartesian(xlim = gg_opts$xlim, ylim = gg_opts$ylim)
    }

    return(p)
  }

  # facet = FALSE -> return list of per-athlete plots
  out_list <- stats::setNames(vector("list", length(athlete_id)), athlete_id)
  for (j in seq_along(athlete_id)) {
    one <- pieces[[j]]
    aid <- athlete_id[j]

    p <- ggplot2::ggplot()

    if (!is.null(one$band_df)) {
      p <- p + ggplot2::geom_ribbon(
        data = one$band_df,
        ggplot2::aes(x = age, ymin = lo, ymax = hi),
        alpha = gg_opts$band_alpha,
        fill  = gg_opts$band_fill
      )
    }

    p <- p + ggplot2::geom_line(
      data = one$plot_df,
      ggplot2::aes(x = age, y = mu),
      size = gg_opts$line_size,
      color = gg_opts$curve_color
    )

    if (isTRUE(show_points)) {
      p <- p + ggplot2::geom_point(
        data = one$obs_df,
        ggplot2::aes(x = age, y = y),
        size = gg_opts$point_size,
        alpha = gg_opts$point_alpha,
        color = gg_opts$point_color
      )
    }

    p <- p + ggplot2::geom_vline(
      xintercept = one$T_hat,
      color = gg_opts$vline_color,
      linetype = gg_opts$vline_linetype
    )

    local_title <- title
    if (is.null(local_title)) local_title <- paste0("Athlete ", aid, " (", curve, ")")
    p <- p + ggplot2::labs(title = local_title, x = age, y = y_label, caption = caption) + gg_opts$theme

    out_list[[j]] <- p
  }

  out_list
}

.peakage_observation_shift <- function(fit, data, delta_mean) {
  frame <- stats::model.frame(fit$observation_terms, data = data,
                              xlev = fit$observation_xlevels, na.action = stats::na.fail)
  design <- stats::model.matrix(fit$observation_terms, data = frame,
                                contrasts.arg = fit$observation_contrasts)
  columns <- fit$design$observation
  if (!all(columns %in% colnames(design)) || any(!is.finite(design))) {
    stop("`df` must match the competition design used in the extended fit.")
  }
  as.vector(design[, columns, drop = FALSE] %*% delta_mean)
}

# Center and standardize the two-piece Student-t used by the extended likelihood.
.peakage_split_t_rng <- function(n, shape, nu) {
  asymmetry <- tanh(shape / 2)
  half_t_mean <- sqrt(nu / pi) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))
  raw_t_variance <- nu / (nu - 2)
  variance <- raw_t_variance + asymmetry^2 * (3 * raw_t_variance - 4 * half_t_mean^2)
  p_left <- stats::plogis(-shape)
  width <- ifelse(stats::runif(n) < p_left, -2 * p_left, 2 * stats::plogis(shape))
  w <- width * abs(stats::rt(n, df = nu))
  (w - 2 * asymmetry * half_t_mean) / sqrt(variance)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
