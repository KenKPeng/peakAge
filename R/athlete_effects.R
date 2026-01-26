#' Athlete-level posterior estimates from a peakage fit
#'
#' Returns posterior point estimates and 95% credible intervals
#' for athlete-specific quantities from Stan generated quantities.
#'
#' @param fit A `peakage_fit` object.
#' @param ids Optional vector of athlete IDs (original IDs, e.g. "A2").
#'   If NULL, returns all athletes.
#' @param point "mean" or "median".
#' @param include_ci Logical; include 95% CI if TRUE.
#' @param include_p Logical; include p if TRUE.
#' @param include_k Logical; include k (or kL/kR) if TRUE.
#'
#' @return A data.frame with one row per athlete.
#' @export
athlete_estimates <- function(fit,
                              ids = NULL,
                              point = c("mean", "median"),
                              include_ci = TRUE,
                              include_p = TRUE,
                              include_k = TRUE) {
  point <- match.arg(point)

  if (!is.list(fit) || is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    stop("`fit` must be a peakage_fit with CmdStanMCMC.")
  }
  if (is.null(fit$athlete_ids)) {
    stop("`fit$athlete_ids` missing. Please refit using updated fit_peakage().")
  }

  curve <- fit$curve %||% "symmetric"
  athlete_ids_all <- fit$athlete_ids
  I <- length(athlete_ids_all)

  # ---- select athletes by original Athlete_id ----
  if (is.null(ids)) {
    keep_i <- seq_len(I)
  } else {
    keep_i <- match(as.character(ids), as.character(athlete_ids_all))
    keep_i <- keep_i[!is.na(keep_i)]
    if (length(keep_i) == 0) stop("No `ids` matched athlete IDs in fit.")
  }

  vars <- if (curve == "symmetric") {
    c("T_hat", "p_hat", "logk_hat")
  } else {
    c("T_hat", "p_hat", "logkL_hat", "logkR_hat")
  }

  # ---- pull draws and compute summaries (robust 95% CI) ----
  dr <- fit$fit$draws(variables = vars, format = "draws_matrix")
  cn <- colnames(dr)

  idx <- .parse_indexed_param(cn)
  ok <- idx$i %in% keep_i & idx$i >= 1 & idx$i <= I & idx$base %in% vars

  dr <- dr[, ok, drop = FALSE]
  base <- idx$base[ok]
  ii   <- idx$i[ok]

  pt <- if (point == "mean") colMeans(dr) else apply(dr, 2, stats::median)
  if (include_ci) {
    lo <- apply(dr, 2, stats::quantile, probs = 0.025, names = FALSE)
    hi <- apply(dr, 2, stats::quantile, probs = 0.975, names = FALSE)
  }

  long <- data.frame(
    i = ii,
    base = base,
    est = as.numeric(pt),
    stringsAsFactors = FALSE
  )
  if (include_ci) {
    long$lo_95 <- as.numeric(lo)
    long$hi_95 <- as.numeric(hi)
  }

  out <- data.frame(
    Athlete_id = athlete_ids_all[keep_i],
    stringsAsFactors = FALSE
  )

  add_family <- function(base_name, prefix) {
    sub <- long[long$base == base_name, , drop = FALSE]
    if (nrow(sub) == 0) return(invisible(NULL))
    sub <- sub[order(sub$i), , drop = FALSE]
    m <- match(keep_i, sub$i)

    out[[paste0(prefix, "_", point)]] <<- sub$est[m]
    if (include_ci) {
      out[[paste0(prefix, "_lo_95")]] <<- sub$lo_95[m]
      out[[paste0(prefix, "_hi_95")]] <<- sub$hi_95[m]
    }
  }

  add_family("T_hat", "T")
  if (include_p) add_family("p_hat", "p")

  if (include_k) {
    if (curve == "symmetric") {
      add_family("logk_hat", "logk_tmp")
      out[[paste0("k_", point)]] <- exp(out[[paste0("logk_tmp_", point)]])
      if (include_ci) {
        out[["k_lo_95"]] <- exp(out[["logk_tmp_lo_95"]])
        out[["k_hi_95"]] <- exp(out[["logk_tmp_hi_95"]])
      }
      out[grep("^logk_tmp_", names(out))] <- NULL
    } else {
      add_family("logkL_hat", "logkL_tmp")
      add_family("logkR_hat", "logkR_tmp")
      out[[paste0("kL_", point)]] <- exp(out[[paste0("logkL_tmp_", point)]])
      out[[paste0("kR_", point)]] <- exp(out[[paste0("logkR_tmp_", point)]])
      if (include_ci) {
        out[["kL_lo_95"]] <- exp(out[["logkL_tmp_lo_95"]])
        out[["kL_hi_95"]] <- exp(out[["logkL_tmp_hi_95"]])
        out[["kR_lo_95"]] <- exp(out[["logkR_tmp_lo_95"]])
        out[["kR_hi_95"]] <- exp(out[["logkR_tmp_hi_95"]])
      }
      out[grep("^logk[LR]_tmp_", names(out))] <- NULL
    }
  }

  rownames(out) <- NULL
  out
}

# ---- helper ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

.parse_indexed_param <- function(x) {
  list(
    base = sub("\\[.*$", "", x),
    i = suppressWarnings(as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", x)))
  )
}
