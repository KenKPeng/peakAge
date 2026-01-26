#' Summarise posterior for a peak-age fit (via CmdStanMCMC$summary)
#'
#' This summary uses cmdstanr's built-in summary table, which is stable and fast.
#' It groups parameters by submodel blocks and optionally renames regression
#' coefficients using the design-matrix column names.
#'
#' Note: prior settings (if stored in the fit object) are kept in the returned object
#' as `res$prior`, but are not printed by default. Please refer to the documentation.
#'
#' @param object A `peakage_fit` object returned by [fit_peakage()].
#' @param probs Numeric vector of probabilities for credible intervals.
#'   Default c(0.05, 0.5, 0.95).
#' @param include_diagnostics Logical; include rhat/ess columns if available.
#' @param include_population Logical; include population/baseline parameters.
#' @param rename_coefs Logical; if TRUE, rename beta rows using design matrix column names.
#' @param format Output format: "data.frame" or "html" (knitr::kable if available).
#' @return An object of class `peakage_summary`.
#' @export
summarise_peakage <- function(object,
                              probs = c(0.05, 0.5, 0.95),
                              include_diagnostics = TRUE,
                              include_population = TRUE,
                              rename_coefs = TRUE,
                              format = c("data.frame", "html")) {
  format <- match.arg(format)

  if (!is.list(object) || is.null(object$fit) || !inherits(object$fit, "CmdStanMCMC")) {
    stop("`object` must be a `peakage_fit` with `object$fit` being a CmdStanMCMC.")
  }

  probs <- as.numeric(probs)
  if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("`probs` must be probabilities strictly between 0 and 1.")
  }
  probs <- sort(unique(probs))

  include_diagnostics <- isTRUE(include_diagnostics)
  include_population  <- isTRUE(include_population)
  rename_coefs        <- isTRUE(rename_coefs)

  curve <- object$curve %||% "symmetric"
  blocks <- .peakage_param_blocks(curve = curve, stan_data = object$stan_data)

  # If user doesn't want population/baseline params, drop them from blocks
  if (!include_population) {
    if (!is.null(blocks$other)) {
      blocks$other <- setdiff(blocks$other, c("p_pop", "logk_pop", "logkL_pop", "logkR_pop", "sigma_y"))
    }
    # optionally hide sigma_T as well when excluding population-level params
    if (!is.null(blocks$AFT_logT)) {
      blocks$AFT_logT <- setdiff(blocks$AFT_logT, "sigma_T")
    }
  }

  out_tables <- list()

  for (bn in names(blocks)) {
    vars <- blocks[[bn]]
    if (length(vars) == 0) next

    tab <- .peakage_cmdstan_summary(
      fit = object$fit,
      variables = vars,
      probs = probs,
      include_diagnostics = include_diagnostics
    )
    if (nrow(tab) == 0) next

    if (rename_coefs) {
      tab$term <- .peakage_pretty_term(
        variable = tab$variable,
        stan_data = object$stan_data
      )
    } else {
      tab$term <- tab$variable
    }

    tab <- .peakage_reorder_summary_cols(tab, probs, include_diagnostics)
    out_tables[[bn]] <- tab
  }

  res <- list(
    curve = curve,
    probs = probs,
    tables = out_tables,
    prior = object$prior %||% NULL,  # keep, but do not print
    call = match.call()
  )
  class(res) <- "peakage_summary"

  if (format == "html") {
    res$html <- .peakage_summary_to_html(res)
  }

  res
}

#' @export
print.peakage_summary <- function(x, ...) {
  cat("peakage_summary\n")
  cat("  curve:", x$curve, "\n")
  cat("  probs:", paste0(x$probs, collapse = ", "), "\n\n")

  if (length(x$tables) == 0) {
    cat("(no tables)\n")
    return(invisible(x))
  }

  for (nm in names(x$tables)) {
    cat("== ", nm, " ==\n", sep = "")
    print(x$tables[[nm]])
    cat("\n")
  }

  invisible(x)
}

# ---- internal helpers ----

`%||%` <- function(a, b) if (!is.null(a)) a else b

.peakage_param_blocks <- function(curve = c("symmetric", "asymmetric"), stan_data = NULL) {
  curve <- match.arg(curve)

  K_T <- if (!is.null(stan_data) && !is.null(stan_data$K_T)) as.integer(stan_data$K_T) else NA_integer_
  K_p <- if (!is.null(stan_data) && !is.null(stan_data$K_p)) as.integer(stan_data$K_p) else NA_integer_
  K_k <- if (!is.null(stan_data) && !is.null(stan_data$K_k)) as.integer(stan_data$K_k) else NA_integer_

  blocks <- list()

  blocks$AFT_logT <- c("beta_T", "sigma_T")

  if (!is.na(K_p) && K_p > 0) {
    blocks$p_submodel <- c("beta_p")
  } else {
    blocks$p_submodel <- character(0)
  }

  if (curve == "symmetric") {
    if (!is.na(K_k) && K_k > 0) {
      blocks$k_submodel <- c("beta_k")
    } else {
      blocks$k_submodel <- character(0)
    }
    blocks$other <- c("p_pop", "logk_pop", "sigma_y")
  } else {
    if (!is.na(K_k) && K_k > 0) {
      blocks$k_submodel <- c("beta_kL", "beta_kR")
    } else {
      blocks$k_submodel <- character(0)
    }
    blocks$other <- c("p_pop", "logkL_pop", "logkR_pop", "sigma_y")
  }

  blocks <- blocks[vapply(blocks, length, integer(1)) > 0]
  blocks
}

.peakage_cmdstan_summary <- function(fit, variables, probs, include_diagnostics) {
  # Prefer passing probs (cmdstanr supports it in newer versions). Fallback if not supported.
  tab <- tryCatch(
    fit$summary(variables = variables, probs = probs),
    error = function(e) NULL
  )
  if (is.null(tab)) {
    tab <- tryCatch(
      fit$summary(variables = variables),
      error = function(e) NULL
    )
  }
  if (is.null(tab)) return(data.frame())

  tab <- as.data.frame(tab)

  if (!include_diagnostics) {
    drop <- intersect(c("rhat", "ess_bulk", "ess_tail"), names(tab))
    if (length(drop) > 0) tab[drop] <- NULL
  }

  tab
}

.peakage_pretty_term <- function(variable, stan_data) {
  v <- variable

  replace_beta <- function(prefix, X, label) {
    if (is.null(X)) return(v)
    cn <- colnames(X)
    if (is.null(cn)) return(v)

    pat <- paste0("^", prefix, "\\[(\\d+)\\]$")
    m <- regexec(pat, v)
    reg <- regmatches(v, m)
    out <- v
    for (i in seq_along(reg)) {
      if (length(reg[[i]]) == 2) {
        idx <- as.integer(reg[[i]][2])
        if (!is.na(idx) && idx >= 1 && idx <= length(cn)) {
          out[i] <- paste0(label, ":", cn[idx])
        }
      }
    }
    out
  }

  v <- replace_beta("beta_T", stan_data$X_T, "AFT")
  v <- replace_beta("beta_p", stan_data$X_p, "p")
  v <- replace_beta("beta_k",  stan_data$X_k, "k")
  v <- replace_beta("beta_kL", stan_data$X_k, "kL")
  v <- replace_beta("beta_kR", stan_data$X_k, "kR")

  v
}

.peakage_reorder_summary_cols <- function(tab, probs, include_diagnostics) {
  qcols <- grep("^q\\d+", names(tab), value = TRUE)

  preferred <- c("term", "variable", "mean", "median", "sd", "mad", qcols)
  if (include_diagnostics) preferred <- c(preferred, "rhat", "ess_bulk", "ess_tail")

  keep <- c(preferred, setdiff(names(tab), preferred))
  tab[, unique(keep), drop = FALSE]
}

.peakage_summary_to_html <- function(x) {
  if (!requireNamespace("knitr", quietly = TRUE)) return(NULL)

  tabs <- x$tables
  if (length(tabs) == 0) return(NULL)

  pieces <- character(0)
  pieces <- c(pieces, sprintf("<h3>peakage summary (curve: %s)</h3>", x$curve))
  pieces <- c(pieces, sprintf("<p>probs: %s</p>", paste0(x$probs, collapse = ", ")))

  for (nm in names(tabs)) {
    pieces <- c(pieces, sprintf("<h4>%s</h4>", nm))
    pieces <- c(pieces, knitr::kable(tabs[[nm]], format = "html"))
  }

  paste(pieces, collapse = "\n")
}

