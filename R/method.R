# ---- small utils ----
`%||%` <- function(a, b) if (!is.null(a)) a else b


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

  # ---- parameter blocks based on stan_data dimensions ----
  blocks <- .peakage_param_blocks(curve = curve, stan_data = object$stan_data)

  # If user doesn't want population/baseline params, drop them from blocks
  if (!include_population) {
    if (!is.null(blocks$other)) {
      blocks$other <- setdiff(blocks$other, c("p_pop", "logk_pop", "logkL_pop", "logkR_pop", "sigma_y"))
    }
    if (!is.null(blocks$AFT_logT)) {
      blocks$AFT_logT <- setdiff(blocks$AFT_logT, "sigma_T")
    }
  }

  # ---- robust filtering: keep only variables that actually exist in the fitted output ----
  avail <- .peakage_available_vars(object$fit)

  out_tables <- list()

  for (bn in names(blocks)) {
    vars <- blocks[[bn]]
    if (length(vars) == 0) next

    # Intersect with available vars. This avoids errors when parameters are absent.
    vars <- vars[vars %in% avail]
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

.peakage_available_vars <- function(fit) {
  # cmdstanr provides $summary() table with a "variable" column.
  tab <- tryCatch(fit$summary(), error = function(e) NULL)
  if (is.null(tab)) return(character(0))
  tab <- as.data.frame(tab)
  if (!("variable" %in% names(tab))) return(character(0))
  unique(tab$variable)
}

.peakage_param_blocks <- function(curve = c("symmetric", "asymmetric"), stan_data = NULL) {
  curve <- match.arg(curve)

  K_T <- if (!is.null(stan_data) && !is.null(stan_data$K_T)) as.integer(stan_data$K_T) else NA_integer_
  K_p <- if (!is.null(stan_data) && !is.null(stan_data$K_p)) as.integer(stan_data$K_p) else NA_integer_
  K_k <- if (!is.null(stan_data) && !is.null(stan_data$K_k)) as.integer(stan_data$K_k) else NA_integer_

  blocks <- list()

  # AFT is always present in our Stan programs
  blocks$AFT_logT <- c("beta_T", "sigma_T")

  # p submodel
  blocks$p_submodel <- if (!is.na(K_p) && K_p > 0) "beta_p" else character(0)

  # k submodel differs by curve type
  if (curve == "symmetric") {
    blocks$k_submodel <- if (!is.na(K_k) && K_k > 0) "beta_k" else character(0)
    blocks$other <- c("p_pop", "logk_pop", "sigma_y")
  } else {
    # asymmetric Stan must define beta_kL and beta_kR if K_k > 0
    blocks$k_submodel <- if (!is.na(K_k) && K_k > 0) c("beta_kL", "beta_kR") else character(0)
    blocks$other <- c("p_pop", "logkL_pop", "logkR_pop", "sigma_y")
  }

  blocks <- blocks[vapply(blocks, length, integer(1)) > 0]
  blocks
}

.peakage_cmdstan_summary <- function(fit, variables, probs, include_diagnostics) {
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

  # Prefer stored names (robust), fallback to colnames(X)
  get_cn <- function(X, nm) {
    if (!is.null(nm) && length(nm) > 0) return(nm)
    cn <- colnames(X)
    if (!is.null(cn)) return(cn)
    NULL
  }

  # Make dummy names more readable:
  #   "sexMale" -> "sex = Male"
  #   "groupB"  -> "group = B"
  #   "x1"      -> "x1" (unchanged)
  #   "(Intercept)" stays as-is
  pretty_col <- function(cn) {
    cn <- as.character(cn)

    out <- cn
    is_int <- cn %in% c("(Intercept)", "Intercept")
    out[is_int] <- "(Intercept)"

    # Try split: leading letters/underscore/dot + rest
    # This works well for model.matrix default naming for factors.
    ok <- !is_int & grepl("^[A-Za-z.][A-Za-z0-9._]*[A-Za-z][0-9A-Za-z._]*$", cn)
    if (any(ok)) {
      m <- regexec("^([A-Za-z.][A-Za-z0-9._]*?)([-]?[A-Za-z0-9].*)$", cn[ok])
      reg <- regmatches(cn[ok], m)
      out_ok <- cn[ok]
      for (i in seq_along(reg)) {
        if (length(reg[[i]]) == 3) {
          var <- reg[[i]][2]
          lvl <- reg[[i]][3]
          # Heuristic: only rewrite if it looks like "var"+"level" (common for factors)
          # Avoid mangling interaction terms like "x1:x2" (those won't match the regex anyway)
          out_ok[i] <- paste0(var, " = ", lvl)
        }
      }
      out[ok] <- out_ok
    }

    out
  }

  replace_beta <- function(prefix, X, label, names_vec = NULL) {
    cn <- get_cn(X, names_vec)
    if (is.null(cn)) return(v)
    cn <- pretty_col(cn)

    pat <- paste0("^", prefix, "\\[(\\d+)\\]$")
    m <- regexec(pat, v)
    reg <- regmatches(v, m)

    out <- v
    for (i in seq_along(reg)) {
      if (length(reg[[i]]) == 2) {
        idx <- as.integer(reg[[i]][2])
        if (!is.na(idx) && idx >= 1 && idx <= length(cn)) {
          out[i] <- paste0(label, ": ", cn[idx])
        }
      }
    }
    out
  }

  v <- replace_beta("beta_T",  stan_data$X_T, "AFT", stan_data$X_T_names %||% NULL)
  v <- replace_beta("beta_p",  stan_data$X_p, "p",   stan_data$X_p_names %||% NULL)
  v <- replace_beta("beta_k",  stan_data$X_k, "k",   stan_data$X_k_names %||% NULL)
  v <- replace_beta("beta_kL", stan_data$X_k, "kL",  stan_data$X_k_names %||% NULL)
  v <- replace_beta("beta_kR", stan_data$X_k, "kR",  stan_data$X_k_names %||% NULL)

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
