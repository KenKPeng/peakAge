#' Build Stan data for symmetric peak-age model (with covariates)
#'
#' Long-format input with one row per competition.
#' Athlete-level covariates are taken from the first row of each athlete after sorting
#' and are assumed constant within athlete (a warning is issued otherwise).
#'
#' This builder is aligned with the Stan model that uses:
#' - AFT: (K_T, X_T) with K_T >= 1
#' - p  : (K_p, X_p) with K_p >= 0
#' - logk: (K_k, X_k) with K_k >= 0
#'
#' NOTE: This builder returns the *shared* data used by both symmetric and asymmetric Stan models.
#' Prior hyperparameters are injected later in `fit_peakage()` based on `curve`.
#'
#' @param df Long-format data with one row per competition.
#' @param peak_formula Athlete-level covariate formula for log peak age (AFT),
#'   e.g. ~ x1 + x2. Default ~ 1.
#' @param p_formula Athlete-level covariate formula for p (peak height) submodel.
#'   Default ~ 1 (treated as "no covariates", i.e. K_p=0).
#' @param k_formula Athlete-level covariate formula for logk (curvature) submodel.
#'   Default ~ 1 (treated as "no covariates", i.e. K_k=0).
#' @param id,age,y Column names in df.
#' @param grainsize grainsize for reduce_sum.
#' @return A list ready for Stan (without prior hyperparameters).
#' @export
build_stan_data_sym <- function(df,
                                peak_formula = ~ 1,
                                p_formula    = ~ 1,
                                k_formula    = ~ 1,
                                id = "Athlete_id",
                                age = "age",
                                y = "y",
                                grainsize = 100) {
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!all(c(id, age, y) %in% names(df))) {
    stop("`df` must contain columns: ", paste(c(id, age, y), collapse = ", "))
  }

  peak_formula <- stats::as.formula(peak_formula)
  p_formula    <- stats::as.formula(p_formula)
  k_formula    <- stats::as.formula(k_formula)

  # vars needed for slicing df (keep minimal)
  vars_needed <- unique(c(
    id, age, y,
    all.vars(peak_formula),
    all.vars(p_formula),
    all.vars(k_formula)
  ))
  vars_needed <- vars_needed[vars_needed %in% names(df)]
  df2 <- df[, vars_needed, drop = FALSE]

  # coerce core columns
  df2[[age]] <- as.numeric(df2[[age]])
  df2[[y]]   <- as.numeric(df2[[y]])

  if (anyNA(df2[[id]]))  stop("Missing values in athlete id column.")
  if (anyNA(df2[[age]])) stop("Missing values in age column.")
  if (anyNA(df2[[y]]))   stop("Missing values in y column.")

  # sort by athlete then age
  df2 <- df2[order(df2[[id]], df2[[age]]), , drop = FALSE]

  ath_ids <- sort(unique(df2[[id]]))
  I <- length(ath_ids)
  N <- nrow(df2)

  id_int <- as.integer(factor(df2[[id]], levels = ath_ids))

  # athlete-level covariates: first row per athlete after sorting
  idx_first <- match(ath_ids, df2[[id]])
  ath_df <- df2[idx_first, , drop = FALSE]

  # warn if covariates vary within athlete
  .warn_if_covariates_vary(
    df_long = df2,
    id_col = id,
    covars = unique(c(all.vars(peak_formula), all.vars(p_formula), all.vars(k_formula)))
  )

  # --- AFT design matrix (keep intercept by default) ---
  X_T <- unname(stats::model.matrix(peak_formula, data = ath_df))
  K_T <- ncol(X_T)
  if (K_T < 1) stop("`peak_formula` produced K_T=0 columns; not allowed for Stan (K_T must be >= 1).")

  # --- p / logk design matrices ---
  # By default, treat ~ 1 as "no covariates" (K=0) to avoid collinearity with p_pop/logk_pop.
  X_p <- .mm_or_empty(p_formula, ath_df)
  X_k <- .mm_or_empty(k_formula, ath_df)
  K_p <- ncol(X_p)
  K_k <- ncol(X_k)

  list(
    N = as.integer(N),
    I = as.integer(I),
    id = as.integer(id_int),
    age = as.numeric(df2[[age]]),
    y = as.numeric(df2[[y]]),
    grainsize = as.integer(grainsize),

    K_T = as.integer(K_T),
    X_T = X_T,

    K_p = as.integer(K_p),
    X_p = X_p,

    K_k = as.integer(K_k),
    X_k = X_k
  )
}

# ---- internal helpers ----

# Return model.matrix for p/logk formulas, but allow "no covariates" => 0-column matrix.
.mm_or_empty <- function(f, data) {
  f <- stats::as.formula(f)

  tt <- stats::terms(f)
  has_terms <- attr(tt, "term.labels")

  # If no non-intercept terms, return empty matrix (I x 0)
  if (length(has_terms) == 0L) {
    I <- nrow(data)
    return(matrix(numeric(0), nrow = I, ncol = 0))
  }

  # Otherwise, drop intercept by default and build matrix
  f0 <- stats::update(f, ~ 0 + .)
  mm <- stats::model.matrix(f0, data = data)
  unname(mm)
}

.warn_if_covariates_vary <- function(df_long, id_col, covars) {
  covars <- unique(covars)
  covars <- covars[covars %in% names(df_long)]
  covars <- setdiff(covars, c(id_col, "age", "y"))
  if (length(covars) == 0) return(invisible(NULL))

  for (v in covars) {
    split_v <- split(df_long[[v]], df_long[[id_col]])
    bad <- vapply(split_v, function(x) {
      x <- x[!is.na(x)]
      length(unique(x)) > 1
    }, logical(1))

    if (any(bad)) {
      warning(
        "Covariate `", v, "` varies within athlete for ",
        sum(bad), " athlete(s). Using the first row per athlete in model.matrix().",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

