# Statistical derivations --------------------------------------------------

#' Execute derivation operations on assay arrays
#'
#' @param arrays Named list of arrays (sample x subject x contrast)
#' @param what Character vector of assay names to derive
#' @param options List of options (e.g., overwrite = TRUE)
#'
#' @return Named list of arrays including derived assays
#' @keywords internal
execute_derive <- function(arrays, what, options = list()) {
  overwrite <- isTRUE(options$overwrite)
  for (target in unique(as.character(what))) {
    if (target %in% names(arrays) && !overwrite) next

    arrays[[target]] <- switch(target,
      var = derive_var(arrays),
      se  = derive_se(arrays),
      t   = derive_t(arrays),
      z   = derive_z(arrays),
      p   = derive_p(arrays),
      stop("Unknown derivation target: ", target, call. = FALSE)
    )
  }
  arrays
}

derive_var <- function(arrays) {
  if (!"se" %in% names(arrays)) {
    stop("Cannot derive var: standard error 'se' not available", call. = FALSE)
  }
  arrays$se^2
}

derive_se <- function(arrays) {
  if (!"var" %in% names(arrays)) {
    stop("Cannot derive se: variance 'var' not available", call. = FALSE)
  }
  sqrt(arrays$var)
}

derive_t <- function(arrays) {
  if (!all(c("beta", "var") %in% names(arrays))) {
    stop("Cannot derive t: require 'beta' and 'var'", call. = FALSE)
  }
  arrays$beta / sqrt(arrays$var)
}

derive_z <- function(arrays) {
  if (all(c("t", "df") %in% names(arrays))) {
    t <- arrays$t
    df <- .broadcast_df(arrays$df, dim(t))
    p_two <- 2 * stats::pt(-abs(t), df)
    return(stats::qnorm(1 - p_two / 2) * sign(t))
  }

  if ("p" %in% names(arrays)) {
    if (!"beta" %in% names(arrays)) {
      stop("Cannot derive signed z from p without effect sign (beta)", call. = FALSE)
    }
    return(stats::qnorm(1 - arrays$p / 2) * sign(arrays$beta))
  }

  # Derive p from F if available, then z if sign information exists
  if (all(c("F", "df1", "df2") %in% names(arrays))) {
    pF <- stats::pf(arrays$F, arrays$df1, arrays$df2, lower.tail = FALSE)
    if (!"beta" %in% names(arrays)) {
      stop("Cannot derive signed z from F without effect sign (beta); use p or Fisher/Lancaster combiners.", call. = FALSE)
    }
    return(stats::qnorm(1 - pF / 2) * sign(arrays$beta))
  }

  stop("Cannot derive z: require {t, df} or {p, beta}", call. = FALSE)
}

derive_p <- function(arrays) {
  if (all(c("t", "df") %in% names(arrays))) {
    return(2 * stats::pt(-abs(arrays$t), arrays$df))
  }
  if ("z" %in% names(arrays)) {
    return(2 * stats::pnorm(-abs(arrays$z)))
  }
  if (all(c("F", "df1", "df2") %in% names(arrays))) {
    return(stats::pf(arrays$F, arrays$df1, arrays$df2, lower.tail = FALSE))
  }
  if (all(c("chi2", "df") %in% names(arrays))) {
    return(stats::pchisq(arrays$chi2, arrays$df, lower.tail = FALSE))
  }
  stop("Cannot derive p: missing required assays", call. = FALSE)
}

# Degrees of freedom aggregation -------------------------------------------

aggregate_df_satterthwaite <- function(M, var_source, df_source) {
  if (is.null(df_source)) {
    stop("Satterthwaite aggregation requires source df", call. = FALSE)
  }

  dims <- dim(var_source)
  n_target <- nrow(M)
  df_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (i in seq_len(n_target)) {
    w <- as.numeric(M[i, ])
    idx <- which(w != 0)
    if (!length(idx)) next
    w <- w[idx]

    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        s2 <- var_source[idx, j, k]
        df <- df_source[idx, j, k]
        numerator <- (sum((w^2) * s2))^2
        denominator <- sum((w^4) * (s2^2) / df)
        df_out[i, j, k] <- if (denominator > 0) numerator / denominator else Inf
      }
    }
  }
  df_out
}

.broadcast_df <- function(df, target_dim) {
  if (length(dim(df)) == 3) {
    return(df)
  }
  if (length(dim(df)) == 2) {
    return(array(df, dim = target_dim))
  }
  if (is.null(dim(df))) {
    return(array(df, dim = target_dim))
  }
  stop("Unsupported df dimensions", call. = FALSE)
}
