# Map execution -----------------------------------------------------------

apply_map_to <- function(node, arrays) {
  combine <- node$combine
  if (!is.null(combine)) {
    combine <- match.arg(combine, c("stouffer", "fisher"))
  }

  map_obj <- node$map
  if (inherits(map_obj, "map_linear")) {
    map <- map_obj$operator
    if (!is.null(map_obj$by_subject)) {
      stop("Subject-specific maps require align(); call align() before map_to().", call. = FALSE)
    }
  } else {
    map <- map_obj
  }

  if (!is.matrix(map) && !inherits(map, "Matrix")) {
    stop("map must be matrix or map_linear", call. = FALSE)
  }
  map <- as.matrix(map)

  has_beta_var <- all(c("beta", "var") %in% names(arrays))

  if (has_beta_var) {
    arrays <- .map_with_effect_scale(arrays, map, node$uncertainty)
  } else {
    arrays <- .map_without_effect_scale(arrays, map, combine)
  }

  list(arrays = arrays, space = node$target_space)
}

.map_with_effect_scale <- function(arrays, map, uncertainty) {
  var_source <- arrays$var
  df_source <- arrays$df

  res <- switch(uncertainty$mode,
    independent = propagate_variance_independent(map, arrays$beta, arrays$var),
    cov_provider = propagate_variance_covariance(map, arrays$beta, arrays$var, uncertainty$cov_provider),
    stop("Unsupported uncertainty mode: ", uncertainty$mode, call. = FALSE)
  )

  arrays$beta <- res$beta
  arrays$var <- res$var

  if ("se" %in% names(arrays)) arrays$se <- sqrt(arrays$var)

  if ("t" %in% names(arrays)) {
    arrays$t <- arrays$beta / sqrt(arrays$var)
    if ("df" %in% names(arrays) && !is.null(df_source) && identical(uncertainty$df_rule, "satterthwaite")) {
      arrays$df <- aggregate_df_satterthwaite(map, var_source, df_source)
    }
  }

  if ("z" %in% names(arrays)) {
    arrays$z <- derive_z(arrays)
  }

  arrays
}

.map_without_effect_scale <- function(arrays, map, combine) {
  if (is.null(combine)) {
    stop("map_to requires beta/var or an explicit combiner when only test statistics are available.", call. = FALSE)
  }

  if (combine == "stouffer") {
    if (!"z" %in% names(arrays)) stop("Stouffer combine requires z-scores", call. = FALSE)
    arrays <- .combine_stouffer(arrays, map)
  } else if (combine == "fisher") {
    if (!"z" %in% names(arrays) && !"p" %in% names(arrays)) {
      stop("Fisher combine requires z-scores or p-values", call. = FALSE)
    }
    arrays <- .combine_fisher(arrays, map)
  }

  arrays
}

.combine_stouffer <- function(arrays, map) {
  map <- as.matrix(map)
  dims <- dim(arrays$z)
  n_target <- nrow(map)
  z_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (i in seq_len(n_target)) {
    w <- map[i, ]
    denom <- sum(w^2)
    if (denom <= 0) next
    idx <- which(w != 0)
    if (!length(idx)) next
    w <- w[idx]
    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        z_out[i, j, k] <- sum(w * arrays$z[idx, j, k]) / sqrt(denom)
      }
    }
  }

  arrays$z <- z_out
  arrays$p <- 2 * stats::pnorm(-abs(z_out))
  arrays
}

.combine_fisher <- function(arrays, map) {
  map <- as.matrix(map)
  dims <- if ("z" %in% names(arrays)) dim(arrays$z) else dim(arrays$p)
  n_target <- nrow(map)
  chi_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  p_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  df_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (i in seq_len(n_target)) {
    w <- map[i, ]
    idx <- which(w != 0)
    if (!length(idx)) next
    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        if ("z" %in% names(arrays)) {
          p_vals <- 2 * stats::pnorm(-abs(arrays$z[idx, j, k]))
        } else {
          p_vals <- arrays$p[idx, j, k]
        }
        chi <- -2 * sum(log(p_vals))
        df_val <- 2 * length(p_vals)
        chi_out[i, j, k] <- chi
        df_out[i, j, k] <- df_val
        p_out[i, j, k] <- stats::pchisq(chi, df_val, lower.tail = FALSE)
      }
    }
  }

  arrays$chi2 <- chi_out
  arrays$df <- df_out
  arrays$p <- p_out
  arrays$z <- stats::qnorm(1 - p_out/2)
  arrays
}
