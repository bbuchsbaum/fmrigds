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
  .check_map_columns(map, arrays, context = "map_to()")

  has_beta_var <- all(c("beta", "var") %in% names(arrays))

  if (has_beta_var) {
    arrays <- .map_with_effect_scale(arrays, map, node$uncertainty)
  } else {
    arrays <- .map_without_effect_scale(arrays, map, combine)
  }

  list(arrays = arrays, space = node$target_space)
}

.check_map_columns <- function(map, arrays, context = "map_to()") {
  n_samples <- if (length(arrays)) dim(arrays[[1L]])[1L] else NA_integer_
  if (!is.na(n_samples) && ncol(map) != n_samples) {
    stop(sprintf(
      paste0(
        "%s: the map has %d columns but the data currently has %d samples. ",
        "The map must be defined on the current sample axis; earlier mask() or ",
        "subset() steps change the number of samples, so build the map for the ",
        "masked/subsetted space (or apply map_to() before masking)."
      ),
      context, ncol(map), n_samples
    ), call. = FALSE)
  }
  invisible(TRUE)
}

# Effect-scale map: beta/var are propagated through the operator, df is
# propagated by Satterthwaite (df_rule = "satterthwaite") or broadcast when it
# is constant per subject/contrast; every other assay is dropped and the
# derived statistics are recomputed (see .reset_after_transform()).
.map_with_effect_scale <- function(arrays, map, uncertainty) {
  var_source <- arrays$var
  df_source <- arrays$df
  had_stats <- any(c("z", "p") %in% names(arrays))

  res <- switch(uncertainty$mode,
    independent = propagate_variance_independent(map, arrays$beta, arrays$var),
    cov_provider = propagate_variance_covariance(map, arrays$beta, arrays$var, uncertainty$cov_provider),
    stop("Unsupported uncertainty mode: ", uncertainty$mode, call. = FALSE)
  )

  out <- list(beta = res$beta, var = res$var)
  if (!is.null(df_source)) {
    out$df <- if (identical(uncertainty$df_rule, "satterthwaite")) {
      aggregate_df_satterthwaite(map, var_source, df_source)
    } else {
      .broadcast_constant_df(df_source, nrow(map))
    }
  }

  .reset_after_transform(out, c("beta", "var", "df"), had_stats = had_stats)
}

.map_without_effect_scale <- function(arrays, map, combine) {
  if (is.null(combine)) {
    stop("map_to requires beta/var or an explicit combiner when only test statistics are available.", call. = FALSE)
  }

  if (combine == "stouffer") {
    if (!"z" %in% names(arrays)) stop("Stouffer combine requires z-scores", call. = FALSE)
    arrays <- .combine_stouffer(arrays, map)
    keep <- c("z", "p")
  } else if (combine == "fisher") {
    if (!"z" %in% names(arrays) && !"p" %in% names(arrays)) {
      stop("Fisher combine requires z-scores or p-values", call. = FALSE)
    }
    arrays <- .combine_fisher(arrays, map)
    keep <- c("chi2", "df", "p", "z")
  }

  # Only the combiner outputs live on the target sample axis; drop the rest.
  arrays[intersect(keep, names(arrays))]
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
