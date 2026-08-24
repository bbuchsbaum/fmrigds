# Model-conditioned robust sensitivity -------------------------------------

.robust_examination_methods <- c(
  "meta:fe", "meta:re", "meta:fe_reg", "meta:re_reg"
)

.robust_feature_status <- c(
  available = 1L,
  insufficient_samples = 2L,
  rank_deficient = 3L,
  leverage_limit = 4L,
  nonconverged = 5L,
  invalid_variance = 6L
)

.robust_fit_failure <- function(status,
                                n_subject,
                                eligible_n = 0L,
                                iterations = 0L) {
  list(
    status = status,
    coefficients = NULL,
    robust_factor = rep(NA_real_, n_subject),
    eligible_n = as.integer(eligible_n),
    iterations = as.integer(iterations),
    downweighted_n = NA_integer_,
    max_robust_leverage = NA_real_,
    objective = NA_real_
  )
}

.huber_factor <- function(residual, tuning_constant) {
  magnitude <- abs(residual)
  out <- rep(1, length(residual))
  downweighted <- magnitude > tuning_constant
  out[downweighted] <- tuning_constant / magnitude[downweighted]
  out
}

.huber_objective <- function(residual, tuning_constant) {
  magnitude <- abs(residual)
  quadratic <- magnitude <= tuning_constant
  value <- numeric(length(residual))
  value[quadratic] <- 0.5 * residual[quadratic]^2
  value[!quadratic] <- tuning_constant * magnitude[!quadratic] -
    0.5 * tuning_constant^2
  sum(value)
}

.robust_weighted_least_squares <- function(X,
                                           y,
                                           weight,
                                           rank_tolerance) {
  if (any(!is.finite(weight)) || any(weight <= 0)) return(NULL)
  root_weight <- sqrt(weight)
  weighted_X <- X * root_weight
  decomposition <- qr(weighted_X, tol = rank_tolerance, LAPACK = FALSE)
  if (decomposition$rank < ncol(X)) return(NULL)
  coefficients <- tryCatch(
    qr.coef(decomposition, y * root_weight),
    error = function(e) NULL
  )
  if (is.null(coefficients) || any(!is.finite(coefficients))) return(NULL)
  inverse <- .diagnostic_inverse(crossprod(weighted_X))
  if (is.null(inverse)) return(NULL)
  leverage <- weight * rowSums((X %*% inverse) * X)
  if (any(!is.finite(leverage))) return(NULL)
  list(
    coefficients = as.numeric(coefficients),
    leverage = leverage
  )
}

# Solve sum_i X_i psi((y_i - X_i beta) / sqrt(v_i)) / sqrt(v_i) = 0.
# The returned robust factors multiply the original inverse-variance weights.
.huber_ivw_fit <- function(y,
                           X,
                           variance,
                           settings,
                           rank_tolerance = sqrt(.Machine$double.eps)) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  variance <- as.numeric(variance)
  n_subject <- length(y)
  if (nrow(X) != n_subject || length(variance) != n_subject) {
    stop("Huber IVW inputs have incompatible subject dimensions.", call. = FALSE)
  }
  if (!ncol(X)) {
    return(.robust_fit_failure("rank_deficient", n_subject))
  }
  valid <- is.finite(y) & is.finite(variance) & variance > 0 &
    rowSums(!is.finite(X)) == 0L
  eligible_n <- sum(valid)
  if (eligible_n < ncol(X) + 1L) {
    return(.robust_fit_failure(
      "insufficient_samples", n_subject, eligible_n = eligible_n
    ))
  }

  yv <- y[valid]
  Xv <- X[valid, , drop = FALSE]
  vv <- variance[valid]
  base_weight <- 1 / vv
  initial <- .robust_weighted_least_squares(
    Xv, yv, base_weight, rank_tolerance
  )
  if (is.null(initial)) {
    return(.robust_fit_failure(
      "rank_deficient", n_subject, eligible_n = eligible_n
    ))
  }
  if (max(initial$leverage) >= settings$leverage_limit) {
    return(.robust_fit_failure(
      "leverage_limit", n_subject, eligible_n = eligible_n
    ))
  }

  coefficients <- initial$coefficients
  converged <- FALSE
  iterations <- 0L
  candidate <- initial
  robust_factor <- rep(1, eligible_n)
  for (iteration in seq_len(settings$max_iterations)) {
    iterations <- as.integer(iteration)
    standardized <- (yv - drop(Xv %*% coefficients)) / sqrt(vv)
    robust_factor <- .huber_factor(
      standardized,
      settings$tuning_constant
    )
    candidate <- .robust_weighted_least_squares(
      Xv,
      yv,
      base_weight * robust_factor,
      rank_tolerance
    )
    if (is.null(candidate)) {
      return(.robust_fit_failure(
        "rank_deficient", n_subject,
        eligible_n = eligible_n,
        iterations = iterations
      ))
    }
    if (max(candidate$leverage) >= settings$leverage_limit) {
      return(.robust_fit_failure(
        "leverage_limit", n_subject,
        eligible_n = eligible_n,
        iterations = iterations
      ))
    }
    scale <- pmax(1, abs(coefficients), abs(candidate$coefficients))
    step <- max(abs(candidate$coefficients - coefficients) / scale)
    coefficients <- candidate$coefficients
    if (is.finite(step) && step <= settings$convergence_tolerance) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    return(.robust_fit_failure(
      "nonconverged", n_subject,
      eligible_n = eligible_n,
      iterations = iterations
    ))
  }

  standardized <- (yv - drop(Xv %*% coefficients)) / sqrt(vv)
  robust_factor <- .huber_factor(standardized, settings$tuning_constant)
  final <- .robust_weighted_least_squares(
    Xv,
    yv,
    base_weight * robust_factor,
    rank_tolerance
  )
  if (is.null(final)) {
    return(.robust_fit_failure(
      "rank_deficient", n_subject,
      eligible_n = eligible_n,
      iterations = iterations
    ))
  }
  if (max(final$leverage) >= settings$leverage_limit) {
    return(.robust_fit_failure(
      "leverage_limit", n_subject,
      eligible_n = eligible_n,
      iterations = iterations
    ))
  }
  fixed_point_scale <- pmax(1, abs(coefficients), abs(final$coefficients))
  fixed_point_step <- max(
    abs(final$coefficients - coefficients) / fixed_point_scale
  )
  if (!is.finite(fixed_point_step) ||
      fixed_point_step > 10 * settings$convergence_tolerance) {
    return(.robust_fit_failure(
      "nonconverged", n_subject,
      eligible_n = eligible_n,
      iterations = iterations
    ))
  }
  coefficients <- final$coefficients
  standardized <- (yv - drop(Xv %*% coefficients)) / sqrt(vv)
  robust_factor <- .huber_factor(standardized, settings$tuning_constant)
  subject_factor <- rep(NA_real_, n_subject)
  subject_factor[valid] <- robust_factor
  list(
    status = "available",
    coefficients = coefficients,
    robust_factor = subject_factor,
    eligible_n = as.integer(eligible_n),
    iterations = iterations,
    downweighted_n = as.integer(sum(
      robust_factor < 1 - sqrt(.Machine$double.eps)
    )),
    max_robust_leverage = max(final$leverage),
    objective = .huber_objective(standardized, settings$tuning_constant)
  )
}

.initialize_robust_sensitivity <- function(scan_context,
                                           model_context,
                                           control) {
  settings <- control$robust %||% NULL
  if (is.null(settings)) return(NULL)
  n_subject <- length(model_context$subjects)
  n_contrast <- length(scan_context$contrasts)
  n_estimand <- nrow(model_context$estimand_matrix)
  supported <- model_context$method %in% .robust_examination_methods &&
    identical(model_context$variance_mode, "measured_first_level")
  reason <- if (supported) {
    NA_character_
  } else {
    paste0(
      "Huber IVW sensitivity requires a measured-variance meta reducer; got '",
      model_context$method, "' with variance mode '",
      model_context$variance_mode, "'."
    )
  }
  matrix_sc <- function(value = 0) matrix(value, n_subject, n_contrast)
  array_ce <- function(value = 0) array(value, c(n_contrast, n_estimand))
  list(
    supported = supported,
    unsupported_reason = reason,
    method = settings$method,
    mode = if (model_context$method %in% c("meta:re", "meta:re_reg")) {
      "huber_ivw_tau2_fixed_full"
    } else {
      "huber_ivw"
    },
    settings = settings,
    subjects = model_context$subjects,
    contrasts = scan_context$contrasts,
    estimands = rownames(model_context$estimand_matrix),
    n_sample = length(scan_context$source_samples),
    maps = list(),
    feature_count = numeric(n_contrast),
    status_counts = matrix(
      0,
      nrow = n_contrast,
      ncol = length(.robust_feature_status),
      dimnames = list(scan_context$contrasts, names(.robust_feature_status))
    ),
    iteration_sum = numeric(n_contrast),
    iteration_max = integer(n_contrast),
    downweighted_sum = numeric(n_contrast),
    subject_factor_sum = matrix_sc(),
    subject_factor_count = matrix_sc(),
    subject_downweighted_count = matrix_sc(),
    subject_factor_min = matrix_sc(NA_real_),
    raw_shift_sum_sq = array_ce(),
    raw_shift_count = array_ce(),
    scaled_shift_sum_sq = array_ce(),
    scaled_shift_count = array_ce(),
    scaled_shift_max = array_ce(NA_real_)
  )
}

.fit_robust_examination_block <- function(beta,
                                          var,
                                          fit,
                                          diagnostic,
                                          model_context,
                                          control) {
  n_subject <- nrow(beta)
  n_sample <- ncol(beta)
  estimands <- model_context$estimand_matrix
  n_estimand <- nrow(estimands)
  X <- model_context$X
  if (is.null(X)) {
    X <- matrix(
      1,
      nrow = n_subject,
      ncol = 1L,
      dimnames = list(model_context$subjects, "pooled_effect")
    )
  }
  random_effects <- model_context$method %in% c("meta:re", "meta:re_reg")
  tau2 <- if (random_effects) as.numeric(fit$tau2) else rep(0, n_sample)
  if (length(tau2) != n_sample) tau2 <- rep(NA_real_, n_sample)
  out <- list(
    status = rep("invalid_variance", n_sample),
    status_code = rep(.robust_feature_status[["invalid_variance"]], n_sample),
    effects = matrix(NA_real_, n_estimand, n_sample),
    delta_effect = matrix(NA_real_, n_estimand, n_sample),
    scaled_delta_effect = matrix(NA_real_, n_estimand, n_sample),
    robust_factor = matrix(NA_real_, n_subject, n_sample),
    iterations = rep(NA_integer_, n_sample),
    eligible_n = rep(NA_integer_, n_sample),
    downweighted_n = rep(NA_integer_, n_sample),
    max_robust_leverage = rep(NA_real_, n_sample),
    objective = rep(NA_real_, n_sample)
  )
  for (b in seq_len(n_sample)) {
    if (!is.finite(tau2[b]) || tau2[b] < 0) next
    variance <- var[, b] + tau2[b]
    robust_fit <- .huber_ivw_fit(
      beta[, b],
      X,
      variance,
      settings = control$robust,
      rank_tolerance = control$tolerance$rank
    )
    out$status[b] <- robust_fit$status
    out$status_code[b] <- .robust_feature_status[[robust_fit$status]]
    out$iterations[b] <- robust_fit$iterations
    out$eligible_n[b] <- robust_fit$eligible_n
    out$downweighted_n[b] <- robust_fit$downweighted_n
    out$max_robust_leverage[b] <- robust_fit$max_robust_leverage
    out$objective[b] <- robust_fit$objective
    out$robust_factor[, b] <- robust_fit$robust_factor
    if (!identical(robust_fit$status, "available")) next
    effect <- drop(estimands %*% robust_fit$coefficients)
    out$effects[, b] <- effect
    out$delta_effect[, b] <- effect - diagnostic$full_effect[, b]
    scalable <- is.finite(diagnostic$full_se[, b]) &
      diagnostic$full_se[, b] > control$tolerance$degeneracy
    out$scaled_delta_effect[scalable, b] <-
      out$delta_effect[scalable, b] / diagnostic$full_se[scalable, b]
  }
  out
}

.accumulate_robust_maps <- function(robust,
                                    values,
                                    ordinal,
                                    contrast_index) {
  for (name in names(values)) {
    block_values <- values[[name]]
    if (length(block_values) != length(ordinal)) {
      stop("Robust map '", name, "' has the wrong block length.", call. = FALSE)
    }
    if (is.null(robust$maps[[name]])) {
      missing_value <- if (is.integer(block_values)) NA_integer_ else NA_real_
      robust$maps[[name]] <- array(
        missing_value,
        c(robust$n_sample, 1L, length(robust$contrasts)),
        dimnames = list(NULL, "meta", robust$contrasts)
      )
    }
    robust$maps[[name]][ordinal, 1L, contrast_index] <- block_values
  }
  robust
}

.accumulate_robust_sensitivity <- function(robust,
                                           beta,
                                           var,
                                           fit,
                                           diagnostic,
                                           block,
                                           contrast_index,
                                           model_context,
                                           control) {
  if (is.null(robust) || !isTRUE(robust$supported)) return(robust)
  result <- .fit_robust_examination_block(
    beta, var, fit, diagnostic, model_context, control
  )
  map_values <- list(
    status_code = as.integer(result$status_code),
    iterations = as.integer(result$iterations),
    eligible_n = as.integer(result$eligible_n),
    downweighted_n = as.integer(result$downweighted_n),
    max_robust_leverage = result$max_robust_leverage,
    objective = result$objective
  )
  for (e in seq_along(robust$estimands)) {
    name <- robust$estimands[e]
    map_values[[paste0("effect:", name)]] <- result$effects[e, ]
    map_values[[paste0("delta_effect:", name)]] <- result$delta_effect[e, ]
    map_values[[paste0("delta_effect_primary_se:", name)]] <-
      result$scaled_delta_effect[e, ]
  }
  robust <- .accumulate_robust_maps(
    robust, map_values, block$ordinal, contrast_index
  )
  robust$feature_count[contrast_index] <-
    robust$feature_count[contrast_index] + ncol(beta)
  counts <- table(factor(result$status, levels = names(.robust_feature_status)))
  robust$status_counts[contrast_index, ] <-
    robust$status_counts[contrast_index, ] + as.numeric(counts)

  available <- result$status == "available"
  if (any(available)) {
    iterations <- result$iterations[available]
    robust$iteration_sum[contrast_index] <-
      robust$iteration_sum[contrast_index] + sum(iterations)
    robust$iteration_max[contrast_index] <- max(
      robust$iteration_max[contrast_index],
      iterations
    )
    robust$downweighted_sum[contrast_index] <-
      robust$downweighted_sum[contrast_index] +
      sum(result$downweighted_n[available])
  }

  for (i in seq_along(robust$subjects)) {
    factors <- result$robust_factor[i, ]
    ok <- is.finite(factors) & available
    if (!any(ok)) next
    values <- factors[ok]
    robust$subject_factor_sum[i, contrast_index] <-
      robust$subject_factor_sum[i, contrast_index] + sum(values)
    robust$subject_factor_count[i, contrast_index] <-
      robust$subject_factor_count[i, contrast_index] + length(values)
    robust$subject_downweighted_count[i, contrast_index] <-
      robust$subject_downweighted_count[i, contrast_index] +
      sum(values < 1 - sqrt(.Machine$double.eps))
    robust$subject_factor_min[i, contrast_index] <- min(
      c(robust$subject_factor_min[i, contrast_index], values),
      na.rm = TRUE
    )
  }

  for (e in seq_along(robust$estimands)) {
    raw <- result$delta_effect[e, ]
    raw_ok <- available & is.finite(raw)
    if (any(raw_ok)) {
      robust$raw_shift_sum_sq[contrast_index, e] <-
        robust$raw_shift_sum_sq[contrast_index, e] + sum(raw[raw_ok]^2)
      robust$raw_shift_count[contrast_index, e] <-
        robust$raw_shift_count[contrast_index, e] + sum(raw_ok)
    }
    scaled <- result$scaled_delta_effect[e, ]
    scaled_ok <- available & is.finite(scaled)
    if (any(scaled_ok)) {
      capped <- pmin(abs(scaled[scaled_ok]), control$geometry$cap)
      robust$scaled_shift_sum_sq[contrast_index, e] <-
        robust$scaled_shift_sum_sq[contrast_index, e] + sum(capped^2)
      robust$scaled_shift_count[contrast_index, e] <-
        robust$scaled_shift_count[contrast_index, e] + sum(scaled_ok)
      robust$scaled_shift_max[contrast_index, e] <- max(
        c(robust$scaled_shift_max[contrast_index, e], abs(scaled[scaled_ok])),
        na.rm = TRUE
      )
    }
  }
  robust
}

.robust_dominant_status <- function(counts) {
  total <- sum(counts)
  if (!total) return("insufficient_samples")
  if (counts[["available"]] == total) return("available")
  if (counts[["available"]] > 0) return("partial")
  names(counts)[which.max(counts)]
}

.finalize_robust_maps <- function(robust, compiled, model_context) {
  if (!isTRUE(robust$supported) || !length(robust$maps)) return(NULL)
  components <- .examination_output_components(compiled)
  new_gds(
    assays = robust$maps,
    space = components$space,
    subjects = "meta",
    contrasts = robust$contrasts,
    row_data = components$row_data,
    metadata = list(
      examination = list(
        scope = "robust_sensitivity",
        categorical_assays = "status_code",
        non_interpolable_assays = "status_code",
        status_lookup = data.frame(
          code = unname(.robust_feature_status),
          status = names(.robust_feature_status),
          stringsAsFactors = FALSE
        ),
        method = robust$method,
        mode = robust$mode,
        model_context_digest = model_context$digest
      )
    )
  )
}

.finalize_robust_sensitivity <- function(robust,
                                         compiled,
                                         model_context,
                                         control) {
  if (is.null(robust)) return(NULL)
  contrast_rows <- lapply(seq_along(robust$contrasts), function(k) {
    counts <- robust$status_counts[k, ]
    available_n <- unname(counts[["available"]])
    data.frame(
      contrast = robust$contrasts[k],
      mode = robust$mode,
      status = if (robust$supported) {
        .robust_dominant_status(counts)
      } else {
        "unsupported_reducer"
      },
      feature_n = as.integer(robust$feature_count[k]),
      available_n = as.integer(available_n),
      availability_fraction = if (robust$feature_count[k] > 0) {
        available_n / robust$feature_count[k]
      } else {
        NA_real_
      },
      insufficient_samples_n = as.integer(counts[["insufficient_samples"]]),
      rank_deficient_n = as.integer(counts[["rank_deficient"]]),
      leverage_limit_n = as.integer(counts[["leverage_limit"]]),
      nonconverged_n = as.integer(counts[["nonconverged"]]),
      invalid_variance_n = as.integer(counts[["invalid_variance"]]),
      mean_iterations = if (available_n > 0) {
        robust$iteration_sum[k] / available_n
      } else {
        NA_real_
      },
      max_iterations = if (available_n > 0) robust$iteration_max[k] else NA_integer_,
      mean_downweighted_n = if (available_n > 0) {
        robust$downweighted_sum[k] / available_n
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  contrast_data <- do.call(rbind, contrast_rows)

  subject_rows <- list()
  index <- 1L
  for (k in seq_along(robust$contrasts)) {
    for (i in seq_along(robust$subjects)) {
      count <- robust$subject_factor_count[i, k]
      subject_rows[[index]] <- data.frame(
        subject = robust$subjects[i],
        contrast = robust$contrasts[k],
        mean_downweight_factor = if (count > 0) {
          robust$subject_factor_sum[i, k] / count
        } else {
          NA_real_
        },
        min_downweight_factor = if (count > 0) {
          robust$subject_factor_min[i, k]
        } else {
          NA_real_
        },
        downweighted_fraction = if (count > 0) {
          robust$subject_downweighted_count[i, k] / count
        } else {
          NA_real_
        },
        eligible_feature_n = as.integer(count),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  subject_data <- do.call(rbind, subject_rows)

  estimand_rows <- list()
  index <- 1L
  for (k in seq_along(robust$contrasts)) {
    for (e in seq_along(robust$estimands)) {
      raw_n <- robust$raw_shift_count[k, e]
      scaled_n <- robust$scaled_shift_count[k, e]
      estimand_rows[[index]] <- data.frame(
        contrast = robust$contrasts[k],
        estimand = robust$estimands[e],
        mode = robust$mode,
        effect_shift_rms = if (raw_n > 0) {
          sqrt(robust$raw_shift_sum_sq[k, e] / raw_n)
        } else {
          NA_real_
        },
        effect_shift_primary_se_energy = if (scaled_n > 0) {
          sqrt(robust$scaled_shift_sum_sq[k, e] / scaled_n)
        } else {
          NA_real_
        },
        max_abs_effect_shift_primary_se = if (scaled_n > 0) {
          robust$scaled_shift_max[k, e]
        } else {
          NA_real_
        },
        eligible_feature_n = as.integer(scaled_n),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  estimand_data <- do.call(rbind, estimand_rows)

  availability <- estimand_data[, c("contrast", "estimand", "mode"), drop = FALSE]
  availability$status <- vapply(availability$contrast, function(contrast) {
    row <- contrast_data[contrast_data$contrast == contrast, , drop = FALSE]
    as.character(row$status[1L])
  }, character(1L))
  availability$reason <- vapply(seq_len(nrow(availability)), function(i) {
    status <- availability$status[i]
    if (identical(status, "available")) return(NA_character_)
    if (!robust$supported) return(robust$unsupported_reason)
    row <- contrast_data[
      contrast_data$contrast == availability$contrast[i],
      ,
      drop = FALSE
    ]
    if (identical(status, "partial")) {
      return(paste0(
        row$available_n[1L], " of ", row$feature_n[1L],
        " robust feature fits are available; inspect status counts and maps."
      ))
    }
    paste0("No robust feature fits are available; feature status: ", status, ".")
  }, character(1L))

  structure(
    list(
      method = robust$method,
      mode = robust$mode,
      interpretation = paste(
        "Alternate Huber M-functional using the primary model design and",
        "variance contract; it does not replace the primary reducer, define",
        "an outlier probability, or alter review status."
      ),
      availability = availability,
      contrast_data = contrast_data,
      estimand_data = estimand_data,
      subject_data = subject_data,
      maps = .finalize_robust_maps(robust, compiled, model_context),
      config = robust$settings,
      provenance = list(
        model_context_digest = model_context$digest,
        variance_mode = model_context$variance_mode,
        tau2_contract = if (grepl("tau2_fixed_full", robust$mode, fixed = TRUE)) {
          "full-data tau2 held fixed"
        } else {
          "not applicable"
        },
        robust_inference = "not provided"
      )
    ),
    class = "gds_robust_sensitivity"
  )
}
