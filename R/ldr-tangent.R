# Local displacement rescue: whole-search tangent screen -----------------

.ldr_tangent_fold_ids <- function(subject_ids, folds) {
  labels <- as.character(subject_ids)
  if (anyNA(labels) || anyDuplicated(labels)) {
    stop("LDR subject identifiers must be unique and non-missing.",
         call. = FALSE)
  }
  key <- vapply(
    labels,
    digest::digest,
    character(1L),
    algo = "xxhash32",
    serialize = FALSE
  )
  sorted <- order(key, labels)
  assigned <- rep(seq_len(as.integer(folds)), length.out = length(labels))
  out <- integer(length(labels))
  out[sorted] <- assigned
  out
}

.ldr_noise_correlation <- function(y_a, y_b, var_a, var_b, indices,
                                   shrinkage = 0.25) {
  shrinkage <- .ldr_scalar(
    shrinkage, "covariance_shrinkage", lower = 0, upper = 1
  )
  indices <- as.integer(indices)
  difference_variance <- var_a[, indices, drop = FALSE] +
    var_b[, indices, drop = FALSE]
  standardized <- (y_a[, indices, drop = FALSE] -
                      y_b[, indices, drop = FALSE]) /
    sqrt(difference_variance)
  if (ncol(standardized) < 3L || any(!is.finite(standardized))) {
    return(NULL)
  }
  standardized <- sweep(
    standardized, 1L, rowMeans(standardized), "-"
  )
  covariance <- tcrossprod(standardized) /
    max(1, ncol(standardized) - 1L)
  scale <- sqrt(pmax(diag(covariance), sqrt(.Machine$double.eps)))
  correlation <- covariance / outer(scale, scale)
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  correlation <- (1 - shrinkage) * correlation +
    shrinkage * diag(nrow(correlation))
  correlation <- (correlation + t(correlation)) / 2
  eig <- eigen(correlation, symmetric = TRUE)
  floor_value <- max(eig$values) * 1e-8
  values <- pmax(eig$values, floor_value)
  correlation <- eig$vectors %*% (values * t(eig$vectors))
  final_scale <- sqrt(pmax(diag(correlation), sqrt(.Machine$double.eps)))
  correlation <- correlation / outer(final_scale, final_scale)
  diag(correlation) <- 1
  list(
    correlation = correlation,
    raw_min_eigenvalue = min(eig$values),
    min_eigenvalue = min(eigen(correlation, symmetric = TRUE,
                               only.values = TRUE)$values),
    shrinkage = shrinkage
  )
}

.ldr_tangent_basis <- function(template, geometry, tolerance = 1e-7) {
  shifted <- .ldr_shift_templates(template, geometry$shift_operators)
  nonzero <- which(geometry$shift_distance_mm > sqrt(.Machine$double.eps))
  if (!length(nonzero)) return(NULL)
  displacement <- geometry$shift_world[nonzero, , drop = FALSE]
  displacement_svd <- svd(displacement)
  if (!length(displacement_svd$d) || max(displacement_svd$d) <= 0) {
    return(NULL)
  }
  spatial_rank <- sum(
    displacement_svd$d > tolerance * max(displacement_svd$d)
  )
  if (!spatial_rank) return(NULL)
  world_basis <- displacement_svd$v[, seq_len(spatial_rank), drop = FALSE]
  coordinates <- displacement %*% world_basis
  difference <- shifted[, nonzero, drop = FALSE] - template
  gram <- crossprod(coordinates)
  derivative <- tryCatch(
    -t(solve(gram, crossprod(coordinates, t(difference)))),
    error = function(e) NULL
  )
  if (is.null(derivative) || any(!is.finite(derivative))) return(NULL)

  # Remove the amplitude direction. This makes an aligned change of magnitude
  # live in g rather than in the derivative subspace.
  derivative <- derivative -
    template %o% as.numeric(crossprod(template, derivative))
  tangent_svd <- svd(derivative)
  if (!length(tangent_svd$d) || max(tangent_svd$d) <= 0) return(NULL)
  tangent_rank <- sum(tangent_svd$d > tolerance * max(tangent_svd$d))
  if (!tangent_rank) return(NULL)
  rotation <- tangent_svd$v[, seq_len(tangent_rank), drop = FALSE]
  list(
    derivative = derivative %*% rotation,
    world_basis = world_basis %*% rotation,
    rank = tangent_rank
  )
}

.ldr_gls_feature_coefficients <- function(y, variance, correlation,
                                          background, taper, template,
                                          derivative) {
  y <- as.numeric(y)
  variance <- as.numeric(variance)
  if (any(!is.finite(y)) || any(!is.finite(variance)) ||
      any(variance <= 0)) {
    return(NULL)
  }
  sd <- sqrt(variance / taper)
  covariance <- correlation * outer(sd, sd)
  precision <- tryCatch(
    chol2inv(chol(covariance)),
    error = function(e) NULL
  )
  if (is.null(precision)) return(NULL)
  design <- cbind(background, template, derivative)
  normal <- crossprod(design, precision %*% design)
  if (qr(normal, tol = 1e-9)$rank < ncol(normal)) return(NULL)
  coefficient <- tryCatch(
    as.numeric(solve(normal, crossprod(design, precision %*% y))),
    error = function(e) NULL
  )
  if (is.null(coefficient) || any(!is.finite(coefficient))) return(NULL)
  amplitude_index <- ncol(background) + 1L
  list(
    amplitude = coefficient[amplitude_index],
    tangent = coefficient[seq.int(amplitude_index + 1L,
                                  length(coefficient))]
  )
}

.ldr_screen_activation_fit <- function(patch, template, train) {
  positive_distance <- patch$geometry$shift_distance_mm[
    patch$geometry$shift_distance_mm > sqrt(.Machine$double.eps)
  ]
  shift_radius <- max(patch$geometry$shift_distance_mm)
  tau <- max(min(positive_distance), shift_radius / 2)
  prior <- .ldr_shift_prior(patch$geometry$shift_distance_mm, tau)
  pilot <- vapply(train, function(i) {
    evidence <- .ldr_subject_evidence(
      patch$subjects[[i]], template, patch$geometry$shift_operators,
      prior, prevalence = 0.95, amplitude_scale = 1,
      posterior = TRUE
    )
    if (is.null(evidence)) return(NA_real_)
    sum(evidence$joint_probability * evidence$amplitude_mean)
  }, numeric(1L))
  amplitude_scale <- stats::median(
    pilot[is.finite(pilot) & pilot > sqrt(.Machine$double.eps)]
  )
  if (!is.finite(amplitude_scale) || amplitude_scale <= 0) {
    amplitude_scale <- 1
  }
  list(
    template = template,
    shift_prior = prior,
    prevalence = 0.95,
    amplitude_scale = amplitude_scale
  )
}

.ldr_directional_activation <- function(subject, fit, operators) {
  positive <- .ldr_subject_evidence(
    subject, fit$template, operators, fit$shift_prior,
    fit$prevalence, fit$amplitude_scale
  )
  negative <- .ldr_subject_evidence(
    subject, -fit$template, operators, fit$shift_prior,
    fit$prevalence, fit$amplitude_scale
  )
  if (is.null(positive) || is.null(negative)) NA_real_ else positive - negative
}

.ldr_tangent_center <- function(y_a, y_b, var_a, var_b, geometry,
                                subject_ids, folds = 5L,
                                covariance_shrinkage = 0.25) {
  patch <- .ldr_prepare_patch(y_a, y_b, var_a, var_b, geometry)
  fold_id <- .ldr_tangent_fold_ids(subject_ids, folds)
  n_subject <- length(subject_ids)
  activation <- displacement <- amplitude_cross <-
    tangent_cross <- rep(NA_real_, n_subject)
  delta_a <- delta_b <- matrix(NA_real_, n_subject, 3L)
  fitted_tangent_a <- fitted_tangent_b <- matrix(
    NA_real_, nrow(y_a), n_subject
  )
  fold_receipts <- vector("list", max(fold_id))

  raw_center <- (y_a[geometry$center, ] / var_a[geometry$center, ] +
                   y_b[geometry$center, ] / var_b[geometry$center, ]) /
    (1 / var_a[geometry$center, ] + 1 / var_b[geometry$center, ])

  for (fold in sort(unique(fold_id))) {
    train <- which(fold_id != fold)
    test <- which(fold_id == fold)
    noise <- .ldr_noise_correlation(
      y_a, y_b, var_a, var_b, train,
      shrinkage = covariance_shrinkage
    )
    if (is.null(noise)) next
    template <- .ldr_initial_template(patch, train)
    tangent <- .ldr_tangent_basis(template, geometry)
    if (is.null(tangent)) next
    activation_fit <- .ldr_screen_activation_fit(patch, template, train)
    for (i in test) {
      coefficient_a <- .ldr_gls_feature_coefficients(
        y_a[, i], var_a[, i], noise$correlation,
        patch$background, geometry$taper, template,
        tangent$derivative
      )
      coefficient_b <- .ldr_gls_feature_coefficients(
        y_b[, i], var_b[, i], noise$correlation,
        patch$background, geometry$taper, template,
        tangent$derivative
      )
      if (is.null(coefficient_a) || is.null(coefficient_b)) next
      activation[i] <- .ldr_directional_activation(
        patch$subjects[[i]], activation_fit, geometry$shift_operators
      )
      fitted_a <- as.numeric(
        tangent$derivative %*% coefficient_a$tangent
      )
      fitted_b <- as.numeric(
        tangent$derivative %*% coefficient_b$tangent
      )
      fitted_tangent_a[, i] <- fitted_a
      fitted_tangent_b[, i] <- fitted_b
      amplitude_cross[i] <- coefficient_a$amplitude *
        coefficient_b$amplitude
      if (coefficient_a$amplitude > sqrt(.Machine$double.eps) &&
          coefficient_b$amplitude > sqrt(.Machine$double.eps)) {
        delta_a[i, ] <- -as.numeric(
          tangent$world_basis %*%
            (coefficient_a$tangent / coefficient_a$amplitude)
        )
        delta_b[i, ] <- -as.numeric(
          tangent$world_basis %*%
            (coefficient_b$tangent / coefficient_b$amplitude)
        )
        bound <- max(geometry$shift_distance_mm)
        norm_a <- sqrt(sum(delta_a[i, ]^2))
        norm_b <- sqrt(sum(delta_b[i, ]^2))
        if (norm_a > bound) delta_a[i, ] <- delta_a[i, ] * bound / norm_a
        if (norm_b > bound) delta_b[i, ] <- delta_b[i, ] * bound / norm_b
      }
    }
    fold_receipts[[fold]] <- list(
      train = subject_ids[train],
      test = subject_ids[test],
      template = template,
      tangent_rank = tangent$rank,
      covariance_shrinkage = noise$shrinkage,
      covariance_raw_min_eigenvalue = noise$raw_min_eigenvalue,
      covariance_min_eigenvalue = noise$min_eigenvalue
    )
  }

  tangent_valid <- is.finite(activation) & is.finite(raw_center) &
    apply(fitted_tangent_a, 2L, function(value) all(is.finite(value))) &
    apply(fitted_tangent_b, 2L, function(value) all(is.finite(value)))
  if (sum(tangent_valid) < 3L) return(NULL)
  centered_tangent_a <- centered_tangent_b <- matrix(
    NA_real_, nrow(y_a), n_subject
  )
  for (fold in sort(unique(fold_id[tangent_valid]))) {
    held_out <- which(tangent_valid & fold_id == fold)
    centered_tangent_a[, held_out] <- sweep(
      fitted_tangent_a[, held_out, drop = FALSE], 1L,
      rowMeans(fitted_tangent_a[, held_out, drop = FALSE]), "-"
    )
    centered_tangent_b[, held_out] <- sweep(
      fitted_tangent_b[, held_out, drop = FALSE], 1L,
      rowMeans(fitted_tangent_b[, held_out, drop = FALSE]), "-"
    )
  }
  tangent_cross[tangent_valid] <- colSums(
    centered_tangent_a[, tangent_valid, drop = FALSE] *
      centered_tangent_b[, tangent_valid, drop = FALSE]
  )
  displacement[tangent_valid] <- tangent_cross[tangent_valid]
  tangent_at_center <- rep(NA_real_, n_subject)
  tangent_at_center[tangent_valid] <- colMeans(rbind(
    centered_tangent_a[geometry$center, tangent_valid, drop = FALSE],
    centered_tangent_b[geometry$center, tangent_valid, drop = FALSE]
  ))
  valid <- tangent_valid & is.finite(displacement) &
    is.finite(tangent_at_center)
  if (sum(valid) < 3L) return(NULL)
  counterfactual <- raw_center - tangent_at_center
  valid_delta <- stats::complete.cases(delta_a) &
    stats::complete.cases(delta_b)
  shift_rms_mm <- NA_real_
  minimum_delta <- max(3L, ceiling(0.25 * sum(tangent_valid)))
  if (sum(valid_delta) >= minimum_delta) {
    centered_a <- sweep(
      delta_a[valid_delta, , drop = FALSE], 2L,
      colMeans(delta_a[valid_delta, , drop = FALSE]), "-"
    )
    centered_b <- sweep(
      delta_b[valid_delta, , drop = FALSE], 2L,
      colMeans(delta_b[valid_delta, , drop = FALSE]), "-"
    )
    reproducible_shift <- mean(rowSums(centered_a * centered_b))
    if (is.finite(reproducible_shift) && reproducible_shift > 0) {
      shift_rms_mm <- sqrt(reproducible_shift)
    }
  }
  tangent_energy <- max(mean(tangent_cross[valid], na.rm = TRUE), 0)
  amplitude_energy <- max(mean(amplitude_cross[valid], na.rm = TRUE), 0)
  energy_total <- tangent_energy + amplitude_energy
  rho_shift <- if (energy_total > sqrt(.Machine$double.eps)) {
    tangent_energy / energy_total
  } else 0

  list(
    activation_t = .ldr_t_statistic(activation[valid]),
    displacement_t = .ldr_t_statistic(displacement[valid]),
    ordinary_t = .ldr_t_statistic(raw_center[valid]),
    counterfactual_t = .ldr_t_statistic(counterfactual[valid]),
    mni_loss = .ldr_t_statistic(counterfactual[valid]) -
      .ldr_t_statistic(raw_center[valid]),
    shift_rms_mm = shift_rms_mm,
    rho_shift = rho_shift,
    activation = activation,
    displacement = displacement,
    delta_a = delta_a,
    delta_b = delta_b,
    raw_center = raw_center,
    counterfactual = counterfactual,
    fold_id = fold_id,
    folds = fold_receipts
  )
}

.ldr_search_configurations <- function(patch_radius_mm, shift_radius_mm) {
  if (!is.numeric(patch_radius_mm) || !length(patch_radius_mm) ||
      any(!is.finite(patch_radius_mm)) || any(patch_radius_mm <= 0)) {
    stop("`patch_radius_mm` must contain positive finite radii.",
         call. = FALSE)
  }
  if (!is.numeric(shift_radius_mm) || !length(shift_radius_mm) ||
      any(!is.finite(shift_radius_mm)) || any(shift_radius_mm <= 0)) {
    stop("`shift_radius_mm` must contain positive finite radii.",
         call. = FALSE)
  }
  configurations <- unique(expand.grid(
    patch_radius_mm = sort(unique(as.numeric(patch_radius_mm))),
    shift_radius_mm = sort(unique(as.numeric(shift_radius_mm))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ))
  configurations <- configurations[
    configurations$patch_radius_mm >=
      2 * configurations$shift_radius_mm, , drop = FALSE
  ]
  if (!nrow(configurations)) {
    stop("No patch/shift configuration satisfies patch >= 2 * shift.",
         call. = FALSE)
  }
  configurations$config_id <- seq_len(nrow(configurations))
  configurations
}

.ldr_search_units <- function(centers, configurations, dims3, affine,
                              mask_idx = NULL) {
  active <- mask_idx %||% seq_len(prod(dims3))
  if (is.null(centers)) {
    centers <- active
  } else {
    if (!is.numeric(centers) || !length(centers) ||
        any(!is.finite(centers)) ||
        any(centers != as.integer(centers))) {
      stop("`centers` must contain full-volume voxel indices.",
           call. = FALSE)
    }
    centers <- unique(as.integer(centers))
    if (any(!centers %in% active)) {
      stop("Every requested center must lie in the active mask.",
           call. = FALSE)
    }
  }
  units <- list()
  failures <- 0L
  for (config_row in seq_len(nrow(configurations))) {
    config <- configurations[config_row, , drop = FALSE]
    for (center in centers) {
      geometry <- tryCatch(
        .ldr_patch_geometry(
          center = center,
          dims3 = dims3,
          affine = affine,
          patch_radius_mm = config$patch_radius_mm,
          shift_radius_mm = config$shift_radius_mm,
          mask_idx = mask_idx
        ),
        error = function(e) NULL
      )
      if (is.null(geometry) ||
          !any(geometry$shift_distance_mm > sqrt(.Machine$double.eps))) {
        failures <- failures + 1L
        next
      }
      units[[length(units) + 1L]] <- list(
        center = as.integer(center),
        config_id = config$config_id,
        patch_radius_mm = config$patch_radius_mm,
        shift_radius_mm = config$shift_radius_mm,
        geometry = geometry
      )
    }
  }
  if (!length(units)) {
    stop("The LDR search contains no center with a complete eligible patch.",
         call. = FALSE)
  }
  list(units = units, requested_centers = centers, failures = failures)
}

.ldr_tangent_scan <- function(beta_a, beta_b, var_a, var_b, units,
                              subject_ids, folds,
                              covariance_shrinkage = 0.25,
                              signs = NULL) {
  if (!is.null(signs)) {
    signs <- as.numeric(signs)
    if (length(signs) != length(subject_ids) ||
        any(!signs %in% c(-1, 1))) {
      stop("Tangent sign flips require one -1/+1 value per subject.",
           call. = FALSE)
    }
    beta_a <- sweep(beta_a, 2L, signs, "*")
    beta_b <- sweep(beta_b, 2L, signs, "*")
  }
  fits <- lapply(units, function(unit) {
    index <- unit$geometry$sample_idx
    .ldr_tangent_center(
      beta_a[index, , drop = FALSE],
      beta_b[index, , drop = FALSE],
      var_a[index, , drop = FALSE],
      var_b[index, , drop = FALSE],
      geometry = unit$geometry,
      subject_ids = subject_ids,
      folds = folds,
      covariance_shrinkage = covariance_shrinkage
    )
  })
  extract <- function(name) vapply(
    fits,
    function(fit) if (is.null(fit)) NA_real_ else fit[[name]],
    numeric(1L)
  )
  subject_matrix <- function(name) {
    vapply(
      fits,
      function(fit) {
        if (is.null(fit)) rep(NA_real_, length(subject_ids)) else fit[[name]]
      },
      numeric(length(subject_ids))
    )
  }
  list(
    fits = fits,
    activation_t = extract("activation_t"),
    displacement_t = extract("displacement_t"),
    ordinary_t = extract("ordinary_t"),
    counterfactual_t = extract("counterfactual_t"),
    mni_loss = extract("mni_loss"),
    shift_rms_mm = extract("shift_rms_mm"),
    rho_shift = extract("rho_shift"),
    activation = subject_matrix("activation"),
    displacement = subject_matrix("displacement")
  )
}

.ldr_maximum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) -Inf else max(x)
}

.ldr_maximum_p <- function(observed, null_maximum) {
  null_maximum <- null_maximum[is.finite(null_maximum)]
  vapply(observed, function(value) {
    if (is.na(value) || !length(null_maximum)) return(NA_real_)
    (1 + sum(null_maximum >= value)) / (length(null_maximum) + 1)
  }, numeric(1L))
}

.ldr_wild_displacement <- function(contributions, multipliers) {
  vapply(seq_len(ncol(contributions)), function(column) {
    value <- contributions[, column]
    valid <- is.finite(value)
    if (sum(valid) < 3L) return(NA_real_)
    centered <- value[valid] - mean(value[valid])
    .ldr_t_statistic(multipliers[valid] * centered)
  }, numeric(1L))
}

.ldr_tangent_calibrate <- function(beta_a, beta_b, var_a, var_b, units,
                                   subject_ids, folds,
                                   covariance_shrinkage, observed,
                                   n_resamples, seed) {
  n_resamples <- as.integer(n_resamples)
  if (length(n_resamples) != 1L || is.na(n_resamples) ||
      n_resamples < 1L) {
    stop("`n_resamples` must be positive for whole-search calibration.",
         call. = FALSE)
  }
  n_unit <- length(units)
  activation_null <- ordinary_null <- displacement_null <-
    matrix(NA_real_, nrow = n_resamples, ncol = n_unit)
  .ldr_with_seed(seed, {
    for (b in seq_len(n_resamples)) {
      multipliers <- sample(c(-1, 1), length(subject_ids), replace = TRUE)
      flipped <- .ldr_tangent_scan(
        beta_a, beta_b, var_a, var_b, units,
        subject_ids = subject_ids,
        folds = folds,
        covariance_shrinkage = covariance_shrinkage,
        signs = multipliers
      )
      activation_null[b, ] <- flipped$activation_t
      ordinary_null[b, ] <- flipped$ordinary_t
      displacement_null[b, ] <- .ldr_wild_displacement(
        observed$displacement, multipliers
      )
    }
  })
  activation_max <- apply(activation_null, 1L, .ldr_maximum)
  displacement_max <- apply(displacement_null, 1L, .ldr_maximum)
  ordinary_max <- apply(ordinary_null, 1L, .ldr_maximum)
  point_p <- function(observed_value, null) {
    vapply(seq_along(observed_value), function(column) {
      .ldr_resample_p(observed_value[column], null[, column])
    }, numeric(1L))
  }
  activation_point <- point_p(observed$activation_t, activation_null)
  displacement_point <- point_p(
    observed$displacement_t, displacement_null
  )
  ordinary_point <- point_p(observed$ordinary_t, ordinary_null)
  activation_fwer <- .ldr_maximum_p(
    observed$activation_t, activation_max
  )
  displacement_fwer <- .ldr_maximum_p(
    observed$displacement_t, displacement_max
  )
  ordinary_fwer <- .ldr_maximum_p(observed$ordinary_t, ordinary_max)
  if (any(activation_fwer + 1e-12 < activation_point, na.rm = TRUE) ||
      any(displacement_fwer + 1e-12 < displacement_point, na.rm = TRUE) ||
      any(ordinary_fwer + 1e-12 < ordinary_point, na.rm = TRUE)) {
    stop("Internal LDR invariant failed: FWER p < pointwise p.",
         call. = FALSE)
  }
  list(
    activation_fwer_p = activation_fwer,
    displacement_fwer_p = displacement_fwer,
    ordinary_fwer_p = ordinary_fwer,
    ldr_fwer_p = pmax(activation_fwer, displacement_fwer),
    pointwise_p = list(
      activation = activation_point,
      displacement = displacement_point,
      ordinary = ordinary_point
    ),
    null_maximum = list(
      activation = activation_max,
      displacement = displacement_max,
      ordinary = ordinary_max
    ),
    n_resamples = n_resamples,
    seed = as.integer(seed),
    minimum_p = 1 / (n_resamples + 1)
  )
}
