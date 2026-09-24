# Local displacement rescue: exact patch reference core -------------------

.ldr_logsumexp <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(-Inf)
  finite <- is.finite(x)
  if (!any(finite)) return(max(x))
  anchor <- max(x[finite])
  anchor + log(sum(exp(x - anchor)))
}

.ldr_halfnormal_log_bf <- function(u, q, scale) {
  if (!is.numeric(scale) || length(scale) != 1L ||
      !is.finite(scale) || scale <= 0) {
    stop("`scale` must be one positive finite value.", call. = FALSE)
  }
  u <- as.numeric(u)
  q <- as.numeric(q)
  if (length(u) != length(q) || any(!is.finite(u)) ||
      any(!is.finite(q)) || any(q < 0)) {
    stop("`u` and `q` must be matching finite vectors with q >= 0.",
         call. = FALSE)
  }
  precision <- q + 1 / scale^2
  standardized <- u / sqrt(precision)
  log(2) - log(scale) - 0.5 * log(precision) +
    u^2 / (2 * precision) +
    stats::pnorm(standardized, log.p = TRUE)
}

.ldr_halfnormal_moments <- function(u, q, scale) {
  log_bf <- .ldr_halfnormal_log_bf(u, q, scale)
  precision <- q + 1 / scale^2
  sd <- 1 / sqrt(precision)
  mean <- u / precision
  standardized <- mean / sd
  log_mills <- stats::dnorm(standardized, log = TRUE) -
    stats::pnorm(standardized, log.p = TRUE)
  mills <- exp(pmin(log_mills, log(.Machine$double.xmax)))
  posterior_mean <- mean + sd * mills
  posterior_second <- mean^2 + sd^2 + mean * sd * mills
  list(
    log_bf = log_bf,
    mean = pmax(posterior_mean, 0),
    second = pmax(posterior_second, 0)
  )
}

.ldr_normal_moments <- function(u, q, scale) {
  if (!is.numeric(scale) || length(scale) != 1L ||
      !is.finite(scale) || scale <= 0) {
    stop("`scale` must be one positive finite value.", call. = FALSE)
  }
  u <- as.numeric(u)
  q <- as.numeric(q)
  if (length(u) != length(q) || any(!is.finite(u)) ||
      any(!is.finite(q)) || any(q < 0)) {
    stop("`u` and `q` must be matching finite vectors with q >= 0.",
         call. = FALSE)
  }
  precision <- q + 1 / scale^2
  mean <- u / precision
  list(
    log_bf = -log(scale) - 0.5 * log(precision) +
      u^2 / (2 * precision),
    mean = mean,
    second = 1 / precision + mean^2
  )
}

.ldr_scalar <- function(x, name, lower = -Inf, upper = Inf,
                        lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", name, "` must be one finite numeric value.", call. = FALSE)
  }
  lower_bad <- if (lower_open) x <= lower else x < lower
  upper_bad <- if (upper_open) x >= upper else x > upper
  if (lower_bad || upper_bad) {
    stop("`", name, "` lies outside its permitted range.", call. = FALSE)
  }
  as.numeric(x)
}

.ldr_mm_offsets <- function(dims3, affine, radius_mm) {
  radius_mm <- .ldr_scalar(radius_mm, "radius_mm", lower = 0)
  dims3 <- as.integer(dims3)
  linear <- affine[seq_len(3L), seq_len(3L), drop = FALSE]
  axis_mm <- sqrt(colSums(linear^2))
  if (any(!is.finite(axis_mm)) || any(axis_mm <= 0)) {
    stop("Voxel affine must have three finite nonzero spatial axes.",
         call. = FALSE)
  }
  bounds <- pmin(
    pmax(dims3 - 1L, 0L),
    as.integer(ceiling(radius_mm / axis_mm))
  )
  values <- lapply(seq_len(3L), function(axis) {
    if (dims3[axis] <= 1L) 0L else seq.int(-bounds[axis], bounds[axis])
  })
  offsets <- as.matrix(expand.grid(
    dx = values[[1L]], dy = values[[2L]], dz = values[[3L]]
  ))
  world <- t(linear %*% t(offsets))
  distance_mm <- sqrt(rowSums(world^2))
  keep <- distance_mm <= radius_mm + sqrt(.Machine$double.eps)
  offsets <- offsets[keep, , drop = FALSE]
  distance_mm <- distance_mm[keep]
  order <- order(distance_mm, offsets[, 3L], offsets[, 2L], offsets[, 1L])
  list(
    voxel = matrix(as.integer(offsets[order, , drop = FALSE]), ncol = 3L,
                   dimnames = list(NULL, c("dx", "dy", "dz"))),
    world = world[keep, , drop = FALSE][order, , drop = FALSE],
    distance_mm = distance_mm[order]
  )
}

.ldr_patch_geometry <- function(center, dims3, affine, patch_radius_mm,
                                shift_radius_mm, mask_idx = NULL) {
  dims3 <- as.integer(dims3)
  n_full <- prod(dims3)
  center <- as.integer(center)
  if (length(center) != 1L || is.na(center) || center < 1L ||
      center > n_full) {
    stop("`center` must identify one full-volume voxel.", call. = FALSE)
  }
  patch <- .ldr_mm_offsets(dims3, affine, patch_radius_mm)
  shifts <- .ldr_mm_offsets(dims3, affine, shift_radius_mm)
  center_coord <- as.integer(arrayInd(center, dims3))
  coords <- sweep(patch$voxel, 2L, center_coord, "+")
  inside <- rowSums(coords < 1L | coords > rep(dims3, each = nrow(coords))) == 0L
  if (!all(inside)) {
    stop("The requested LDR patch crosses the image boundary.", call. = FALSE)
  }
  full_idx <- coords[, 1L] +
    (coords[, 2L] - 1L) * dims3[1L] +
    (coords[, 3L] - 1L) * dims3[1L] * dims3[2L]
  active <- mask_idx %||% seq_len(n_full)
  sample_idx <- match(full_idx, active)
  if (anyNA(sample_idx)) {
    stop("The requested LDR patch contains voxels outside the active mask.",
         call. = FALSE)
  }

  offset_key <- function(x) paste(x[, 1L], x[, 2L], x[, 3L], sep = ":")
  keys <- offset_key(patch$voxel)
  shift_index <- vapply(seq_len(nrow(shifts$voxel)), function(j) {
    source <- sweep(patch$voxel, 2L, shifts$voxel[j, ], "-")
    match(offset_key(source), keys)
  }, integer(nrow(patch$voxel)))
  if (is.null(dim(shift_index))) {
    shift_index <- matrix(shift_index, ncol = 1L)
  }
  shift_operators <- lapply(seq_len(ncol(shift_index)), function(j) {
    operator <- matrix(0, nrow = nrow(shift_index), ncol = nrow(shift_index))
    valid <- which(!is.na(shift_index[, j]))
    operator[cbind(valid, shift_index[valid, j])] <- 1
    operator
  })
  taper <- if (patch_radius_mm <= 0) {
    rep(1, length(patch$distance_mm))
  } else {
    0.5 * (1 + cos(pi * patch$distance_mm / patch_radius_mm))
  }
  taper <- pmax(taper, 1e-6)
  center_in_patch <- which(rowSums(abs(patch$voxel)) == 0L)
  # The v0.1 orientation constraint is deliberately narrow: only the anchor
  # voxel is forced nonnegative. Neighbouring negative lobes are part of the
  # signed feature and must not be rectified away.
  core <- center_in_patch
  list(
    full_idx = as.integer(full_idx),
    sample_idx = as.integer(sample_idx),
    offsets = patch$voxel,
    coordinates_mm = patch$world,
    distance_mm = patch$distance_mm,
    taper = taper,
    center = center_in_patch,
    core = core,
    shift_offsets = shifts$voxel,
    shift_world = shifts$world,
    shift_distance_mm = shifts$distance_mm,
    shift_index = shift_index,
    shift_operators = shift_operators
  )
}

.ldr_background <- function(coordinates_mm) {
  coordinates_mm <- as.matrix(coordinates_mm)
  out <- matrix(1, nrow = nrow(coordinates_mm), ncol = 1L,
                dimnames = list(NULL, "intercept"))
  for (axis in seq_len(ncol(coordinates_mm))) {
    value <- coordinates_mm[, axis]
    if (stats::sd(value) > sqrt(.Machine$double.eps)) {
      value <- (value - mean(value)) / stats::sd(value)
      out <- cbind(out, value)
      colnames(out)[ncol(out)] <- c("x", "y", "z")[axis]
    }
  }
  out
}

.ldr_projection <- function(y, variance, background, taper) {
  ok <- is.finite(y) & is.finite(variance) & variance > 0 &
    is.finite(taper) & taper > 0
  if (sum(ok) <= ncol(background)) {
    return(NULL)
  }
  weight <- numeric(length(y))
  weight[ok] <- taper[ok] / variance[ok]
  y_observed <- y
  y_observed[!ok] <- 0
  bwb <- crossprod(background, weight * background)
  inverse <- tryCatch(
    solve(bwb),
    error = function(e) NULL
  )
  if (is.null(inverse)) return(NULL)
  projection <- diag(length(y)) -
    background %*% inverse %*% crossprod(background, diag(weight))
  residual <- as.numeric(projection %*% y_observed)
  list(weight = weight, projection = projection, residual = residual)
}

.ldr_prepare_patch <- function(y_a, y_b, var_a, var_b, geometry) {
  y_a <- as.matrix(y_a)
  y_b <- as.matrix(y_b)
  var_a <- as.matrix(var_a)
  var_b <- as.matrix(var_b)
  if (!identical(dim(y_a), dim(y_b)) || !identical(dim(y_a), dim(var_a)) ||
      !identical(dim(y_a), dim(var_b))) {
    stop("Paired LDR patch effects and variances must have matching dimensions.",
         call. = FALSE)
  }
  if (nrow(y_a) != length(geometry$sample_idx)) {
    stop("LDR patch data do not match the patch geometry.", call. = FALSE)
  }
  background <- .ldr_background(geometry$coordinates_mm)
  subjects <- vector("list", ncol(y_a))
  for (i in seq_len(ncol(y_a))) {
    a <- .ldr_projection(y_a[, i], var_a[, i], background, geometry$taper)
    b <- .ldr_projection(y_b[, i], var_b[, i], background, geometry$taper)
    subjects[[i]] <- if (is.null(a) || is.null(b)) NULL else list(a = a, b = b)
  }
  list(
    y_a = y_a,
    y_b = y_b,
    var_a = var_a,
    var_b = var_b,
    background = background,
    geometry = geometry,
    subjects = subjects
  )
}

.ldr_shift_templates <- function(template, operators) {
  vapply(operators, function(operator) as.numeric(operator %*% template),
         numeric(length(template)))
}

.ldr_subject_sufficient <- function(subject, shifted_templates) {
  if (is.null(subject)) return(NULL)
  calculate <- function(split) {
    projected <- split$projection %*% shifted_templates
    list(
      u = colSums(projected * (split$weight * split$residual)),
      q = colSums(projected * (split$weight * projected))
    )
  }
  a <- calculate(subject$a)
  b <- calculate(subject$b)
  list(u = a$u + b$u, q = pmax(a$q + b$q, 0))
}

.ldr_shift_prior <- function(distance_mm, tau_mm) {
  tau_mm <- .ldr_scalar(tau_mm, "tau_mm", lower = 0)
  distance_mm <- as.numeric(distance_mm)
  if (tau_mm == 0) {
    probability <- as.numeric(distance_mm == min(distance_mm))
    probability <- probability / sum(probability)
    return(probability)
  }
  log_probability <- -0.5 * (distance_mm / tau_mm)^2
  exp(log_probability - .ldr_logsumexp(log_probability))
}

.ldr_subject_evidence <- function(subject, template, operators, shift_prior,
                                  prevalence, amplitude_scale,
                                  amplitude_family = "positive",
                                  posterior = FALSE) {
  if (is.null(subject)) return(NULL)
  shifted <- .ldr_shift_templates(template, operators)
  sufficient <- .ldr_subject_sufficient(subject, shifted)
  moments <- switch(
    amplitude_family,
    positive = .ldr_halfnormal_moments(
      sufficient$u, sufficient$q, amplitude_scale
    ),
    signed = .ldr_normal_moments(
      sufficient$u, sufficient$q, amplitude_scale
    ),
    stop("Unknown LDR amplitude family.", call. = FALSE)
  )
  log_active_shift <- log(shift_prior) + moments$log_bf
  log_active <- .ldr_logsumexp(log_active_shift)
  log_null_weight <- log1p(-prevalence)
  log_active_weight <- log(prevalence) + log_active
  log_mixture <- .ldr_logsumexp(c(log_null_weight, log_active_weight))
  if (!posterior) return(log_mixture)
  active_probability <- exp(log_active_weight - log_mixture)
  shift_probability <- exp(log_active_shift - log_active)
  list(
    log_bf = log_mixture,
    active_probability = active_probability,
    shift_probability = shift_probability,
    joint_probability = active_probability * shift_probability,
    amplitude_mean = moments$mean,
    amplitude_second = moments$second,
    shifted_templates = shifted,
    u = sufficient$u,
    q = sufficient$q
  )
}

.ldr_initial_template <- function(patch, indices) {
  numerator <- numeric(nrow(patch$y_a))
  denominator <- numeric(nrow(patch$y_a))
  for (i in indices) {
    subject <- patch$subjects[[i]]
    if (is.null(subject)) next
    for (split in subject) {
      numerator <- numerator + split$weight * split$residual
      denominator <- denominator + split$weight
    }
  }
  template <- numerator / denominator
  template[!is.finite(template)] <- 0
  template[patch$geometry$core] <- pmax(template[patch$geometry$core], 0)
  if (sum(template[patch$geometry$core]) <= sqrt(.Machine$double.eps)) {
    template[] <- 0
    template[patch$geometry$center] <- 1
  }
  norm <- sqrt(sum(template^2))
  template / norm
}

.ldr_update_template <- function(patch, indices, fit, ridge) {
  n_voxel <- length(fit$template)
  normal <- diag(ridge, n_voxel)
  rhs <- numeric(n_voxel)
  active_total <- 0
  second_total <- 0
  for (i in indices) {
    evidence <- .ldr_subject_evidence(
      patch$subjects[[i]], fit$template, patch$geometry$shift_operators,
      fit$shift_prior, fit$prevalence, fit$amplitude_scale,
      amplitude_family = fit$amplitude_family,
      posterior = TRUE
    )
    if (is.null(evidence)) next
    active_total <- active_total + evidence$active_probability
    second_total <- second_total +
      sum(evidence$joint_probability * evidence$amplitude_second)
    for (d in seq_along(patch$geometry$shift_operators)) {
      probability <- evidence$joint_probability[d]
      if (!is.finite(probability) || probability <= 1e-10) next
      operator <- patch$geometry$shift_operators[[d]]
      for (split in patch$subjects[[i]]) {
        design <- split$projection %*% operator
        normal <- normal + probability * evidence$amplitude_second[d] *
          crossprod(design, split$weight * design)
        rhs <- rhs + probability * evidence$amplitude_mean[d] *
          as.numeric(crossprod(design, split$weight * split$residual))
      }
    }
  }
  template <- tryCatch(
    as.numeric(solve(normal, rhs)),
    error = function(e) fit$template
  )
  template[!is.finite(template)] <- 0
  template[patch$geometry$core] <- pmax(template[patch$geometry$core], 0)
  norm <- sqrt(sum(template^2))
  if (!is.finite(norm) || norm <= sqrt(.Machine$double.eps)) {
    template <- fit$template
    norm <- 1
  } else {
    template <- template / norm
  }
  valid_n <- sum(vapply(indices, function(i) !is.null(patch$subjects[[i]]),
                        logical(1L)))
  prevalence <- if (valid_n) active_total / valid_n else fit$prevalence
  prevalence <- min(0.995, max(0.01, prevalence))
  amplitude_scale <- if (active_total > sqrt(.Machine$double.eps)) {
    sqrt(second_total / active_total) * norm
  } else fit$amplitude_scale
  if (!is.finite(amplitude_scale) || amplitude_scale <= 0) {
    amplitude_scale <- fit$amplitude_scale
  }
  list(
    template = template,
    prevalence = prevalence,
    amplitude_scale = amplitude_scale
  )
}

.ldr_model_score <- function(patch, indices, fit) {
  scores <- vapply(indices, function(i) {
    value <- .ldr_subject_evidence(
      patch$subjects[[i]], fit$template, patch$geometry$shift_operators,
      fit$shift_prior, fit$prevalence, fit$amplitude_scale,
      amplitude_family = fit$amplitude_family
    )
    value %||% NA_real_
  }, numeric(1L))
  sum(scores, na.rm = TRUE)
}

.ldr_fit_tau <- function(patch, indices, tau_mm, iterations = 20L,
                         ridge = 1e-4, tolerance = 1e-5,
                         amplitude_family = "positive") {
  template <- .ldr_initial_template(patch, indices)
  projections <- vapply(indices, function(i) {
    evidence <- .ldr_subject_evidence(
      patch$subjects[[i]], template, patch$geometry$shift_operators,
      .ldr_shift_prior(patch$geometry$shift_distance_mm, tau_mm),
      prevalence = 0.75, amplitude_scale = 1,
      amplitude_family = amplitude_family,
      posterior = TRUE
    )
    if (is.null(evidence)) return(NA_real_)
    max(evidence$amplitude_mean, na.rm = TRUE)
  }, numeric(1L))
  amplitude_scale <- stats::median(projections[is.finite(projections) &
                                                projections > 0])
  if (!is.finite(amplitude_scale)) amplitude_scale <- 1
  fit <- list(
    template = template,
    prevalence = 0.75,
    amplitude_scale = amplitude_scale,
    amplitude_family = amplitude_family,
    tau_mm = tau_mm,
    shift_prior = .ldr_shift_prior(
      patch$geometry$shift_distance_mm, tau_mm
    )
  )
  converged <- FALSE
  previous <- -Inf
  for (iteration in seq_len(as.integer(iterations))) {
    updated <- .ldr_update_template(patch, indices, fit, ridge)
    fit$template <- updated$template
    fit$prevalence <- updated$prevalence
    fit$amplitude_scale <- updated$amplitude_scale
    score <- .ldr_model_score(patch, indices, fit)
    if (is.finite(previous) && abs(score - previous) <=
        tolerance * (1 + abs(previous))) {
      converged <- TRUE
      break
    }
    previous <- score
  }
  fit$score <- .ldr_model_score(patch, indices, fit)
  fit$converged <- converged
  fit$iterations <- iteration
  fit$shift_rms_mm <- sqrt(sum(
    fit$shift_prior * patch$geometry$shift_distance_mm^2
  ))
  fit
}

.ldr_fit_model <- function(patch, indices, tau_grid_mm, iterations = 20L,
                           ridge = 1e-4, tolerance = 1e-5,
                           amplitude_family = "positive") {
  fits <- lapply(tau_grid_mm, function(tau) {
    .ldr_fit_tau(
      patch, indices, tau,
      iterations = iterations,
      ridge = ridge,
      tolerance = tolerance,
      amplitude_family = amplitude_family
    )
  })
  scores <- vapply(fits, `[[`, numeric(1L), "score")
  fits[[which.max(scores)]]
}

.ldr_fold_ids <- function(subject_ids, folds) {
  n <- length(subject_ids)
  folds <- min(as.integer(folds), n)
  if (folds < 2L) stop("LDR cross-fitting requires at least two folds.",
                       call. = FALSE)
  sorted <- order(as.character(subject_ids))
  assigned <- rep(seq_len(folds), length.out = n)
  out <- integer(n)
  out[sorted] <- assigned
  out
}

.ldr_t_statistic <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  scale <- stats::sd(x)
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
    return(if (mean(x) > 0) Inf else if (mean(x) < 0) -Inf else 0)
  }
  mean(x) / (scale / sqrt(length(x)))
}

.ldr_crossfit <- function(patch, subject_ids, tau_grid_mm, folds = 5L,
                          iterations = 20L, ridge = 1e-4,
                          tolerance = 1e-5) {
  fold_id <- .ldr_fold_ids(subject_ids, folds)
  n_subject <- length(subject_ids)
  aligned_log_bf <- heterogeneity_log_bf <- jitter_log_bf <-
    rep(NA_real_, n_subject)
  observed_signal <- aligned_signal <- rep(NA_real_, n_subject)
  fold_receipts <- vector("list", max(fold_id))
  for (fold in sort(unique(fold_id))) {
    train <- which(fold_id != fold)
    test <- which(fold_id == fold)
    aligned <- .ldr_fit_model(
      patch, train, tau_grid_mm = 0,
      iterations = iterations, ridge = ridge, tolerance = tolerance,
      amplitude_family = "positive"
    )
    heterogeneity <- .ldr_fit_model(
      patch, train, tau_grid_mm = 0,
      iterations = iterations, ridge = ridge, tolerance = tolerance,
      amplitude_family = "signed"
    )
    jitter <- .ldr_fit_model(
      patch, train, tau_grid_mm = tau_grid_mm,
      iterations = iterations, ridge = ridge, tolerance = tolerance
    )
    for (i in test) {
      aligned_log_bf[i] <- .ldr_subject_evidence(
        patch$subjects[[i]], aligned$template,
        patch$geometry$shift_operators, aligned$shift_prior,
        aligned$prevalence, aligned$amplitude_scale
      ) %||% NA_real_
      heterogeneity_log_bf[i] <- .ldr_subject_evidence(
        patch$subjects[[i]], heterogeneity$template,
        patch$geometry$shift_operators, heterogeneity$shift_prior,
        heterogeneity$prevalence, heterogeneity$amplitude_scale,
        amplitude_family = "signed"
      ) %||% NA_real_
      evidence <- .ldr_subject_evidence(
        patch$subjects[[i]], jitter$template,
        patch$geometry$shift_operators, jitter$shift_prior,
        jitter$prevalence, jitter$amplitude_scale,
        amplitude_family = "positive",
        posterior = TRUE
      )
      if (is.null(evidence)) next
      jitter_log_bf[i] <- evidence$log_bf
      signal_weight <- evidence$joint_probability * evidence$amplitude_mean
      observed_signal[i] <- sum(
        signal_weight * evidence$shifted_templates[patch$geometry$center, ]
      )
      aligned_signal[i] <- sum(signal_weight) *
        jitter$template[patch$geometry$center]
    }
    fold_receipts[[fold]] <- list(
      train = subject_ids[train],
      test = subject_ids[test],
      aligned = aligned,
      heterogeneity = heterogeneity,
      jitter = jitter
    )
  }
  valid <- is.finite(aligned_log_bf) & is.finite(heterogeneity_log_bf) &
    is.finite(jitter_log_bf)
  if (sum(valid) < 3L) {
    stop("Too few subjects produced finite held-out LDR scores.", call. = FALSE)
  }
  raw_center <- (patch$y_a[patch$geometry$center, ] /
                   patch$var_a[patch$geometry$center, ] +
                   patch$y_b[patch$geometry$center, ] /
                   patch$var_b[patch$geometry$center, ]) /
    (1 / patch$var_a[patch$geometry$center, ] +
       1 / patch$var_b[patch$geometry$center, ])
  counterfactual <- raw_center + aligned_signal - observed_signal
  selected_rms <- vapply(fold_receipts, function(receipt) {
    receipt$jitter$shift_rms_mm
  }, numeric(1L))
  selected_prevalence <- vapply(fold_receipts, function(receipt) {
    receipt$jitter$prevalence
  }, numeric(1L))
  aligned_difference <- jitter_log_bf[valid] - aligned_log_bf[valid]
  heterogeneity_difference <-
    jitter_log_bf[valid] - heterogeneity_log_bf[valid]
  displacement_components <- c(
    aligned = 2 * sum(aligned_difference),
    heterogeneity = 2 * sum(heterogeneity_difference)
  )
  displacement_t_components <- c(
    aligned = .ldr_t_statistic(aligned_difference),
    heterogeneity = .ldr_t_statistic(heterogeneity_difference)
  )
  list(
    activation_score = 2 * sum(jitter_log_bf[valid]),
    displacement_score = min(displacement_components),
    displacement_components = displacement_components,
    activation_t = .ldr_t_statistic(jitter_log_bf[valid]),
    displacement_t = min(displacement_t_components),
    displacement_t_components = displacement_t_components,
    ordinary_t = .ldr_t_statistic(raw_center),
    counterfactual_t = .ldr_t_statistic(counterfactual),
    mni_loss = .ldr_t_statistic(counterfactual) -
      .ldr_t_statistic(raw_center),
    shift_rms_mm = mean(selected_rms),
    prevalence = mean(selected_prevalence),
    aligned_log_bf = aligned_log_bf,
    heterogeneity_log_bf = heterogeneity_log_bf,
    jitter_log_bf = jitter_log_bf,
    observed_signal = observed_signal,
    aligned_signal = aligned_signal,
    fold_id = fold_id,
    folds = fold_receipts
  )
}
