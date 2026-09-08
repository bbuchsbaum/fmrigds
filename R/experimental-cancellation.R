# Experimental cancellation and spatial-dispersion diagnostics ------------

.experimental_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                 lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", name, "` must be one finite numeric value.", call. = FALSE)
  }
  lower_bad <- if (lower_open) x <= lower else x < lower
  upper_bad <- if (upper_open) x >= upper else x > upper
  if (lower_bad || upper_bad) {
    interval <- paste0(
      if (lower_open) "(" else "[", lower, ", ", upper,
      if (upper_open) ")" else "]"
    )
    stop("`", name, "` must lie in ", interval, ".", call. = FALSE)
  }
  as.numeric(x)
}

.experimental_variance <- function(x) {
  values <- assays(x)
  if (!is.null(values$var)) return(values$var)
  if (!is.null(values$se)) return(values$se^2)
  stop("Experimental cancellation requires genuine `var` or `se` assays.", call. = FALSE)
}

.experimental_validate_input <- function(x) {
  if (!inherits(x, "gds")) {
    stop("`x` must be a realised GDS object.", call. = FALSE)
  }
  if (!inherits(space(x), "space_voxel")) {
    stop("Experimental spatial rescue requires a voxel-space GDS input.", call. = FALSE)
  }
  beta <- assay(x, "beta")
  if (is.null(beta)) {
    stop("Experimental cancellation requires a `beta` assay.", call. = FALSE)
  }
  variance <- .experimental_variance(x)
  if (isTRUE(metadata(x)$synthetic_var) ||
      isTRUE(attr(variance, "synthetic_unit_variance"))) {
    stop(
      "Experimental cancellation requires genuine first-level `var` or `se`; ",
      "synthetic unit variance is not valid.",
      call. = FALSE
    )
  }
  if (!identical(dim(beta), dim(variance))) {
    stop("Beta and variance dimensions must match.", call. = FALSE)
  }
  invisible(TRUE)
}

.experimental_probability_kernel <- function(beta, variance, delta,
                                             equivalence, prevalence,
                                             min_subjects) {
  fit <- core_meta_re_dl_kernel(
    beta = beta,
    var = variance,
    opts = list(min_subjects = min_subjects, alternative = "two.sided")
  )
  tau <- sqrt(pmax(fit$tau2, 0))
  prevalence_quantile <- stats::qnorm(prevalence)

  # Under theta ~ N(mu, tau^2), both tail prevalences exceed `prevalence`
  # exactly when mu lies in [sign_lower, sign_upper]. Intersect this with
  # practical equivalence for mu and integrate the plug-in normal posterior.
  sign_lower <- delta + tau * prevalence_quantile
  sign_upper <- -delta - tau * prevalence_quantile
  lower <- pmax(sign_lower, -equivalence)
  upper <- pmin(sign_upper, equivalence)
  se <- fit$se_g

  probability <- rep(NA_real_, length(lower))
  valid <- is.finite(tau) & is.finite(fit$beta_g) &
    is.finite(se) & se > 0 & is.finite(lower) & is.finite(upper)
  probability[valid & tau <= 0] <- 0
  possible <- valid & tau > 0 & upper > lower
  probability[valid & tau > 0 & !possible] <- 0
  probability[possible] <-
    stats::pnorm((upper[possible] - fit$beta_g[possible]) / se[possible]) -
    stats::pnorm((lower[possible] - fit$beta_g[possible]) / se[possible])
  probability <- pmin(1, pmax(0, probability))
  probability[fit$n_eff < min_subjects] <- NA_real_
  probability
}

.experimental_meaningful_probability <- function(beta, variance, delta) {
  ok <- is.finite(beta) & is.finite(variance) & variance > 0
  se <- sqrt(ifelse(ok, variance, NA_real_))
  probability <- stats::pnorm((-delta - beta) / se) +
    stats::pnorm((beta - delta) / se)
  probability[!ok] <- NA_real_
  probability
}

.experimental_local_support <- function(probability, spec, radius,
                                        support_probability, min_subjects) {
  offsets <- .experimental_offset_grid(radius, spec$dim)
  source_maps <- lapply(seq_len(nrow(offsets)), function(j) {
    .experimental_shift_source(
      spec$mask_idx, offsets[j, ], spec$dim
    )
  })
  local_cutoff <- 1 - (1 - support_probability) / nrow(offsets)
  n_valid <- integer(length(spec$mask_idx))
  n_detected <- integer(length(spec$mask_idx))
  for (subject in seq_len(ncol(probability))) {
    full <- rep.int(NA_real_, spec$n_full)
    full[spec$mask_idx] <- probability[, subject]
    local <- rep.int(NA_real_, length(spec$mask_idx))
    for (source in source_maps) {
      candidate <- rep.int(NA_real_, length(source))
      valid <- is.finite(source)
      candidate[valid] <- full[source[valid]]
      local <- pmax(local, candidate, na.rm = TRUE)
    }
    valid <- is.finite(local)
    n_valid <- n_valid + valid
    n_detected <- n_detected + (valid & local >= local_cutoff)
  }
  support <- n_detected / n_valid
  support[n_valid < min_subjects | !is.finite(support)] <- NA_real_
  support
}

.experimental_voxel_spec <- function(x) {
  sp <- space(x)
  dims3 <- as.integer(sp$dim)
  n_full <- prod(dims3)
  mask_idx <- sp$mask_idx %||% seq_len(n_full)
  if (length(mask_idx) != dim(assay(x, "beta"))[1L]) {
    stop("Voxel mask length does not match the GDS sample dimension.", call. = FALSE)
  }
  list(dim = dims3, n_full = n_full, mask_idx = as.integer(mask_idx))
}

.experimental_unpack_matrix <- function(x, contrast, spec) {
  out <- matrix(NA_real_, nrow = spec$n_full, ncol = dim(x)[2L])
  values <- matrix(
    x[, , contrast], nrow = dim(x)[1L], ncol = dim(x)[2L]
  )
  out[spec$mask_idx, ] <- values
  out
}

.experimental_pack_vector <- function(x, spec) as.numeric(x[spec$mask_idx])

.experimental_offset_grid <- function(radius, dims3 = NULL) {
  values <- seq.int(-radius, radius)
  axis_values <- if (is.null(dims3)) {
    rep(list(values), 3L)
  } else {
    lapply(as.integer(dims3), function(size) if (size <= 1L) 0L else values)
  }
  offsets <- as.matrix(expand.grid(
    dx = axis_values[[1L]], dy = axis_values[[2L]], dz = axis_values[[3L]]
  ))
  distance2 <- rowSums(offsets^2)
  offsets[order(distance2, offsets[, 3L], offsets[, 2L], offsets[, 1L]), , drop = FALSE]
}

.experimental_shift_source <- function(target_idx, offset, dims3) {
  coords <- arrayInd(target_idx, dims3)
  source <- sweep(coords, 2L, as.integer(offset), "+")
  valid <- rowSums(source < 1L | source > rep(dims3, each = nrow(source))) == 0L
  source_idx <- rep.int(NA_integer_, length(target_idx))
  source_idx[valid] <- source[valid, 1L] +
    (source[valid, 2L] - 1L) * dims3[1L] +
    (source[valid, 3L] - 1L) * dims3[1L] * dims3[2L]
  source_idx
}

.experimental_dilate <- function(idx, dims3, radius) {
  if (!length(idx) || radius <= 0L) return(as.integer(idx))
  offsets <- .experimental_offset_grid(radius, dims3)
  expanded <- lapply(seq_len(nrow(offsets)), function(j) {
    .experimental_shift_source(idx, offsets[j, ], dims3)
  })
  sort(unique(as.integer(stats::na.omit(unlist(expanded, use.names = FALSE)))))
}

.experimental_connected_components <- function(active, dims3) {
  labels <- integer(length(active))
  active_idx <- which(active)
  component <- 0L
  plane <- dims3[1L] * dims3[2L]
  queue <- integer(length(active_idx))

  for (start in active_idx) {
    if (labels[start] != 0L) next
    component <- component + 1L
    queue[1L] <- start
    head <- 1L
    tail <- 1L
    labels[start] <- component

    while (head <= tail) {
      current <- queue[head]
      head <- head + 1L
      zero <- current - 1L
      x <- zero %% dims3[1L] + 1L
      y <- (zero %/% dims3[1L]) %% dims3[2L] + 1L
      z <- zero %/% plane + 1L
      neighbours <- c(
        if (x > 1L) current - 1L else NA_integer_,
        if (x < dims3[1L]) current + 1L else NA_integer_,
        if (y > 1L) current - dims3[1L] else NA_integer_,
        if (y < dims3[2L]) current + dims3[1L] else NA_integer_,
        if (z > 1L) current - plane else NA_integer_,
        if (z < dims3[3L]) current + plane else NA_integer_
      )
      neighbours <- neighbours[is.finite(neighbours)]
      add <- neighbours[active[neighbours] & labels[neighbours] == 0L]
      if (length(add)) {
        labels[add] <- component
        queue[tail + seq_along(add)] <- add
        tail <- tail + length(add)
      }
    }
  }
  labels
}

.experimental_regions <- function(regions, active, spec, merge_radius) {
  eligible <- logical(spec$n_full)
  eligible[spec$mask_idx] <- TRUE
  if (is.null(regions)) {
    seed <- which(active & eligible)
    expanded <- logical(spec$n_full)
    expanded[.experimental_dilate(seed, spec$dim, merge_radius)] <- TRUE
    return(.experimental_connected_components(expanded & eligible, spec$dim))
  }
  if (!(is.atomic(regions) &&
        length(regions) %in% c(length(spec$mask_idx), spec$n_full))) {
    stop("`regions` must have one value per GDS sample or full voxel.", call. = FALSE)
  }
  full <- rep.int(NA, spec$n_full)
  if (length(regions) == spec$n_full) {
    full[] <- regions
  } else {
    full[spec$mask_idx] <- regions
  }
  included <- !is.na(full) & full != 0 & eligible
  labels <- integer(spec$n_full)
  labels[included] <- as.integer(as.factor(full[included]))
  labels
}

.experimental_weighted_cor <- function(x, y, weight, min_pairs = 5L) {
  ok <- is.finite(x) & is.finite(y) & is.finite(weight) & weight > 0
  if (sum(ok) < min_pairs) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  weight <- weight[ok]
  weight <- weight / sum(weight)
  x <- x - sum(weight * x)
  y <- y - sum(weight * y)
  denominator <- sqrt(sum(weight * x^2) * sum(weight * y^2))
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    return(NA_real_)
  }
  sum(weight * x * y) / denominator
}

.experimental_loo_template <- function(beta, variance, subject, target_idx,
                                       eps = 1e-12) {
  others <- setdiff(seq_len(ncol(beta)), subject)
  values <- beta[target_idx, others, drop = FALSE]
  variances <- variance[target_idx, others, drop = FALSE]
  ok <- is.finite(values) & is.finite(variances) & variances > 0
  weight <- matrix(0, nrow(values), ncol(values))
  weight[ok] <- 1 / pmax(variances[ok], eps)
  values[!ok] <- 0
  sum_weight <- rowSums(weight)
  template <- rowSums(weight * values) / sum_weight
  template_variance <- 1 / sum_weight
  template[!is.finite(template) | sum_weight <= 0] <- NA_real_
  template_variance[!is.finite(template_variance) | sum_weight <= 0] <- NA_real_
  list(beta = template, var = template_variance)
}

.experimental_shift_values <- function(x, subject, target_idx, offset, dims3) {
  source <- .experimental_shift_source(target_idx, offset, dims3)
  out <- rep.int(NA_real_, length(source))
  valid <- is.finite(source)
  out[valid] <- x[source[valid], subject]
  out
}

.experimental_subject_shift <- function(beta, variance, subject, patch_idx,
                                        dims3, offsets, delta,
                                        translation_penalty, mass_retention,
                                        min_correlation_gain,
                                        min_aligned_correlation) {
  template <- .experimental_loo_template(beta, variance, subject, patch_idx)
  correlation <- rep.int(NA_real_, nrow(offsets))
  mass <- rep.int(NA_real_, nrow(offsets))
  distance2 <- rowSums(offsets^2)
  min_pairs <- min(5L, length(patch_idx))

  for (j in seq_len(nrow(offsets))) {
    values <- .experimental_shift_values(
      beta, subject, patch_idx, offsets[j, ], dims3
    )
    value_variance <- .experimental_shift_values(
      variance, subject, patch_idx, offsets[j, ], dims3
    )
    weight <- 1 / (value_variance + template$var)
    correlation[j] <- .experimental_weighted_cor(
      values, template$beta, weight, min_pairs = min_pairs
    )
    mass[j] <- sum(pmax(abs(values) - delta, 0), na.rm = TRUE)
  }

  zero <- which(distance2 == 0)[1L]
  raw_correlation <- correlation[zero]
  score <- correlation - translation_penalty * distance2
  score[!is.finite(score)] <- -Inf
  best <- if (any(is.finite(correlation))) which.max(score) else zero
  aligned_correlation <- correlation[best]
  retention <- mass[best] / mass[zero]
  valid <- is.finite(raw_correlation) & is.finite(aligned_correlation) &
    is.finite(retention) & mass[zero] > sqrt(.Machine$double.eps)
  nonzero <- distance2[best] > 0
  accepted <- valid & nonzero &
    aligned_correlation - raw_correlation > min_correlation_gain &
    aligned_correlation >= min_aligned_correlation &
    retention >= mass_retention[1L] & retention <= mass_retention[2L]

  rescue <- if (valid) 0 else NA_real_
  selected <- c(0L, 0L, 0L)
  if (accepted) {
    rescue <- (aligned_correlation - raw_correlation) /
      pmax(1 - raw_correlation, sqrt(.Machine$double.eps))
    rescue <- pmin(1, pmax(0, rescue))
    selected <- as.integer(offsets[best, ])
  } else if (valid) {
    aligned_correlation <- raw_correlation
    retention <- 1
  }

  list(
    rescue = rescue,
    valid = valid,
    shifted = accepted,
    shift = selected,
    raw_correlation = raw_correlation,
    aligned_correlation = aligned_correlation,
    retention = retention
  )
}

.experimental_region_rescue <- function(beta, variance, patch_idx, dims3,
                                        offsets, delta, min_subjects,
                                        translation_penalty, mass_retention,
                                        min_correlation_gain,
                                        min_aligned_correlation) {
  if (length(patch_idx) < 3L) {
    return(list(
      rescue = NA_real_, n_valid = 0L, n_shifted = 0L,
      shift_prevalence = NA_real_, mean_raw_correlation = NA_real_,
      mean_aligned_correlation = NA_real_, shift_counts = integer()
    ))
  }
  fits <- lapply(seq_len(ncol(beta)), function(subject) {
    .experimental_subject_shift(
      beta, variance, subject, patch_idx, dims3, offsets, delta,
      translation_penalty, mass_retention, min_correlation_gain,
      min_aligned_correlation
    )
  })
  valid <- vapply(fits, function(x) isTRUE(x$valid), logical(1L))
  n_valid <- sum(valid)
  if (n_valid < min_subjects) {
    return(list(
      rescue = NA_real_, n_valid = n_valid, n_shifted = 0L,
      shift_prevalence = NA_real_, mean_raw_correlation = NA_real_,
      mean_aligned_correlation = NA_real_, shift_counts = integer()
    ))
  }

  rescue <- vapply(fits[valid], `[[`, numeric(1L), "rescue")
  shifted <- vapply(fits[valid], `[[`, logical(1L), "shifted")
  raw_correlation <- vapply(
    fits[valid], `[[`, numeric(1L), "raw_correlation"
  )
  aligned_correlation <- vapply(
    fits[valid], `[[`, numeric(1L), "aligned_correlation"
  )
  shifts <- do.call(rbind, lapply(fits[valid], `[[`, "shift"))
  shift_key <- apply(shifts, 1L, paste, collapse = ",")

  list(
    rescue = mean(rescue),
    n_valid = n_valid,
    n_shifted = sum(shifted),
    shift_prevalence = mean(shifted),
    mean_raw_correlation = mean(raw_correlation),
    mean_aligned_correlation = mean(aligned_correlation),
    shift_counts = table(shift_key)
  )
}

#' Experimental sign-cancellation and spatial-dispersion maps
#'
#' Estimates two complementary diagnostics from one set of subject-level maps.
#' `cancellation_probability` is the plug-in posterior probability, under a
#' normal random-effects model, that meaningful positive and negative population
#' prevalences both exceed `prevalence` while the population mean lies inside
#' `[-equivalence, equivalence]`. `shift_rescue_fraction` asks a different,
#' region-level question: what fraction of each subject's mismatch to an
#' independently formed leave-one-subject-out template is removed by a small
#' local translation? The region summary is projected back to its voxels.
#'
#' @section Experimental:
#' This interface and its statistical model may change. Neither output is
#' corrected whole-brain inference. Spatial rescue describes between-subject
#' localisation variability; it does not test movement within a subject. A
#' subject is excluded from its own template, but its shift and rescue are
#' estimated from the same subject-level map. Treat the result as a descriptive
#' diagnostic unless it is confirmed in independent data.
#'
#' @section Automatic regions:
#' When `regions` is omitted, a subject supports a voxel when the
#' measurement-model probability of an effect outside `[-delta, delta]`
#' somewhere within `shift_radius` reaches the search-adjusted threshold implied
#' by `support_probability`. Voxels enter the spatial diagnostic when at least
#' `prevalence` of subjects support them; nearby support is joined into
#' six-connected local regions. This rule is independent of
#' `cancellation_probability`, so a same-signed activation that varies in
#' location can be evaluated.
#'
#' @param x A realised subject-level [`gds`] with aligned `beta` and genuine
#'   `var` or `se` assays in voxel space. One map per subject is sufficient.
#' @param delta Positive smallest meaningful subject-level effect, in beta units.
#' @param equivalence Positive practical-equivalence bound for the population
#'   mean. Defaults to `delta`.
#' @param prevalence Required positive and negative population prevalence for
#'   `cancellation_probability`; also the minimum fraction of subjects supporting
#'   an automatic local region. Must be strictly between 0 and 0.5.
#' @param support_probability Nominal subject-level measurement-model
#'   probability of an effect outside `[-delta, delta]` used for automatic
#'   region support. The per-location threshold is Bonferroni-adjusted over the
#'   searched offsets. Defaults to 0.95.
#' @param shift_radius Maximum integer translation in voxels along each axis.
#' @param patch_radius Integer dilation around each local region used to estimate
#'   subject-to-template pattern correlation.
#' @param regions Optional region labels with one value per active sample or full
#'   voxel. Zero and `NA` are excluded. When omitted, automatic local regions are
#'   formed from meaningful-effect support.
#' @param min_subjects Minimum number of valid subjects for either diagnostic.
#' @param translation_penalty Nonnegative penalty per squared voxel shifted.
#' @param mass_retention Length-two accepted interval for retained meaningful
#'   local signal mass after translation.
#' @param min_correlation_gain Minimum correlation improvement required to accept
#'   a nonzero subject displacement.
#' @param min_aligned_correlation Minimum signed correlation with the
#'   leave-one-subject-out template after alignment. This prevents an inverted
#'   pattern from appearing rescued merely because it was shifted away.
#'
#' @return A group-level [`gds`] with exactly two assays:
#'   `cancellation_probability` and `shift_rescue_fraction`. Spatial rescue is
#'   `NA` where no supported local region or too few valid subject patterns are
#'   available; an evaluated region with no accepted displacement returns zero.
#' @export
experimental_cancellation <- function(x,
                                      delta,
                                      equivalence = delta,
                                      prevalence = 0.2,
                                      support_probability = 0.95,
                                      shift_radius = 1L,
                                      patch_radius = 2L,
                                      regions = NULL,
                                      min_subjects = 8L,
                                      translation_penalty = 0.01,
                                      mass_retention = c(0.75, 1.25),
                                      min_correlation_gain = 0.05,
                                      min_aligned_correlation = 0.5) {
  .experimental_validate_input(x)
  delta <- .experimental_scalar(delta, "delta", lower = 0, lower_open = TRUE)
  equivalence <- .experimental_scalar(
    equivalence, "equivalence", lower = 0, lower_open = TRUE
  )
  prevalence <- .experimental_scalar(
    prevalence, "prevalence", lower = 0, upper = 0.5,
    lower_open = TRUE, upper_open = TRUE
  )
  support_probability <- .experimental_scalar(
    support_probability, "support_probability", lower = 0.5, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  translation_penalty <- .experimental_scalar(
    translation_penalty, "translation_penalty", lower = 0
  )
  min_correlation_gain <- .experimental_scalar(
    min_correlation_gain, "min_correlation_gain", lower = 0
  )
  min_aligned_correlation <- .experimental_scalar(
    min_aligned_correlation, "min_aligned_correlation", lower = -1, upper = 1
  )
  shift_radius <- as.integer(shift_radius)
  patch_radius <- as.integer(patch_radius)
  min_subjects <- as.integer(min_subjects)
  if (length(shift_radius) != 1L || is.na(shift_radius) || shift_radius < 0L) {
    stop("`shift_radius` must be one nonnegative integer.", call. = FALSE)
  }
  if (length(patch_radius) != 1L || is.na(patch_radius) || patch_radius < 0L) {
    stop("`patch_radius` must be one nonnegative integer.", call. = FALSE)
  }
  if (length(min_subjects) != 1L || is.na(min_subjects) || min_subjects < 3L) {
    stop("`min_subjects` must be one integer of at least 3.", call. = FALSE)
  }
  if (!is.numeric(mass_retention) || length(mass_retention) != 2L ||
      any(!is.finite(mass_retention)) || mass_retention[1L] <= 0 ||
      mass_retention[2L] < mass_retention[1L]) {
    stop("`mass_retention` must be an increasing positive length-two interval.", call. = FALSE)
  }

  beta <- assay(x, "beta")
  variance <- .experimental_variance(x)
  dims <- dim(beta)
  probability <- array(NA_real_, c(dims[1L], 1L, dims[3L]))
  rescue <- array(NA_real_, dim(probability))
  spec <- .experimental_voxel_spec(x)
  offsets <- .experimental_offset_grid(shift_radius, spec$dim)
  region_receipts <- vector("list", dims[3L])

  for (k in seq_len(dims[3L])) {
    beta_k <- matrix(beta[, , k], nrow = dims[1L], ncol = dims[2L])
    variance_k <- matrix(
      variance[, , k], nrow = dims[1L], ncol = dims[2L]
    )
    probability[, 1L, k] <- .experimental_probability_kernel(
      beta = t(beta_k),
      variance = t(variance_k),
      delta = delta,
      equivalence = equivalence,
      prevalence = prevalence,
      min_subjects = min_subjects
    )

    meaningful_probability <- .experimental_meaningful_probability(
      beta_k, variance_k, delta
    )
    support <- .experimental_local_support(
      meaningful_probability, spec, shift_radius, support_probability,
      min_subjects
    )
    support_full <- rep.int(NA_real_, spec$n_full)
    support_full[spec$mask_idx] <- support
    active <- is.finite(support_full) & support_full >= prevalence
    labels <- .experimental_regions(regions, active, spec, shift_radius)
    beta_full <- .experimental_unpack_matrix(beta, k, spec)
    variance_full <- .experimental_unpack_matrix(variance, k, spec)
    rescue_full <- rep.int(NA_real_, spec$n_full)
    contrast_receipt <- list()

    for (label in sort(unique(labels[labels > 0L]))) {
      region_idx <- which(labels == label)
      if (!any(active[region_idx])) next
      patch_idx <- .experimental_dilate(region_idx, spec$dim, patch_radius)
      patch_idx <- intersect(patch_idx, spec$mask_idx)
      region_fit <- .experimental_region_rescue(
        beta = beta_full,
        variance = variance_full,
        patch_idx = patch_idx,
        dims3 = spec$dim,
        offsets = offsets,
        delta = delta,
        min_subjects = min_subjects,
        translation_penalty = translation_penalty,
        mass_retention = mass_retention,
        min_correlation_gain = min_correlation_gain,
        min_aligned_correlation = min_aligned_correlation
      )
      rescue_full[region_idx] <- region_fit$rescue
      contrast_receipt[[as.character(label)]] <- c(
        list(n_voxel = length(region_idx), n_patch = length(patch_idx)),
        region_fit
      )
    }
    rescue[, 1L, k] <- .experimental_pack_vector(rescue_full, spec)
    region_receipts[[k]] <- contrast_receipt
  }

  dimnames(probability) <- list(NULL, "meta", contrasts(x))
  dimnames(rescue) <- dimnames(probability)
  metadata_out <- utils::modifyList(
    metadata(x),
    list(
      experimental = list(
        method = "cancellation_spatial_dispersion_v2",
        status = "experimental",
        cancellation_estimand = "opposing_prevalence_with_equivalent_mean",
        shift_estimand = "mean_fraction_of_loo_pattern_mismatch_removed",
        shift_scope = "between_subject",
        validation = "leave_one_subject_out_template_descriptive",
        delta = delta,
        equivalence = equivalence,
        prevalence = prevalence,
        support_probability = support_probability,
        shift_radius = shift_radius,
        patch_radius = patch_radius,
        min_subjects = min_subjects,
        region_receipts = region_receipts
      )
    )
  )
  metadata_out <- add_provenance_node(
    metadata_out,
    "experimental_cancellation",
    list(
      delta = delta,
      equivalence = equivalence,
      prevalence = prevalence,
      support_probability = support_probability,
      shift_radius = shift_radius,
      patch_radius = patch_radius
    )
  )

  new_gds(
    assays = list(
      cancellation_probability = probability,
      shift_rescue_fraction = rescue
    ),
    space = space(x),
    subjects = "meta",
    contrasts = contrasts(x),
    row_data = row_data(x),
    metadata = metadata_out
  )
}
