# Experimental cancellation and spatial-rescue diagnostics -----------------

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

.experimental_validate_splits <- function(split_a, split_b) {
  if (!inherits(split_a, "gds") || !inherits(split_b, "gds")) {
    stop("`split_a` and `split_b` must be realised GDS objects.", call. = FALSE)
  }
  if (!inherits(space(split_a), "space_voxel") ||
      !inherits(space(split_b), "space_voxel")) {
    stop("Experimental spatial rescue requires voxel-space GDS inputs.", call. = FALSE)
  }
  assert_compatible_spaces(split_a, split_b)
  if (!isTRUE(all.equal(space(split_a)$affine, space(split_b)$affine,
                        tolerance = 1e-8))) {
    stop("Split voxel affines must match.", call. = FALSE)
  }
  if (!identical(space(split_a)$mask_idx, space(split_b)$mask_idx)) {
    stop("Split voxel masks and sample order must match.", call. = FALSE)
  }
  if (!identical(subjects(split_a), subjects(split_b))) {
    stop("Split subjects and their order must match.", call. = FALSE)
  }
  if (!identical(contrasts(split_a), contrasts(split_b))) {
    stop("Split contrasts and their order must match.", call. = FALSE)
  }
  beta_a <- assay(split_a, "beta")
  beta_b <- assay(split_b, "beta")
  if (is.null(beta_a) || is.null(beta_b)) {
    stop("Experimental cancellation requires a `beta` assay in both splits.", call. = FALSE)
  }
  if (!identical(dim(beta_a), dim(beta_b))) {
    stop("Split beta assays must have identical dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

.experimental_combine_splits <- function(beta_a, var_a, beta_b, var_b,
                                         eps = 1e-12) {
  ok_a <- is.finite(beta_a) & is.finite(var_a) & var_a > 0
  ok_b <- is.finite(beta_b) & is.finite(var_b) & var_b > 0
  both <- ok_a & ok_b
  wa <- 1 / pmax(var_a, eps)
  wb <- 1 / pmax(var_b, eps)
  sw <- wa + wb
  beta <- (wa * beta_a + wb * beta_b) / sw
  variance <- 1 / sw
  beta[!both] <- NA_real_
  variance[!both] <- NA_real_
  list(beta = beta, var = variance)
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
  # the practical-equivalence interval for the population mean, then
  # integrate the plug-in normal posterior for mu.
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
  out[spec$mask_idx, ] <- x[, , contrast, drop = TRUE]
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

.experimental_signed_cor <- function(x, y, min_pairs = 5L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_pairs) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  x <- x - mean(x)
  y <- y - mean(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    return(NA_real_)
  }
  sum(x * y) / denominator
}

.experimental_estimate_shifts <- function(train, patch_idx, dims3, offsets,
                                          translation_penalty) {
  n_subject <- ncol(train)
  shifts <- matrix(0L, nrow = n_subject, ncol = 3L)
  source_maps <- lapply(seq_len(nrow(offsets)), function(j) {
    .experimental_shift_source(patch_idx, offsets[j, ], dims3)
  })

  for (i in seq_len(n_subject)) {
    others <- setdiff(seq_len(n_subject), i)
    template <- rowMeans(train[patch_idx, others, drop = FALSE], na.rm = TRUE)
    template[!is.finite(template)] <- NA_real_
    score <- rep(-Inf, nrow(offsets))
    for (j in seq_len(nrow(offsets))) {
      source <- source_maps[[j]]
      values <- rep.int(NA_real_, length(source))
      valid <- is.finite(source)
      values[valid] <- train[source[valid], i]
      correlation <- .experimental_signed_cor(values, template)
      if (is.finite(correlation)) {
        score[j] <- correlation - translation_penalty * sum(offsets[j, ]^2)
      }
    }
    if (any(is.finite(score))) shifts[i, ] <- offsets[which.max(score), ]
  }
  shifts
}

.experimental_apply_shifts <- function(test, target_idx, shifts, dims3) {
  out <- matrix(NA_real_, nrow = length(target_idx), ncol = ncol(test))
  for (i in seq_len(ncol(test))) {
    source <- .experimental_shift_source(target_idx, shifts[i, ], dims3)
    valid <- is.finite(source)
    out[valid, i] <- test[source[valid], i]
  }
  out
}

.experimental_gap <- function(x, delta, min_subjects) {
  meaningful <- sign(x) * pmax(abs(x) - delta, 0)
  meaningful[!is.finite(x)] <- NA_real_
  n <- rowSums(is.finite(meaningful))
  unsigned <- rowMeans(abs(meaningful), na.rm = TRUE)
  signed <- abs(rowMeans(meaningful, na.rm = TRUE))
  gap <- unsigned - signed
  gap[n < min_subjects | !is.finite(gap)] <- NA_real_
  pmax(gap, 0)
}

.experimental_patch_gain <- function(train, test, aligned_patch, patch_idx,
                                     shifts, dims3) {
  gain <- rep(NA_real_, ncol(test))
  aligned_correlation <- rep(NA_real_, ncol(test))
  for (i in seq_len(ncol(test))) {
    others <- setdiff(seq_len(ncol(test)), i)
    template <- rowMeans(train[patch_idx, others, drop = FALSE], na.rm = TRUE)
    raw <- test[patch_idx, i]
    shifted <- aligned_patch[, i]
    raw_cor <- .experimental_signed_cor(raw, template)
    shifted_cor <- .experimental_signed_cor(shifted, template)
    if (is.finite(raw_cor) && is.finite(shifted_cor)) {
      gain[i] <- shifted_cor - raw_cor
      aligned_correlation[i] <- shifted_cor
    }
  }
  list(
    gain = mean(gain, na.rm = TRUE),
    aligned_correlation = mean(aligned_correlation, na.rm = TRUE)
  )
}

.experimental_fold_rescue <- function(train, test, region_idx, eval_idx,
                                      patch_idx, dims3, offsets, delta,
                                      min_subjects, translation_penalty,
                                      mass_retention, min_correlation_gain,
                                      min_aligned_correlation) {
  shifts <- .experimental_estimate_shifts(
    train, patch_idx, dims3, offsets, translation_penalty
  )
  aligned_eval <- .experimental_apply_shifts(test, eval_idx, shifts, dims3)
  raw_gap <- .experimental_gap(test[eval_idx, , drop = FALSE], delta, min_subjects)
  aligned_gap <- .experimental_gap(aligned_eval, delta, min_subjects)
  rescue <- 1 - aligned_gap / raw_gap
  rescue[!is.finite(rescue) | raw_gap <= sqrt(.Machine$double.eps)] <- NA_real_
  rescue <- pmin(1, pmax(0, rescue))

  aligned_patch <- .experimental_apply_shifts(test, patch_idx, shifts, dims3)
  raw_mass <- sum(rowMeans(abs(test[patch_idx, , drop = FALSE]), na.rm = TRUE),
                  na.rm = TRUE)
  aligned_mass <- sum(rowMeans(abs(aligned_patch), na.rm = TRUE), na.rm = TRUE)
  retention <- aligned_mass / raw_mass
  correlation <- .experimental_patch_gain(
    train, test, aligned_patch, patch_idx, shifts, dims3
  )
  valid_region <- is.finite(retention) & retention >= mass_retention[1L] &
    retention <= mass_retention[2L] & is.finite(correlation$gain) &
    correlation$gain > min_correlation_gain &
    is.finite(correlation$aligned_correlation) &
    correlation$aligned_correlation >= min_aligned_correlation
  if (!valid_region) rescue[] <- NA_real_

  list(
    rescue = rescue,
    valid = valid_region,
    retention = retention,
    correlation_gain = correlation$gain,
    aligned_correlation = correlation$aligned_correlation,
    shifts = shifts
  )
}

.experimental_regions <- function(regions, probability_full, candidate_threshold,
                                  spec) {
  candidate <- is.finite(probability_full) & probability_full >= candidate_threshold
  if (is.null(regions)) {
    return(.experimental_connected_components(candidate, spec$dim))
  }
  if (!(is.atomic(regions) && length(regions) %in% c(length(spec$mask_idx), spec$n_full))) {
    stop("`regions` must have one value per GDS sample or full voxel.", call. = FALSE)
  }
  full <- rep.int(NA, spec$n_full)
  if (length(regions) == spec$n_full) {
    full[] <- regions
  } else {
    full[spec$mask_idx] <- regions
  }
  labels <- as.integer(as.factor(full))
  labels[is.na(full) | full == 0] <- 0L
  labels
}

#' Experimental cancellation and fast spatial-rescue maps
#'
#' Estimates two deliberately narrow voxelwise diagnostics from independent
#' subject-level split maps. `cancellation_probability` is the plug-in posterior
#' probability, under a normal random-effects model, that meaningful positive
#' and negative population prevalences both exceed `prevalence` while the mean
#' lies inside `[-equivalence, equivalence]`. `shift_rescue_fraction` estimates
#' the fraction of thresholded cancellation removed by a small local translation
#' learned on one split and evaluated on the other, averaged over both split
#' directions.
#'
#' @section Experimental:
#' This interface and its statistical model may change. Spatial rescue is a
#' diagnostic, not corrected whole-brain inference. It is returned only where
#' cancellation exceeds `candidate_threshold` and both held-out folds preserve
#' local unsigned mass and improve signed patch correlation.
#'
#' @param split_a,split_b Realised subject-level [`gds`] objects from independent
#'   runs or balanced data splits. Both require aligned `beta` and genuine `var`
#'   or `se` assays, identical voxel space, subjects, and contrasts.
#' @param delta Positive smallest meaningful subject-level effect, in beta units.
#' @param equivalence Positive practical-equivalence bound for the population
#'   mean. Defaults to `delta`.
#' @param prevalence Required positive and negative population prevalence;
#'   strictly between 0 and 0.5.
#' @param shift_radius Maximum integer translation in voxels along each axis.
#' @param patch_radius Integer dilation around each candidate region used to
#'   estimate signed patch correlation.
#' @param candidate_threshold Minimum cancellation probability at which spatial
#'   rescue is evaluated.
#' @param regions Optional region labels with one value per active sample or full
#'   voxel. Zero and `NA` are excluded. When omitted, six-connected candidate
#'   components are used.
#' @param min_subjects Minimum number of valid subjects at a voxel.
#' @param translation_penalty Nonnegative penalty per squared voxel shifted.
#' @param mass_retention Length-two accepted interval for held-out unsigned patch
#'   mass after translation.
#' @param min_correlation_gain Required held-out signed patch-correlation gain.
#' @param min_aligned_correlation Minimum mean held-out signed correlation after
#'   alignment. This guards against apparent rescue caused by shifting an
#'   inverted pattern away from a strongly negative match.
#'
#' @return A group-level [`gds`] with exactly two assays:
#'   `cancellation_probability` and `shift_rescue_fraction`.
#' @export
experimental_cancellation <- function(split_a,
                                      split_b,
                                      delta,
                                      equivalence = delta,
                                      prevalence = 0.2,
                                      shift_radius = 1L,
                                      patch_radius = 2L,
                                      candidate_threshold = 0.5,
                                      regions = NULL,
                                      min_subjects = 8L,
                                      translation_penalty = 0.01,
                                      mass_retention = c(0.75, 1.25),
                                      min_correlation_gain = 0,
                                      min_aligned_correlation = 0.5) {
  .experimental_validate_splits(split_a, split_b)
  delta <- .experimental_scalar(delta, "delta", lower = 0, lower_open = TRUE)
  equivalence <- .experimental_scalar(
    equivalence, "equivalence", lower = 0, lower_open = TRUE
  )
  prevalence <- .experimental_scalar(
    prevalence, "prevalence", lower = 0, upper = 0.5,
    lower_open = TRUE, upper_open = TRUE
  )
  candidate_threshold <- .experimental_scalar(
    candidate_threshold, "candidate_threshold", lower = 0, upper = 1
  )
  translation_penalty <- .experimental_scalar(
    translation_penalty, "translation_penalty", lower = 0
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
  if (length(min_correlation_gain) != 1L || !is.numeric(min_correlation_gain) ||
      !is.finite(min_correlation_gain)) {
    stop("`min_correlation_gain` must be one finite numeric value.", call. = FALSE)
  }
  min_aligned_correlation <- .experimental_scalar(
    min_aligned_correlation, "min_aligned_correlation", lower = -1, upper = 1
  )

  beta_a <- assay(split_a, "beta")
  beta_b <- assay(split_b, "beta")
  var_a <- .experimental_variance(split_a)
  var_b <- .experimental_variance(split_b)
  if (!identical(dim(beta_a), dim(var_a)) || !identical(dim(beta_b), dim(var_b))) {
    stop("Beta and variance dimensions must match within each split.", call. = FALSE)
  }
  combined <- .experimental_combine_splits(beta_a, var_a, beta_b, var_b)
  dims <- dim(beta_a)
  probability <- array(NA_real_, c(dims[1L], 1L, dims[3L]))
  rescue <- array(NA_real_, dim(probability))
  spec <- .experimental_voxel_spec(split_a)
  offsets <- .experimental_offset_grid(shift_radius, spec$dim)
  shift_receipts <- vector("list", dims[3L])

  for (k in seq_len(dims[3L])) {
    probability[, 1L, k] <- .experimental_probability_kernel(
      beta = t(combined$beta[, , k, drop = TRUE]),
      variance = t(combined$var[, , k, drop = TRUE]),
      delta = delta,
      equivalence = equivalence,
      prevalence = prevalence,
      min_subjects = min_subjects
    )
    probability_full <- rep.int(NA_real_, spec$n_full)
    probability_full[spec$mask_idx] <- probability[, 1L, k]
    labels <- .experimental_regions(
      regions, probability_full, candidate_threshold, spec
    )
    candidate <- is.finite(probability_full) &
      probability_full >= candidate_threshold & labels > 0L
    beta_a_full <- .experimental_unpack_matrix(beta_a, k, spec)
    beta_b_full <- .experimental_unpack_matrix(beta_b, k, spec)
    rescue_full <- rep.int(NA_real_, spec$n_full)
    contrast_receipt <- list()

    for (label in sort(unique(labels[labels > 0L]))) {
      region_idx <- which(labels == label)
      eval_idx <- region_idx[candidate[region_idx]]
      if (!length(eval_idx)) next
      patch_idx <- .experimental_dilate(region_idx, spec$dim, patch_radius)
      fold_ab <- .experimental_fold_rescue(
        train = beta_a_full,
        test = beta_b_full,
        region_idx = region_idx,
        eval_idx = eval_idx,
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
      fold_ba <- .experimental_fold_rescue(
        train = beta_b_full,
        test = beta_a_full,
        region_idx = region_idx,
        eval_idx = eval_idx,
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
      both_valid <- is.finite(fold_ab$rescue) & is.finite(fold_ba$rescue)
      rescued <- rep.int(NA_real_, length(eval_idx))
      rescued[both_valid] <-
        (fold_ab$rescue[both_valid] + fold_ba$rescue[both_valid]) / 2
      rescue_full[eval_idx] <- rescued
      contrast_receipt[[as.character(label)]] <- list(
        n_voxel = length(region_idx),
        n_candidate = length(eval_idx),
        fold_ab_valid = fold_ab$valid,
        fold_ba_valid = fold_ba$valid,
        fold_ab_mass_retention = fold_ab$retention,
        fold_ba_mass_retention = fold_ba$retention,
        fold_ab_correlation_gain = fold_ab$correlation_gain,
        fold_ba_correlation_gain = fold_ba$correlation_gain,
        fold_ab_aligned_correlation = fold_ab$aligned_correlation,
        fold_ba_aligned_correlation = fold_ba$aligned_correlation
      )
    }
    rescue[, 1L, k] <- .experimental_pack_vector(rescue_full, spec)
    shift_receipts[[k]] <- contrast_receipt
  }

  dimnames(probability) <- list(NULL, "meta", contrasts(split_a))
  dimnames(rescue) <- dimnames(probability)
  metadata_out <- utils::modifyList(
    gds_metadata(),
    list(
      experimental = list(
        method = "cancellation_shift_rescue_v1",
        status = "experimental",
        delta = delta,
        equivalence = equivalence,
        prevalence = prevalence,
        shift_radius = shift_radius,
        patch_radius = patch_radius,
        candidate_threshold = candidate_threshold,
        min_subjects = min_subjects,
        cross_fit = "A_to_B_and_B_to_A",
        shift_receipts = shift_receipts
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
      shift_radius = shift_radius,
      patch_radius = patch_radius,
      candidate_threshold = candidate_threshold
    )
  )

  new_gds(
    assays = list(
      cancellation_probability = probability,
      shift_rescue_fraction = rescue
    ),
    space = space(split_a),
    subjects = "meta",
    contrasts = contrasts(split_a),
    row_data = row_data(split_a),
    metadata = metadata_out
  )
}
