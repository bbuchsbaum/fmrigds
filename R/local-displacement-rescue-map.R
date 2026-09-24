# Local displacement rescue: corrected search-map interface --------------

.ldr_matrix_contrast <- function(x, assay_name, contrast_index) {
  value <- assay(x, assay_name)[, , contrast_index, drop = FALSE]
  matrix(value, nrow = dim(value)[1L], ncol = dim(value)[2L])
}

.ldr_best_search_units <- function(units, scan, calibration = NULL) {
  center <- vapply(units, `[[`, integer(1L), "center")
  groups <- split(seq_along(units), center)
  vapply(groups, function(indices) {
    if (is.null(calibration)) {
      score <- pmin(
        scan$activation_t[indices], scan$displacement_t[indices]
      )
      score[!is.finite(score)] <- -Inf
      return(indices[which.max(score)])
    }
    p_value <- calibration$ldr_fwer_p[indices]
    score <- pmin(
      scan$activation_t[indices], scan$displacement_t[indices]
    )
    order_value <- order(
      ifelse(is.na(p_value), Inf, p_value),
      -ifelse(is.finite(score), score, -Inf),
      indices
    )
    indices[order_value[1L]]
  }, integer(1L))
}

.ldr_exact_search_refit <- function(beta_a, beta_b, var_a, var_b, unit,
                                    subject_ids, folds, iterations, ridge,
                                    tolerance) {
  geometry <- unit$geometry
  patch <- .ldr_prepare_patch(
    beta_a[geometry$sample_idx, , drop = FALSE],
    beta_b[geometry$sample_idx, , drop = FALSE],
    var_a[geometry$sample_idx, , drop = FALSE],
    var_b[geometry$sample_idx, , drop = FALSE],
    geometry
  )
  .ldr_crossfit(
    patch,
    subject_ids = subject_ids,
    tau_grid_mm = .ldr_default_tau_grid(unit$shift_radius_mm),
    folds = folds,
    iterations = iterations,
    ridge = ridge,
    tolerance = tolerance
  )
}

#' Whole-search Local Displacement Rescue map
#'
#' Screens every eligible center in a prespecified voxel search with a
#' cross-fitted spatial-derivative tangent model, corrects activation,
#' displacement, and ordinary fixed-coordinate evidence over the complete
#' center-by-scale family, and passes screen discoveries through the exact
#' discrete-shift LDR reference model.
#'
#' @section Scientific contract:
#' Local noise correlation is estimated from standardized paired-split
#' differences in training subjects and used only for held-out generalized
#' least-squares scoring. Activation and ordinary-MNI maxima are calibrated by
#' rerunning the complete tangent search under paired whole-subject sign flips.
#' Displacement maxima use a subject-level Rademacher wild bootstrap of centered
#' held-out cross-split tangent products, with one multiplier shared across the
#' search for each subject. `ldr_fwer_p` is the maximum of corrected activation
#' and displacement p-values.
#'
#' A nonzero `ldr_flag` additionally requires a nonsignificant corrected
#' ordinary-MNI test, positive MNI loss, and exact soft-shift improvement over
#' both the aligned-positive and aligned signed-amplitude competitors. Exact
#' refitting can remove a tangent-screen discovery but cannot create one.
#'
#' @section Experimental:
#' This is a translation-only voxel-space procedure requiring independent
#' paired split maps. With `centers = NULL`, the multiplicity family includes
#' every active center having a complete patch at each valid requested scale.
#' Supplying `centers` produces search-set, not whole-brain, correction. With
#' `n_resamples = 0`, corrected p-values are missing, no exact refits are run,
#' and `ldr_flag` is zero. Use at least 999 resamples for confirmatory work.
#' The tangent screen uses split-difference spatial covariance; the exact
#' confirmation currently retains the v0.1 diagonal likelihood.
#'
#' @param split_a,split_b Realised independent subject-level split maps with
#'   matching voxel space, subjects, contrasts, `beta`, and genuine `var` or
#'   `se` assays.
#' @param centers Optional full-volume voxel indices defining the search. `NULL`
#'   searches every active center with a complete patch.
#' @param contrast Optional contrast name or index.
#' @param patch_radius_mm One or more positive physical patch radii.
#' @param shift_radius_mm One or more positive maximum translation radii. Every
#'   searched configuration must satisfy `patch >= 2 * shift`.
#' @param folds Number of deterministic subject cross-fitting folds.
#' @param min_subjects Minimum number of paired subjects.
#' @param covariance_shrinkage Fixed shrinkage of split-difference local noise
#'   correlation toward identity, between zero and one.
#' @param iterations,ridge,tolerance Exact-refit optimization controls.
#' @param n_resamples Number of whole-search calibration resamples. Zero returns
#'   uncalibrated tangent diagnostics only.
#' @param seed Integer calibration seed. The caller's random state is restored.
#' @param alpha Familywise level used for candidate admission and `ldr_flag`.
#' @param min_prevalence Minimum exact-model activation prevalence required for
#'   a rescue flag. This separates a broadly shared displaced feature from a
#'   small active subgroup.
#' @param max_refits Maximum exact confirmations. `Inf` confirms every admitted
#'   candidate; a finite cap is conservative and is recorded as truncation.
#'
#' @return A group-level [`gds`] containing the corrected LDR map and component
#'   maps. Inferential values are `NA` outside analyzable centers; `ldr_flag` is
#'   zero wherever the complete rescue conjunction is not established.
#' @export
local_displacement_rescue_map <- function(split_a,
                                          split_b,
                                          centers = NULL,
                                          contrast = NULL,
                                          patch_radius_mm = 12,
                                          shift_radius_mm = 4,
                                          folds = 5L,
                                          min_subjects = 12L,
                                          covariance_shrinkage = 0.25,
                                          iterations = 20L,
                                          ridge = 1e-4,
                                          tolerance = 1e-5,
                                          n_resamples = 0L,
                                          seed = 1L,
                                          alpha = 0.05,
                                          min_prevalence = 0.7,
                                          max_refits = Inf) {
  .ldr_validate_splits(split_a, split_b)
  folds <- as.integer(folds)
  min_subjects <- as.integer(min_subjects)
  iterations <- as.integer(iterations)
  n_resamples <- as.integer(n_resamples)
  seed <- as.integer(seed)
  if (length(folds) != 1L || is.na(folds) || folds < 2L) {
    stop("`folds` must be one integer of at least 2.", call. = FALSE)
  }
  if (length(min_subjects) != 1L || is.na(min_subjects) ||
      min_subjects < 6L) {
    stop("`min_subjects` must be one integer of at least 6.",
         call. = FALSE)
  }
  if (length(iterations) != 1L || is.na(iterations) || iterations < 1L) {
    stop("`iterations` must be one positive integer.", call. = FALSE)
  }
  if (length(n_resamples) != 1L || is.na(n_resamples) ||
      n_resamples < 0L) {
    stop("`n_resamples` must be one nonnegative integer.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }
  alpha <- .ldr_scalar(
    alpha, "alpha", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  min_prevalence <- .ldr_scalar(
    min_prevalence, "min_prevalence", lower = 0, upper = 1,
    lower_open = TRUE
  )
  covariance_shrinkage <- .ldr_scalar(
    covariance_shrinkage, "covariance_shrinkage", lower = 0, upper = 1
  )
  ridge <- .ldr_scalar(ridge, "ridge", lower = 0, lower_open = TRUE)
  tolerance <- .ldr_scalar(
    tolerance, "tolerance", lower = 0, lower_open = TRUE
  )
  if (!is.numeric(max_refits) || length(max_refits) != 1L ||
      is.na(max_refits) || max_refits < 0 ||
      (!is.infinite(max_refits) && max_refits != as.integer(max_refits))) {
    stop("`max_refits` must be a nonnegative integer or `Inf`.",
         call. = FALSE)
  }
  subject_ids <- subjects(split_a)
  if (length(subject_ids) < min_subjects) {
    stop("LDR has fewer subjects than `min_subjects`.", call. = FALSE)
  }
  if (folds > length(subject_ids) - 3L) {
    stop("`folds` leaves fewer than three training subjects.",
         call. = FALSE)
  }

  contrast_index <- .ldr_resolve_contrast(split_a, contrast)
  configurations <- .ldr_search_configurations(
    patch_radius_mm, shift_radius_mm
  )
  sp <- space(split_a)
  search <- .ldr_search_units(
    centers = centers,
    configurations = configurations,
    dims3 = sp$dim,
    affine = sp$affine,
    mask_idx = sp$mask_idx
  )
  units <- search$units
  beta_a <- .ldr_matrix_contrast(split_a, "beta", contrast_index)
  beta_b <- .ldr_matrix_contrast(split_b, "beta", contrast_index)
  variance_name_a <- if (!is.null(assay(split_a, "var"))) "var" else "se"
  variance_name_b <- if (!is.null(assay(split_b, "var"))) "var" else "se"
  var_a <- .ldr_matrix_contrast(split_a, variance_name_a, contrast_index)
  var_b <- .ldr_matrix_contrast(split_b, variance_name_b, contrast_index)
  if (identical(variance_name_a, "se")) var_a <- var_a^2
  if (identical(variance_name_b, "se")) var_b <- var_b^2

  screen <- .ldr_tangent_scan(
    beta_a, beta_b, var_a, var_b, units,
    subject_ids = subject_ids,
    folds = folds,
    covariance_shrinkage = covariance_shrinkage
  )
  analyzable <- is.finite(screen$activation_t) &
    is.finite(screen$displacement_t) & is.finite(screen$ordinary_t)
  if (!any(analyzable)) {
    stop("No eligible LDR center produced a full-rank tangent screen.",
         call. = FALSE)
  }
  calibration <- if (n_resamples > 0L) {
    .ldr_tangent_calibrate(
      beta_a, beta_b, var_a, var_b, units,
      subject_ids = subject_ids,
      folds = folds,
      covariance_shrinkage = covariance_shrinkage,
      observed = screen,
      n_resamples = n_resamples,
      seed = seed
    )
  } else NULL
  best <- .ldr_best_search_units(units, screen, calibration)
  best_centers <- vapply(units[best], `[[`, integer(1L), "center")

  admitted <- integer()
  if (!is.null(calibration)) {
    admitted <- best[
      is.finite(calibration$ldr_fwer_p[best]) &
        calibration$ldr_fwer_p[best] <= alpha &
        screen$activation_t[best] > 0 &
        screen$displacement_t[best] > 0
    ]
  }
  if (length(admitted)) {
    rank_order <- order(
      calibration$ldr_fwer_p[admitted],
      -pmin(screen$activation_t[admitted],
            screen$displacement_t[admitted]),
      vapply(units[admitted], `[[`, integer(1L), "center")
    )
    admitted <- admitted[rank_order]
  }
  admitted_total <- length(admitted)
  if (is.finite(max_refits) && length(admitted) > max_refits) {
    admitted <- admitted[seq_len(as.integer(max_refits))]
  }

  exact <- vector("list", length(admitted))
  names(exact) <- as.character(vapply(
    units[admitted], `[[`, integer(1L), "center"
  ))
  for (candidate in seq_along(admitted)) {
    unit_index <- admitted[candidate]
    exact[[candidate]] <- tryCatch(
      .ldr_exact_search_refit(
        beta_a, beta_b, var_a, var_b, units[[unit_index]],
        subject_ids = subject_ids,
        folds = folds,
        iterations = iterations,
        ridge = ridge,
        tolerance = tolerance
      ),
      error = function(e) structure(
        list(message = conditionMessage(e)), class = "ldr_refit_error"
      )
    )
  }

  n_sample <- nrow(beta_a)
  active_index <- sp$mask_idx %||% seq_len(prod(sp$dim))
  output <- function(default = NA_real_) rep(default, n_sample)
  ldr_fwer_p <- activation_fwer_p <- displacement_fwer_p <-
    ordinary_fwer_p <- output()
  activation_t <- displacement_t <- ordinary_t <- tangent_mni_loss <-
    tangent_shift_rms_mm <- rho_shift <- output()
  exact_mni_loss <- exact_shift_rms_mm <- output()
  ldr_flag <- output(0)
  best_sample <- match(best_centers, active_index)
  activation_t[best_sample] <- screen$activation_t[best]
  displacement_t[best_sample] <- screen$displacement_t[best]
  ordinary_t[best_sample] <- screen$ordinary_t[best]
  tangent_mni_loss[best_sample] <- screen$mni_loss[best]
  tangent_shift_rms_mm[best_sample] <- screen$shift_rms_mm[best]
  rho_shift[best_sample] <- screen$rho_shift[best]
  if (!is.null(calibration)) {
    ldr_fwer_p[best_sample] <- calibration$ldr_fwer_p[best]
    activation_fwer_p[best_sample] <-
      calibration$activation_fwer_p[best]
    displacement_fwer_p[best_sample] <-
      calibration$displacement_fwer_p[best]
    ordinary_fwer_p[best_sample] <- calibration$ordinary_fwer_p[best]
  }

  exact_receipts <- vector("list", length(exact))
  names(exact_receipts) <- names(exact)
  for (candidate in seq_along(exact)) {
    unit_index <- admitted[candidate]
    unit <- units[[unit_index]]
    sample_index <- match(unit$center, active_index)
    fit <- exact[[candidate]]
    if (inherits(fit, "ldr_refit_error")) {
      exact_receipts[[candidate]] <- list(
        center = unit$center,
        confirmed = FALSE,
        error = fit$message
      )
      next
    }
    confirmed <- is.finite(fit$activation_t) && fit$activation_t > 0 &&
      all(is.finite(fit$displacement_t_components)) &&
      all(fit$displacement_t_components > 0) &&
      is.finite(fit$mni_loss) && fit$mni_loss > 0
    confirmed <- confirmed && is.finite(fit$prevalence) &&
      fit$prevalence >= min_prevalence
    exact_mni_loss[sample_index] <- fit$mni_loss
    exact_shift_rms_mm[sample_index] <- fit$shift_rms_mm
    if (confirmed && is.finite(ordinary_fwer_p[sample_index]) &&
        ordinary_fwer_p[sample_index] > alpha) {
      ldr_flag[sample_index] <- -log10(ldr_fwer_p[sample_index])
    }
    exact_receipts[[candidate]] <- list(
      center = unit$center,
      config_id = unit$config_id,
      patch_radius_mm = unit$patch_radius_mm,
      shift_radius_mm = unit$shift_radius_mm,
      confirmed = confirmed,
      activation_t = fit$activation_t,
      displacement_t_components = fit$displacement_t_components,
      mni_loss = fit$mni_loss,
      shift_rms_mm = fit$shift_rms_mm,
      prevalence = fit$prevalence
    )
  }

  contrast_label <- contrasts(split_a)[contrast_index]
  make_assay <- function(value) {
    out <- array(value, c(n_sample, 1L, 1L))
    dimnames(out) <- list(NULL, "meta", contrast_label)
    out
  }
  output_assays <- list(
    ldr_fwer_p = make_assay(ldr_fwer_p),
    ldr_flag = make_assay(ldr_flag),
    ldr_activation_fwer_p = make_assay(activation_fwer_p),
    ldr_displacement_fwer_p = make_assay(displacement_fwer_p),
    ordinary_fwer_p = make_assay(ordinary_fwer_p),
    tangent_activation_t = make_assay(activation_t),
    tangent_displacement_t = make_assay(displacement_t),
    ordinary_t = make_assay(ordinary_t),
    tangent_mni_loss = make_assay(tangent_mni_loss),
    tangent_shift_rms_mm = make_assay(tangent_shift_rms_mm),
    rho_shift = make_assay(rho_shift),
    mni_loss = make_assay(exact_mni_loss),
    shift_rms_mm = make_assay(exact_shift_rms_mm)
  )

  unit_table <- data.frame(
    center = vapply(units, `[[`, integer(1L), "center"),
    config_id = vapply(units, `[[`, integer(1L), "config_id"),
    patch_radius_mm = vapply(
      units, `[[`, numeric(1L), "patch_radius_mm"
    ),
    shift_radius_mm = vapply(
      units, `[[`, numeric(1L), "shift_radius_mm"
    ),
    analyzable = analyzable
  )
  merged_lineages <- .merge_provenance_lineages(
    metadata(split_a), metadata(split_b)
  )
  ldr_metadata <- list(
    status = "experimental-whole-search",
    contract = "ldr-v0.2",
    calibrated = !is.null(calibration),
    calibration_scope = if (is.null(centers)) {
      "all-eligible-active-centers"
    } else "prespecified-search-set",
    contrast = contrast_label,
    configurations = configurations,
    requested_center_count = length(search$requested_centers),
    ineligible_unit_count = search$failures,
    searched_units = unit_table,
    best_unit = unit_table[best, , drop = FALSE],
    folds = folds,
    min_subjects = min_subjects,
    covariance_model = "training-split-difference-shrunk-correlation",
    covariance_shrinkage = covariance_shrinkage,
    tangent_estimand = "cross-split-derivative-energy",
    alpha = alpha,
    min_prevalence = min_prevalence,
    calibration = if (is.null(calibration)) {
      list(n_resamples = 0L, seed = seed)
    } else list(
      n_resamples = calibration$n_resamples,
      seed = calibration$seed,
      minimum_p = calibration$minimum_p,
      null_maximum = calibration$null_maximum,
      pointwise_p = lapply(
        calibration$pointwise_p,
        function(value) value[best]
      )
    ),
    exact_refits = exact_receipts,
    admitted_candidates = admitted_total,
    completed_refits = length(admitted),
    refits_truncated = admitted_total > length(admitted),
    max_refits = max_refits,
    flag_definition = paste(
      "screen activation FWER <= alpha AND screen displacement FWER <= alpha",
      "AND exact soft-shift confirmation AND exact mni_loss > 0",
      "AND exact prevalence >= min_prevalence",
      "AND ordinary MNI FWER > alpha"
    )
  )
  metadata_out <- .merge_gds_metadata(
    gds_metadata(),
    list(
      ldr = ldr_metadata,
      synthetic_var = FALSE,
      sample_labels_synthetic = isTRUE(
        metadata(split_a)$sample_labels_synthetic
      ) || isTRUE(metadata(split_b)$sample_labels_synthetic)
    )
  )
  metadata_out$provenance <- merged_lineages$provenance
  metadata_out <- add_provenance_node(
    metadata_out,
    "local_displacement_rescue_map",
    list(
      centers = if (is.null(centers)) NULL else as.integer(centers),
      contrast = contrast_label,
      patch_radius_mm = as.numeric(patch_radius_mm),
      shift_radius_mm = as.numeric(shift_radius_mm),
      folds = folds,
      covariance_shrinkage = covariance_shrinkage,
      n_resamples = n_resamples,
      seed = seed,
      alpha = alpha,
      min_prevalence = min_prevalence,
      max_refits = max_refits
    ),
    inputs = unique(c(
      merged_lineages$parents_a,
      merged_lineages$parents_b
    ))
  )

  new_gds(
    assays = output_assays,
    space = space(split_a),
    subjects = "meta",
    contrasts = contrast_label,
    row_data = row_data(split_a),
    metadata = metadata_out
  )
}
