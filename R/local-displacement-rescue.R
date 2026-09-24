# Local displacement rescue: experimental ROI interface -------------------

.ldr_variance <- function(x) {
  values <- assays(x)
  if (!is.null(values$var)) return(values$var)
  if (!is.null(values$se)) return(values$se^2)
  stop("LDR requires genuine `var` or `se` assays.", call. = FALSE)
}

.ldr_validate_splits <- function(split_a, split_b) {
  if (!inherits(split_a, "gds") || !inherits(split_b, "gds")) {
    stop("`split_a` and `split_b` must be realised GDS objects.",
         call. = FALSE)
  }
  if (!inherits(space(split_a), "space_voxel") ||
      !inherits(space(split_b), "space_voxel")) {
    stop("LDR requires voxel-space GDS inputs.", call. = FALSE)
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
    stop("LDR requires a `beta` assay in both splits.", call. = FALSE)
  }
  var_a <- .ldr_variance(split_a)
  var_b <- .ldr_variance(split_b)
  if (!identical(dim(beta_a), dim(beta_b)) ||
      !identical(dim(beta_a), dim(var_a)) ||
      !identical(dim(beta_b), dim(var_b))) {
    stop("Paired LDR beta and variance assays must have matching dimensions.",
         call. = FALSE)
  }
  synthetic <- isTRUE(metadata(split_a)$synthetic_var) ||
    isTRUE(metadata(split_b)$synthetic_var)
  if (synthetic) {
    stop("LDR does not accept synthetic variance placeholders.", call. = FALSE)
  }
  invisible(TRUE)
}

.ldr_resolve_contrast <- function(x, contrast) {
  labels <- contrasts(x)
  if (is.null(contrast)) {
    if (length(labels) != 1L) {
      stop("`contrast` is required when the input has multiple contrasts.",
           call. = FALSE)
    }
    return(1L)
  }
  if (is.character(contrast) && length(contrast) == 1L &&
      !is.na(contrast)) {
    index <- match(contrast, labels)
  } else if (is.numeric(contrast) && length(contrast) == 1L &&
             is.finite(contrast) && contrast == as.integer(contrast)) {
    index <- as.integer(contrast)
  } else {
    stop("`contrast` must be one contrast name or index.", call. = FALSE)
  }
  if (is.na(index) || index < 1L || index > length(labels)) {
    stop("`contrast` does not identify an input contrast.", call. = FALSE)
  }
  index
}

.ldr_default_tau_grid <- function(shift_radius_mm) {
  sort(unique(shift_radius_mm * c(0.25, 0.5, 1)))
}

#' Experimental local displacement rescue at one voxel
#'
#' Fits an exact, paired-split reference model for a compact signed feature at
#' one prespecified voxel. The function compares a positive aligned feature, a
#' positive randomly displaced feature, and an aligned signed-amplitude
#' heterogeneity model. Subject shifts are integrated rather than selected in
#' held-out data, and every reported score is evaluated by subject
#' cross-fitting.
#'
#' @section Scientific contract:
#' The activation component compares the displaced positive-feature model with
#' no local feature. The displacement component is the smaller held-out improvement
#' over (1) the aligned positive-feature model and (2) the aligned signed-
#' amplitude heterogeneity model. Requiring both comparisons prevents stable
#' polarity mixtures from being attributed to spatial displacement.
#' `mni_loss` is the change in the ordinary one-sample t statistic after a
#' cross-fitted posterior correction that sets displacement to zero while
#' retaining each subject's observed center value. `shift_rms_mm` is the mean
#' RMS scale of the shift prior selected inside the training folds.
#'
#' @section Experimental:
#' This is the exact ROI reference implementation, not whole-brain inference.
#' When `n_resamples > 0`, `ldr_p` is an ROI-level intersection-union p-value:
#' activation is calibrated by paired whole-subject sign flips, while
#' displacement and MNI loss are calibrated by parametric bootstraps from both
#' aligned competitor models. The complete cross-fitting and scale selection
#' procedure is rerun in every resample. No spatial search, FWER correction, or
#' LDR flag is returned. The current noise likelihood uses the supplied diagonal
#' `var`/`se` precision after removing a constant and available linear spatial
#' trends. Spatial-covariance estimation, the tangent screen, and whole-brain
#' maximum-statistic inference remain separate future stages.
#'
#' @param split_a,split_b Realised, independent subject-level split maps with
#'   matching voxel space, subjects, contrasts, `beta`, and genuine `var` or
#'   `se` assays.
#' @param center One-based full-volume voxel index anchoring the positive core.
#'   For packed inputs this is still the index in the full voxel grid.
#' @param contrast Optional contrast name or index. Required for multi-contrast
#'   inputs.
#' @param patch_radius_mm Positive physical radius of the fitted patch. It must
#'   be at least twice `shift_radius_mm`.
#' @param shift_radius_mm Positive maximum translation radius in millimetres.
#' @param tau_grid_mm Positive candidate displacement scales. Selection occurs
#'   only in training folds.
#' @param folds Number of deterministic subject cross-fitting folds.
#' @param min_subjects Minimum number of paired subjects.
#' @param iterations Maximum template EM iterations per model and fold.
#' @param ridge Positive template normal-equation ridge.
#' @param tolerance Positive relative convergence tolerance.
#' @param n_resamples Number of resamples for ROI calibration. Zero skips
#'   calibration and returns `NA` for `ldr_p`. Use at least 199 for exploratory
#'   work and more for small tail probabilities.
#' @param seed Integer seed used only inside calibration. The caller's random
#'   number state is restored.
#'
#' @return A group-level [`gds`] with three assays, finite only at `center`:
#'   `ldr_p`, `mni_loss`, and `shift_rms_mm`. `ldr_p` is `NA` when calibration
#'   is skipped. Supporting scores and model quantities are stored under
#'   `metadata(result)$ldr`.
#' @export
local_displacement_rescue <- function(split_a,
                                      split_b,
                                      center,
                                      contrast = NULL,
                                      patch_radius_mm = 12,
                                      shift_radius_mm = 4,
                                      tau_grid_mm = NULL,
                                      folds = 5L,
                                      min_subjects = 12L,
                                      iterations = 20L,
                                      ridge = 1e-4,
                                      tolerance = 1e-5,
                                      n_resamples = 0L,
                                      seed = 1L) {
  .ldr_validate_splits(split_a, split_b)
  patch_radius_mm <- .ldr_scalar(
    patch_radius_mm, "patch_radius_mm", lower = 0, lower_open = TRUE
  )
  shift_radius_mm <- .ldr_scalar(
    shift_radius_mm, "shift_radius_mm", lower = 0, lower_open = TRUE
  )
  if (patch_radius_mm < 2 * shift_radius_mm) {
    stop("`patch_radius_mm` must be at least twice `shift_radius_mm`.",
         call. = FALSE)
  }
  folds <- as.integer(folds)
  min_subjects <- as.integer(min_subjects)
  iterations <- as.integer(iterations)
  if (length(folds) != 1L || is.na(folds) || folds < 2L) {
    stop("`folds` must be one integer of at least 2.", call. = FALSE)
  }
  if (length(min_subjects) != 1L || is.na(min_subjects) ||
      min_subjects < 6L) {
    stop("`min_subjects` must be one integer of at least 6.", call. = FALSE)
  }
  if (length(iterations) != 1L || is.na(iterations) || iterations < 1L) {
    stop("`iterations` must be one positive integer.", call. = FALSE)
  }
  ridge <- .ldr_scalar(ridge, "ridge", lower = 0, lower_open = TRUE)
  tolerance <- .ldr_scalar(
    tolerance, "tolerance", lower = 0, lower_open = TRUE
  )
  n_resamples <- as.integer(n_resamples)
  seed <- as.integer(seed)
  if (length(n_resamples) != 1L || is.na(n_resamples) || n_resamples < 0L) {
    stop("`n_resamples` must be one nonnegative integer.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }
  if (length(subjects(split_a)) < min_subjects) {
    stop("LDR has fewer subjects than `min_subjects`.", call. = FALSE)
  }
  if (folds > length(subjects(split_a)) - 3L) {
    stop("`folds` leaves fewer than three training subjects.", call. = FALSE)
  }
  if (is.null(tau_grid_mm)) {
    tau_grid_mm <- .ldr_default_tau_grid(shift_radius_mm)
  }
  if (!is.numeric(tau_grid_mm) || !length(tau_grid_mm) ||
      any(!is.finite(tau_grid_mm)) || any(tau_grid_mm <= 0) ||
      any(tau_grid_mm > shift_radius_mm)) {
    stop("`tau_grid_mm` must contain positive scales no larger than `shift_radius_mm`.",
         call. = FALSE)
  }
  tau_grid_mm <- sort(unique(as.numeric(tau_grid_mm)))
  contrast_index <- .ldr_resolve_contrast(split_a, contrast)
  sp <- space(split_a)
  geometry <- .ldr_patch_geometry(
    center = center,
    dims3 = sp$dim,
    affine = sp$affine,
    patch_radius_mm = patch_radius_mm,
    shift_radius_mm = shift_radius_mm,
    mask_idx = sp$mask_idx
  )
  if (!any(geometry$shift_distance_mm > 0)) {
    stop("`shift_radius_mm` permits no nonzero translation in this voxel space.",
         call. = FALSE)
  }
  beta_a <- assay(split_a, "beta")
  beta_b <- assay(split_b, "beta")
  var_a <- .ldr_variance(split_a)
  var_b <- .ldr_variance(split_b)
  patch <- .ldr_prepare_patch(
    beta_a[geometry$sample_idx, , contrast_index, drop = TRUE],
    beta_b[geometry$sample_idx, , contrast_index, drop = TRUE],
    var_a[geometry$sample_idx, , contrast_index, drop = TRUE],
    var_b[geometry$sample_idx, , contrast_index, drop = TRUE],
    geometry
  )
  fit <- .ldr_crossfit(
    patch,
    subject_ids = subjects(split_a),
    tau_grid_mm = tau_grid_mm,
    folds = folds,
    iterations = iterations,
    ridge = ridge,
    tolerance = tolerance
  )
  calibration <- if (n_resamples > 0L) {
    .ldr_calibrate(
      patch,
      subject_ids = subjects(split_a),
      observed = fit,
      tau_grid_mm = tau_grid_mm,
      folds = folds,
      iterations = iterations,
      ridge = ridge,
      tolerance = tolerance,
      n_resamples = n_resamples,
      seed = seed
    )
  } else {
    list(
      p_ldr = NA_real_,
      p_activation = NA_real_,
      p_displacement = NA_real_,
      p_loss = NA_real_,
      component_p = NULL,
      n_resamples = 0L,
      seed = seed,
      minimum_p = NA_real_
    )
  }

  n_sample <- dim(beta_a)[1L]
  active_index <- sp$mask_idx %||% seq_len(prod(sp$dim))
  center_sample <- match(as.integer(center), active_index)
  make_assay <- function(value) {
    out <- array(NA_real_, c(n_sample, 1L, 1L))
    out[center_sample, 1L, 1L] <- value
    dimnames(out) <- list(NULL, "meta", contrasts(split_a)[contrast_index])
    out
  }
  output_assays <- list(
    ldr_p = make_assay(calibration$p_ldr),
    mni_loss = make_assay(fit$mni_loss),
    shift_rms_mm = make_assay(fit$shift_rms_mm)
  )

  merged_lineages <- .merge_provenance_lineages(
    metadata(split_a), metadata(split_b)
  )
  ldr_metadata <- list(
    status = "experimental-reference",
    calibrated = n_resamples > 0L,
    calibration_scope = if (n_resamples > 0L) "prespecified-ROI" else "none",
    scope = "one-center-one-contrast",
    center_full_index = as.integer(center),
    contrast = contrasts(split_a)[contrast_index],
    patch_radius_mm = patch_radius_mm,
    shift_radius_mm = shift_radius_mm,
    tau_grid_mm = tau_grid_mm,
    folds = folds,
    min_subjects = min_subjects,
    iterations = iterations,
    noise_model = "supplied-diagonal-precision",
    comparators = c("no-feature", "aligned-positive", "aligned-signed-amplitude"),
    ordinary_t = fit$ordinary_t,
    ordinary_p_greater = stats::pt(
      fit$ordinary_t,
      df = length(subjects(split_a)) - 1L,
      lower.tail = FALSE
    ),
    counterfactual_t = fit$counterfactual_t,
    prevalence = fit$prevalence,
    activation_score = fit$activation_score,
    displacement_score = fit$displacement_score,
    displacement_components = fit$displacement_components,
    displacement_t_components = fit$displacement_t_components,
    activation_t = fit$activation_t,
    displacement_t = fit$displacement_t,
    calibration = calibration[c(
      "p_ldr", "p_activation", "p_displacement", "p_loss",
      "component_p", "n_resamples", "seed", "minimum_p"
    )],
    fold_id = stats::setNames(fit$fold_id, subjects(split_a)),
    fold_receipts = lapply(fit$folds, function(receipt) {
      list(
        train = receipt$train,
        test = receipt$test,
        aligned = list(
          prevalence = receipt$aligned$prevalence,
          amplitude_scale = receipt$aligned$amplitude_scale,
          converged = receipt$aligned$converged,
          iterations = receipt$aligned$iterations
        ),
        heterogeneity = list(
          prevalence = receipt$heterogeneity$prevalence,
          amplitude_scale = receipt$heterogeneity$amplitude_scale,
          converged = receipt$heterogeneity$converged,
          iterations = receipt$heterogeneity$iterations
        ),
        jitter = list(
          prevalence = receipt$jitter$prevalence,
          amplitude_scale = receipt$jitter$amplitude_scale,
          tau_mm = receipt$jitter$tau_mm,
          shift_rms_mm = receipt$jitter$shift_rms_mm,
          converged = receipt$jitter$converged,
          iterations = receipt$jitter$iterations
        )
      )
    })
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
    "local_displacement_rescue",
    list(
      center = as.integer(center),
      contrast = contrasts(split_a)[contrast_index],
      patch_radius_mm = patch_radius_mm,
      shift_radius_mm = shift_radius_mm,
      tau_grid_mm = tau_grid_mm,
      folds = folds,
      n_resamples = n_resamples,
      seed = seed
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
    contrasts = contrasts(split_a)[contrast_index],
    row_data = row_data(split_a),
    metadata = metadata_out
  )
}
