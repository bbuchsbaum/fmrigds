# =============================================================================
# neuroim2 Interoperability - Import and Export
# =============================================================================
#
# Provides functions to:
# - Import: Convert neuroim2 NeuroVol/NeuroVec objects to GDS
# - Export: Convert GDS objects back to neuroim2 NeuroVol/NeuroVec
# - Split GDS by grouping variables

# =============================================================================
# IMPORT: neuroim2 -> GDS
# =============================================================================

# -----------------------------------------------------------------------------
# as_gds.NeuroVol - Import a single 3D volume
# -----------------------------------------------------------------------------

#' Convert a NeuroVol to GDS
#'
#' Import a single `neuroim2::NeuroVol` (3D brain volume) as a GDS object.
#' The volume becomes a single subject with a single contrast.
#'
#' @param x A `neuroim2::NeuroVol` object.
#' @param assay_name Name for the assay (default: `"beta"`).
#' @param subject Subject identifier (default: `"subject1"`).
#' @param contrast Contrast identifier (default: `"contrast1"`).
#' @param mask Optional mask. Can be:
#'   - `NULL`: use all non-zero voxels as mask
#'   - A `NeuroVol`: use as binary mask
#'   - A logical array: use directly as mask
#'   - `"none"`: include all voxels (no masking)
#' @param ... Additional arguments passed to [new_gds()].
#'
#' @return A `gds` object with voxel space.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vol <- neuroim2::read_vol("beta.nii.gz")
#' gds_obj <- as_gds(vol, subject = "sub-01", contrast = "faces")
#' }
as_gds.NeuroVol <- function(x,
                            assay_name = "beta",
                            subject = "subject1",
                            contrast = "contrast1",
                            mask = NULL,
                            ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required", call. = FALSE)
  }

  arr <- as.array(x)
  vdim <- dim(arr)
  nspace <- neuroim2::space(x)

  # Handle masking
  mask_result <- .resolve_neuroim_mask(arr, mask, vdim)
  mask_idx <- mask_result$mask_idx
  storage <- mask_result$storage

  # Extract data (pack if masked)
  if (storage == "packed") {
    data_vec <- arr[mask_idx]
  } else {
    data_vec <- as.vector(arr)
  }

  # Build GDS arrays: [samples x subjects x contrasts]
  n_samples <- length(data_vec)
  beta_arr <- array(data_vec, dim = c(n_samples, 1L, 1L))

 # Create a dummy variance array (required by GDS)
  var_arr <- array(1, dim = c(n_samples, 1L, 1L))

  # Build space
  affine <- neuroim2::trans(nspace)
  sp <- space_voxel(
    dim = vdim,
    affine = affine,
    mask_idx = if (storage == "packed") mask_idx else NULL,
    storage = storage
  )

  new_gds(
    assays = setNames(list(beta_arr, var_arr), c(assay_name, "var")),
    space = sp,
    subjects = as.character(subject),
    contrasts = as.character(contrast),
    ...
  )
}

#' @rdname as_gds.NeuroVol
#' @export
as_gds.DenseNeuroVol <- as_gds.NeuroVol

# -----------------------------------------------------------------------------
# as_gds.NeuroVec - Import a 4D volume
# -----------------------------------------------------------------------------

#' Convert a NeuroVec to GDS
#'
#' Import a `neuroim2::NeuroVec` (4D brain volume) as a GDS object.
#' The 4th dimension can represent subjects OR contrasts.
#'
#' @param x A `neuroim2::NeuroVec` object.
#' @param assay_name Name for the assay (default: `"beta"`).
#' @param along What does the 4th dimension represent? One of:
#'   - `"subject"`: Each volume is a different subject (single contrast)
#'   - `"contrast"`: Each volume is a different contrast (single subject)
#' @param subjects Character vector of subject IDs. Length must match the 4th
#'   dimension when `along = "subject"`, or length 1 when `along = "contrast"`.
#' @param contrasts Character vector of contrast names. Length must match the 4th
#'   dimension when `along = "contrast"`, or length 1 when `along = "subject"`.
#' @param mask Optional mask. See `as_gds.NeuroVol` for options.
#' @param ... Additional arguments passed to [new_gds()].
#'
#' @return A `gds` object with voxel space.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # 4D volume where each volume is a subject
#' vec4d <- neuroim2::read_vec("all_subjects.nii.gz")
#' gds_obj <- as_gds(vec4d, along = "subject",
#'                   subjects = c("sub-01", "sub-02", "sub-03"))
#'
#' # 4D volume where each volume is a contrast
#' gds_obj <- as_gds(vec4d, along = "contrast",
#'                   contrasts = c("faces", "places", "objects"))
#' }
as_gds.NeuroVec <- function(x,
                            assay_name = "beta",
                            along = c("subject", "contrast"),
                            subjects = NULL,
                            contrasts = NULL,
                            mask = NULL,
                            ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required", call. = FALSE)
  }

  along <- match.arg(along)

  arr <- as.array(x)
  dims <- dim(arr)
  if (length(dims) != 4L) {
    stop("NeuroVec must be 4-dimensional", call. = FALSE
    )
  }

  vdim <- dims[1:3]
  n_vols <- dims[4]
  nspace <- neuroim2::space(x)

  # Infer/validate subjects and contrasts
  if (along == "subject") {
    subjects <- subjects %||% paste0("subject", seq_len(n_vols))
    contrasts <- contrasts %||% "contrast1"
    if (length(subjects) != n_vols) {
      stop("Length of 'subjects' (", length(subjects), ") must equal dim[4] (",
           n_vols, ")", call. = FALSE)
    }
    if (length(contrasts) != 1L) {
      stop("When along='subject', 'contrasts' must have length 1", call. = FALSE)
    }
    n_subj <- n_vols
    n_con <- 1L
  } else {
    subjects <- subjects %||% "subject1"
    contrasts <- contrasts %||% paste0("contrast", seq_len(n_vols))
    if (length(contrasts) != n_vols) {
      stop("Length of 'contrasts' (", length(contrasts), ") must equal dim[4] (",
           n_vols, ")", call. = FALSE)
    }
    if (length(subjects) != 1L) {
      stop("When along='contrast', 'subjects' must have length 1", call. = FALSE)
    }
    n_subj <- 1L
    n_con <- n_vols
  }

  # Handle masking - use first volume to determine mask
  first_vol <- arr[, , , 1]
  mask_result <- .resolve_neuroim_mask(first_vol, mask, vdim)
  mask_idx <- mask_result$mask_idx
  storage <- mask_result$storage

  n_samples <- if (storage == "packed") length(mask_idx) else prod(vdim)

  # Extract and reshape data
  # GDS format: [samples x subjects x contrasts]
  beta_arr <- array(NA_real_, dim = c(n_samples, n_subj, n_con))

  for (i in seq_len(n_vols)) {
    vol_data <- arr[, , , i]
    if (storage == "packed") {
      vec <- vol_data[mask_idx]
    } else {
      vec <- as.vector(vol_data)
    }

    if (along == "subject") {
      beta_arr[, i, 1] <- vec
    } else {
      beta_arr[, 1, i] <- vec
    }
  }

  # Create dummy variance
  var_arr <- array(1, dim = c(n_samples, n_subj, n_con))

  # Build space
  affine <- neuroim2::trans(nspace)
  sp <- space_voxel(
    dim = vdim,
    affine = affine,
    mask_idx = if (storage == "packed") mask_idx else NULL,
    storage = storage
  )

  new_gds(
    assays = setNames(list(beta_arr, var_arr), c(assay_name, "var")),
    space = sp,
    subjects = as.character(subjects),
    contrasts = as.character(contrasts),
    ...
  )
}

#' @rdname as_gds.NeuroVec
#' @export
as_gds.DenseNeuroVec <- as_gds.NeuroVec

# -----------------------------------------------------------------------------
# gds_from_neurovols - Import paired lists of NeuroVols (single contrast)
# -----------------------------------------------------------------------------

#' Create GDS from lists of NeuroVol objects
#'
#' Import matched lists of NeuroVol objects (e.g., beta and variance maps)
#' into a GDS. This handles the common case of multiple subjects with a
#' single contrast. For multiple contrasts, use [gds_from_neurovol_nested()].
#'
#' @param beta A named list of `NeuroVol` objects for beta/effect estimates.
#'   Names become subject IDs.
#' @param var Optional named list of `NeuroVol` objects for variance estimates.
#'   Must have same names as `beta`.
#' @param se Optional named list of `NeuroVol` objects for standard error.
#'   Must have same names as `beta`. Provide either `var` or `se`, not both.
#' @param contrast Contrast name (default: `"contrast1"`).
#' @param mask Optional mask. See `as_gds.NeuroVol` for options.
#' @param col_data Optional data frame of subject-level covariates.
#'   Row names must match subject IDs (names of `beta`).
#' @param ... Additional arguments passed to [new_gds()].
#'
#' @return A `gds` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load beta and variance maps for each subject
#' beta_vols <- list(
#'   `sub-01` = neuroim2::read_vol("sub-01_beta.nii.gz"),
#'   `sub-02` = neuroim2::read_vol("sub-02_beta.nii.gz")
#' )
#' var_vols <- list(
#'   `sub-01` = neuroim2::read_vol("sub-01_var.nii.gz"),
#'   `sub-02` = neuroim2::read_vol("sub-02_var.nii.gz")
#' )
#'
#' gds_obj <- gds_from_neurovols(beta_vols, var = var_vols)
#'
#' # With subject covariates
#' covars <- data.frame(age = c(25, 30), group = c("A", "B"),
#'                      row.names = c("sub-01", "sub-02"))
#' gds_obj <- gds_from_neurovols(beta_vols, var = var_vols, col_data = covars)
#' }
gds_from_neurovols <- function(beta,
                                var = NULL,
                                se = NULL,
                                contrast = "contrast1",
                                mask = NULL,
                                col_data = NULL,
                                ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required", call. = FALSE)
  }

  # Validate inputs
  if (!is.list(beta) || !length(beta)) {
    stop("'beta' must be a non-empty list of NeuroVol objects", call. = FALSE)
  }
  if (is.null(names(beta)) || any(!nzchar(names(beta)))) {
    stop("'beta' must be a named list (names become subject IDs)", call. = FALSE)
  }

  subjects <- names(beta)
  n_subj <- length(subjects)

  # Check all are NeuroVol
  if (!all(vapply(beta, inherits, logical(1L), "NeuroVol"))) {
    stop("All elements of 'beta' must be NeuroVol objects", call. = FALSE)
  }

  # Validate var/se
  has_var <- !is.null(var)
  has_se <- !is.null(se)
  if (has_var && has_se) {
    stop("Provide either 'var' or 'se', not both", call. = FALSE)
  }
  if (!has_var && !has_se) {
    .warn_synthetic_variance(once = FALSE)
  }

  if (has_var) {
    if (!setequal(names(var), subjects)) {
      stop("Names of 'var' must match names of 'beta'", call. = FALSE)
    }
    var <- var[subjects]  # Reorder to match
  }
  if (has_se) {
    if (!setequal(names(se), subjects)) {
      stop("Names of 'se' must match names of 'beta'", call. = FALSE)
    }
    se <- se[subjects]
  }

  # Get dimensions from first volume
  first_vol <- beta[[1]]
  vdim <- dim(as.array(first_vol))
  nspace <- neuroim2::space(first_vol)

  # Handle masking
  mask_result <- .resolve_neuroim_mask(as.array(first_vol), mask, vdim)
  mask_idx <- mask_result$mask_idx
  storage <- mask_result$storage

  n_samples <- if (storage == "packed") length(mask_idx) else prod(vdim)

  # Build arrays
  beta_arr <- array(NA_real_, dim = c(n_samples, n_subj, 1L))
  var_arr <- array(NA_real_, dim = c(n_samples, n_subj, 1L))

  for (j in seq_len(n_subj)) {
    subj <- subjects[j]

    # Beta
    vol_data <- as.array(beta[[subj]])
    if (storage == "packed") {
      beta_arr[, j, 1] <- vol_data[mask_idx]
    } else {
      beta_arr[, j, 1] <- as.vector(vol_data)
    }

    # Variance
    if (has_var) {
      vol_data <- as.array(var[[subj]])
      if (storage == "packed") {
        var_arr[, j, 1] <- vol_data[mask_idx]
      } else {
        var_arr[, j, 1] <- as.vector(vol_data)
      }
    } else if (has_se) {
      vol_data <- as.array(se[[subj]])
      if (storage == "packed") {
        var_arr[, j, 1] <- vol_data[mask_idx]^2
      } else {
        var_arr[, j, 1] <- as.vector(vol_data)^2
      }
    } else {
      var_arr[, j, 1] <- 1
    }
  }

  # Build space
  affine <- neuroim2::trans(nspace)
  sp <- space_voxel(
    dim = vdim,
    affine = affine,
    mask_idx = if (storage == "packed") mask_idx else NULL,
    storage = storage
  )

  # Tag the synthetic unit-variance placeholder so variance-weighted reducers
  # can refuse it (issue #5), both at the array level (consumption-site
  # backstop) and in metadata (reduce()-verb guard).
  synthetic_var <- !has_var && !has_se
  if (synthetic_var) attr(var_arr, "synthetic_unit_variance") <- TRUE

  dots <- list(...)
  meta <- dots$metadata %||% list()
  if (synthetic_var) meta$synthetic_var <- TRUE
  dots$metadata <- NULL

  do.call(new_gds, c(list(
    assays = list(beta = beta_arr, var = var_arr),
    space = sp,
    subjects = subjects,
    contrasts = as.character(contrast),
    col_data = col_data,
    metadata = meta
  ), dots))
}

# -----------------------------------------------------------------------------
# gds_from_neurovol_nested - Full subjects x contrasts import
# -----------------------------------------------------------------------------

#' Create GDS from nested structure of NeuroVol objects
#'
#' The most flexible import function for creating GDS objects from neuroim2
#' data. Supports multiple subjects and multiple contrasts, with proper
#' variance/SE handling.
#'
#' @param beta A nested list structure for beta/effect estimates. Can be:
#'   - Nested list: `list(subject1 = list(con1 = vol, con2 = vol), ...)`
#'   - List of 4D NeuroVec: `list(subject1 = vec4d, subject2 = vec4d, ...)`
#'     where each NeuroVec has contrasts along the 4th dimension
#' @param var Optional matching structure for variance estimates.
#' @param se Optional matching structure for standard errors (alternative to var).
#' @param contrasts Character vector of contrast names. Required when using
#'   NeuroVec inputs; inferred from inner list names otherwise.
#' @param mask Optional mask. See `as_gds.NeuroVol` for options.
#' @param col_data Optional data frame of subject-level covariates.
#' @param metadata Optional list of additional metadata to attach.
#' @param ... Additional arguments passed to [new_gds()].
#'
#' @return A `gds` object with dimensions `samples x subjects x contrasts`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Nested list structure (subjects x contrasts)
#' beta <- list(
#'   `sub-01` = list(
#'     faces = neuroim2::read_vol("sub-01_faces_beta.nii.gz"),
#'     places = neuroim2::read_vol("sub-01_places_beta.nii.gz")
#'   ),
#'   `sub-02` = list(
#'     faces = neuroim2::read_vol("sub-02_faces_beta.nii.gz"),
#'     places = neuroim2::read_vol("sub-02_places_beta.nii.gz")
#'   )
#' )
#' var <- list(
#'   `sub-01` = list(
#'     faces = neuroim2::read_vol("sub-01_faces_var.nii.gz"),
#'     places = neuroim2::read_vol("sub-01_places_var.nii.gz")
#'   ),
#'   `sub-02` = list(
#'     faces = neuroim2::read_vol("sub-02_faces_var.nii.gz"),
#'     places = neuroim2::read_vol("sub-02_places_var.nii.gz")
#'   )
#' )
#'
#' gds_obj <- gds_from_neurovol_nested(beta, var = var)
#'
#' # Example 2: List of 4D NeuroVecs (each subject has multi-contrast 4D file)
#' beta_vecs <- list(
#'   `sub-01` = neuroim2::read_vec("sub-01_allcons_beta.nii.gz"),
#'   `sub-02` = neuroim2::read_vec("sub-02_allcons_beta.nii.gz")
#' )
#' gds_obj <- gds_from_neurovol_nested(beta_vecs,
#'                                      contrasts = c("faces", "places", "objects"))
#' }
gds_from_neurovol_nested <- function(beta,
                                  var = NULL,
                                  se = NULL,
                                  contrasts = NULL,
                                  mask = NULL,
                                  col_data = NULL,
                                  metadata = list(),
                                  ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required", call. = FALSE)
  }

  # Validate beta structure
  if (!is.list(beta) || !length(beta)) {
    stop("'beta' must be a non-empty named list", call. = FALSE)
  }
  if (is.null(names(beta)) || any(!nzchar(names(beta)))) {
    stop("'beta' must be a named list (names become subject IDs)", call. = FALSE)
  }

  subjects <- names(beta)
  n_subj <- length(subjects)

  # Detect input type: nested list of NeuroVol or list of NeuroVec
  first_elem <- beta[[1]]
  is_neurovec_input <- inherits(first_elem, "NeuroVec")
  is_nested_list <- is.list(first_elem) && !inherits(first_elem, "NeuroVol")

  if (!is_neurovec_input && !is_nested_list && !inherits(first_elem, "NeuroVol")) {
    stop("'beta' elements must be NeuroVol, NeuroVec, or nested lists", call. = FALSE)
  }

  # Handle single-contrast case (flat list of NeuroVol)
  if (inherits(first_elem, "NeuroVol") && !is_nested_list) {
    return(gds_from_neurovols(beta, var = var, se = se,
                               contrast = contrasts %||% "contrast1",
                               mask = mask, col_data = col_data, ...))
  }

  # Determine contrasts
  if (is_neurovec_input) {
    if (is.null(contrasts)) {
      # Infer from 4th dimension
      n_con <- dim(as.array(first_elem))[4]
      contrasts <- paste0("contrast", seq_len(n_con))
    }
    n_con <- length(contrasts)
  } else {
    # Nested list - get contrasts from inner list names
    if (is.null(contrasts)) {
      contrasts <- names(first_elem)
      if (is.null(contrasts)) {
        stop("Inner lists must be named (contrast names) or provide 'contrasts' argument",
             call. = FALSE)
      }
    }
    n_con <- length(contrasts)
  }

  # Validate var/se
  has_var <- !is.null(var)
  has_se <- !is.null(se)
  if (has_var && has_se) {
    stop("Provide either 'var' or 'se', not both", call. = FALSE)
  }
  synthetic_var <- !has_var && !has_se
  if (synthetic_var) {
    .warn_synthetic_variance(once = FALSE)
    metadata$synthetic_var <- TRUE
  }

  # Get dimensions from first volume
  if (is_neurovec_input) {
    first_vol_arr <- as.array(first_elem)[, , , 1]
  } else {
    first_vol_arr <- as.array(first_elem[[1]])
  }
  vdim <- dim(first_vol_arr)

  # Get space from first element
  if (is_neurovec_input) {
    nspace <- neuroim2::space(first_elem)
  } else {
    nspace <- neuroim2::space(first_elem[[1]])
  }

  # Handle masking
  mask_result <- .resolve_neuroim_mask(first_vol_arr, mask, vdim)
  mask_idx <- mask_result$mask_idx
  storage <- mask_result$storage

  n_samples <- if (storage == "packed") length(mask_idx) else prod(vdim)

  # Build arrays: [samples x subjects x contrasts]
  beta_arr <- array(NA_real_, dim = c(n_samples, n_subj, n_con))
  var_arr <- array(NA_real_, dim = c(n_samples, n_subj, n_con))

  # Extract data
  for (j in seq_len(n_subj)) {
    subj <- subjects[j]
    beta_subj <- beta[[subj]]
    var_subj <- if (has_var) var[[subj]] else NULL
    se_subj <- if (has_se) se[[subj]] else NULL

    for (k in seq_len(n_con)) {
      con <- contrasts[k]

      # Get beta volume data
      if (is_neurovec_input) {
        vol_data <- as.array(beta_subj)[, , , k]
      } else {
        vol_data <- as.array(beta_subj[[con]])
      }

      # Pack if needed
      if (storage == "packed") {
        beta_arr[, j, k] <- vol_data[mask_idx]
      } else {
        beta_arr[, j, k] <- as.vector(vol_data)
      }

      # Get variance/SE
      if (has_var) {
        if (is_neurovec_input) {
          vol_data <- as.array(var_subj)[, , , k]
        } else {
          vol_data <- as.array(var_subj[[con]])
        }
        if (storage == "packed") {
          var_arr[, j, k] <- vol_data[mask_idx]
        } else {
          var_arr[, j, k] <- as.vector(vol_data)
        }
      } else if (has_se) {
        if (is_neurovec_input) {
          vol_data <- as.array(se_subj)[, , , k]
        } else {
          vol_data <- as.array(se_subj[[con]])
        }
        if (storage == "packed") {
          var_arr[, j, k] <- vol_data[mask_idx]^2
        } else {
          var_arr[, j, k] <- as.vector(vol_data)^2
        }
      } else {
        var_arr[, j, k] <- 1
      }
    }
  }

  # Build space
  affine <- neuroim2::trans(nspace)
  sp <- space_voxel(
    dim = vdim,
    affine = affine,
    mask_idx = if (storage == "packed") mask_idx else NULL,
    storage = storage
  )

  if (synthetic_var) attr(var_arr, "synthetic_unit_variance") <- TRUE

  new_gds(
    assays = list(beta = beta_arr, var = var_arr),
    space = sp,
    subjects = subjects,
    contrasts = contrasts,
    col_data = col_data,
    metadata = metadata,
    ...
  )
}

# -----------------------------------------------------------------------------
# Helper: Resolve mask for neuroim2 import
# -----------------------------------------------------------------------------

.resolve_neuroim_mask <- function(arr, mask, vdim) {
  if (is.null(mask)) {
    # Default: mask to non-zero voxels
    mask_idx <- which(arr != 0)
    if (length(mask_idx) == 0) {
      # All zeros - use all voxels
      return(list(mask_idx = NULL, storage = "dense"))
    }
    return(list(mask_idx = mask_idx, storage = "packed"))
  }

  if (is.character(mask) && mask == "none") {
    return(list(mask_idx = NULL, storage = "dense"))
  }

  if (inherits(mask, "NeuroVol")) {
    mask_arr <- as.array(mask) > 0
    mask_idx <- which(mask_arr)
    return(list(mask_idx = mask_idx, storage = "packed"))
  }

  if (is.logical(mask) && is.array(mask)) {
    if (!identical(dim(mask), vdim)) {
      stop("Mask dimensions must match volume dimensions", call. = FALSE)
    }
    mask_idx <- which(mask)
    return(list(mask_idx = mask_idx, storage = "packed"))
  }

  stop("Invalid mask specification", call. = FALSE)
}

# =============================================================================
# EXPORT: GDS -> neuroim2
# =============================================================================

# -----------------------------------------------------------------------------
# as_neurovol_list - Extract assay as list of NeuroVol objects
# -----------------------------------------------------------------------------

#' Convert GDS assay to list of NeuroVol objects
#'
#' Extract a single assay from a GDS object and return it as a list of
#' `neuroim2::NeuroVol` objects, one per subject/contrast combination.
#'
#' @param x A GDS object with a voxel space.
#' @param assay Name of the assay to extract (default: `"beta"`).
#' @param by How to organize the output list. One of:
#'   - `"subject"`: list indexed by subject, each element is a NeuroVol
#'     (only works when there's a single contrast)
#'   - `"contrast"`: list indexed by contrast, each element is a NeuroVol
#'     (only works when there's a single subject)
#'
#' @param drop_dim If `TRUE` (default), drop singleton dimensions when
#'   extracting. If `FALSE`, always return a nested list structure.
#'
#' @return A named list of `neuroim2::NeuroVol` objects. The structure depends
#'   on the `by` parameter and the dimensions of the data:
#'   - Single contrast, `by = "subject"`: list of NeuroVol named by subject
#'   - Single subject, `by = "contrast"`: list of NeuroVol named by contrast
#'   - Multiple of both: nested list `[[subject]][[contrast]]`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract beta maps as list of NeuroVols (one per subject)
#' vols <- as_neurovol_list(gds, assay = "beta")
#' vols$sub01  # NeuroVol for subject "sub01"
#'
#' # Extract by contrast
#' vols_by_con <- as_neurovol_list(gds, assay = "t", by = "contrast")
#' }
as_neurovol_list <- function(x,
                              assay = "beta",
                              by = c("subject", "contrast"),
                              drop_dim = TRUE) {
  stopifnot(inherits(x, "gds"))

  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required for as_neurovol_list()", call. = FALSE)
  }

  sp <- space(x)
  if (!inherits(sp, "space_voxel"))
    stop("as_neurovol_list() requires a GDS with voxel space", call. = FALSE)

  by <- match.arg(by)

  arr <- assay(x, assay)
  if (is.null(arr)) {
    stop("Assay '", assay, "' not found in GDS", call. = FALSE)
  }

  dims <- dim(arr)  # [samples, subjects, contrasts]
  subj_ids <- subjects(x)
  con_ids <- contrasts(x)

  # Build neuroim2 space from GDS space
  neuro_space <- .gds_space_to_neuroim(sp)

  # Helper to unpack a single vector to a NeuroVol
  unpack_vol <- function(vec) {
    full <- numeric(prod(sp$dim))
    if (!is.null(sp$mask_idx) && identical(sp$storage, "packed")) {
      full[sp$mask_idx] <- vec
    } else {
      full <- vec
    }
    vol_arr <- array(full, dim = sp$dim)
    neuroim2::NeuroVol(vol_arr, neuro_space)
  }

  # Build the list structure
  n_subj <- dims[2]
  n_con <- dims[3]

  if (n_subj == 1 && n_con == 1) {
    # Single subject, single contrast: return a single NeuroVol in a list
    result <- list(unpack_vol(arr[, 1, 1]))
    names(result) <- subj_ids[1]
    return(result)

  }

  if (by == "subject") {
    if (n_con == 1 || !drop_dim) {
      # One contrast or nested: list by subject
      result <- lapply(seq_len(n_subj), function(j) {
        if (n_con == 1 && drop_dim) {
          unpack_vol(arr[, j, 1])
        } else {
          # Nested: list of contrasts per subject
          setNames(
            lapply(seq_len(n_con), function(k) unpack_vol(arr[, j, k])),
            con_ids
          )
        }
      })
      names(result) <- subj_ids
    } else {
      stop("Multiple contrasts present. Use by = 'contrast' or set drop_dim = FALSE",
           call. = FALSE)
    }
  } else {
    # by == "contrast"
    if (n_subj == 1 || !drop_dim) {
      result <- lapply(seq_len(n_con), function(k) {
        if (n_subj == 1 && drop_dim) {
          unpack_vol(arr[, 1, k])
        } else {
          # Nested: list of subjects per contrast
          setNames(
            lapply(seq_len(n_subj), function(j) unpack_vol(arr[, j, k])),
            subj_ids
          )
        }
      })
      names(result) <- con_ids
    } else {
      stop("Multiple subjects present. Use by = 'subject' or set drop_dim = FALSE",
           call. = FALSE)
    }
  }

  result
}

# -----------------------------------------------------------------------------
# as_neurovec - Extract assay as 4D NeuroVec
# -----------------------------------------------------------------------------

#' Convert GDS assay to a 4D NeuroVec
#'
#' Extract a single assay from a GDS object and return it as a
#' `neuroim2::NeuroVec` (4D neuroimaging array), stacking along the 4th
#' dimension.
#'
#' @param x A GDS object with a voxel space.
#' @param assay Name of the assay to extract (default: `"beta"`).
#' @param along Which dimension to stack along the 4th dimension. One of:
#'   - `"subject"`: Stack subjects (collapses contrasts if multiple).
#'   - `"contrast"`: Stack contrasts (collapses subjects if multiple).
#'   - `"both"`: Stack all subject x contrast combinations.
#'
#' @param subset_subjects Optional character vector of subject IDs to include.
#' @param subset_contrasts Optional character vector of contrast names to include
#'
#' @return A `neuroim2::NeuroVec` object with dimensions `[x, y, z, n]` where
#'   `n` depends on the `along` parameter.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Stack all subjects into a 4D volume
#' vec4d <- as_neurovec(gds, assay = "beta", along = "subject")
#'
#' # Stack specific contrasts
#' vec4d <- as_neurovec(gds, assay = "t", along = "contrast",
#'                      subset_contrasts = c("faces", "places"))
#' }
as_neurovec <- function(x,
                        assay = "beta",
                        along = c("subject", "contrast", "both"),
                        subset_subjects = NULL,
                        subset_contrasts = NULL) {
  stopifnot(inherits(x, "gds"))

  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required for as_neurovec()", call. = FALSE)
  }

  sp <- space(x)
  if (!inherits(sp, "space_voxel")) {
    stop("as_neurovec() requires a GDS with voxel space", call. = FALSE)
  }

  along <- match.arg(along)

  arr <- assay(x, assay)
  if (is.null(arr)) {
    stop("Assay '", assay, "' not found in GDS", call. = FALSE)
  }

  dims <- dim(arr)  # [samples, subjects, contrasts]
  subj_ids <- subjects(x)
  con_ids <- contrasts(x)

  # Apply subsetting
  subj_idx <- if (!is.null(subset_subjects)) {
    match(subset_subjects, subj_ids)
  } else {
    seq_along(subj_ids)
  }
  con_idx <- if (!is.null(subset_contrasts)) {
    match(subset_contrasts, con_ids)
  } else {
    seq_along(con_ids)
  }

  if (any(is.na(subj_idx))) {
    stop("Unknown subjects: ",
         paste(subset_subjects[is.na(subj_idx)], collapse = ", "),
         call. = FALSE)
  }
  if (any(is.na(con_idx))) {
    stop("Unknown contrasts: ",
         paste(subset_contrasts[is.na(con_idx)], collapse = ", "),
         call. = FALSE)
  }

  arr <- arr[, subj_idx, con_idx, drop = FALSE]
  subj_ids <- subj_ids[subj_idx]
  con_ids <- con_ids[con_idx]
  dims <- dim(arr)

  # Build neuroim2 space
  neuro_space <- .gds_space_to_neuroim(sp)

  # Determine stacking
  n_subj <- dims[2]
  n_con <- dims[3]
  n_vox <- prod(sp$dim)

  if (along == "subject") {
    if (n_con > 1) {
      warning("Multiple contrasts present; using first contrast only. ",
              "Use along = 'both' to include all.", call. = FALSE)
    }
    n_vols <- n_subj
    vol_labels <- subj_ids
    vec_data <- arr[, , 1, drop = FALSE]  # [samples, subjects, 1]
    dim(vec_data) <- c(dims[1], n_subj)
  } else if (along == "contrast") {
    if (n_subj > 1) {
      warning("Multiple subjects present; using first subject only. ",
              "Use along = 'both' to include all.", call. = FALSE)
    }
    n_vols <- n_con
    vol_labels <- con_ids
    vec_data <- arr[, 1, , drop = FALSE]  # [samples, 1, contrasts]
    dim(vec_data) <- c(dims[1], n_con)
  } else {
    # along == "both": flatten subject x contrast
    n_vols <- n_subj * n_con
    vol_labels <- as.vector(outer(subj_ids, con_ids, paste, sep = "_"))
    vec_data <- matrix(arr, nrow = dims[1], ncol = n_vols)
  }

  # Unpack to full voxel grid
  full_data <- array(0, dim = c(sp$dim, n_vols))

  for (i in seq_len(n_vols)) {
    vec <- vec_data[, i]
    if (!is.null(sp$mask_idx) && identical(sp$storage, "packed")) {
      full_vec <- numeric(n_vox)
      full_vec[sp$mask_idx] <- vec
    } else {
      full_vec <- vec
    }
    full_data[, , , i] <- array(full_vec, dim = sp$dim)
  }

  # Create 4D NeuroSpace
  neuro_space_4d <- neuroim2::NeuroSpace(
    dim = c(sp$dim, n_vols),
    spacing = neuroim2::spacing(neuro_space),
    origin = neuroim2::origin(neuro_space),
    trans = sp$affine
  )

  neuroim2::NeuroVec(full_data, neuro_space_4d)
}

# -----------------------------------------------------------------------------
# split.gds - Split GDS by grouping variable
# -----------------------------------------------------------------------------

#' Split a GDS object by a grouping variable
#'
#' Divide a GDS object into a list of GDS objects based on a grouping factor
#' from `col_data` (subject-level metadata).
#'
#' @param x A GDS object.
#' @param f A factor or character vector for splitting, OR a column name
#'   (as a string) from `col_data(x)`.
#' @param drop Logical; if `TRUE`, drop unused factor levels. Default `TRUE`.
#' @param ... Additional arguments (ignored).
#'
#' @return A named list of GDS objects, one per level of `f`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Split by group column in col_data
#' gds_by_group <- split(gds, "group")
#' gds_by_group$patients  # GDS with only patient subjects
#' gds_by_group$controls  # GDS with only control subjects
#'
#' # Or provide a factor directly
#' groups <- factor(c("A", "B", "A", "B"))
#' gds_split <- split(gds, groups)
#' }
split.gds <- function(x, f, drop = TRUE, ...) {

  stopifnot(inherits(x, "gds"))

  cd <- col_data(x)
  subj_ids <- subjects(x)

  # If f is a column name, extract from col_data
  if (is.character(f) && length(f) == 1L && f %in% names(cd)) {
    f <- cd[[f]]
  }

  # Validate f
  if (length(f) != length(subj_ids)) {
    stop("Length of 'f' (", length(f), ") must equal number of subjects (",
         length(subj_ids), ")", call. = FALSE)
  }

  f <- as.factor(f)
  if (drop) f <- droplevels(f)

  levels_f <- levels(f)

  result <- lapply(levels_f, function(lev) {
    idx <- which(f == lev)
    .subset_gds_subjects(x, idx)
  })

  names(result) <- levels_f
  result
}

# -----------------------------------------------------------------------------
# Helper: subset GDS by subject indices
# -----------------------------------------------------------------------------

.subset_gds_subjects <- function(x, subj_idx) {
  # Subset all assays along subject dimension (dim 2)
  new_assays <- lapply(assays(x), function(a) {
    a[, subj_idx, , drop = FALSE]
  })

  new_subjects <- subjects(x)[subj_idx]

  # Subset col_data
  cd <- col_data(x)
  new_col_data <- if (!is.null(cd) && nrow(cd) > 0) {
    cd[subj_idx, , drop = FALSE]
  } else {
    NULL
  }

  new_gds(
    assays = new_assays,
    space = space(x),
    subjects = new_subjects,
    contrasts = contrasts(x),
    col_data = new_col_data,
    row_data = row_data(x),
    metadata = metadata(x)
  )
}

# -----------------------------------------------------------------------------
# Helper: Convert GDS space to neuroim2 NeuroSpace
# -----------------------------------------------------------------------------

.gds_space_to_neuroim <- function(sp) {
  if (!inherits(sp, "space_voxel")) {
    stop("Only voxel spaces can be converted to neuroim2 NeuroSpace", call. = FALSE)
  }

  # Extract spacing from affine (diagonal elements of upper-left 3x3)
  # This is a simplification; full affine may have rotations
  spacing <- abs(c(sp$affine[1, 1], sp$affine[2, 2], sp$affine[3, 3]))

  # Origin is the translation column
  origin <- sp$affine[1:3, 4]

  neuroim2::NeuroSpace(
    dim = sp$dim,
    spacing = spacing,
    origin = origin,
    trans = sp$affine
  )
}

# -----------------------------------------------------------------------------
# Convenience: extract_group - Get NeuroVols for a specific group
# -----------------------------------------------------------------------------

#' Extract NeuroVol list for a specific group
#'
#' Convenience function that combines [split.gds()] and [as_neurovol_list()]
#' to extract neuroimaging volumes for subjects in a specific group.
#'
#' @param x A GDS object.
#' @param group_var Column name in `col_data(x)` containing group labels.
#' @param group_level The specific group level to extract.
#' @param assay Name of the assay to extract (default: `"beta"`).
#' @param ... Additional arguments passed to [as_neurovol_list()].
#'
#' @return A named list of `neuroim2::NeuroVol` objects for subjects in the
#'   specified group.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get beta maps for all patients
#' patient_vols <- extract_group(gds, "diagnosis", "patient", assay = "beta")
#'
#' # Get t-maps for controls
#' control_ts <- extract_group(gds, "diagnosis", "control", assay = "t")
#' }
extract_group <- function(x, group_var, group_level, assay = "beta", ...) {
  stopifnot(inherits(x, "gds"))

  cd <- col_data(x)
  if (!group_var %in% names(cd)) {
    stop("Column '", group_var, "' not found in col_data", call. = FALSE)
  }

  grp_col <- cd[[group_var]]
  if (!group_level %in% grp_col) {
    stop("Group level '", group_level, "' not found in column '", group_var, "'",
         call. = FALSE)
  }

  # Subset to the group
  subj_idx <- which(grp_col == group_level)
  gds_subset <- .subset_gds_subjects(x, subj_idx)

  # Convert to NeuroVol list
  as_neurovol_list(gds_subset, assay = assay, ...)
}
