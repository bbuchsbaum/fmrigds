#' Register NIfTI Adapter
#'
#' Registers the NIfTI adapter for reading neuroimaging data in NIfTI format.
#' Uses the neuroim2 package for reading NIfTI files.
#'
#' @details
#' The NIfTI adapter supports:
#' - 3D and 4D NIfTI files (.nii, .nii.gz)
#' - Optional mask specification (separate NIfTI file with 1=included, 0=excluded)
#' - Automatic spatial information extraction (affine transformation)
#' - BIDS-compatible file naming
#'
#' To provide a mask, pass it to gds() via the mask parameter:
#' \code{gds(source, mask = "path/to/mask.nii")}
#'
#' @return NULL (invisibly)
#' @keywords internal
register_nifti_adapter <- function() {
  register_adapter(
    name = "nifti",
    detect = .nifti_detect,
    open = .nifti_open,
    probe = .nifti_probe,
    read = .nifti_read,
    close = .nifti_close
  )
}

.nifti_detect <- function(source) {
  if (!is.character(source)) return(FALSE)
  files <- .nifti_normalise_source(source)
  if (!length(files)) return(FALSE)
  if (all(file.exists(files))) return(0.9)
  FALSE
}

.nifti_normalise_source <- function(source) {
  if (length(source) == 1L && dir.exists(source)) {
    list.files(source, pattern = "\\.nii(\\.gz)?$", full.names = TRUE)
  } else {
    source
  }
}

.nifti_open <- function(source, mode = "r", ...) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required for NIfTI adapter", call. = FALSE)
  }
  files <- .nifti_normalise_source(source)
  if (!length(files)) stop("No NIfTI files found", call. = FALSE)

  list(files = files)
}

.nifti_probe <- function(handle, ...) {
  # Read header using neuroim2
  meta <- neuroim2::read_header(handle$files[1])

  # Get dimensions from metadata
  dim_all <- meta@dims

  if (length(dim_all) == 3L) {
    spatial_dim <- dim_all
    n_contrasts <- 1L
  } else if (length(dim_all) == 4L) {
    spatial_dim <- dim_all[1:3]
    n_contrasts <- dim_all[4]
  } else {
    stop("Unsupported NIfTI dimensionality: ", length(dim_all), call. = FALSE)
  }

  # Load the first volume to get spatial information
  # For 3D files, read directly; for 4D, extract first volume
  first_vol <- if (n_contrasts == 1L) {
    neuroim2::read_vol(handle$files[1], index = 1)
  } else {
    # For 4D files, use sub_vector to extract first volume
    vec <- neuroim2::read_vec(handle$files[1])
    neuroim2::sub_vector(vec, 1)
  }

  # Load mask from separate file if provided via ... arguments
  # Otherwise create a full mask (all voxels included)
  dots <- list(...)
  mask <- dots$mask

  if (!is.null(mask)) {
    # Mask provided as file path or NeuroVol object
    mask_vol <- if (is.character(mask)) {
      neuroim2::read_vol(mask)
    } else if (inherits(mask, "NeuroVol")) {
      mask
    } else {
      stop("mask must be a file path or NeuroVol object", call. = FALSE)
    }

    # Convert to logical bitmap (1 = included, 0 = excluded)
    mask_bitmap <- as.array(mask_vol) > 0
  } else {
    # No mask provided: include all voxels
    mask_bitmap <- array(TRUE, dim = spatial_dim)
  }

  mask_idx <- which(as.vector(mask_bitmap))

  subjects <- .nifti_subject_ids(handle$files)
  contrasts <- if (n_contrasts == 1L) {
    "contrast1"
  } else {
    paste0("contrast", seq_len(n_contrasts))
  }

  # Extract affine transformation from neuroim2 NeuroSpace
  nspace <- neuroim2::space(first_vol)
  affine <- neuroim2::trans(nspace)

  space <- space_voxel(
    dim = spatial_dim,
    affine = affine,
    mask_bitmap = mask_bitmap,
    storage = "packed"
  )

  list(
    assays = c("beta"),
    dims = c(sample = length(mask_idx), subject = length(handle$files), contrast = n_contrasts),
    subjects = subjects,
    contrasts = contrasts,
    space = space,
    maps = list(),
    metadata = list(schema_version = "0.1.0", source_files = handle$files),
    columns = list(effect_cols = NULL, subject_col = NULL, sample_col = NULL, contrast_col = NULL),
    mask_idx = mask_idx,
    spatial_dim = spatial_dim,
    n_contrasts = n_contrasts
  )
}

.nifti_subject_ids <- function(files) {
  ids <- tools::file_path_sans_ext(basename(files))
  sub("\\.nii$", "", ids)
}

.nifti_read <- function(handle,
                        assays,
                        block = NULL,
                        mask_idx = NULL,
                        spatial_dim = NULL,
                        n_contrasts = NULL,
                        ...) {
  if (!identical(assays, "beta")) stop("NIfTI adapter currently supports only 'beta' assay", call. = FALSE)
  if (is.null(mask_idx) || is.null(spatial_dim) || is.null(n_contrasts)) {
    stop("Probe metadata must include mask_idx, spatial_dim, n_contrasts", call. = FALSE)
  }

  samples_total <- length(mask_idx)
  sample_idx <- if (!is.null(block) && !is.null(block$sample)) block$sample else seq_len(samples_total)
  n_subjects <- length(handle$files)

  arr <- array(NA_real_, dim = c(length(sample_idx), n_subjects, n_contrasts))

  for (j in seq_along(handle$files)) {
    # Read using neuroim2
    if (n_contrasts == 1L) {
      # 3D file: use read_vol
      img <- neuroim2::read_vol(handle$files[j])
      vec <- .nifti_extract(img, mask_idx, n_contrasts)
    } else {
      # 4D file: use read_vec
      img <- neuroim2::read_vec(handle$files[j])
      vec <- .nifti_extract(img, mask_idx, n_contrasts)
    }
    arr[, j, ] <- vec[sample_idx, , drop = FALSE]
  }

  list(beta = arr)
}

.nifti_extract <- function(img, mask_idx, n_contrasts) {
  # Convert neuroim2 object to array
  # neuroim2::read_vol returns NeuroVol, neuroim2::read_vec returns NeuroVec
  arr <- as.array(img)

  # Get spatial dimensions (first 3 dims)
  dims <- dim(arr)
  spatial_dims <- if (length(dims) == 3L) dims else dims[1:3]
  n_voxels <- prod(spatial_dims)

  # Reshape to [voxels × contrasts] matrix
  vals <- if (n_contrasts == 1L) {
    # 3D volume: flatten to column vector
    matrix(as.vector(arr), nrow = n_voxels, ncol = 1L)
  } else {
    # 4D volume: flatten spatial dims, keep 4th dim as columns
    # Array is [x, y, z, contrasts], we want [voxels, contrasts]
    matrix(as.vector(arr), nrow = n_voxels, ncol = n_contrasts)
  }

  # Extract only the masked voxels
  vals[mask_idx, , drop = FALSE]
}

.nifti_close <- function(handle) {
  invisible(NULL)
}
