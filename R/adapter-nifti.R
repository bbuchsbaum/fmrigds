#' Register NIfTI Adapter
#'
#' Registers the NIfTI adapter for reading neuroimaging data in NIfTI format.
#' Uses the neuroim2 package for reading NIfTI files.
#'
#' @details
#' The NIfTI adapter supports:
#' - 3D, 4D, and 5D+ NIfTI files (.nii, .nii.gz)
#'   - For 5D+, non-spatial axes (dims 4+) are flattened into the GDS contrast axis
#' - Optional mask specification (separate NIfTI file with 1=included, 0=excluded)
#' - Automatic spatial information extraction (affine transformation)
#' - BIDS-compatible file naming
#'
#' Source specification options:
#' - Character vector or directory: pass one or more `.nii`/`.nii.gz` paths, or a
#'   directory. The adapter will attempt to classify files into beta vs se sets by
#'   filename patterns: `beta|cope|effect` for beta, and `(se|stderr|sterr|sigma|std(err)?)`
#'   occurring with separators like `_`, `.`, or `-`.
#' - Named list: pass `list(beta = <paths/dir>, se = <paths/dir>)` to explicitly
#'   provide the sets. See [nifti_source()].
#'
#' Beta-only NIfTI sources are accepted for raw map workflows. When a realised
#' GDS is materialised from beta maps without standard-error images, the adapter
#' creates a synthetic `var` assay filled with 1 and emits a warning. That
#' variance is only a placeholder; use unweighted reducers such as
#' `method = "ols:voxelwise"` for group analyses of subject-level stat maps
#' rather than fixed/random-effects meta-analysis.
#'
#' To provide a mask, pass it to gds() via the mask parameter:
#' `gds(source, mask = "path/to/mask.nii")`
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
  # Accept character vector (files/dir) or a named list with beta/se entries
  if (is.data.frame(source)) return(FALSE)
  if (is.character(source)) {
    files <- .nifti_normalise_source(source)
    if (!length(files)) return(FALSE)
    if (!all(file.exists(files))) return(FALSE)
    if (!all(grepl("\\.nii(\\.gz)?$", files, ignore.case = TRUE))) return(FALSE)
    return(0.9)
  }
  if (is.list(source)) {
    beta_files <- if (!is.null(source$beta)) .nifti_normalise_source(source$beta) else NULL
    se_files <- if (!is.null(source$se)) .nifti_normalise_source(source$se) else NULL
    beta_ok <- !is.null(beta_files) &&
      length(beta_files) > 0L &&
      all(file.exists(beta_files)) &&
      all(grepl("\\.nii(\\.gz)?$", beta_files, ignore.case = TRUE))
    se_ok <- !is.null(se_files) &&
      length(se_files) > 0L &&
      all(file.exists(se_files)) &&
      all(grepl("\\.nii(\\.gz)?$", se_files, ignore.case = TRUE))
    if (beta_ok || se_ok) return(0.8)
  }
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
  if (is.character(source)) {
    files <- .nifti_normalise_source(source)
    if (!length(files)) stop("No NIfTI files found", call. = FALSE)
    # Try to classify by filename patterns when both beta and se present
    bpat <- "beta|cope|effect"
    spat <- "(^|[_.-])(se|stderr|sterr|sigma|std(err)?)([_.-]|$)"
    base <- basename(files)
    beta_files <- files[grepl(bpat, base, ignore.case = TRUE)]
    se_files <- files[grepl(spat, base, ignore.case = TRUE)]
    if (length(beta_files) && length(se_files)) {
      return(.nifti_align_file_sets(beta_files, se_files))
    }
    # Fallback: treat all as beta if unclassified
    return(list(files_beta = files, files_se = NULL))
  }
  if (is.list(source)) {
    files_beta <- if (!is.null(source$beta)) .nifti_normalise_source(source$beta) else NULL
    files_se <- if (!is.null(source$se)) .nifti_normalise_source(source$se) else NULL
    if (is.null(files_beta) && is.null(files_se)) stop("No NIfTI files found in list", call. = FALSE)
    return(.nifti_align_file_sets(files_beta, files_se))
  }
  stop("Unsupported NIfTI source", call. = FALSE)
}

.nifti_probe <- function(handle, ...) {
  # Read header using neuroim2
  files1 <- handle$files_beta %||% handle$files_se
  meta <- neuroim2::read_header(files1[1])

  # Get dimensions from metadata
  dim_all <- as.integer(meta@dims)
  if (length(dim_all) < 3L) {
    stop("NIfTI dimensionality must be at least 3, got ", length(dim_all), call. = FALSE)
  }
  spatial_dim <- dim_all[1:3]
  n_contrasts <- .nifti_contrast_count(dim_all)

  # Read image object (NeuroVol/NeuroVec/NeuroHyperVec depending on source ndim)
  first_img <- .nifti_read_image(files1[1], dim_all = dim_all)

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

  subjects <- .nifti_subject_key(files1)
  contrasts <- if (n_contrasts == 1L) {
    "contrast1"
  } else {
    paste0("contrast", seq_len(n_contrasts))
  }

  # Extract affine transformation from neuroim2 NeuroSpace
  nspace <- neuroim2::space(first_img)
  affine <- neuroim2::trans(nspace)

  space <- space_voxel(
    dim = spatial_dim,
    affine = affine,
    mask_bitmap = mask_bitmap,
    storage = "packed"
  )

  assays_avail <- c()
  if (!is.null(handle$files_beta)) assays_avail <- c(assays_avail, "beta")
  if (!is.null(handle$files_se)) assays_avail <- c(assays_avail, "se")
  if (!is.null(handle$files_beta) && is.null(handle$files_se)) assays_avail <- c(assays_avail, "var")
  out <- list(
    assays = assays_avail,
    dims = gds_dims(sample = length(mask_idx), subject = length(files1), contrast = n_contrasts),
    subjects = subjects,
    contrasts = contrasts,
    space = space,
    maps = list(),
    metadata = list(
      schema_version = "0.1.0",
      source_files = c(handle$files_beta %||% character(), handle$files_se %||% character()),
      # Beta-only sources have no real uncertainty; the `var` assay is a
      # synthetic unit-variance placeholder. Flagged so variance-weighted
      # reducers can refuse it (see reduce()).
      synthetic_var = is.null(handle$files_se) && !is.null(handle$files_beta)
    ),
    columns = list(effect_cols = NULL, subject_col = NULL, sample_col = NULL, contrast_col = NULL),
    mask_idx = mask_idx,
    spatial_dim = spatial_dim,
    n_contrasts = n_contrasts,
    source_ndim = length(dim_all)
  )
  probe_contract(out)
}

.nifti_subject_ids <- function(files) {
  ids <- tools::file_path_sans_ext(basename(files))
  sub("\\.nii$", "", ids)
}

.nifti_subject_key <- function(files) {
  ids <- .nifti_subject_ids(files)
  ids <- sub(
    "([_.-](beta|cope|effect|se|stderr|sterr|sigma|std(err)?))+$",
    "",
    ids,
    ignore.case = TRUE,
    perl = TRUE
  )
  sub("[_.-]+$", "", ids, perl = TRUE)
}

.nifti_align_file_sets <- function(files_beta, files_se) {
  if (is.null(files_beta) || is.null(files_se)) {
    return(list(files_beta = files_beta, files_se = files_se))
  }

  beta_key <- .nifti_subject_key(files_beta)
  se_key <- .nifti_subject_key(files_se)
  if (anyDuplicated(beta_key) || anyDuplicated(se_key)) {
    stop("NIfTI beta/se files must not contain duplicate subject identifiers", call. = FALSE)
  }

  ord_beta <- order(beta_key)
  ord_se <- order(se_key)
  beta_key <- beta_key[ord_beta]
  se_key <- se_key[ord_se]
  if (!identical(beta_key, se_key)) {
    stop("NIfTI beta and se files must refer to the same subjects", call. = FALSE)
  }

  list(files_beta = files_beta[ord_beta], files_se = files_se[ord_se])
}

.nifti_read <- function(handle,
                        assays,
                        block = NULL,
                        mask_idx = NULL,
                        spatial_dim = NULL,
                        n_contrasts = NULL,
                        ...) {
  # Support reading "beta" and/or "se"
  if (is.null(mask_idx) || is.null(spatial_dim) || is.null(n_contrasts)) {
    stop("Probe metadata must include mask_idx, spatial_dim, n_contrasts", call. = FALSE)
  }

  samples_total <- length(mask_idx)
  sample_idx <- if (!is.null(block) && !is.null(block$sample)) block$sample else seq_len(samples_total)
  n_subjects <- length(handle$files_beta %||% handle$files_se)
  out <- list()
  read_stack <- function(file_vec, assay_name) {
    if (length(file_vec) != n_subjects) {
      stop("NIfTI assay file count mismatch for ", assay_name, call. = FALSE)
    }
    arr <- array(NA_real_, dim = c(length(sample_idx), n_subjects, n_contrasts))
    for (j in seq_along(file_vec)) {
      img <- .nifti_read_image(file_vec[j])
      vec <- .nifti_extract(img, mask_idx, n_contrasts)
      arr[, j, ] <- vec[sample_idx, , drop = FALSE]
    }
    arr
  }
  for (nm in assays) {
    if (identical(nm, "beta") && !is.null(handle$files_beta)) out$beta <- read_stack(handle$files_beta, "beta")
    if (identical(nm, "se") && !is.null(handle$files_se)) out$se <- read_stack(handle$files_se, "se")
    if (identical(nm, "var") && !is.null(handle$files_beta) && is.null(handle$files_se)) {
      warning("No variance or SE provided; using unit variance", call. = FALSE)
      out$var <- array(1, dim = c(length(sample_idx), n_subjects, n_contrasts))
    }
  }
  out
}

.nifti_extract <- function(img, mask_idx, n_contrasts) {
  # Prefer direct dense conversion; some NeuroHyperVec variants do not implement
  # base as.array() and require explicit voxel-slice reconstruction.
  arr <- tryCatch(as.array(img), error = function(e) NULL)
  vals <- if (!is.null(arr)) {
    dims <- dim(arr)
    if (length(dims) < 3L) {
      stop("NIfTI image must have at least 3 dimensions", call. = FALSE)
    }
    n_voxels <- prod(dims[1:3])
    contrast_count <- .nifti_contrast_count(dims)
    if (!identical(as.integer(contrast_count), as.integer(n_contrasts))) {
      stop(
        "Inconsistent contrast count for NIfTI image. Expected ", n_contrasts,
        " but found ", contrast_count, ".", call. = FALSE
      )
    }
    matrix(as.vector(arr), nrow = n_voxels, ncol = n_contrasts)
  } else if (inherits(img, "NeuroHyperVec")) {
    .nifti_hyper_to_matrix(img, n_contrasts)
  } else {
    stop("Unable to coerce NIfTI image object to an array", call. = FALSE)
  }

  vals[mask_idx, , drop = FALSE]
}

.nifti_contrast_count <- function(dims) {
  dims <- as.integer(dims)
  if (length(dims) <= 3L) return(1L)
  as.integer(prod(dims[-(1:3)]))
}

.nifti_read_image <- function(file, dim_all = NULL) {
  dim_all <- if (is.null(dim_all)) as.integer(neuroim2::read_header(file)@dims) else as.integer(dim_all)
  ndim <- length(dim_all)

  if (ndim == 3L) {
    return(neuroim2::read_vol(file))
  }
  if (ndim == 4L) {
    return(neuroim2::read_vec(file))
  }
  if (ndim >= 5L) {
    ns <- asNamespace("neuroim2")
    if (exists("read_hyper_vec", envir = ns, mode = "function", inherits = FALSE)) {
      return(get("read_hyper_vec", envir = ns, inherits = FALSE)(file))
    }
    stop(
      "Reading 5D+ NIfTI requires neuroim2::read_hyper_vec(). ",
      "Please update the neuroim2 package.", call. = FALSE
    )
  }

  stop("Unsupported NIfTI dimensionality: ", ndim, call. = FALSE)
}

.nifti_hyper_to_matrix <- function(img, n_contrasts) {
  dims <- as.integer(dim(img))
  if (length(dims) < 4L) {
    stop("NeuroHyperVec must have at least 4 dimensions", call. = FALSE)
  }

  n_voxels <- prod(dims[1:3])
  contrast_dims <- dims[-(1:3)]
  contrast_count <- as.integer(prod(contrast_dims))
  if (!identical(as.integer(contrast_count), as.integer(n_contrasts))) {
    stop(
      "Inconsistent contrast count for NeuroHyperVec. Expected ", n_contrasts,
      " but found ", contrast_count, ".", call. = FALSE
    )
  }

  idx_grid <- expand.grid(lapply(contrast_dims, seq_len), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vals <- matrix(NA_real_, nrow = n_voxels, ncol = nrow(idx_grid))
  for (i in seq_len(nrow(idx_grid))) {
    tail_idx <- lapply(idx_grid[i, , drop = TRUE], as.integer)
    # Programmatic equivalent of img[,,,i4,i5,..., drop = TRUE]
    idx_txt <- paste(unlist(tail_idx, use.names = FALSE), collapse = ",")
    slice_expr <- parse(text = paste0("img[,,,", idx_txt, ", drop = TRUE]"))[[1]]
    slice <- eval(slice_expr, envir = list(img = img), enclos = baseenv())
    vals[, i] <- as.vector(slice)
  }
  vals
}

.nifti_close <- function(handle) {
  invisible(NULL)
}

#' Construct a NIfTI source specification
#'
#' Build an explicit NIfTI source list for use with `gds(..., format = "nifti")`.
#' If `beta` is provided without `se`, materialising the GDS creates a unit
#' variance placeholder with a warning. This supports raw beta/stat-map
#' workflows; the synthetic variance should not be interpreted as subject-level
#' uncertainty, and variance-weighted reducers will refuse such a GDS.
#'
#' When both `beta` and `se` are supplied, files are paired by subject key: the
#' file basename (extension removed) with a **trailing** statistic token
#' stripped, matching `([_.-](beta|cope|effect|se|stderr|sterr|sigma|std(err)?))+$`
#' (case insensitive). For example `sub-01_beta.nii.gz` and `sub-01_se.nii.gz`
#' both reduce to the key `sub-01`. Filenames where the statistic token is not
#' last (e.g. `sub-01_stat-beta_sm-2.nii.gz`) will fail to pair; rename so the
#' token is trailing, or pass `beta`/`se` as already-aligned vectors in matching
#' subject order. For beta/stat-only maps, prefer [gds_from_scalar_maps()], which
#' takes an explicit `subject` vector and avoids filename heuristics entirely.
#'
#' @param beta Character vector of NIfTI paths or a directory containing beta/effect images
#'             (optional, can be `NULL` if only `se` is provided)
#' @param se   Character vector of NIfTI paths or a directory containing standard error images
#'             (optional, can be `NULL` if only `beta` is provided)
#'
#' @return A list suitable for `gds(source = <returned>, format = "nifti")`
#' @export
#' @examples
#' src <- nifti_source(beta = c("sub-01_beta.nii.gz", "sub-02_beta.nii.gz"),
#'                     se   = c("sub-01_se.nii.gz",   "sub-02_se.nii.gz"))
#' str(src)
nifti_source <- function(beta = NULL, se = NULL) {
  if (is.null(beta) && is.null(se)) {
    stop("Provide at least one of `beta` or `se`", call. = FALSE)
  }
  list(beta = beta, se = se)
}
