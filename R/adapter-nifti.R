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
#'   provide the sets. See [nifti_source()]. Supplying explicit `subjects`
#'   labels pairs beta/se vectors positionally instead of by filename key.
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
    close = .nifti_close,
    capabilities = list(
      sample_blocks = TRUE,
      subject_blocks = FALSE,
      contrast_blocks = FALSE,
      persistent_handle = TRUE,
      preferred_axis = "sample",
      cheap_revisit = FALSE,
      requires_staging = TRUE
    )
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
    if (!is.null(source$files)) {
      files <- .nifti_normalise_source(source$files)
      files_ok <- length(files) > 0L &&
        all(file.exists(files)) &&
        all(grepl("\\.nii(\\.gz)?$", files, ignore.case = TRUE))
      if (files_ok) return(0.8)
    }
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
  dots <- list(...)
  axis <- .nifti_axis_metadata(source, dots)
  if (is.character(source)) {
    files <- .nifti_normalise_source(source)
    if (!length(files)) stop("No NIfTI files found", call. = FALSE)
    return(.nifti_open_file_vector(files, subjects = axis$subjects, contrasts = axis$contrasts))
  }
  if (is.list(source)) {
    if (!is.null(source$files)) {
      files <- .nifti_normalise_source(source$files)
      if (!length(files)) stop("No NIfTI files found in list", call. = FALSE)
      return(.nifti_open_file_vector(files, subjects = axis$subjects, contrasts = axis$contrasts))
    }
    files_beta <- if (!is.null(source$beta)) .nifti_normalise_source(source$beta) else NULL
    files_se <- if (!is.null(source$se)) .nifti_normalise_source(source$se) else NULL
    if (is.null(files_beta) && is.null(files_se)) stop("No NIfTI files found in list", call. = FALSE)
    out <- .nifti_align_file_sets(files_beta, files_se, subjects = axis$subjects)
    out$contrasts <- axis$contrasts
    return(out)
  }
  stop("Unsupported NIfTI source", call. = FALSE)
}

.nifti_attach_source_metadata <- function(source, dots = list()) {
  axis <- .nifti_axis_metadata(source, dots)
  if (is.null(axis$subjects) && is.null(axis$contrasts)) return(source)

  if (is.character(source)) {
    source <- list(files = source)
  } else if (!is.list(source)) {
    return(source)
  }

  source$subject <- NULL
  source$contrast <- NULL
  if (!is.null(axis$subjects)) source$subjects <- axis$subjects
  if (!is.null(axis$contrasts)) source$contrasts <- axis$contrasts
  source
}

.nifti_axis_metadata <- function(source, dots = list()) {
  list(
    subjects = .nifti_axis_value(source, dots, "subject", "subjects"),
    contrasts = .nifti_axis_value(source, dots, "contrast", "contrasts")
  )
}

.nifti_axis_value <- function(source, dots, singular, plural) {
  values <- list()
  if (is.list(source)) {
    if (!is.null(source[[singular]])) values[[length(values) + 1L]] <- source[[singular]]
    if (!is.null(source[[plural]])) values[[length(values) + 1L]] <- source[[plural]]
  }
  if (!is.null(dots[[singular]])) values[[length(values) + 1L]] <- dots[[singular]]
  if (!is.null(dots[[plural]])) values[[length(values) + 1L]] <- dots[[plural]]
  if (!length(values)) return(NULL)

  values <- lapply(values, as.character)
  ref <- values[[1L]]
  for (val in values[-1L]) {
    if (!identical(ref, val)) {
      stop("Conflicting `", singular, "`/`", plural, "` values supplied.", call. = FALSE)
    }
  }
  ref
}

.nifti_open_file_vector <- function(files, subjects = NULL, contrasts = NULL) {
  bpat <- "beta|cope|effect"
  spat <- "(^|[_.-])(se|stderr|sterr|sigma|std(err)?)([_.-]|$)"
  base <- basename(files)
  beta_files <- files[grepl(bpat, base, ignore.case = TRUE)]
  se_files <- files[grepl(spat, base, ignore.case = TRUE)]
  if (length(beta_files) && length(se_files)) {
    out <- .nifti_align_file_sets(beta_files, se_files, subjects = subjects)
  } else {
    out <- list(files_beta = files, files_se = NULL)
    if (!is.null(subjects)) {
      out$subjects <- .nifti_validate_subjects(subjects, length(files))
    }
  }
  out$contrasts <- contrasts
  out
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

  subjects <- handle$subjects %||% .nifti_subject_key(files1)
  contrasts <- .nifti_contrast_labels(handle$contrasts, n_contrasts)

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
      synthetic_var = is.null(handle$files_se) && !is.null(handle$files_beta),
      sample_labels_synthetic = TRUE
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
  vapply(ids, .nifti_one_subject_key, character(1), USE.NAMES = FALSE)
}

# Derive a subject key from one filename stem. BIDS-named maps (e.g. from
# fmrireg::write_results) carry the statistic in a `desc-` infix and end in a
# suffix such as `_bold`, so the legacy "strip a trailing stat token" rule never
# fires and beta/se files for the same subject get different keys. When a BIDS
# `sub-<label>` entity is present, key on it; otherwise fall back to the legacy
# trailing-token strip (preserving behaviour for non-BIDS inputs).
.nifti_one_subject_key <- function(id) {
  m <- regmatches(id, regexpr("sub-[A-Za-z0-9]+", id, perl = TRUE))
  if (length(m) == 1L && nzchar(m)) {
    return(m)
  }
  stripped <- sub(
    "([_.-](beta|cope|effect|se|stderr|sterr|sigma|std(err)?))+$",
    "",
    id,
    ignore.case = TRUE,
    perl = TRUE
  )
  sub("[_.-]+$", "", stripped, perl = TRUE)
}

# Filename "context" used to confirm a beta file and the se file it was paired
# with (by subject) really are the same map. It is the basename with the file
# extension and the statistic indicator removed -- whether that indicator is a
# BIDS `desc-`/`stat-` entity value (e.g. desc-beta, stat-se) or a legacy
# trailing token (e.g. _beta, .se). Two files that differ ONLY in the statistic
# therefore share a context; files differing in task/run/session/space do not.
.nifti_pair_context <- function(ids) {
  b <- basename(ids)
  b <- sub("\\.nii(\\.gz)?$", "", b, ignore.case = TRUE, perl = TRUE)
  stat <- "(beta|cope|effect|se|stderr|sterr|sigma|std(err)?)"
  b <- gsub(paste0("(?i)[_-](desc|stat)-", stat, "(?=([_.-]|$))"), "", b, perl = TRUE)
  b <- sub(paste0("(?i)([_.-]", stat, ")+$"), "", b, perl = TRUE)
  gsub("[_.-]+$", "", b, perl = TRUE)
}

.nifti_align_file_sets <- function(files_beta, files_se, subjects = NULL) {
  if (is.null(files_beta) || is.null(files_se)) {
    out <- list(files_beta = files_beta, files_se = files_se)
    if (!is.null(subjects)) {
      out$subjects <- .nifti_validate_subjects(subjects, length(files_beta %||% files_se))
    }
    return(out)
  }

  if (!is.null(subjects)) {
    if (length(files_beta) != length(files_se)) {
      stop("NIfTI beta and se file vectors must have the same length when explicit subjects are supplied", call. = FALSE)
    }
    return(list(
      files_beta = files_beta,
      files_se = files_se,
      subjects = .nifti_validate_subjects(subjects, length(files_beta))
    ))
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

  files_beta <- files_beta[ord_beta]
  files_se <- files_se[ord_se]

  # Subject-keyed pairing alone is unsafe: two files sharing `sub-<label>` but
  # differing in task/run/session/space would be paired silently. Require the
  # non-statistic filename context to match for each paired subject; surface a
  # mismatch so the user can fix the inputs or pass explicit `subjects=`.
  ctx_mismatch <- .nifti_pair_context(files_beta) != .nifti_pair_context(files_se)
  if (any(ctx_mismatch)) {
    i <- which(ctx_mismatch)[1L]
    stop(
      "NIfTI beta and se files share subject '", beta_key[i],
      "' but differ in non-statistic filename entities (e.g. task/run/session/space):\n  beta: ",
      basename(files_beta[i]), "\n  se:   ", basename(files_se[i]),
      "\nEnsure beta/se maps correspond, or pass explicit `subject=`/`subjects=` to override pairing.",
      call. = FALSE
    )
  }

  list(files_beta = files_beta, files_se = files_se, subjects = beta_key)
}

.nifti_validate_subjects <- function(subjects, n_subjects) {
  subjects <- as.character(subjects)
  if (length(subjects) != n_subjects) {
    stop("Length of `subjects` must equal number of NIfTI subject files.", call. = FALSE)
  }
  if (anyNA(subjects) || any(!nzchar(subjects))) {
    stop("`subjects` must contain non-missing, non-empty labels.", call. = FALSE)
  }
  if (anyDuplicated(subjects)) {
    stop("`subjects` must not contain duplicates.", call. = FALSE)
  }
  subjects
}

.nifti_contrast_labels <- function(contrasts, n_contrasts) {
  if (is.null(contrasts)) {
    if (n_contrasts == 1L) return("contrast1")
    return(paste0("contrast", seq_len(n_contrasts)))
  }
  contrasts <- as.character(contrasts)
  if (length(contrasts) != n_contrasts) {
    stop("Length of `contrasts` must equal the NIfTI contrast count.", call. = FALSE)
  }
  if (anyNA(contrasts) || any(!nzchar(contrasts))) {
    stop("`contrasts` must contain non-missing, non-empty labels.", call. = FALSE)
  }
  if (anyDuplicated(contrasts)) {
    stop("`contrasts` must not contain duplicates.", call. = FALSE)
  }
  contrasts
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
      .warn_synthetic_variance(once = TRUE)
      out$var <- array(1, dim = c(length(sample_idx), n_subjects, n_contrasts))
      attr(out$var, "synthetic_unit_variance") <- TRUE
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
#' When both `beta` and `se` are supplied, files are paired by subject key unless
#' explicit `subject`/`subjects` labels are supplied. The fallback key is the file
#' basename (extension removed) with a **trailing** statistic token stripped,
#' matching `([_.-](beta|cope|effect|se|stderr|sterr|sigma|std(err)?))+$`
#' (case insensitive). For example `sub-01_beta.nii.gz` and `sub-01_se.nii.gz`
#' both reduce to the key `sub-01`. Filenames where the statistic token is not
#' last (e.g. `sub-01_stat-beta_sm-2.nii.gz`) will fail to pair unless explicit
#' subject labels are provided. For beta/stat-only maps, prefer
#' [gds_from_scalar_maps()], which takes an explicit `subject` vector and avoids
#' filename heuristics entirely.
#'
#' @param beta Character vector of NIfTI paths or a directory containing beta/effect images
#'             (optional, can be `NULL` if only `se` is provided)
#' @param se   Character vector of NIfTI paths or a directory containing standard error images
#'             (optional, can be `NULL` if only `beta` is provided)
#' @param subject,subjects Optional subject labels. When both `beta` and `se`
#'   are supplied, these labels declare that the vectors are already aligned in
#'   matching subject order and bypass filename-key pairing.
#' @param contrast,contrasts Optional contrast labels. Length must match the
#'   number of contrasts in each NIfTI image (one for ordinary 3D maps).
#'
#' @return A list suitable for `gds(source = <returned>, format = "nifti")`
#' @export
#' @examples
#' src <- nifti_source(beta = c("sub-01_beta.nii.gz", "sub-02_beta.nii.gz"),
#'                     se   = c("sub-01_se.nii.gz",   "sub-02_se.nii.gz"))
#' str(src)
nifti_source <- function(beta = NULL,
                         se = NULL,
                         subject = NULL,
                         subjects = NULL,
                         contrast = NULL,
                         contrasts = NULL) {
  if (is.null(beta) && is.null(se)) {
    stop("Provide at least one of `beta` or `se`", call. = FALSE)
  }
  source <- list(beta = beta, se = se)
  .nifti_attach_source_metadata(
    source,
    list(subject = subject, subjects = subjects, contrast = contrast, contrasts = contrasts)
  )
}
