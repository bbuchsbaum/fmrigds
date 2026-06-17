.gds_tidy_frame <- function(gds, stats, drop_na = FALSE, include_col_data = FALSE) {
  stats <- intersect(stats, names(assays(gds)))
  if (!length(stats)) {
    stop("No matching assays found for export.", call. = FALSE)
  }

  dims <- dim(assays(gds)[[stats[1]]])
  samples <- seq_len(dims[1])
  subjects <- gds$subjects
  contrasts <- gds$contrasts

  sample_labels <- samples
  if (inherits(gds$space, "space_parcels") && !is.null(gds$space$labels)) {
    sample_labels <- gds$space$labels
  } else if (inherits(gds$space, "space_sample_labels") && !is.null(gds$space$labels)) {
    sample_labels <- gds$space$labels
  } else if (inherits(gds$space, "space_voxel") && !is.null(gds$space$mask_idx)) {
    sample_labels <- gds$space$mask_idx
  }
  if (length(sample_labels) != length(samples)) {
    if (length(sample_labels) >= length(samples)) {
      sample_labels <- sample_labels[seq_len(length(samples))]
    } else {
      sample_labels <- samples
    }
  }

  grid <- expand.grid(
    sample = sample_labels,
    subject = subjects,
    contrast = contrasts,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  for (name in stats) {
    grid[[name]] <- as.vector(assay(gds, name))
  }

  if (include_col_data && !is.null(gds$col_data) && nrow(gds$col_data)) {
    cd <- gds$col_data
    cd_subjects <- rownames(cd) %||% gds$subjects
    for (nm in colnames(cd)) {
      column <- rep(NA, nrow(grid))
      idx <- match(grid$subject, cd_subjects)
      column[!is.na(idx)] <- unlist(cd[idx[!is.na(idx)], nm], use.names = FALSE)
      grid[[nm]] <- column
    }
  }

  if (drop_na) {
    keep <- rowSums(!is.na(grid[, stats, drop = FALSE])) > 0
    grid <- grid[keep, , drop = FALSE]
  }
  grid
}

.write_gds_csv <- function(gds, path, options = list()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  stats <- options$stats %||% names(assays(gds))
  drop_na <- isTRUE(options$drop_na)
  include_col_data <- isTRUE(options$include_col_data)
  df <- .gds_tidy_frame(gds, stats, drop_na = drop_na, include_col_data = include_col_data)
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

.write_gds_parquet <- function(gds, path, options = list()) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("The 'arrow' package is required for Parquet export.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  stats <- options$stats %||% names(assays(gds))
  drop_na <- isTRUE(options$drop_na)
  include_col_data <- isTRUE(options$include_col_data)
  df <- .gds_tidy_frame(gds, stats, drop_na = drop_na, include_col_data = include_col_data)
  arrow::write_parquet(df, sink = path)
  invisible(path)
}

# Write a voxel-space numeric array to NIfTI, preserving the spatial affine.
#
# The GDS voxel `affine` originates from neuroim2 (`neuroim2::trans()`), so for
# 3D images we round-trip through neuroim2: this reproduces the affine exactly
# and writes a valid `scl_slope = 1`. For higher-dimensional arrays we fall back
# to RNifti, stamping the affine into the sform/qform and pixdim so the geometry
# is still preserved (RNifti writes `scl_slope = 0`, which is spec-valid "no
# scaling").
.write_voxel_nifti <- function(arr, affine, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  d <- dim(arr) %||% length(arr)
  # Drop singleton non-spatial axes (e.g. a collapsed subject axis after reduce)
  if (length(d) > 3L) {
    keep <- c(rep(TRUE, 3L), d[-(1:3)] != 1L)
    arr <- array(as.numeric(arr), dim = d[keep])
    d <- dim(arr)
  }
  affine_ok <- is.matrix(affine) && all(dim(affine) == 4L)
  if (affine_ok && length(d) == 3L && requireNamespace("neuroim2", quietly = TRUE)) {
    sp <- neuroim2::NeuroSpace(dim = as.integer(d), trans = affine)
    vol <- neuroim2::NeuroVol(array(as.numeric(arr), dim = as.integer(d)), sp)
    neuroim2::write_vol(vol, path)
    return(invisible(path))
  }
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("The 'RNifti' or 'neuroim2' package is required for NIfTI export.", call. = FALSE)
  }
  RNifti::writeNifti(.stamp_nifti_xform(arr, if (affine_ok) affine else NULL), path)
  invisible(path)
}

# Stamp a 4x4 voxel->world affine onto an RNifti image (fallback path for
# non-3D arrays). Sets pixdim from the affine column norms and both sform and
# qform with a generic ALIGNED_ANAT code so viewers honour the geometry.
#
# Two subtleties make the naive version silently wrong for >3D outputs:
#   1. `RNifti::pixdim(vol) <- vox` does NOT persist on an image built from a
#      bare R array, so the qform (quaternion + pixdim) ends up encoding unit
#      scale even though the sform is correct. RNifti's own xform() prefers the
#      qform, so the affine reads back as identity-scaled. We instead stamp
#      pixdim as an attribute on the array *before* asNifti(), where it sticks.
#   2. The qform can only represent an orthogonal rotation * scale (with an
#      optional handedness flip). For a sheared / non-orthogonal affine the
#      quaternion encoding is lossy but assignment still "succeeds", so the
#      tryCatch never fires. We only stamp the qform when the direction cosines
#      are orthogonal; otherwise we leave it unset (code 0) and let readers fall
#      back to the exact sform.
.stamp_nifti_xform <- function(arr, affine) {
  if (!(is.matrix(affine) && all(dim(affine) == 4L))) {
    return(RNifti::asNifti(arr))
  }
  rot <- affine[1:3, 1:3, drop = FALSE]
  vox <- sqrt(colSums(rot^2))
  vox[!is.finite(vox) | vox == 0] <- 1
  d <- dim(arr) %||% length(arr)
  attr(arr, "pixdim") <- c(vox, rep(1, max(0L, length(d) - 3L)))
  vol <- RNifti::asNifti(arr)
  xf <- structure(affine, code = 2L)
  RNifti::sform(vol) <- xf
  cosines <- sweep(rot, 2L, vox, "/")
  if (max(abs(crossprod(cosines) - diag(3))) < 1e-4) {
    tryCatch({ RNifti::qform(vol) <- xf }, error = function(e) NULL)
  }
  vol
}

.write_gds_nifti <- function(gds, path, options = list()) {
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("The 'RNifti' package is required for NIfTI export.", call. = FALSE)
  }
  space <- gds$space
  if (!inherits(space, "space_voxel")) {
    stop("NIfTI export requires a voxel space.", call. = FALSE)
  }
  stat <- options$stat %||% if ("beta" %in% names(assays(gds))) "beta" else names(assays(gds))[1]
  arr <- assay(gds, stat)
  dims <- dim(arr)
  vox_dim <- space$dim
  n_vox <- prod(vox_dim)
  if (dims[1] != length(space$mask_idx %||% seq_len(n_vox))) {
    stop("Sample axis does not match voxel space dimensions.", call. = FALSE)
  }

  out <- array(0, dim = c(vox_dim, dims[2], dims[3]))
  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      vec <- arr[, j, k]
      if (identical(space$storage, "packed") && !is.null(space$mask_idx)) {
        full <- numeric(n_vox)
        full[space$mask_idx] <- vec
        vec <- full
      }
      out[,,, j, k] <- array(vec, dim = vox_dim)
    }
  }

  .write_voxel_nifti(out, space$affine, path)
}

.write_export <- function(gds, fmt, path, options) {
  switch(fmt,
    h5 = write_gds_h5(gds, path, options),
    csv = .write_gds_csv(gds, path, options),
    parquet = .write_gds_parquet(gds, path, options),
    nifti = .write_gds_nifti(gds, path, options),
    stop("Unsupported write_out format: ", fmt, call. = FALSE)
  )
}
