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
  } else if (inherits(gds$space, "space_voxel") && !is.null(gds$space$mask_idx)) {
    sample_labels <- gds$space$mask_idx
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

  vol <- RNifti::asNifti(out)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(vol, path)
  invisible(path)
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
