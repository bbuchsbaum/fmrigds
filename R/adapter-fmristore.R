# fmristore adapter (Path A + Path B preference) ------------------------------ # nocov start

# Detection priority scores
.FMRI_SCORE_GDS_NATIVE <- 1.0   # Perfect match: /gds group
.FMRI_SCORE_FMRS_LEGACY <- 0.95 # Legacy fmristore layouts (labeled, latent)

#' Register fmristore Adapter
#'
#' Registers the fmristore adapter for reading fmristore HDF5 layouts.
#' Supports LabeledVolumeSet, LatentNeuroVec, and native /gds layouts.
#'
#' Details:
#' - Path B: Native /gds group (delegates to h5 adapter)
#' - Path A: LabeledVolumeSet layout (/header, /mask, /labels, /data)
#' - Path A: LatentNeuroVec layout (/basis, /scans)
#' - Multi-file subject composition via vector of file paths
#'
#' @return NULL (invisibly)
#' @keywords internal
register_fmristore_adapter <- function() {
  register_adapter(
    name = "fmristore",
    detect = .fmri_detect,
    open   = .fmri_open,
    probe  = .fmri_probe,
    read   = .fmri_read,
    close  = .fmri_close
  )
}

.is_h5_file <- function(x) {
  is.character(x) && length(x) == 1L && tolower(tools::file_ext(x)) %in% c("h5", "hdf5")
}

.fmri_detect <- function(source) {
  # Support single file or vector of files
  paths <- if (is.character(source)) source else return(FALSE)
  if (!length(paths)) return(FALSE)
  if (!all(file.exists(paths))) return(FALSE)

  # Quick check first file
  p <- paths[[1]]
  if (!.is_h5_file(p)) return(FALSE)
  if (!requireNamespace("hdf5r", quietly = TRUE)) return(FALSE)
  h5 <- hdf5r::H5File$new(p, mode = "r")
  on.exit(h5$close(), add = TRUE)

  # Path B: native /gds
  .exists <- function(h, p) tryCatch({ obj <- h$open(p); obj$close(); TRUE }, error = function(e) FALSE)
  if (.exists(h5, "/gds")) return(.FMRI_SCORE_GDS_NATIVE)
  # Path A: fmristore-like labeled volume layout
  if (.exists(h5, "/header") && .exists(h5, "/mask") && .exists(h5, "/labels") && .exists(h5, "/data")) {
    return(.FMRI_SCORE_FMRS_LEGACY)
  }
  # Path A: Parcellated layout (cluster map + scans/data)
  if (.exists(h5, "/cluster_map") && (.exists(h5, "/scans") || .exists(h5, "/data"))) {
    return(.FMRI_SCORE_FMRS_LEGACY)
  }
  # Path A: latent layout
  if (.exists(h5, "/basis/basis_matrix") && .exists(h5, "/scans")) {
    return(.FMRI_SCORE_FMRS_LEGACY)
  }
  FALSE
}

.fmri_open <- function(source, mode = "r", ...) {
  stopifnot(identical(mode, "r"))
  paths <- if (is.character(source)) source else stop("Unsupported source for fmristore adapter", call. = FALSE)
  pre_read_stats <- lapply(paths, .source_stat)
  # Open first file for probe; keep all paths for multi-subject composition
  h5 <- hdf5r::H5File$new(paths[[1]], mode = "r")
  list(paths = paths, h5 = h5, pre_read_stats = pre_read_stats)
}

.fmri_close <- function(handle) {
  if (!is.null(handle$h5)) handle$h5$close_all()
  invisible(NULL)
}

.read_string_vec <- function(h5, path) {
  ds <- h5$open(path)
  on.exit(ds$close(), add = TRUE)
  as.character(ds$read())
}

.read_uint8_array <- function(h5, path) {
  ds <- h5$open(path)
  on.exit(ds$close(), add = TRUE)
  as.array(ds$read())
}

.h5_read_numvec <- function(h5, path) {
  ds <- h5$open(path)
  on.exit(ds$close(), add = TRUE)
  as.numeric(ds$read())
}

.sanitize_label <- function(label) {
  # Match fmristore's label sanitization: replace non-alphanumeric with underscore
  gsub("[^A-Za-z0-9_.]", "_", as.character(label))
}

## h5 helpers come from R/h5-utils.R: .h5_safe_exists(), .h5_first_dataset_in_group()

.fmri_probe <- function(handle, temporal_policy = c("as_is","mean","design"), contrast_matrix = NULL, contrast_names = NULL, source_identity = "sha256", ...) {
  h5 <- handle$h5
  source_identity <- .normalize_source_identity_policy(source_identity)
  source_records <- lapply(seq_along(handle$paths), function(i) list(
    path = handle$paths[[i]],
    role = "fmristore-container",
    ordinal = i,
    pair = i,
    kind = "file",
    expected_stat = handle$pre_read_stats[[i]]
  ))
  source_metadata <- function(extra = list()) {
    utils::modifyList(
      list(schema_version = "0.2.0", source_files = handle$paths),
      extra
    )
  }
  finish_probe <- function(out) {
    out$source_identity_records <- source_records
    out$source_identity_policy <- source_identity
    probe_contract(out)
  }
  # Delegate to h5 adapter if /gds present
  if (.h5_safe_exists(h5, "/gds")) {
    # Use existing h5 adapter for probe
    h5_close <- FALSE
    adapter <- get_adapter("h5")
    # Build a pseudo handle compatible with h5 adapter
    out <- adapter$probe(
      list(
        file = h5,
        path = handle$paths[[1]],
        pre_read_stat = handle$pre_read_stats[[1L]]
      ),
      source_identity = source_identity
    )
    return(probe_contract(out))
  }

  # Path A: Parcellated (cluster/ROI) layout
  if (.h5_safe_exists(h5, "/cluster_map") && (.h5_safe_exists(h5, "/scans") || .h5_safe_exists(h5, "/data"))) {
    temporal_policy <- match.arg(temporal_policy)
    # Read cluster map and compute parcel labels
    cmap_ds <- h5$open("/cluster_map"); on.exit(cmap_ds$close(), add = TRUE)
    cluster_map <- as.array(cmap_ds$read())
    cluster_ids <- sort(unique(as.integer(cluster_map[cluster_map > 0])))
    n_parcels <- length(cluster_ids)
    parcel_labels <- as.character(cluster_ids)
    # Optional metadata labels
    if (.h5_safe_exists(h5, "/cluster_metadata/names")) {
      nm <- .read_string_vec(h5, "/cluster_metadata/names")
      if (length(nm) == n_parcels) parcel_labels <- nm
    }
    # Subjects
    subjects <- NULL
    if (.h5_safe_exists(h5, "/scans")) {
      scans_grp <- h5$open("/scans"); on.exit(scans_grp$close(), add = TRUE)
      scl <- scans_grp$ls(); subjects <- if (is.data.frame(scl)) as.character(scl$name) else as.character(scl)
    } else {
      subjects <- if (length(handle$paths) > 1) basename(handle$paths) else basename(handle$paths[[1]])
    }
    # Determine time dimension from first subject or /data
    get_first_mat_dim <- function() {
      if (.h5_safe_exists(h5, "/scans")) {
        for (u in subjects) {
          grp_path <- paste0("/scans/", u)
          ds_path <- .h5_first_dataset_in_group(h5, grp_path)
          if (is.na(ds_path)) next
          ds <- h5$open(ds_path); on.exit(ds$close(), add = TRUE)
          dm_try <- tryCatch(as.integer(ds$dims), error = function(e) NULL)
          if (!is.null(dm_try) && length(dm_try) == 2L) return(dm_try)
        }
      }
      if (.h5_safe_exists(h5, "/data")) {
        ds <- h5$open("/data"); on.exit(ds$close(), add = TRUE)
        return(as.integer(ds$dims))
      }
      stop("Parcellated: no 2D cluster x time matrix found under /scans or /data", call. = FALSE)
    }
    dm <- tryCatch(get_first_mat_dim(), error = function(e) NULL)
    # Infer orientation: prefer rows=clusters; if not found, fallback to 1 contrast
    if (is.null(dm)) {
      time_len <- 1L
    } else if (length(dm) == 2L && dm[1] == n_parcels) {
      time_len <- dm[2]
    } else if (length(dm) == 2L && dm[2] == n_parcels) {
      time_len <- dm[1]
    } else {
      time_len <- 1L
    }
    # Temporal policy determines contrasts
    n_contrasts <- switch(temporal_policy,
      mean = 1L,
      design = if (!is.null(contrast_matrix)) ncol(contrast_matrix) else time_len,
      as_is = time_len
    )
    # Contrast names
    contrasts <- NULL
    if (temporal_policy == "design" && !is.null(contrast_names)) {
      contrasts <- as.character(contrast_names)
    } else if (.h5_safe_exists(h5, "/labels")) {
      # reuse labels if present to name time axis
      lbl <- .read_string_vec(h5, "/labels")
      contrasts <- if (length(lbl) >= n_contrasts) lbl[seq_len(n_contrasts)] else lbl
    }
    if (is.null(contrasts)) contrasts <- paste0("t", seq_len(n_contrasts))

    # Build space
    # Also compute membership map vox->parcels
    mask_bitmap <- cluster_map > 0
    mask_idx <- which(as.vector(mask_bitmap))
    # Build sparse membership matrix (rows = parcels, cols = masked voxels)
    # Compute cluster IDs in mask order
    cl_mask <- as.integer(as.vector(cluster_map)[mask_idx])
    # Map parcel label cluster_ids to positions
    id_to_row <- match(cl_mask, cluster_ids)
    # Prepare sparse i (row), j (col) with 1s, later normalize by row counts
    j <- seq_along(mask_idx)
    i <- id_to_row
    keep <- which(!is.na(i))
    i <- i[keep]; j <- j[keep]
    M <- Matrix::sparseMatrix(i = i, j = j, x = 1.0, dims = c(n_parcels, length(mask_idx)))
    rs <- Matrix::rowSums(M)
    rs[rs == 0] <- 1
    M <- Matrix::Diagonal(x = 1/as.numeric(rs)) %*% M

    space_par <- space_parcels(labels = parcel_labels)
    # Optional voxel ref space
    # Use cluster_map dimensions to ensure consistency with mask_bitmap
    vox_dim <- as.integer(dim(cluster_map))
    if (length(vox_dim) < 3L) {
      stop("cluster_map must be at least 3D to define voxel space", call. = FALSE)
    }
    if (length(vox_dim) > 3L) vox_dim <- vox_dim[1:3]
    affine <- diag(4)
    space_vox <- space_voxel(dim = vox_dim, affine = affine, mask_bitmap = mask_bitmap, mask_idx = mask_idx, storage = "packed")

    maps <- list(
      vox_to_parcels_mean = list(type = "linear", from = space_vox, to = space_par, operator = M, uncertainty = list(mode = "independent"))
    )

    # Detect assays: assume beta present; try to detect var/se groups
    assays <- c("beta")
    if (.h5_safe_exists(h5, "/scans_var") || .h5_safe_exists(h5, "/var")) assays <- c(assays, "var")
    if (.h5_safe_exists(h5, "/scans_se")  || .h5_safe_exists(h5, "/se"))  assays <- c(assays, "se")

    out <- list(
      assays = unique(assays),
      dims = gds_dims(sample = n_parcels, subject = length(subjects), contrast = n_contrasts),
      subjects = as.character(subjects),
      contrasts = as.character(contrasts),
      space = space_par,
      maps = maps,
      metadata = source_metadata(list(
        cluster_ids = cluster_ids,
        temporal_policy = temporal_policy
      )),
      columns = list()
    )
    return(finish_probe(out))
  }

  # Path A: LatentNeuroVec-like layout
  if (.h5_safe_exists(h5, "/basis/basis_matrix") && .h5_safe_exists(h5, "/scans")) {
    # Read basis matrix and determine k and V
    bds <- h5$open("/basis/basis_matrix")
    on.exit(bds$close(), add = TRUE)
    bmat <- bds$read()
    dm <- dim(bmat)
    # Expect [k, V]; if [V, k], transpose
    if (length(dm) != 2L) stop("/basis/basis_matrix must be 2D", call. = FALSE)
    if (dm[1] < dm[2]) {
      k <- as.integer(dm[1])
      projector <- bmat
    } else {
      # transpose to [k, V]
      projector <- t(bmat)
      k <- nrow(projector)
    }
    basis_name <- tryCatch(as.character(h5$open("/basis/basis_method")$read()), error = function(e) NULL)
    # Optional voxel reference
    voxel_space <- NULL
    if (.h5_safe_exists(h5, "/mask") || .h5_safe_exists(h5, "/header/dim")) {
      mask_bitmap <- if (.h5_safe_exists(h5, "/mask")) (.read_uint8_array(h5, "/mask") > 0) else NULL
      dim_hdr <- tryCatch(as.integer(h5$open("/header/dim")$read()), error = function(e) NULL)
      vox_dim <- if (!is.null(dim_hdr) && length(dim_hdr) >= 4) as.integer(dim_hdr[2:4]) else if (!is.null(mask_bitmap)) dim(mask_bitmap) else NULL
      if (!is.null(vox_dim)) {
        affine <- diag(4)
        voxel_space <- space_voxel(dim = vox_dim, affine = affine, mask_bitmap = mask_bitmap, storage = if (!is.null(mask_bitmap)) "packed" else "dense")
      }
    }
    sp <- space_basis(k = k, basis_name = basis_name, projector = projector, voxel_space = voxel_space)
    # Subjects from /scans children
    scans_grp <- h5$open("/scans")
    on.exit(scans_grp$close(), add = TRUE)
    scl <- scans_grp$ls()
    subj <- if (is.data.frame(scl)) scl$name else scl
    # Determine contrasts by reading first embedding shape
    emb_path <- paste0("/scans/", subj[[1]], "/embedding")
    e1 <- h5$open(emb_path)
    on.exit(e1$close(), add = TRUE)
    em <- e1$read()
    ed <- dim(em)
    n_contrasts <- if (is.null(ed)) {
      1L
    } else if (length(ed) == 1L) {
      1L
    } else if (ed[1] == k) {
      ed[2]
    } else {
      ed[1]
    }
    # Construct contrasts names
    contrasts <- paste0("c", seq_len(n_contrasts))
    out <- list(
      assays = c("beta"),
      dims = gds_dims(sample = k, subject = length(subj), contrast = length(contrasts)),
      subjects = as.character(subj),
      contrasts = contrasts,
      space = sp,
      maps = list(),
      metadata = source_metadata(),
      columns = list()
    )
    return(finish_probe(out))
  }

  # Path A: LabeledVolumeSet-like layout
  labels <- .read_string_vec(h5, "/labels")
  mask <- .read_uint8_array(h5, "/mask")
  if (length(dim(mask)) != 3L) stop("/mask must be a 3D array", call. = FALSE)
  mask_bitmap <- mask > 0
  mask_idx <- which(as.vector(mask_bitmap))

  # Header affine/dims (best-effort)
  dim_hdr <- tryCatch(as.integer(h5$open("/header/dim")$read()), error = function(e) NULL)
  # dim may be c(nDim, X, Y, Z, T, ...). Use the first 3 spatial dims when available
  vox_dim <- if (!is.null(dim_hdr) && length(dim_hdr) >= 4) as.integer(dim_hdr[2:4]) else as.integer(dim(mask))
  affine <- diag(4)
  # Try srow_x/y/z if present
  if (.h5_safe_exists(h5, "/header/srow_x") && .h5_safe_exists(h5, "/header/srow_y") && .h5_safe_exists(h5, "/header/srow_z")) {
    ax <- .h5_read_numvec(h5, "/header/srow_x")
    ay <- .h5_read_numvec(h5, "/header/srow_y")
    az <- .h5_read_numvec(h5, "/header/srow_z")
    if (length(ax) == 4L && length(ay) == 4L && length(az) == 4L) {
      affine[1, ] <- ax
      affine[2, ] <- ay
      affine[3, ] <- az
    } else {
      warning("srow vectors have unexpected length; using identity affine", call. = FALSE)
    }
  } else {
    # Try quaternion if available and neuroim2 present
    if (.h5_safe_exists(h5, "/header/quatern_b") && .h5_safe_exists(h5, "/header/qoffset_x") &&
        requireNamespace("neuroim2", quietly = TRUE)) {
      qb <- .h5_read_numvec(h5, "/header/quatern_b")
      qc <- .h5_read_numvec(h5, "/header/quatern_c")
      qd <- .h5_read_numvec(h5, "/header/quatern_d")
      qx <- .h5_read_numvec(h5, "/header/qoffset_x")
      qy <- .h5_read_numvec(h5, "/header/qoffset_y")
      qz <- .h5_read_numvec(h5, "/header/qoffset_z")
      pixdim <- tryCatch(.h5_read_numvec(h5, "/header/pixdim"), error = function(e) NULL)
      qfac <- if (.h5_safe_exists(h5, "/header/qfac")) .h5_read_numvec(h5, "/header/qfac") else 1.0
      if (!is.null(pixdim) && length(pixdim) >= 4L) {
        affine <- neuroim2::quaternToMatrix(
          quat = c(qb, qc, qd),
          origin = c(qx, qy, qz),
          stepSize = pixdim[2:4],
          qfac = qfac
        )
      }
    }
  }

  space <- space_voxel(dim = vox_dim, affine = affine, mask_bitmap = mask_bitmap, mask_idx = mask_idx, storage = "packed")

  # Detect assays by scanning /data children
  data_grp <- h5$open("/data")
  on.exit(data_grp$close(), add = TRUE)
  entries <- data_grp$ls()
  dnames <- if (is.data.frame(entries)) entries$name else entries
  # Build canon assay list based on suffixes
  base_labels <- labels
  assays <- character(0)
  # Always support beta
  assays <- unique(c(assays, "beta"))
  if (any(paste0(base_labels, "_var") %in% dnames)) assays <- unique(c(assays, "var"))
  if (any(paste0(base_labels, "_stderr") %in% dnames)) assays <- unique(c(assays, "se"))

  dims <- gds_dims(sample = length(mask_idx), subject = length(handle$paths), contrast = length(labels))
  subjects <- basename(handle$paths)

  out <- list(
    assays = assays,
    dims = dims,
    subjects = as.character(subjects),
    contrasts = as.character(labels),
    space = space,
    maps = list(),
    metadata = source_metadata(),
    columns = list()
  )
  finish_probe(out)
}

.fmri_read <- function(handle,
                       assays,
                       block = NULL,
                       temporal_policy = c("as_is","mean","design"),
                       contrast_matrix = NULL,
                       contrast_names = NULL,
                       ...) {
  h5 <- handle$h5
  # Delegate to h5 adapter if /gds present
  if (h5$exists("/gds")) {
    adapter <- get_adapter("h5")
    return(adapter$read(list(file = h5, path = handle$paths[[1]]), assays = assays, block = block, ...))
  }

  # Path A: Parcellated (cluster/ROI)
  if (.h5_safe_exists(h5, "/cluster_map") && (.h5_safe_exists(h5, "/scans") || .h5_safe_exists(h5, "/data"))) {
    temporal_policy <- match.arg(temporal_policy)
    # cluster ids and labels
    cmap_ds <- h5$open("/cluster_map"); on.exit(cmap_ds$close(), add = TRUE)
    cluster_map <- as.array(cmap_ds$read())
    cluster_ids <- sort(unique(as.integer(cluster_map[cluster_map > 0])))
    parcel_labels <- as.character(cluster_ids)
    # Determine subjects
    subjects <- NULL
    if (.h5_safe_exists(h5, "/scans")) {
      scans_grp <- h5$open("/scans"); on.exit(scans_grp$close(), add = TRUE)
      scl <- scans_grp$ls(); subjects <- if (is.data.frame(scl)) as.character(scl$name) else as.character(scl)
    } else {
      subjects <- if (length(handle$paths) > 1) basename(handle$paths) else basename(handle$paths[[1]])
    }
    n_s <- length(subjects)

    # Helper to resolve dataset path for an assay
    assay_path_for <- function(h5file, subj, assay) {
      if (assay == "beta") {
        if (!is.null(subj) && .h5_safe_exists(h5file, paste0("/scans/", subj))) {
          ds_path <- .h5_first_dataset_in_group(h5file, paste0("/scans/", subj))
          if (!is.na(ds_path)) return(ds_path)
        }
        return("/data")
      }
      if (assay == "var") {
        if (!is.null(subj) && .h5_safe_exists(h5file, paste0("/scans_var/", subj))) {
          ds_path <- .h5_first_dataset_in_group(h5file, paste0("/scans_var/", subj))
          if (!is.na(ds_path)) return(ds_path)
        }
        if (.h5_safe_exists(h5file, "/var")) return("/var")
      }
      if (assay == "se") {
        if (!is.null(subj) && .h5_safe_exists(h5file, paste0("/scans_se/", subj))) {
          ds_path <- .h5_first_dataset_in_group(h5file, paste0("/scans_se/", subj))
          if (!is.na(ds_path)) return(ds_path)
        }
        if (.h5_safe_exists(h5file, "/se")) return("/se")
      }
      return(NA_character_)
    }

    # Read subject matrix and orient as [parcels x time]
    read_subject_matrix <- function(h5file, subj = NULL, assay = "beta") {
      path <- assay_path_for(h5file, subj, assay)
      if (is.na(path)) return(NULL)
      ds <- h5file$open(path); on.exit(ds$close(), add = TRUE)
      mat <- ds$read()
      dm <- dim(mat)
      if (length(dm) != 2L) stop("Parcellated data matrix must be 2D", call. = FALSE)
      if (dm[1] == length(cluster_ids)) {
        mat <- as.matrix(mat)
      } else if (dm[2] == length(cluster_ids)) {
        mat <- t(as.matrix(mat))
      } else {
        mat <- as.matrix(mat) # best effort
      }
      mat
    }

    # Derive time length and contrasts count
    m0 <- read_subject_matrix(h5, if (.h5_safe_exists(h5, "/scans")) subjects[[1]] else NULL, assay = "beta")
    time_len <- ncol(m0)
    samples_idx <- if (!is.null(block) && !is.null(block$sample)) {
      bs <- block$sample; if (is.logical(bs)) which(bs) else as.integer(bs)
    } else seq_along(parcel_labels)
    n_i <- length(samples_idx)

    # Setup contrast mapping
    if (temporal_policy == "as_is") {
      time_idx <- if (!is.null(block) && !is.null(block$contrast)) {
        bc <- block$contrast; if (is.logical(bc)) which(bc) else as.integer(bc)
      } else seq_len(time_len)
      n_c <- length(time_idx)
    } else if (temporal_policy == "mean") {
      time_idx <- NULL; n_c <- 1L
    } else { # design
      C <- contrast_matrix
      if (is.null(C)) {
        # fall back to identity
        C <- diag(time_len)
      }
      n_c <- ncol(C)
      time_idx <- NULL
    }

    out <- list()
    for (assay_name in assays) {
      arr <- array(NA_real_, dim = c(n_i, n_s, n_c))
      for (s in seq_len(n_s)) {
        h5s <- h5
        if (!.h5_safe_exists(h5s, "/scans") && s != 1 && length(handle$paths) >= s) {
          h5s <- hdf5r::H5File$new(handle$paths[[s]], mode = "r"); on.exit(h5s$close(), add = TRUE)
        }
        mat <- read_subject_matrix(h5s, if (.h5_safe_exists(h5, "/scans")) subjects[[s]] else NULL, assay = assay_name)
        if (is.null(mat)) {
          # No dataset for this assay; leave NA
          next
        }
        # Apply assay-specific companion if available (var/se not implemented separately yet)
        # Build selected rows
        rows <- samples_idx
        if (temporal_policy == "as_is") {
          vals <- mat[rows, time_idx, drop = FALSE]
        } else if (temporal_policy == "mean") {
          if (assay_name == "beta") {
            vals <- rowMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
            vals <- matrix(vals, nrow = length(rows), ncol = 1)
          } else {
            # conservative: do not fabricate aggregated uncertainty/statistics
            vals <- matrix(NA_real_, nrow = length(rows), ncol = 1)
          }
        } else {
          C <- contrast_matrix; if (is.null(C)) C <- diag(ncol(mat))
          if (assay_name == "beta") {
            vals <- mat[rows, , drop = FALSE] %*% C
          } else if (assay_name == "var") {
            # independence across timepoints: Var(C^T y) = (C^2)^T Var(y)
            C2 <- C^2
            vals <- mat[rows, , drop = FALSE] %*% C2
          } else if (assay_name == "se") {
            # derive later from var if present
            vals <- matrix(NA_real_, nrow = length(rows), ncol = ncol(C))
          } else {
            vals <- matrix(NA_real_, nrow = length(rows), ncol = ncol(C))
          }
        }
        arr[, s, ] <- vals
      }
      out[[assay_name]] <- arr
    }
    # If only se is present and var requested later, let derive handle it
    return(out)
  }

  # Path A: LatentNeuroVec-like
  if (.h5_safe_exists(h5, "/basis/basis_matrix") && .h5_safe_exists(h5, "/scans")) {
    scans_grp <- h5$open("/scans")
    on.exit(scans_grp$close(), add = TRUE)
    scl <- scans_grp$ls()
    subj <- if (is.data.frame(scl)) scl$name else scl
    # Read first embedding to get shapes
    e1 <- h5$open(paste0("/scans/", subj[[1]], "/embedding"))
    on.exit(e1$close(), add = TRUE)
    em <- e1$read()
    k <- if (is.null(dim(em))) length(em) else if (nrow(em) >= ncol(em)) nrow(em) else ncol(em)
    C <- if (is.null(dim(em))) 1L else if (nrow(em) >= ncol(em)) ncol(em) else nrow(em)
    samples_idx <- if (!is.null(block) && !is.null(block$sample)) {
      bs <- block$sample
      if (is.logical(bs)) which(bs) else as.integer(bs)
    } else seq_len(k)
    n_i <- length(samples_idx)
    out <- list()
    if ("beta" %in% assays) {
      arr <- array(NA_real_, dim = c(n_i, length(subj), C))
      for (s in seq_along(subj)) {
        emb_ds <- h5$open(paste0("/scans/", subj[[s]], "/embedding"))
        on.exit(emb_ds$close(), add = TRUE)
        emb <- emb_ds$read()
        # Coerce to [k, C]
        if (is.null(dim(emb))) {
          mat <- matrix(emb, nrow = k, ncol = 1)
        } else if (nrow(emb) == k) {
          mat <- emb
        } else if (ncol(emb) == k) {
          mat <- t(emb)
        } else {
          stop("Embedding shape incompatible with k", call. = FALSE)
        }
        arr[, s, ] <- mat[samples_idx, , drop = FALSE]
      }
      out$beta <- arr
      # provide placeholder se to satisfy GDS minimal constraint
      se_arr <- array(NA_real_, dim = dim(arr))
      out$se <- se_arr
    }
    return(out)
  }

  # Path A: LabeledVolumeSet-like
  labels <- .read_string_vec(h5, "/labels")
  mask <- .read_uint8_array(h5, "/mask")
  mask_idx <- which(as.vector(mask > 0))
  samples_idx <- if (!is.null(block) && !is.null(block$sample)) {
    bs <- block$sample
    if (is.logical(bs)) which(bs) else as.integer(bs)
  } else seq_along(mask_idx)
  # Determine selected subjects and contrasts
  n_s_total <- length(handle$paths)
  subj_idx <- if (!is.null(block) && !is.null(block$subject)) {
    bsj <- block$subject
    if (is.logical(bsj)) which(bsj) else as.integer(bsj)
  } else seq_len(n_s_total)
  ctr_idx <- if (!is.null(block) && !is.null(block$contrast)) {
    bc <- block$contrast
    if (is.logical(bc)) which(bc) else as.integer(bc)
  } else seq_along(labels)
  n_s <- length(subj_idx)
  n_c <- length(ctr_idx)
  n_i <- length(samples_idx)

  read_label_vec <- function(h5file, label, suffix = "") {
    safe_label <- .sanitize_label(label)
    path <- paste0("/data/", safe_label, suffix)
    if (!h5file$exists(path)) return(NULL)
    ds <- h5file$open(path)
    on.exit(ds$close(), add = TRUE)
    vec <- as.numeric(ds$read())
    if (length(vec) != length(mask_idx)) stop("Data length does not match mask", call. = FALSE)
    vec[samples_idx]
  }

  out <- list()
  for (assay_name in assays) {
    arr <- array(NA_real_, dim = c(n_i, n_s, n_c))
    for (si in seq_along(subj_idx)) {
      s <- subj_idx[[si]]
      # Multi-subject: open each file when reading
      h5s <- if (s == 1) h5 else hdf5r::H5File$new(handle$paths[[s]], mode = "r")
      if (s != 1) on.exit(h5s$close(), add = TRUE)
      for (ki in seq_along(ctr_idx)) {
        k <- ctr_idx[[ki]]
        lab <- labels[[k]]
        suffix <- switch(assay_name, beta = "", var = "_var", se = "_stderr", stop("Unsupported assay", call. = FALSE))
        vec <- read_label_vec(h5s, lab, suffix)
        if (!is.null(vec)) arr[, si, ki] <- vec
      }
    }
    out[[assay_name]] <- arr
  }
  out
}

# nocov end
