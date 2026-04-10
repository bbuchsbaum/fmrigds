#' Register the neurotabs (NFTab) adapter
#'
#' Makes \code{nftab} objects from the \pkg{neurotabs} package a first-class
#' data source for \code{\link{gds}()}.  Called automatically during
#' \code{.onLoad()} when \pkg{neurotabs} is installed; can also be called
#' manually.
#'
#' @return Invisibly, the adapter entry.
#' @export
register_nftab_adapter <- function() {
  register_adapter(
    name    = "nftab",
    detect  = .nftab_detect,
    open    = .nftab_open,
    probe   = .nftab_probe,
    read    = .nftab_read,
    close   = .nftab_close
  )
}

# ---------------------------------------------------------------------------
# detect
# ---------------------------------------------------------------------------
.nftab_detect <- function(source) {
  if (inherits(source, "nftab")) return(0.95)
  0
}

# ---------------------------------------------------------------------------
# open
# ---------------------------------------------------------------------------
.nftab_open <- function(source, mode = "r", ...) {
  if (!identical(mode, "r")) stop("nftab adapter is read-only", call. = FALSE)
  if (!requireNamespace("neurotabs", quietly = TRUE)) {
    stop("The 'neurotabs' package is required for the nftab adapter", call. = FALSE)
  }
  dots <- list(...)
  list(
    nftab        = source,
    feature      = dots$feature,
    features     = dots$features,
    subject_axis = dots$subject_axis,
    contrast_axis = dots$contrast_axis
  )
}

# ---------------------------------------------------------------------------
# probe
# ---------------------------------------------------------------------------
.nftab_probe <- function(handle, ...) {
  dots <- list(...)
  tab  <- handle$nftab

  # --- Resolve feature → assay mapping ---
  feature_map <- .nftab_resolve_features(handle, dots)

  # --- Discover subject & contrast axes ---
  manifest    <- tab$manifest
  design      <- neurotabs::nf_design(tab)
  subject_col <- .nftab_find_axis(manifest, "subject",
                                  handle$subject_axis %||% dots$subject_axis)
  contrast_col <- .nftab_find_axis(manifest, c("contrast", "condition"),
                                   handle$contrast_axis %||% dots$contrast_axis)

  subjects  <- as.character(unique(design[[subject_col]]))
  contrasts <- .nftab_contrast_labels(design, contrast_col, subject_col)

  # --- Determine spatial extent from first feature ---
  first_feat   <- feature_map[[1L]]
  logical_schema <- neurotabs::nf_feature_schema(tab, first_feat)
  n_samples    <- .nftab_sample_count(logical_schema, tab, feature_name = first_feat)
  space        <- .nftab_build_space(logical_schema, tab, feature_name = first_feat)

  dims <- gds_dims(
    sample   = n_samples,
    subject  = length(subjects),
    contrast = length(contrasts)
  )

  # --- Build col_data / row_data ---
  col_data <- .nftab_extract_col_data(design, subjects, subject_col)

  out <- list(
    assays   = names(feature_map),
    dims     = dims,
    subjects = subjects,
    contrasts = contrasts,
    space    = space,
    maps     = list(),
    metadata = list(
      schema_version = "0.1.0",
      source         = "nftab",
      dataset_id     = manifest$dataset_id %||% "<unknown>"
    ),
    columns  = list(
      effect_cols  = feature_map,
      subject_col  = subject_col,
      contrast_col = contrast_col
    ),
    col_data = col_data
  )

  probe_contract(out)
}

# ---------------------------------------------------------------------------
# read
# ---------------------------------------------------------------------------
.nftab_read <- function(handle, assays, block = NULL, ...) {
  dots    <- list(...)
  tab     <- handle$nftab
  # compute() passes adapter_columns as effect_cols/subject_col/contrast_col
  fmap    <- dots$effect_cols
  sub_col <- dots$subject_col
  con_col <- dots$contrast_col

  design   <- neurotabs::nf_design(tab)
  subjects <- as.character(unique(design[[sub_col]]))
  contrasts <- .nftab_contrast_labels(design, con_col, sub_col)

  # Apply block subsetting to subjects/contrasts
  if (!is.null(block$subject)) {
    idx <- if (is.logical(block$subject)) which(block$subject) else as.integer(block$subject)
    subjects <- subjects[idx]
  }
  if (!is.null(block$contrast)) {
    idx <- if (is.logical(block$contrast)) which(block$contrast) else as.integer(block$contrast)
    contrasts <- contrasts[idx]
  }

  out <- lapply(assays, function(assay_name) {
    feat_name <- fmap[[assay_name]]
    .nftab_read_feature(tab, feat_name, design, subjects, contrasts,
                        sub_col, con_col, block$sample)
  })
  names(out) <- assays
  out
}

# ---------------------------------------------------------------------------
# close
# ---------------------------------------------------------------------------
.nftab_close <- function(handle) {
  invisible(NULL)
}

# ===========================================================================
# Internal helpers
# ===========================================================================

#' Resolve the feature → assay name mapping from user arguments
#' @keywords internal
.nftab_resolve_features <- function(handle, dots) {
  feature  <- handle$feature %||% dots$feature
  features <- handle$features %||% dots$features

  if (!is.null(features)) {
    # Named vector: c(beta = "statmap", var = "statmap_var")
    if (is.null(names(features))) {
      names(features) <- features
    }
    return(as.list(features))
  }

  tab <- handle$nftab
  avail <- neurotabs::nf_feature_names(tab)

  if (!is.null(feature)) {
    return(list(beta = feature))
  }

  # Default: use first feature as "beta"
  if (length(avail) == 0L) {
    stop("nftab has no features", call. = FALSE)
  }
  list(beta = tab$manifest$primary_feature %||% avail[1L])
}

#' Find an observation axis column by semantic role or user override
#' @keywords internal
.nftab_find_axis <- function(manifest, roles, user_override = NULL) {
  if (!is.null(user_override)) {
    if (!user_override %in% manifest$observation_axes) {
      stop(
        sprintf("User-specified axis '%s' is not in observation_axes: %s",
                user_override, paste(manifest$observation_axes, collapse = ", ")),
        call. = FALSE
      )
    }
    return(user_override)
  }

  # Search by semantic role
  obs_cols <- manifest$observation_columns
  for (role in roles) {
    for (col_name in names(obs_cols)) {
      if (identical(obs_cols[[col_name]]$semantic_role, role)) {
        return(col_name)
      }
    }
  }

  # Fallback: look for columns named like the roles
  axes <- manifest$observation_axes
  for (role in roles) {
    matches <- grep(role, axes, ignore.case = TRUE, value = TRUE)
    if (length(matches)) return(matches[1L])
  }

  # Last resort for subject: first axis
  if ("subject" %in% roles) return(axes[1L])

  # For contrast: synthesize a single "all" contrast

  NULL
}

#' Build contrast labels from the design table
#' @keywords internal
.nftab_contrast_labels <- function(design, contrast_col, subject_col) {
  if (is.null(contrast_col)) {
    # No contrast axis — every subject has one observation → single contrast
    return("all")
  }
  as.character(unique(design[[contrast_col]]))
}

#' Compute the number of spatial samples from the logical schema
#' @keywords internal
.nftab_sample_count <- function(logical_schema, tab, feature_name = NULL) {
  kind <- logical_schema$kind
  shape <- logical_schema$shape

  if (kind == "volume" && !is.null(shape) && length(shape) == 3L) {
    return(prod(shape))
  }
  if (kind %in% c("vector", "matrix", "surface", "array") &&
      !is.null(shape) && length(shape) >= 1L) {
    return(prod(shape))
  }
  # Fallback: resolve first row and check length
  feat_name <- feature_name %||% tab$manifest$primary_feature %||%
    names(tab$manifest$features)[1L]
  first <- neurotabs::nf_resolve(tab, 1L, feat_name)
  length(as.vector(first))
}

#' Map neurotabs support/logical schema to a GDS space object
#' @keywords internal
.nftab_build_space <- function(logical_schema, tab, feature_name = NULL) {
  kind <- logical_schema$kind
  shape <- logical_schema$shape
  support_ref <- logical_schema$support_ref

  # Try to get support info

  support <- NULL
  if (!is.null(support_ref) && !is.null(tab$manifest$supports)) {
    support <- tryCatch(
      neurotabs::nf_support_info(tab, support_ref),
      error = function(e) NULL
    )
  }

  if (!is.null(support)) {
    stype <- support$support_type

    if (identical(stype, "volume")) {
      # Volume: need dim and affine
      vdim <- if (!is.null(shape) && length(shape) == 3L) as.integer(shape) else c(1L, 1L, 1L)
      affine <- if (!is.null(support$metadata$affine)) {
        matrix(as.numeric(support$metadata$affine), 4L, 4L)
      } else {
        diag(4)
      }
      return(.nftab_annotate_space(
        space_voxel(dim = vdim, affine = affine),
        support
      ))
    }

    if (identical(stype, "parcel")) {
      n <- support$n_parcels %||% shape[1L] %||% 1L
      labels <- .nftab_parcel_labels(support, tab, n)
      return(.nftab_annotate_space(
        space_parcels(labels = labels),
        support
      ))
    }

    if (identical(stype, "surface")) {
      n <- if (!is.null(shape)) prod(shape) else 1L
      return(.nftab_surface_space(support, n))
    }
  }

  # No support or generic: use sample labels
  n <- if (!is.null(shape)) prod(shape) else 1L
  # Try to get labels from axis domains
  labels <- .nftab_sample_labels(
    logical_schema,
    tab,
    n,
    feature_name = feature_name
  )
  space_sample_labels(labels = labels)
}

.nftab_annotate_space <- function(space, support) {
  if (is.null(support)) return(space)

  space$support_type <- support$support_type %||% NULL
  space$support_id <- support$support_id %||% NULL
  space$reference_space <- support$space %||% NULL
  space$grid_id <- support$grid_id %||% NULL
  space$affine_id <- support$affine_id %||% NULL
  if (is.null(space$template_id)) {
    space$template_id <- support$template %||% NULL
  }
  space$mesh_id <- support$mesh_id %||% NULL
  space$topology_id <- support$topology_id %||% NULL
  if (is.null(space$hemi)) {
    hemi <- support$hemisphere %||% NULL
    if (!is.null(hemi)) {
      space$hemi <- .nftab_surface_hemi(hemi)
    }
  }
  space
}

.nftab_surface_space <- function(support, n) {
  .nftab_annotate_space(
    structure(
      list(
        type = "surface",
        vertices = NULL,
        faces = NULL,
        hemi = .nftab_surface_hemi(support$hemisphere %||% NA_character_),
        template_id = support$template %||% NULL,
        labels = as.character(seq_len(n)),
        n_vertices = as.integer(n)
      ),
      class = c("Space", "space_surface", "gds_space")
    ),
    support
  )
}

.nftab_surface_hemi <- function(x) {
  hemi <- tolower(as.character(x))
  switch(
    hemi,
    left = "L",
    l = "L",
    right = "R",
    r = "R",
    both = "LR",
    lr = "LR",
    midline = "LR",
    x
  )
}

#' Extract parcel labels from support info or generate defaults
#' @keywords internal
.nftab_parcel_labels <- function(support, tab, n) {
  # Try parcel_map TSV first
  if (!is.null(support$parcel_map) && !is.null(tab$.root)) {
    map_path <- file.path(tab$.root, support$parcel_map)
    if (file.exists(map_path)) {
      map_df <- tryCatch(
        utils::read.delim(map_path, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (!is.null(map_df) && "label" %in% names(map_df)) {
        return(as.character(map_df$label))
      }
    }
  }
  paste0("P", seq_len(n))
}

#' Get sample-level labels from axis domains or generate defaults
#' @keywords internal
.nftab_sample_labels <- function(logical_schema, tab, n, feature_name = NULL) {
  # Try axis_domains on first axis
  if (!is.null(logical_schema$axis_domains)) {
    first_domain <- logical_schema$axis_domains[[1L]]
    if (!is.null(first_domain) && !is.null(first_domain$labels)) {
      labels <- tryCatch({
        feat_name <- feature_name %||% tab$manifest$primary_feature %||%
          names(tab$manifest$features)[1L]
        axis_name <- logical_schema$axes[1L]
        lbl <- neurotabs::nf_axis_labels(tab, feat_name, axis_name)
        as.character(lbl$label)
      }, error = function(e) NULL)
      if (!is.null(labels) && length(labels) == n) return(labels)
    }
  }
  as.character(seq_len(n))
}

#' Extract subject-level covariates from the design table
#' @keywords internal
.nftab_extract_col_data <- function(design, subjects, subject_col) {
  # Keep only columns that have one unique value per subject
  non_axis <- setdiff(names(design), subject_col)
  col_data <- data.frame(row.names = subjects)

  for (col in non_axis) {
    vals <- tapply(design[[col]], design[[subject_col]], function(v) {
      uv <- unique(v)
      if (length(uv) == 1L) uv else NA
    }, simplify = TRUE)
    vals <- as.vector(vals[subjects])
    # Only include if all subjects have a unique value
    if (!anyNA(vals)) {
      col_data[[col]] <- vals
    }
  }

  if (ncol(col_data) == 0L) return(NULL)
  col_data
}

#' Read and reshape a single feature into a 3D array
#' @keywords internal
.nftab_read_feature <- function(tab, feat_name, design, subjects, contrasts,
                                subject_col, contrast_col, sample_block) {
  logical_schema <- neurotabs::nf_feature_schema(tab, feat_name)
  n_samples <- .nftab_sample_count(logical_schema, tab, feature_name = feat_name)
  ns <- length(subjects)
  nc <- length(contrasts)

  arr <- array(NA_real_, dim = c(n_samples, ns, nc))

  # Build row lookup: for each (subject, contrast) → row index in the nftab
  for (si in seq_along(subjects)) {
    for (ci in seq_along(contrasts)) {
      subj <- subjects[si]
      if (is.null(contrast_col)) {
        row_mask <- design[[subject_col]] == subj
      } else {
        row_mask <- design[[subject_col]] == subj &
                    design[[contrast_col]] == contrasts[ci]
      }
      row_idx <- which(row_mask)
      if (length(row_idx) == 0L) next
      # Use first matching row if multiple (shouldn't happen with unique axes)
      row_idx <- row_idx[1L]

      value <- neurotabs::nf_resolve(tab, row_idx, feat_name, as_array = TRUE)
      arr[, si, ci] <- as.vector(value)[seq_len(n_samples)]
    }
  }

  # Apply sample block subsetting
  if (!is.null(sample_block)) {
    idx <- if (is.logical(sample_block)) which(sample_block) else as.integer(sample_block)
    arr <- arr[idx, , , drop = FALSE]
  }

  arr
}
