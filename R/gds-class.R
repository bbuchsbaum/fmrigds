#' Create a Group Data Set (GDS)
#'
#' @param assays Named list of 3-D arrays (sample x subject x contrast)
#' @param space A space object describing the sample axis
#' @param subjects Character vector of subject identifiers
#' @param contrasts Character vector of contrast identifiers
#' @param col_data Optional data frame keyed by subjects
#' @param row_data Optional data frame keyed by samples
#' @param metadata Optional list merged with [`gds_metadata()`]
#'
#' @return An object of class `c("gds", "group_data")`
#' @name gds
#' @export
new_gds <- function(assays,
                    space,
                    subjects,
                    contrasts,
                    col_data = NULL,
                    row_data = NULL,
                    metadata = list()) {
  assays <- .validate_gds_assays(assays)
  dims <- dim(assays[[1]])

  .validate_gds_dims(assays, subjects, contrasts)
  .validate_gds_space(space, dims[1L])

  col_data <- .normalise_col_data(col_data, subjects)
  row_data <- .normalise_row_data(row_data, dims[1L])

  metadata <- .merge_gds_metadata(gds_metadata(), metadata)

  structure(
    list(
      assays = assays,
      space = space,
      subjects = subjects,
      contrasts = contrasts,
      col_data = col_data,
      row_data = row_data,
      metadata = metadata
    ),
    class = c("gds", "group_data")
  )
}

# -------------------------------------------------------------------------
# Metadata & provenance ----------------------------------------------------

#' Construct default metadata for a GDS object
#'
#' @param schema_version Schema version string
#' @param units Named list of assay units
#' @param provenance Versioned provenance structure containing portable source
#'   entities, addressable activities, explicit input edges, and graph receipts.
#'   Legacy list-only graphs remain readable but are marked incomplete.
#' @param software Software metadata (package name, version, R version)
#' @param alignment Optional alignment metadata
#' @param map_families Named list of registered map families
#' @param mask_info Optional mask metadata
#' @param contrast_info Optional contrast metadata
#' @param design_mats Optional design matrices metadata
#' @param notes Optional user notes
#' @param created Timestamp for when the metadata was created
#'
#' @return Metadata list
#' @name gds_metadata
#' @export
gds_metadata <- function(schema_version = "0.2.0",
                         units = list(),
                         provenance = NULL,
                         software = list(
                           package = "fmrigds",
                           version = .pkg_version(),
                           R_version = getRversion()
                         ),
                         alignment = NULL,
                         map_families = list(),
                         mask_info = NULL,
                         contrast_info = NULL,
                         design_mats = NULL,
                         notes = NULL,
                         created = Sys.time()) {
  list(
    schema_version = schema_version,
    units = units,
    provenance = .normalize_provenance(provenance),
    software = software,
    alignment = alignment,
    map_families = map_families,
    mask_info = mask_info,
    contrast_info = contrast_info,
    design_mats = design_mats,
    notes = notes,
    created = created
  )
}

#' Create a provenance node
#'
#' @param op_name Operation name
#' @param params Named list of parameters
#' @param inputs Vector/list of parent hashes/ids
#' @param timestamp Timestamp of operation
#' @param software Software metadata
#' @param hash Optional pre-computed hash
#'
#' @return A provenance node list
#' @export
provenance_node <- function(op_name,
                            params,
                            inputs = list(),
                            timestamp = Sys.time(),
                            software = list(
                              package = "fmrigds",
                              version = .pkg_version(),
                              RVersion = as.character(getRversion())
                            ),
                            hash = NULL) {
  descriptor <- .activity_descriptor(op_name, params, inputs, software)
  id <- paste0("activity-", substr(.portable_sha256(list(
    descriptor = descriptor,
    occurrence = .next_portable_occurrence("activity"),
    timestamp = .iso_utc(timestamp)
  ))$value, 1L, 20L))
  node <- .new_activity_record(op_name, params, inputs, id, timestamp, software)
  if (!is.null(hash)) node$legacyHash <- as.character(hash)
  node
}

#' Append a provenance node to metadata
#'
#' @param metadata Metadata list from [`gds_metadata()`]
#' @param op_name Operation name
#' @param params Named list of parameters
#' @param inputs Parent node identifiers
#' @param id Optional unique activity occurrence identifier. Plan execution uses
#'   the plan node identifier; other callers normally let fmrigds mint one.
#' @param timestamp Timestamp of the activity occurrence.
#' @param software Portable software identity stored with the activity.
#'
#' @return Updated metadata list
#' @export
add_provenance_node <- function(metadata,
                                op_name,
                                params,
                                inputs = NULL,
                                id = NULL,
                                timestamp = Sys.time(),
                                software = list(
                                  package = "fmrigds",
                                  version = .pkg_version(),
                                  RVersion = as.character(getRversion())
                                )) {
  provenance <- .normalize_provenance(metadata$provenance %||% NULL)
  if (identical(provenance$status, "legacy-incomplete") && length(provenance$graph)) {
    stop(
      "Cannot append to legacy provenance without explicit addressable endpoints; ",
      "retain it as incomplete or start a new lineage receipt.",
      call. = FALSE
    )
  }
  if (is.null(inputs)) {
    if (length(provenance$graph)) {
      stop(
        "`inputs` is required when appending to a non-empty provenance graph.",
        call. = FALSE
      )
    }
    inputs <- character()
  }
  inputs <- as.character(unlist(inputs, use.names = FALSE))
  descriptor <- .activity_descriptor(op_name, params, inputs, software)
  id <- id %||% .next_activity_id(provenance, descriptor)
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("Provenance activity `id` must be one non-empty string.", call. = FALSE)
  }
  existing_ids <- c(
    vapply(provenance$entities %||% list(), `[[`, character(1L), "id"),
    vapply(provenance$graph %||% list(), `[[`, character(1L), "id")
  )
  if (id %in% existing_ids) {
    stop("Duplicate provenance id: ", id, call. = FALSE)
  }
  missing <- setdiff(inputs, existing_ids)
  if (length(missing)) {
    stop("Unknown provenance input endpoint(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  node <- .new_activity_record(op_name, params, inputs, id, timestamp, software)
  provenance$graph <- c(provenance$graph, list(node))

  log_entry <- sprintf(
    "[%s] %s(%s)",
    node$timestamp,
    op_name,
    paste0(
      names(params),
      "=",
      vapply(params, function(value) {
        tryCatch(.canonical_portable_json(value), error = function(e) paste(class(value), collapse = "/"))
      }, character(1L)),
      collapse = ", "
    )
  )
  provenance$log <- c(provenance$log, log_entry)
  provenance <- .refresh_provenance(provenance)
  problems <- .provenance_graph_problems(provenance, check_declared = TRUE)
  if (length(problems)) stop(paste(problems, collapse = "; "), call. = FALSE)
  metadata$provenance <- provenance

  metadata
}

# -------------------------------------------------------------------------
# Accessors ----------------------------------------------------------------

#' Extract assays from a GDS object
#'
#' @param x A GDS object
#'
#' @return Named list of 3D arrays
#' @export
assays <- function(x) UseMethod("assays")

#' @export
assays.gds <- function(x) x$assays

#' Extract a single assay from a GDS object
#'
#' @param x A GDS object
#' @param name Assay name (default: "beta")
#' @param ... Additional arguments
#'
#' @return A 3D array
#' @export
assay <- function(x, name = "beta", ...) UseMethod("assay")

#' @export
assay.gds <- function(x, name = "beta", ...) x$assays[[name]]

#' List the assays present on a GDS (or that a reducer would produce)
#'
#' Discoverability helper for finding assay names to pass to [assay()],
#' [write_nifti_assays()], or [write_out()] without having to run a job and
#' inspect `names(assays(fit))`. Works on a realised [`gds`], a lazy
#' [`gds_plan`]/[`gds_source`] (reporting the *input* assays the source probed),
#' or---when `reducer` is supplied---reports a reducer's declared output stems
#' without computing anything.
#'
#' @param x A realised [`gds`], a [`gds_plan`], or a [`gds_source`].
#' @param reducer Optional reducer name (e.g. `"meta:re"`, `"ols:voxelwise"`, or
#'   an alias like `"random"`). When supplied, `x` is ignored and the reducer's
#'   declared `provides` stems are returned. Regression/LMM reducers expand
#'   `coef`/`se_coef`/`z_coef`/`t_coef`/`p_coef` into per-term
#'   `coef:<term>` etc. at
#'   compute time, so realised names depend on the model formula.
#' @param info Logical; if `TRUE` (default) return a data frame with `assay`,
#'   `role`, and `units` columns (from [assay_info()] where registered, else
#'   `NA`); if `FALSE` return a plain character vector of names.
#'
#' @return A data frame (or a character vector when `info = FALSE`).
#' @seealso [assays()], [assay()], [assay_info()], [list_reducers()],
#'   [list_posthoc()], [write_nifti_assays()]
#' @export
#' @examples
#' \dontrun{
#' fit <- one_sample(g) |> compute()
#' list_assays(fit)                       # names + roles of computed assays
#' list_assays(reducer = "meta:re")       # what random-effects would produce
#' }
list_assays <- function(x, reducer = NULL, info = TRUE) {
  if (!is.null(reducer)) {
    red <- get_reducer(.normalize_reducer_name(reducer))
    if (is.null(red)) {
      stop("Unknown reducer: ", reducer, ". See list_reducers().", call. = FALSE)
    }
    nm <- as.character(red$provides %||% character(0))
  } else if (inherits(x, "gds")) {
    nm <- names(assays(x))
  } else if (inherits(x, "gds_source")) {
    nm <- as.character(x$probe$assays %||% character(0))
  } else if (inherits(x, "gds_plan")) {
    nm <- as.character(x$source$probe$assays %||% character(0))
  } else {
    stop("`x` must be a gds, gds_plan, or gds_source (or supply `reducer`).", call. = FALSE)
  }
  nm <- unique(nm)
  if (!isTRUE(info)) {
    return(nm)
  }
  meta <- lapply(nm, assay_info)
  data.frame(
    assay = nm,
    role = vapply(meta, function(m) if (is.null(m)) NA_character_ else m$role, character(1)),
    units = vapply(meta, function(m) {
      u <- if (is.null(m)) NA_character_ else m$units %||% NA_character_
      if (length(u) != 1L) NA_character_ else as.character(u)
    }, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Extract space descriptor from a GDS object
#'
#' For non-GDS objects (e.g. a \pkg{neuroim2} `NeuroVol`/`NeuroVec`) this
#' generic falls back to [neuroim2::space()] when that package is available, so
#' that attaching \pkg{fmrigds} does not mask `neuroim2::space()` on the search
#' path.
#'
#' @param x A GDS object, or any object with a `space` method (such as a
#'   \pkg{neuroim2} image).
#' @param ... Additional arguments passed to methods.
#'
#' @return A space object
#' @export
space <- function(x, ...) UseMethod("space")

#' @export
space.gds <- function(x, ...) x$space

#' @export
space.gds_plan <- function(x, ...) x$source$probe$space

#' @export
space.default <- function(x, ...) {
  if (requireNamespace("neuroim2", quietly = TRUE)) {
    return(neuroim2::space(x, ...))
  }
  stop(
    "No applicable 'space()' method for an object of class ",
    paste(deparse(class(x)), collapse = ""),
    ". Install 'neuroim2' to use space() on neuroimaging objects.",
    call. = FALSE
  )
}

#' Extract subject identifiers from a GDS object
#'
#' @param x A GDS object
#'
#' @return Character vector of subject IDs
#' @export
subjects <- function(x) UseMethod("subjects")

#' @export
subjects.gds <- function(x) x$subjects

#' @export
subjects.gds_plan <- function(x) {
  # Prefer subjects from the bound source probe, fall back to plan meta
  x$source$probe$subjects %||% x$meta$subjects
}

#' Extract contrast identifiers from a GDS object
#'
#' @param x A GDS object
#'
#' @return Character vector of contrast names
#' @export
contrasts <- function(x) UseMethod("contrasts")

#' @export
contrasts.gds <- function(x) x$contrasts

#' @export
contrasts.gds_plan <- function(x) {
  x$source$probe$contrasts %||% (x$meta$contrasts %||% character())
}

#' Extract column (subject) metadata from a GDS object
#'
#' @param x A GDS object
#'
#' @return Data frame with subject-level metadata
#' @export
col_data <- function(x) UseMethod("col_data")

#' @export
col_data.gds <- function(x) x$col_data

#' @export
col_data.gds_plan <- function(x) {
  x$meta$col_data %||% x$source$probe$col_data
}

#' Extract row (sample) metadata from a GDS object
#'
#' @param x A GDS object
#'
#' @return Data frame with sample-level metadata
#' @export
row_data <- function(x) UseMethod("row_data")

#' @export
row_data.gds <- function(x) x$row_data

#' @export
row_data.gds_plan <- function(x) {
  x$meta$row_data %||% x$source$probe$row_data
}

#' Extract contrast-level metadata from a GDS object
#'
#' @param x A GDS object or plan
#'
#' @return Data frame with one row per contrast/repeated-measure level
#' @export
contrast_data <- function(x) UseMethod("contrast_data")

#' @export
contrast_data.gds <- function(x) {
  info <- x$metadata$contrast_info %||% NULL
  if (is.null(info)) return(NULL)
  info$data %||% NULL
}

#' @export
contrast_data.gds_plan <- function(x) {
  x$meta$contrast_data %||% x$source$probe$contrast_data
}

#' Extract sample labels from a GDS object
#'
#' @param x A GDS object or plan
#'
#' @return Character vector of sample labels
#' @export
sample_labels <- function(x) UseMethod("sample_labels")

#' @export
sample_labels.gds <- function(x) {
  .sample_labels_from_row_data(row_data(x), space(x), metadata(x))
}

#' @export
sample_labels.gds_plan <- function(x) {
  .sample_labels_from_row_data(row_data(x), space(x), x$source$probe$metadata %||% list())
}

#' Extract sample-group metadata from a GDS object
#'
#' @param x A GDS object or plan
#' @param vars Candidate column names searched in `row_data(x)`
#'
#' @return A vector of group labels or `NULL` when unavailable
#' @export
sample_groups <- function(x, vars = c("feature_group", "spatial_group", "group", "parcel", "parcel_id")) UseMethod("sample_groups")

#' @export
sample_groups.gds <- function(x, vars = c("feature_group", "spatial_group", "group", "parcel", "parcel_id")) {
  .sample_groups_from_row_data(row_data(x), vars = vars)
}

#' @export
sample_groups.gds_plan <- sample_groups.gds

#' Extract metadata from a GDS object
#'
#' @param x A GDS object
#'
#' @return Metadata list
#' @export
metadata <- function(x) UseMethod("metadata")

#' @export
metadata.gds <- function(x) x$metadata

# -------------------------------------------------------------------------
# Helpers ------------------------------------------------------------------

.pkg_version <- function() {
  out <- try(utils::packageVersion("fmrigds"), silent = TRUE)
  if (inherits(out, "try-error")) return("0.0.0")
  as.character(out)
}

.validate_gds_assays <- function(assays) {
  if (!is.list(assays) || !length(assays)) {
    stop("`assays` must be a non-empty named list of arrays", call. = FALSE)
  }
  if (is.null(names(assays)) || any(!nzchar(names(assays)))) {
    stop("`assays` must be a named list", call. = FALSE)
  }

  # Upcast any 2-D arrays to 3-D along contrast dimension
  assays <- lapply(assays, function(a) {
    if (is.array(a) && length(dim(a)) == 2L) {
      a <- array(a, dim = c(nrow(a), ncol(a), 1L))
    }
    a
  })
  dims <- lapply(assays, dim)
  first_dim <- dims[[1L]]
  if (!all(vapply(dims, identical, logical(1L), first_dim))) {
    stop("All assays must share identical dimensions", call. = FALSE)
  }

  if ("beta" %in% names(assays)) {
    has_var <- "var" %in% names(assays)
    has_se <- "se" %in% names(assays)
    if (!has_var && !has_se) {
      stop("If 'beta' is provided, either 'var' or 'se' must also be present", call. = FALSE)
    }
  }

  assays
}

.validate_gds_dims <- function(assays, subjects, contrasts) {
  dims <- dim(assays[[1L]])
  if (length(subjects) != dims[2L]) {
    stop(sprintf("Length of `subjects` (%d) must equal dim[2] (%d)", length(subjects), dims[2L]), call. = FALSE)
  }
  if (length(contrasts) != dims[3L]) {
    stop(sprintf("Length of `contrasts` (%d) must equal dim[3] (%d)", length(contrasts), dims[3L]), call. = FALSE)
  }
  invisible(NULL)
}

.validate_gds_space <- function(space, n_samples) {
  if (!inherits(space, "gds_space")) {
    stop("`space` must inherit from class 'gds_space'", call. = FALSE)
  }
  if (!is.null(space$storage) && identical(space$storage, "packed")) {
    mask_len <- length(space$mask_idx)
    if (!isTRUE(mask_len == n_samples)) {
      stop("Packed spaces must have mask_idx length equal to sample dimension", call. = FALSE)
    }
  } else {
    if (!is.null(space$dim) && prod(space$dim) != n_samples) {
      stop("Space dimensions do not align with sample count", call. = FALSE)
    }
  }
  invisible(NULL)
}

.normalise_col_data <- function(col_data, subjects) {
  if (is.null(col_data)) {
    return(data.frame(subject = subjects, row.names = subjects, check.names = FALSE))
  }
  if (!is.data.frame(col_data)) {
    stop("`col_data` must be a data.frame", call. = FALSE)
  }
  if (!all(subjects %in% rownames(col_data))) {
    stop("`col_data` must contain rownames matching `subjects`", call. = FALSE)
  }
  col_data[subjects, , drop = FALSE]
}

.normalise_row_data <- function(row_data, n_samples) {
  if (is.null(row_data)) {
    return(data.frame(row.names = as.character(seq_len(n_samples)), check.names = FALSE))
  }
  if (!is.data.frame(row_data)) {
    stop("`row_data` must be a data.frame", call. = FALSE)
  }
  if (nrow(row_data) != n_samples) {
    stop("`row_data` must have one row per sample", call. = FALSE)
  }
  row_data
}

.normalise_contrast_data <- function(contrast_data, contrasts) {
  if (is.null(contrast_data)) {
    return(NULL)
  }
  if (!is.data.frame(contrast_data)) {
    stop("`contrast_data` must be a data.frame", call. = FALSE)
  }
  if (!length(contrasts)) {
    if (nrow(contrast_data)) {
      stop("`contrast_data` supplied but no contrasts are available", call. = FALSE)
    }
    return(contrast_data)
  }
  rn <- rownames(contrast_data)
  if (!is.null(rn) && !anyNA(rn) && !anyDuplicated(rn)) {
    missing <- setdiff(contrasts, rn)
    if (length(missing)) {
      stop("`contrast_data` must contain rownames matching `contrasts`", call. = FALSE)
    }
    return(contrast_data[contrasts, , drop = FALSE])
  }
  if (nrow(contrast_data) != length(contrasts)) {
    stop("`contrast_data` must have one row per contrast", call. = FALSE)
  }
  rownames(contrast_data) <- contrasts
  contrast_data
}

.is_positional_sample_labels <- function(x, n) {
  if (is.null(x) || is.null(n) || !is.finite(n) || length(x) != n) {
    return(FALSE)
  }
  identical(as.character(x), as.character(seq_len(n)))
}

.sample_labels_from_row_data <- function(row_data, space = NULL, metadata = list()) {
  n_samples <- if (!is.null(row_data)) {
    nrow(row_data)
  } else if (!is.null(space$labels)) {
    length(space$labels)
  } else if (!is.null(space$mask_idx)) {
    length(space$mask_idx)
  } else {
    NA_integer_
  }
  synthetic <- isTRUE(metadata$sample_labels_synthetic)

  if (!is.null(row_data) && nrow(row_data) > 0L) {
    for (nm in c("label", "roi", "parcel")) {
      if (nm %in% names(row_data)) {
        return(as.character(row_data[[nm]]))
      }
    }
  }

  if (!is.null(space)) {
    if (!is.null(space$labels)) {
      if (synthetic && .is_positional_sample_labels(space$labels, n_samples)) {
        return(NULL)
      }
      return(as.character(space$labels))
    }
    if (!is.null(space$mask_idx)) {
      if (synthetic) {
        return(NULL)
      }
      return(as.character(space$mask_idx))
    }
  }

  if (!is.null(row_data) && nrow(row_data) > 0L) {
    if ("sample" %in% names(row_data)) {
      if (.is_positional_sample_labels(row_data[["sample"]], nrow(row_data))) {
        return(NULL)
      }
      return(as.character(row_data[["sample"]]))
    }
    rn <- rownames(row_data)
    if (!is.null(rn) && !anyNA(rn)) {
      if (.is_positional_sample_labels(rn, nrow(row_data))) {
        return(NULL)
      }
      return(as.character(rn))
    }
  }

  NULL
}

.sample_groups_from_row_data <- function(row_data, vars = c("feature_group", "spatial_group", "group", "parcel", "parcel_id")) {
  if (is.null(row_data) || !nrow(row_data)) {
    return(NULL)
  }

  match_name <- intersect(vars, names(row_data))
  if (!length(match_name)) {
    return(NULL)
  }

  row_data[[match_name[[1L]]]]
}
