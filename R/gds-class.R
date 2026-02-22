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

  metadata <- utils::modifyList(gds_metadata(), metadata)

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
#' @param provenance Provenance structure (`graph`, `log`, `digest`)
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
gds_metadata <- function(schema_version = "0.1.0",
                         units = list(),
                         provenance = list(graph = list(), log = character(), digest = NULL),
                         software = list(
                           package = "gdsfmri",
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
    provenance = provenance,
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
                            software = list(package = "gdsfmri", version = .pkg_version()),
                            hash = NULL) {
  if (is.null(hash)) {
    hash <- digest::digest(list(op_name = op_name, params = params, inputs = inputs))
  }

  structure(
    list(
      op = op_name,
      params = params,
      inputs = inputs,
      timestamp = timestamp,
      software = software,
      hash = hash
    ),
    class = "gds_provenance_node"
  )
}

#' Append a provenance node to metadata
#'
#' @param metadata Metadata list from [`gds_metadata()`]
#' @param op_name Operation name
#' @param params Named list of parameters
#' @param inputs Parent node identifiers
#'
#' @return Updated metadata list
#' @export
add_provenance_node <- function(metadata,
                                op_name,
                                params,
                                inputs = list()) {
  if (is.null(metadata$provenance)) {
    metadata$provenance <- list(graph = list(), log = character(), digest = NULL)
  }

  node <- provenance_node(op_name, params, inputs)
  metadata$provenance$graph <- c(metadata$provenance$graph, list(node))

  log_entry <- sprintf(
    "[%s] %s(%s)",
    format(node$timestamp, "%Y-%m-%d %H:%M:%S"),
    op_name,
    paste0(names(params), "=", vapply(params, format, character(1L)), collapse = ", ")
  )
  metadata$provenance$log <- c(metadata$provenance$log, log_entry)

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

#' Extract space descriptor from a GDS object
#'
#' @param x A GDS object
#'
#' @return A space object
#' @export
space <- function(x) UseMethod("space")

#' @export
space.gds <- function(x) x$space

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

#' Extract row (sample) metadata from a GDS object
#'
#' @param x A GDS object
#'
#' @return Data frame with sample-level metadata
#' @export
row_data <- function(x) UseMethod("row_data")

#' @export
row_data.gds <- function(x) x$row_data

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
  out <- try(utils::packageVersion("gdsfmri"), silent = TRUE)
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
    return(data.frame(sample = seq_len(n_samples), row.names = seq_len(n_samples), check.names = FALSE))
  }
  if (!is.data.frame(row_data)) {
    stop("`row_data` must be a data.frame", call. = FALSE)
  }
  if (nrow(row_data) != n_samples) {
    stop("`row_data` must have one row per sample", call. = FALSE)
  }
  row_data
}
