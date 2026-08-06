# Exception-safe staging for non-block-local prefixes ----------------------

#' Classify the portability of a GDS space
#'
#' @param space A GDS space object.
#'
#' @return One of `"portable"`, `"portable_with_external_reference"`, or
#'   `"nonportable_function"`.
#' @export
space_portability <- function(space) {
  if (inherits(space, c("space_voxel", "space_parcels", "space_sample_labels", "space_surface"))) {
    return("portable")
  }
  if (inherits(space, "space_basis")) {
    if (is.function(space$projector)) return("nonportable_function")
    if (is.null(space$projector) && !is.null(space$basis_name)) {
      return("portable_with_external_reference")
    }
    return("portable")
  }
  "nonportable_function"
}

.source_fingerprint <- function(source) {
  specification <- source$source
  if (is.character(specification) && length(specification) &&
      all(file.exists(specification))) {
    info <- file.info(specification)
    entries <- data.frame(
      path = normalizePath(specification, mustWork = TRUE),
      size = as.numeric(info$size),
      mtime = as.numeric(info$mtime),
      stringsAsFactors = FALSE
    )
    return(list(
      kind = "files",
      entries = entries,
      digest = digest::digest(entries, algo = "xxhash64")
    ))
  }
  value <- list(kind = "source_object", source_hash = source$hash)
  value$digest <- digest::digest(value, algo = "xxhash64")
  value
}

.stage_examination_plan <- function(compiled, control) {
  if (!inherits(compiled, "gds_examination_plan")) {
    stop("compiled must be a gds_examination_plan.", call. = FALSE)
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Staged group examination requires the optional hdf5r package.", call. = FALSE)
  }
  stage_dir <- control$staging$tempdir %||% tempdir()
  if (!dir.exists(stage_dir)) {
    stop("Staging tempdir does not exist: ", stage_dir, call. = FALSE)
  }
  stage_path <- tempfile(
    pattern = "fmrigds-group-examination-",
    tmpdir = stage_dir,
    fileext = ".h5"
  )
  started <- Sys.time()
  prefix_plan <- compiled$plan
  prefix_plan$nodes <- compiled$prefix
  realized <- compute(prefix_plan)
  portability <- space_portability(space(realized))
  h5_portable <- inherits(
    space(realized),
    c("space_voxel", "space_parcels", "space_sample_labels")
  )
  if (!h5_portable) {
    stop(
      "The staged prefix produced a ", portability,
      " space of class ", paste(class(space(realized)), collapse = "/"),
      ", but native HDF5 staging currently supports voxel, parcel, and sample-label spaces.",
      call. = FALSE
    )
  }
  write_gds_h5(realized, stage_path)
  stage_size <- as.numeric(file.info(stage_path)$size)
  staged_plan <- gds(stage_path, format = "h5")
  staged_plan$source$probe$metadata <- utils::modifyList(
    staged_plan$source$probe$metadata %||% list(),
    metadata(realized) %||% list()
  )
  staged_plan$nodes <- c(list(compiled$reducer), compiled$conclusion_tail)
  staged_plan <- .ensure_plan_node_ids(staged_plan)
  staged <- compile_examination_plan(staged_plan)
  staged$source_plan_digest <- compiled$source_plan_digest
  staged$origin_source_hash <- compiled$origin_source_hash
  staged$source_fingerprint <- compiled$source_fingerprint
  staged$source_fingerprint_digest <- compiled$source_fingerprint_digest
  staged$origin_node_ids <- compiled$node_ids
  staged$origin_discarded_writes <- compiled$discarded_writes
  record <- list(
    strategy = "stage",
    source_adapter = compiled$plan$source$adapter,
    staged_adapter = "h5",
    reasons = compiled$scan_reasons,
    portability = portability,
    adapter_reads = 1L,
    bytes_read = sum(vapply(assays(realized), utils::object.size, numeric(1))),
    stage_size_bytes = stage_size,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    cleanup_policy = "always",
    cleanup_succeeded = NA
  )
  list(compiled = staged, path = stage_path, record = record)
}

.cleanup_examination_stage <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(TRUE)
  isTRUE(unlink(path) == 0L) && !file.exists(path)
}
