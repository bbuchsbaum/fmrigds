register_memory_adapter <- function() {
  register_adapter(
    name = "memory",
    detect = .memory_detect,
    open = .memory_open,
    probe = .memory_probe,
    read = .memory_read,
    close = .memory_close,
    capabilities = list(
      sample_blocks = TRUE,
      subject_blocks = TRUE,
      contrast_blocks = TRUE,
      persistent_handle = TRUE,
      preferred_axis = "sample",
      cheap_revisit = TRUE
    )
  )
}

.memory_detect <- function(source) {
  # Only used when detect_adapter(source) is called; explicit use via as_plan.gds sets adapter="memory".
  if (inherits(source, "gds")) return(1.0)
  if (is.list(source) && length(source) && !is.null(names(source))) {
    # crude heuristic for list of arrays
    dims <- try(dim(source[[1L]]), silent = TRUE)
    if (!inherits(dims, "try-error") && length(dims) == 3L) return(0.9)
  }
  0
}

.memory_open <- function(source, mode = "r", ...) {
  if (!identical(mode, "r")) stop("memory adapter is read-only", call. = FALSE)
  if (inherits(source, "gds")) {
    return(list(gds = source))
  }
  if (is.list(source)) {
    return(list(arrays = source))
  }
  stop("Unsupported source for memory adapter", call. = FALSE)
}

.memory_probe <- function(handle, ...) {
  if (!is.null(handle$gds)) {
    x <- handle$gds
    arrs <- x$assays
    dims <- dim(arrs[[1L]])
    out <- list(
      assays = names(arrs),
      dims = gds_dims(sample = dims[1L], subject = dims[2L], contrast = dims[3L]),
      subjects = x$subjects,
      contrasts = x$contrasts,
      space = x$space,
      maps = x$metadata$map_families %||% list(),
      metadata = x$metadata,
      columns = list(),
      col_data = x$col_data
    )
    return(probe_contract(out))
  }
  if (!is.null(handle$arrays)) {
    arrs <- handle$arrays
    dims <- dim(arrs[[1L]])
    out <- list(
      assays = names(arrs),
      dims = gds_dims(sample = dims[1L], subject = dims[2L], contrast = dims[3L]),
      subjects = as.character(seq_len(dims[2L])),
      contrasts = paste0("contrast", seq_len(dims[3L])),
      space = space_sample_labels(labels = as.character(seq_len(dims[1L]))),
      maps = list(),
      metadata = list(schema_version = "0.1.0", source = "<memory>", sample_labels_synthetic = TRUE),
      columns = list()
    )
    return(probe_contract(out))
  }
  stop("memory adapter: invalid handle for probe", call. = FALSE)
}

.memory_read <- function(handle,
                         assays,
                         block = NULL,
                         ...) {
  arrs <- if (!is.null(handle$gds)) handle$gds$assays else handle$arrays
  if (is.null(arrs)) stop("memory adapter: missing arrays", call. = FALSE)
  first <- arrs[[1L]]
  dims <- dim(first)
  # Build indices
  samples_idx <- .normalize_block_index(block$sample %||% NULL, dims[1L], "sample")
  subject_idx <- .normalize_block_index(block$subject %||% NULL, dims[2L], "subject")
  contrast_idx <- .normalize_block_index(block$contrast %||% NULL, dims[3L], "contrast")

  out <- lapply(assays, function(name) {
    a <- arrs[[name]]
    if (is.null(a)) stop("Assay not found in memory source: ", name, call. = FALSE)
    a[samples_idx, subject_idx, contrast_idx, drop = FALSE]
  })
  names(out) <- assays
  out
}

.memory_close <- function(handle) {
  invisible(NULL)
}
