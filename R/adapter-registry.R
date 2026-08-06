# Adapter registry ---------------------------------------------------------

.adapter_registry <- new.env(parent = emptyenv())

#' Register a storage adapter
#'
#' @param name Adapter name
#' @param detect Function (source) -> score in \code{[0,1]} or FALSE
#' @param open Function (source, ...) -> handle
#' @param probe Function (handle, ...) -> list(metadata)
#' @param read Function (handle, assays, block = NULL, ...) -> named list of arrays
#' @param close Function (handle) -> NULL
#' @param capabilities Optional scan-capability declaration. Supported fields
#'   are `sample_blocks`, `subject_blocks`, `contrast_blocks`,
#'   `persistent_handle`, `preferred_axis`, `cheap_revisit`, and
#'   `requires_staging`. Omitted fields default conservatively.
#' @param scan Optional adapter-specific scan callback.
#' @param ... Additional elements stored with adapter entry
#'
#' @return Invisibly, the adapter entry
#' @export
register_adapter <- function(name,
                             detect,
                             open,
                             probe,
                             read,
                             close,
                             capabilities = NULL,
                             scan = NULL,
                             ...) {
  capabilities <- .normalize_adapter_capabilities(capabilities)
  if (!is.null(scan) && !is.function(scan)) {
    stop("`scan` must be NULL or a function.", call. = FALSE)
  }
  .adapter_registry[[name]] <- list(
    name = name,
    detect = detect,
    open = open,
    probe = probe,
    read = read,
    close = close,
    capabilities = capabilities,
    scan = scan,
    ...
  )
  invisible(.adapter_registry[[name]])
}

.adapter_capability_defaults <- function() {
  list(
    sample_blocks = FALSE,
    subject_blocks = FALSE,
    contrast_blocks = FALSE,
    persistent_handle = FALSE,
    preferred_axis = NA_character_,
    cheap_revisit = FALSE,
    requires_staging = FALSE
  )
}

.normalize_adapter_capabilities <- function(capabilities = NULL) {
  defaults <- .adapter_capability_defaults()
  if (is.null(capabilities)) return(defaults)
  if (!is.list(capabilities)) {
    stop("Adapter capabilities must be a list.", call. = FALSE)
  }
  unknown <- setdiff(names(capabilities), names(defaults))
  if (length(unknown)) {
    stop("Unknown adapter capabilities: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- utils::modifyList(defaults, capabilities)
  logical_fields <- setdiff(names(defaults), "preferred_axis")
  for (field in logical_fields) {
    if (!is.logical(out[[field]]) || length(out[[field]]) != 1L || is.na(out[[field]])) {
      stop("Adapter capability ", field, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.na(out$preferred_axis) &&
      (!is.character(out$preferred_axis) || length(out$preferred_axis) != 1L ||
       !out$preferred_axis %in% c("sample", "subject", "contrast"))) {
    stop("Adapter preferred_axis must be sample, subject, contrast, or NA.", call. = FALSE)
  }
  out
}

.adapter_capabilities <- function(adapter) {
  if (is.character(adapter)) adapter <- get_adapter(adapter)
  .normalize_adapter_capabilities(adapter$capabilities %||% NULL)
}

.normalize_block_index <- function(index, n, axis) {
  if (is.null(index)) return(seq_len(n))
  if (is.logical(index)) {
    if (length(index) != n) {
      stop("Logical ", axis, " block index must have length ", n, ".", call. = FALSE)
    }
    index <- which(index)
  }
  if (!is.numeric(index) || anyNA(index) || any(!is.finite(index)) ||
      any(index != as.integer(index)) || any(index < 1L) || any(index > n)) {
    stop("Invalid ", axis, " block index.", call. = FALSE)
  }
  as.integer(index)
}

#' Retrieve a registered adapter
#'
#' @param name Adapter name
#'
#' @return Adapter list
#' @export
get_adapter <- function(name) {
  adapter <- .adapter_registry[[name]]
  if (is.null(adapter)) stop("Adapter not found: ", name, call. = FALSE)
  adapter
}

#' Detect the best adapter for a source
#'
#' @param source Source specification (path, list, etc.)
#' @param prefer Optional adapter name to prefer when available
#'
#' @return Adapter name
#' @export
detect_adapter <- function(source, prefer = NULL) {
  names <- ls(.adapter_registry)
  if (!length(names)) stop("No adapters registered", call. = FALSE)

  scores <- vapply(names, function(name) {
    adapter <- .adapter_registry[[name]]
    score <- tryCatch(adapter$detect(source), error = function(e) 0)
    if (isFALSE(score) || is.null(score)) 0 else as.numeric(score)
  }, numeric(1))

  if (!is.null(prefer) && prefer %in% names(scores) && scores[prefer] > 0) {
    return(prefer)
  }

  best <- names[which.max(scores)]
  if (scores[best] <= 0) stop("No adapter detected for source", call. = FALSE)
  best
}
