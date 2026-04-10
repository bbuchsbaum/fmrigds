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
                             ...) {
  .adapter_registry[[name]] <- list(
    name = name,
    detect = detect,
    open = open,
    probe = probe,
    read = read,
    close = close,
    ...
  )
  invisible(.adapter_registry[[name]])
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
