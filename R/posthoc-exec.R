# Post-hoc execution ----------------------------------------------------------

apply_posthoc <- function(node, arrays) {
  name <- node$method
  ph <- get_posthoc(name)
  if (is.null(ph)) stop("Unknown post-hoc method: ", name, call. = FALSE)
  req <- ph$requires %||% character(0)
  if (length(req)) {
    arrays <- .ensure_required_arrays(arrays, req)
  }
  res <- ph$fun(arrays, node$options %||% list())
  if (!is.list(res) || is.null(names(res))) {
    stop("Post-hoc function must return a named list of arrays", call. = FALSE)
  }
  for (nm in names(res)) {
    arrays[[nm]] <- res[[nm]]
  }
  list(arrays = arrays)
}

