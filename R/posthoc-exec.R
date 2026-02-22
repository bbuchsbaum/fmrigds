# Post-hoc execution ----------------------------------------------------------

apply_posthoc <- function(node, arrays, context = list()) {
  name <- node$method
  ph <- get_posthoc(name)
  if (is.null(ph)) stop("Unknown post-hoc method: ", name, call. = FALSE)
  req <- ph$requires %||% character(0)
  if (length(req)) {
    arrays <- .ensure_required_arrays(arrays, req)
  }
  opts <- utils::modifyList(context %||% list(), node$options %||% list())
  raw_res <- ph$fun(arrays, opts)
  attachments <- list()

  is_structured <- is.list(raw_res) &&
    !is.null(raw_res$arrays) &&
    is.list(raw_res$arrays)

  if (is_structured) {
    res <- raw_res$arrays
    attachments <- raw_res$attachments %||% list()
    if (!is.list(attachments)) {
      stop("Post-hoc attachments must be a list", call. = FALSE)
    }
  } else {
    res <- raw_res
  }

  if (!is.list(res) || is.null(names(res))) {
    stop("Post-hoc function must return a named list of arrays", call. = FALSE)
  }
  for (nm in names(res)) {
    arrays[[nm]] <- res[[nm]]
  }
  list(arrays = arrays, attachments = attachments)
}
