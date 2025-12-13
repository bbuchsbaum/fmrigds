# Adapter probe contract validator (internal)

#' Validate and normalize adapter probe results
#' @keywords internal
probe_contract <- function(x) {
  # Basic structure
  if (!is.list(x)) stop("Adapter probe must return a list", call. = FALSE)
  required <- c("assays", "dims", "subjects", "contrasts", "space", "maps", "metadata", "columns")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(sprintf("Probe missing fields: %s. How to fix: ensure adapter$probe returns all required fields.", paste(missing, collapse = ", ")), call. = FALSE)
  }

  # Types and shapes
  if (!is.character(x$assays) || !length(x$assays)) {
    stop("Probe field 'assays' must be a non-empty character vector", call. = FALSE)
  }
  if (!inherits(x$dims, "gds_dims") || length(x$dims) != 3L) {
    stop("Probe field 'dims' must be a gds_dims(sample, subject, contrast)", call. = FALSE)
  }
  if (!is.character(x$subjects)) stop("Probe field 'subjects' must be character", call. = FALSE)
  if (!is.character(x$contrasts)) stop("Probe field 'contrasts' must be character", call. = FALSE)
  if (!inherits(x$space, "gds_space")) stop("Probe field 'space' must inherit 'gds_space'", call. = FALSE)
  if (!is.list(x$maps)) stop("Probe field 'maps' must be a list (possibly empty)", call. = FALSE)
  if (!is.list(x$metadata)) stop("Probe field 'metadata' must be a list", call. = FALSE)
  if (!is.list(x$columns)) stop("Probe field 'columns' must be a list", call. = FALSE)
  if (!is.null(x$col_data) && !is.data.frame(x$col_data)) {
    stop("Probe optional field 'col_data' must be a data.frame or NULL", call. = FALSE)
  }

  # Length consistency
  if (length(x$subjects) != x$dims[2L]) {
    stop(sprintf("Probe subjects length (%d) must equal dims[2] (%d)", length(x$subjects), x$dims[2L]), call. = FALSE)
  }
  if (length(x$contrasts) != x$dims[3L]) {
    stop(sprintf("Probe contrasts length (%d) must equal dims[3] (%d)", length(x$contrasts), x$dims[3L]), call. = FALSE)
  }
  x
}

