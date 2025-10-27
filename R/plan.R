#' Construct a lazy GDS plan
#'
#' @param source A [`gds_source`] object
#' @param nodes List of operation nodes
#' @param meta Metadata list (planner bookkeeping)
#'
#' @return A plan object (`gds_plan`)
#' @export
gds_plan <- function(source, nodes = list(), meta = list()) {
  if (!inherits(source, "gds_source")) {
    stop("`source` must be a gds_source", call. = FALSE)
  }
  structure(list(source = source, nodes = nodes, meta = meta), class = "gds_plan")
}

#' Define a plan source (adapter binding)
#'
#' @param adapter Adapter identifier string
#' @param source_spec Source specification (paths, handles, etc.)
#' @param probe_result Optional probe metadata
#'
#' @return Source descriptor
#' @export
gds_source <- function(adapter, source_spec, probe_result = NULL) {
  if (!is.character(adapter) || length(adapter) != 1L) {
    stop("`adapter` must be a scalar character", call. = FALSE)
  }
  structure(
    list(
      adapter = adapter,
      source = source_spec,
      probe = probe_result,
      hash = digest::digest(list(adapter = adapter, source = source_spec))
    ),
    class = "gds_source"
  )
}

# Operation node constructors ----------------------------------------------

#' @export
op_subset_axis <- function(sample = NULL, subject = NULL, contrast = NULL) {
  list(op = "subset_axis", sample = sample, subject = subject, contrast = contrast)
}

#' @export
op_derive <- function(what, options = list()) {
  list(op = "derive", what = as.character(what), options = options)
}

#' @export
op_align_to_group <- function(family) {
  list(op = "align_to_group", family = family)
}

#' @export
op_mask_policy <- function(policy) {
  list(op = "mask_policy", policy = policy)
}

#' @export
op_map <- function(target_space, map, uncertainty, combine = NULL) {
  list(op = "map", target_space = target_space, map = map, uncertainty = uncertainty, combine = combine)
}

#' @export
op_reduce <- function(method, weights, by) {
  list(op = "reduce", method = method, weights = weights, by = by)
}

#' @export
op_write <- function(path, format, options = list()) {
  list(op = "write", path = path, format = format, options = options)
}

# Utilities ----------------------------------------------------------------

#' Ensure an object is a plan
#'
#' @param x Plan, source, or realized GDS
#'
#' @return A plan object
#' @export
as_plan <- function(x) {
  if (inherits(x, "gds_plan")) return(x)
  if (inherits(x, "gds_source")) return(gds_plan(x))
  stop("Cannot convert object to plan", call. = FALSE)
}

#' Append an operation to a plan
#'
#' @param plan A `gds_plan`
#' @param node Operation node (list)
#'
#' @return Updated plan
#' @export
add_op <- function(plan, node) {
  if (!inherits(plan, "gds_plan")) {
    stop("`plan` must be a gds_plan", call. = FALSE)
  }
  plan$nodes <- c(plan$nodes, list(node))
  plan
}

# Lazy verbs ---------------------------------------------------------------

#' @export
subset.gds_plan <- function(x, sample = NULL, subject = NULL, contrast = NULL, ...) {
  add_op(x, op_subset_axis(sample, subject, contrast))
}
