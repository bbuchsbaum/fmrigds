#' Construct a lazy GDS plan
#'
#' @param source A [`gds_source`] object
#' @param nodes List of operation nodes
#' @param meta Metadata list (planner bookkeeping)
#'
#' @return A plan object (`gds_plan`)
#' @name gds_plan
#' @export
gds_plan <- function(source, nodes = list(), meta = list()) {
  if (!inherits(source, "gds_source")) {
    stop("`source` must be a gds_source", call. = FALSE)
  }
  plan <- structure(list(source = source, nodes = nodes, meta = meta), class = "gds_plan")
  .ensure_plan_node_ids(plan)
}

#' Define a plan source (adapter binding)
#'
#' @param adapter Adapter identifier string
#' @param source_spec Source specification (paths, handles, etc.)
#' @param probe_result Optional probe metadata
#'
#' @return Source descriptor
#' @name gds_source
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

#' Create a subset operation node
#'
#' @param sample Sample indices or labels
#' @param subject Subject indices or IDs
#' @param contrast Contrast indices or names
#'
#' @return Operation node list
#' @export
op_subset_axis <- function(sample = NULL, subject = NULL, contrast = NULL) {
  list(op = "subset_axis", sample = sample, subject = subject, contrast = contrast)
}

#' Create a derive operation node
#'
#' @param what Character vector of statistics to derive (e.g., "t", "z", "p")
#' @param options Optional list of derivation options
#'
#' @return Operation node list
#' @export
op_derive <- function(what, options = list()) {
  list(op = "derive", what = as.character(what), options = options)
}

#' Create an align-to-group operation node
#'
#' @param family MapFamily object for alignment
#' @param family_name Name of registered map family
#'
#' @return Operation node list
#' @export
op_align_to_group <- function(family = NULL, family_name = NULL) {
  if (!is.null(family) && is.null(family_name)) {
    family_name <- family$name %||% NULL
  }
  list(op = "align_to_group", family = family, family_name = family_name)
}

#' Create a mask policy operation node
#'
#' @param policy MaskPolicy object
#'
#' @return Operation node list
#' @export
op_mask_policy <- function(policy) {
  list(op = "mask_policy", policy = policy)
}

#' Create a map operation node
#'
#' @param target_space Target space descriptor
#' @param map Mapping operator (matrix or function)
#' @param uncertainty UncertaintyRule for variance propagation
#' @param combine Optional combine method
#'
#' @return Operation node list
#' @export
op_map <- function(target_space, map, uncertainty, combine = NULL) {
  list(op = "map", target_space = target_space, map = map, uncertainty = uncertainty, combine = combine)
}

#' Create a reduce operation node
#'
#' @param method Reduction method (e.g., "fixed", "random")
#' @param weights Optional weight vector or method
#' @param by Optional grouping variable
#' @param options Optional reduction options
#' @param formula Optional model formula for meta-regression
#'
#' @return Operation node list
#' @export
op_reduce <- function(method, weights, by, options = list(), formula = NULL) {
  list(op = "reduce", method = method, weights = weights, by = by, options = options, formula = formula)
}

#' Create a write operation node
#'
#' @param path Output file path
#' @param format Output format (e.g., "h5", "csv", "parquet")
#' @param options Optional format-specific options
#'
#' @return Operation node list
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
#' @name as_plan
#' @export
as_plan <- function(x) UseMethod("as_plan")

#' @export
as_plan.gds_plan <- function(x) x

#' @export
as_plan.gds_source <- function(x) gds_plan(x)

#' @export
as_plan.default <- function(x) {
  stop("Cannot convert object to plan", call. = FALSE)
}

#' @export
as_plan.gds <- function(x) {
  stopifnot(inherits(x, "gds"))
  # Build probe_result directly from realized GDS
  arrs <- assays(x)
  dims <- dim(arrs[[1L]])
  probe_result <- list(
    assays = names(arrs),
    dims = gds_dims(sample = dims[1L], subject = dims[2L], contrast = dims[3L]),
    subjects = x$subjects,
    contrasts = x$contrasts,
    space = x$space,
    maps = x$metadata$map_families %||% list(),
    metadata = x$metadata,
    columns = list(),
    col_data = x$col_data,
    row_data = x$row_data,
    contrast_data = contrast_data(x)
  )
  src <- gds_source("memory", x, probe_result)
  plan <- gds_plan(source = src)
  plan$meta$adapter_columns <- list()
  plan$meta$map_families <- probe_result$maps %||% list()
  plan$meta$subjects <- probe_result$subjects
  if (!is.null(x$col_data)) plan$meta$col_data <- x$col_data
  if (!is.null(x$row_data)) plan$meta$row_data <- x$row_data
  if (!is.null(probe_result$contrast_data)) plan$meta$contrast_data <- probe_result$contrast_data
  # Legacy-friendly alias for tests expecting plan$metadata$dims
  plan$metadata <- list(dims = probe_result$dims)
  plan
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
  plan <- .ensure_plan_node_ids(plan)
  if (is.null(node$node_id)) {
    node$node_id <- .make_plan_node_id(plan, node, length(plan$nodes) + 1L)
  }
  existing <- vapply(plan$nodes, function(x) x$node_id %||% "", character(1))
  if (node$node_id %in% existing) {
    stop("Plan node IDs must be unique: ", node$node_id, call. = FALSE)
  }
  plan$nodes <- c(plan$nodes, list(node))
  plan
}

.node_identity_payload <- function(node) {
  node[names(node) != "node_id"]
}

.make_plan_node_id <- function(plan, node, index) {
  hash <- digest::digest(
    list(
      source = plan$source$hash,
      index = as.integer(index),
      node = .node_identity_payload(node)
    ),
    algo = "xxhash64"
  )
  sprintf("node-%04d-%s", as.integer(index), substr(hash, 1L, 12L))
}

.ensure_plan_node_ids <- function(plan) {
  if (!inherits(plan, "gds_plan") || !length(plan$nodes)) return(plan)
  for (i in seq_along(plan$nodes)) {
    if (is.null(plan$nodes[[i]]$node_id) || !nzchar(plan$nodes[[i]]$node_id)) {
      plan$nodes[[i]]$node_id <- .make_plan_node_id(plan, plan$nodes[[i]], i)
    }
  }
  ids <- vapply(plan$nodes, function(x) x$node_id, character(1))
  if (anyDuplicated(ids)) {
    stop("Plan node IDs must be unique.", call. = FALSE)
  }
  plan
}

# Lazy verbs ---------------------------------------------------------------

#' @export
subset.gds_plan <- function(x, sample = NULL, subject = NULL, contrast = NULL, ...) {
  add_op(x, op_subset_axis(sample, subject, contrast))
}
