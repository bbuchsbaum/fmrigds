#' Add a space transformation to a plan
#'
#' @param x Plan, source, or realised GDS
#' @param target_space Target space descriptor or identifier
#' @param map Mapping object (`map_linear` or matrix)
#' @param uncertainty [`UncertaintyRule`] instance controlling variance propagation
#' @param combine Optional evidence combiner when effect scale absent
#'
#' @return Updated plan
#' @export
map_to <- function(x,
                   target_space,
                   map,
                   uncertainty = UncertaintyRule("independent"),
                   combine = NULL) {
  plan <- as_plan(x)
  add_op(plan, op_map(target_space, map, uncertainty, combine))
}

#' @export
map_to_eager <- function(x, ...) {
  compute(map_to(x, ...))
}
