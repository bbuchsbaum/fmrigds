#' Coerce to a lazy plan (alias)
#'
#' Alias for [as_plan()] to improve ergonomics when converting realized GDS
#' objects or sources into a `gds_plan`.
#'
#' @param x A `gds`, `gds_source`, or `gds_plan`
#'
#' @return A `gds_plan`
#' @export
plan <- function(x) {
  as_plan(x)
}

