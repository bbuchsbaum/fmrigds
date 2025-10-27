#' Record a derivation operation in a plan
#'
#' @param x Plan or source
#' @param what Character vector of assays to derive
#' @param options Optional named list of options
#'
#' @return Updated plan
#' @export
derive <- function(x, what = c("var", "se", "t", "z"), options = list()) {
  plan <- as_plan(x)
  what <- unique(as.character(what))
  add_op(plan, op_derive(what, options))
}
