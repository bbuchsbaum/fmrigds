#' Add a post-hoc operation to a plan
#'
#' @param x A plan or object coercible via [as_plan()]
#' @param method Post-hoc method name (e.g., "fdr:bh", "fdr:by")
#' @param options Optional list of method-specific options
#'
#' @return Updated `gds_plan`
#' @seealso [posthoc-methods] for the built-in methods and the assays each adds.
#' @export
#' @examples
#' df <- data.frame(
#'   sample = rep(c("ROI_1","ROI_2"), each = 2),
#'   subject = rep(c("s1","s2"), times = 2),
#'   contrast = "c1",
#'   beta = c(0.5,0.6,0.7,0.8),
#'   var = c(0.04,0.05,0.06,0.07)
#' )
#' tmp <- tempfile(fileext = ".csv"); utils::write.csv(df, tmp, row.names = FALSE)
#' plan <- gds(tmp) |> reduce(method = "fixed") |> posthoc("fdr:bh")
#' # g <- compute(plan)  # returns q-values
posthoc <- function(x, method, options = list()) {
  plan <- as_plan(x)
  add_op(plan, list(op = "posthoc", method = method, options = options))
}
