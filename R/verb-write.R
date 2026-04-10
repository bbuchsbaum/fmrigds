#' Declare an output target for a plan
#'
#' Adds a lazy write operation to the plan. Data are only written when
#' [`compute()`] executes.
#'
#' @param x A plan or object coercible via [as_plan()]
#' @param path Output file path
#' @param format Output format (`"h5"`, `"csv"`, `"parquet"`, `"nifti"`)
#' @param options Optional list of format-specific options
#'
#' @return Updated plan with a write node appended
#' @export
write_out <- function(x,
                      path,
                      format = c("h5", "csv", "parquet", "nifti"),
                      options = list()) {
  plan <- as_plan(x)
  if (!is.list(options)) {
    stop("`options` must be a list", call. = FALSE)
  }
  if (missing(path) || !is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty string", call. = FALSE)
  }
  format <- match.arg(format)
  add_op(plan, op_write(path = path, format = format, options = options))
}
