#' Reduce across subjects (meta-analysis)
#'
#' @param x Plan, source, or realised GDS
#' @param method Reduction method. Built-ins: "fixed", "random", "stouffer", "fisher";
#'   or any registered reducer id like "meta:fe", "meta:re", "combine:stouffer", etc.
#' @param weights Weighting scheme (`"1/var"`, `"n_eff"`, `"equal"`, `"custom"`)
#' @param by Grouping variable (e.g., `"contrast"`)
#' @param formula Optional model formula for meta-regression; if provided and
#'   `x` is a realised GDS with `col_data`, an `X` matrix is built automatically
#'   and passed via `options$X`.
#' @param data Optional data frame for building `X` when `x` is not a realised GDS.
#' @param options Additional options, e.g., `tau2 = NULL`, `custom_weights`
#'
#' Note: Current grouping is per-contrast. Values other than "contrast" are
#' not interpreted yet; future versions may extend grouping to other
#' factors in `row_data`/`col_data`.
#' @return Updated plan
#' @export
reduce <- function(x,
                   method = c("fixed", "random", "stouffer", "fisher"),
                   weights = c("1/var", "n_eff", "equal", "custom"),
                   by = "contrast",
                   options = list(),
                   formula = NULL,
                   data = NULL) {
  # Allow registry-backed methods in addition to legacy names
  method <- as.character(method)[1]
  # Support alias: ivw -> 1/var
  w_in <- as.character(weights)[1]
  if (identical(w_in, "ivw")) w_in <- "1/var"
  weights <- match.arg(w_in, c("1/var", "n_eff", "equal", "custom"))
  # Auto-build X from formula and col_data/data if provided
  if (!is.null(formula)) {
    f <- if (inherits(formula, "formula")) formula else stats::as.formula(formula)
    cd <- NULL
    if (inherits(x, "gds")) cd <- x$col_data
    if (is.null(cd) && !is.null(data)) cd <- data
    if (!is.null(cd)) {
      options$X <- stats::model.matrix(f, data = cd)
      options$formula <- deparse(f)
    } else {
      options$formula <- deparse(f)
    }
  }
  plan <- as_plan(x)
  add_op(plan, op_reduce(method, weights, by = by, options = options, formula = options$formula %||% NULL))
}

#' Eagerly reduce across subjects and compute immediately
#'
#' @param x Plan, source, or realized GDS
#' @param ... Arguments passed to reduce(), then compute()
#'
#' @return A realized GDS object
#' @export
reduce_eager <- function(x, ...) {
  compute(reduce(x, ...))
}
