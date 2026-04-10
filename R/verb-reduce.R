#' Reduce across subjects (meta-analysis)
#'
#' @param x Plan, source, or realised GDS
#' @param method Reduction method. Built-ins include `"fixed"`, `"random"`,
#'   `"stouffer"`, and `"fisher"`. Registry-backed reducers include
#'   meta-analytic regression reducers such as `"meta:fe_reg"` and the
#'   restricted repeated-measures Gaussian LMM reducers `"lmm:ri"` and
#'   `"lmm:ri_slope1"`.
#' @param weights Weighting scheme (`"1/var"`, `"n_eff"`, `"equal"`, `"custom"`)
#' @param by Grouping variable (e.g., `"contrast"`)
#' @param formula Optional model formula. For meta-regression reducers the design
#'   is built from subject-level `col_data`. For `"lmm:ri"` and
#'   `"lmm:ri_slope1"`, the formula specifies fixed effects at the
#'   subject-by-repeat level and must not contain `lmer`-style random-effects
#'   syntax.
#' @param data Optional data frame for building `X` when `x` is not a realised GDS.
#' @param options Additional reducer options. Common examples include
#'   `tau2 = NULL`, `custom_weights`, or for restricted LMM reducers
#'   `list(fit = "REML", theta_mode = "pooled")` plus
#'   `list(slope = "time", covariance = "diag")` for `"lmm:ri_slope1"`.
#'
#' Note: Current grouping is per-contrast. Values other than "contrast" are
#' not interpreted yet; future versions may extend grouping to other
#' factors in `row_data`/`col_data`.
#'
#' Restricted repeated-measures LMM contract:
#' - `"lmm:ri"` fits a Gaussian random-intercept model.
#' - `"lmm:ri_slope1"` fits a Gaussian random intercept plus one within-subject slope.
#' - Repeated-measure metadata must live on the contrast axis via [with_contrast_data()]
#'   or tabular ingestion with `gds(..., contrast_data_cols = ...)`.
#' - `options$theta_mode` can be `"pooled"` or `"voxelwise"`.
#' - The supported family is intentionally narrow: one grouping factor only,
#'   Gaussian responses only, and no general random-effects formula grammar.
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
  if (!is.list(options)) {
    stop("`options` must be a list", call. = FALSE)
  }
  plan <- as_plan(x)
  is_observation_level <- grepl("^lmm:", method)
  # Auto-build X from formula and col_data/data if provided
  if (!is.null(formula) && !is_observation_level) {
    f <- if (inherits(formula, "formula")) formula else stats::as.formula(formula)
    cd <- data
    if (is.null(cd) && inherits(x, "gds")) cd <- x$col_data
    if (is.null(cd)) cd <- col_data(plan)
    if (!is.null(cd)) {
      subjects <- if (inherits(x, "gds")) x$subjects else .plan_subjects_for_model_matrix(plan)
      if (!is.null(subjects) && length(subjects)) {
        if (!is.null(rownames(cd)) && length(rownames(cd))) {
          cd <- .align_col_data_for_subjects(cd, subjects, warn_extra = FALSE, context = "reduce")
        } else if ("subject" %in% names(cd)) {
          idx <- match(subjects, cd$subject)
          if (anyNA(idx)) stop("`data` is missing subjects required by the plan", call. = FALSE)
          cd <- cd[idx, , drop = FALSE]
        } else if (nrow(cd) != length(subjects)) {
          stop("`data` must have rownames, a `subject` column, or one row per subject in plan order", call. = FALSE)
        }
        if ("subject" %in% names(cd)) cd$subject <- NULL
        cd <- data.frame(subject = subjects, cd, check.names = FALSE)
      }
      options$X <- stats::model.matrix(f, data = cd)
    }
    options$formula <- deparse(f)
  } else if (!is.null(formula)) {
    f <- if (inherits(formula, "formula")) formula else stats::as.formula(formula)
    options$formula <- deparse(f)
  }
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
