#' Reduce across subjects (meta-analysis)
#'
#' @param x A [`gds_plan`], [`gds_source`], or realised [`gds`]
#' @param method Reduction method. Built-ins include `"fixed"`, `"random"`,
#'   `"stouffer"`, and `"fisher"`. Registry-backed reducers include
#'   meta-analytic regression reducers such as `"meta:fe_reg"` and the
#'   restricted repeated-measures Gaussian LMM reducers `"lmm:ri"` and
#'   `"lmm:ri_slope1"`, plus permutation reducers `"perm:onesample"` and
#'   `"perm:twosample"`.
#' @param weights Weighting scheme (`"1/var"`, `"n_eff"`, `"equal"`, `"custom"`)
#' @param by Grouping variable (e.g., `"contrast"`)
#' @param formula Optional model formula. For meta-regression reducers the design
#'   is built from subject-level `col_data`. For `"lmm:ri"` and
#'   `"lmm:ri_slope1"`, the formula specifies fixed effects at the
#'   subject-by-repeat level and must not contain `lmer`-style random-effects
#'   syntax.
#' @param data Optional data frame for building `X` when `x` is not a realised GDS.
#'
#' @details
#' `reduce()` always works on a lazy plan internally, so `x` is whichever stage
#' of the fmrigds workflow you already have:
#'
#' - Start from files or another external source with [gds()]. That returns a
#'   [`gds_plan`] you can pipe directly into `reduce()`.
#' - Start from an in-memory result with a realised [`gds`] returned by
#'   [compute()] or created directly with [new_gds()]. `reduce()` will convert
#'   it with [as_plan()] for you.
#' - You can also pass a low-level [`gds_source`] created by [gds_source()], but
#'   most users do not need this because [gds()] creates the source binding
#'   automatically.
#'
#' If you want that conversion to be explicit, use [as_plan()] or its alias
#' [plan()] before calling `reduce()`.
#'
#' For worked examples, see `vignette("fmrigds")` for the basic source -> plan
#' -> compute workflow and `vignette("as-plan-and-spatial-fdr")` for chaining
#' verbs on realised GDS objects.
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
#'
#' Permutation reducer contract:
#' - `"perm:onesample"` performs an unweighted one-sample sign-flip t test
#'   for each sample and contrast.
#' - `"perm:twosample"` performs an unweighted two-sample label-permutation
#'   t test. It can infer the tested two-level group column from `formula`
#'   and `col_data`, or use `options$group` directly.
#' - Common options include `n_perm`, `seed`, `alternative`, and for
#'   `"perm:twosample"` `variance = "welch"` or `"pooled"`.
#' - Outputs include `t_g`, parametric `p_g`, permutation `p_perm`, and
#'   max-|t| family-wise `p_fwer`.
#' @section Output assays:
#' Reducers collapse the subject axis and write group-level assays (named on the
#' realised GDS, retrievable with [assays()]/[assay()]):
#' - `"ols:voxelwise"`: per-term `coef:<term>`, `se_coef:<term>`,
#'   `t_coef:<term>`, `p_coef:<term>`, plus `sigma2`, `df_res`, `n_obs`
#'   (see [reducer-ols-voxelwise]).
#' - `"fixed"`/`"meta:fe"` and `"random"`/`"meta:re"`: `beta_g`, `var_g`,
#'   `se_g`, `z_g`, `p_g`, `Q`, `I2`, `n_eff` (random-effects also adds
#'   `tau2`). `n_eff` is the number of finite subject-level inputs used for
#'   each sample and contrast.
#' - `"meta:fe_reg"`/`"meta:re_reg"`: per-term `coef:<term>`/`se_coef:<term>`
#'   plus `Q`, `df_res` (random-effects also adds `tau2`).
#' - `"stouffer"`/`"combine:stouffer"`: `z_g`, `p_g`.
#' - `"fisher"`/`"combine:fisher"` and `"combine:lancaster"`: `chi2`, `df`,
#'   `p_g`.
#' - `"perm:onesample"`/`"perm:twosample"`: `beta_g`, `se_g`, `t_g`, `df`,
#'   `p_g`, `p_perm`, `p_fwer`.
#'
#' @return Updated plan
#' @seealso [gds()], [compute()], [as_plan()], [plan()], [new_gds()], [gds_source()]
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

  # Guard: variance-weighted reducers on a synthetic unit-variance placeholder
  # (beta/stat maps ingested without standard errors) would yield meaningless
  # group standard errors. Refuse rather than silently mislead.
  synthetic_var <- isTRUE(if (inherits(x, "gds")) x$metadata$synthetic_var else NULL) ||
    isTRUE(tryCatch(metadata(plan)$synthetic_var, error = function(e) NULL))
  if (synthetic_var) {
    red <- get_reducer(.normalize_reducer_name(method))
    needs_var <- (!is.null(red) && "var" %in% (red$requires %||% character())) ||
      weights %in% c("1/var", "n_eff")
    if (needs_var) {
      stop(sprintf(
        paste0(
          "Reducer '%s' is variance-weighted, but this GDS carries only a synthetic ",
          "unit-variance placeholder (beta/stat maps ingested without standard errors), ",
          "so the group standard errors would be meaningless. Use an unweighted reducer ",
          "such as method = \"ols:voxelwise\" (e.g. one_sample()/group_ols()), or supply ",
          "real standard errors via nifti_source(se = ...)."
        ), method), call. = FALSE)
    }
  }

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
