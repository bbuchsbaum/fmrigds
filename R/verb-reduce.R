#' Reduce across subjects (meta-analysis)
#'
#' @param x A [`gds_plan`], [`gds_source`], or realised [`gds`]
#' @param method Reduction method. Built-ins include `"fixed"`, `"random"`,
#'   `"stouffer"`, and `"fisher"`. Registry-backed reducers include
#'   meta-analytic regression reducers such as `"meta:fe_reg"` and the
#'   restricted repeated-measures Gaussian LMM reducers `"lmm:ri"`,
#'   `"lmm:ri_slope1"`, `"lmm:ri_knownvar"`, and
#'   `"lmm:ri_slope1_knownvar"`, plus permutation reducers
#'   `"perm:onesample"` and `"perm:twosample"`.
#' @param weights Weighting scheme (`"1/var"`, `"n_eff"`, `"equal"`, `"custom"`).
#'   When omitted, permutation reducers default to `"equal"`; other reducers
#'   default to `"1/var"`.
#' @param by Grouping variable (e.g., `"contrast"`)
#' @param formula Optional model formula. For meta-regression reducers the design
#'   is built from subject-level `col_data`. For repeated-measures LMM reducers,
#'   the formula specifies fixed effects at the
#'   subject-by-repeat level and must not contain `lmer`-style random-effects
#'   syntax.
#' @param data Optional subject-keyed data frame to attach as `col_data` before
#'   execution-time construction of the model matrix.
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
#' - The corresponding `*_knownvar` methods add the known diagonal sampling
#'   variance from `var` plus an estimated residual heterogeneity component.
#' - Repeated-measure metadata must live on the contrast axis via [with_contrast_data()]
#'   or tabular ingestion with `gds(..., contrast_data_cols = ...)`.
#' - Unweighted methods allow `theta_mode = "pooled"` or `"voxelwise"`;
#'   known-variance methods currently require `"voxelwise"`.
#' - The supported family is intentionally narrow: one grouping factor only,
#'   Gaussian responses only, and no general random-effects formula grammar.
#'
#' Permutation reducer contract:
#' - `"perm:onesample"` performs a one-sample sign-flip t test for each sample
#'   and contrast. `weights = "equal"` preserves the ordinary one-sample t
#'   statistic. Fixed positive `"1/var"`, `"n_eff"`, or `"custom"` weights use
#'   the reliability-weighted mean and variance, with Kish effective sample size
#'   for the standard error and degrees of freedom. Weights are held fixed across
#'   sign flips; multiplying every weight by a positive constant does not change
#'   the result. Inverse-variance weighting requires a genuine `var` or `se`
#'   assay, `n_eff` weighting requires an `n_eff` assay, and custom weighting
#'   requires `options$custom_weights` as a matching 3-D array or one value per
#'   subject.
#' - `"perm:twosample"` performs an unweighted two-sample label-permutation
#'   t test. It can infer the tested two-level group column from `formula`
#'   and `col_data`, or use `options$group` directly. It currently accepts only
#'   `weights = "equal"`.
#' - Common options include `n_perm`, `seed`, `alternative`, and for
#'   `"perm:twosample"` `variance = "welch"` or `"pooled"`.
#' - Outputs include `t_g`, parametric `p_g`, permutation `p_perm`, and
#'   tail-matched max-statistic family-wise `p_fwer`.
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
#' - `"meta:fe_reg"`/`"meta:re_reg"`: per-term `coef:<term>`,
#'   `se_coef:<term>`, normal-approximation `z_coef:<term>` and
#'   `p_coef:<term>`, plus `Q`, `df_res` (random-effects also adds `tau2`).
#'   `"meta:re_reg"`
#'   uses the DerSimonian--Laird meta-regression denominator
#'   \eqn{\mathrm{tr}\{W - WX(X^T W X)^{-1}X^T W\}} and excludes rows with
#'   non-finite effects/design values or non-finite, nonpositive variances.
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
  weights_missing <- missing(weights)
  # Allow registry-backed methods in addition to legacy names
  method <- as.character(method)[1]
  if (weights_missing && method %in% c("perm:onesample", "perm:twosample")) {
    weights <- "equal"
  }
  # Support alias: ivw -> 1/var
  w_in <- as.character(weights)[1]
  if (identical(w_in, "ivw")) w_in <- "1/var"
  weights <- match.arg(w_in, c("1/var", "n_eff", "equal", "custom"))
  if (identical(method, "perm:twosample") && !identical(weights, "equal")) {
    stop("perm:twosample currently supports only weights = \"equal\"", call. = FALSE)
  }
  if (!is.list(options)) {
    stop("`options` must be a list", call. = FALSE)
  }
  plan <- as_plan(x)
  reducer <- get_reducer(.normalize_reducer_name(method))

  # Guard: variance-weighted reducers on a synthetic unit-variance placeholder
  # (beta/stat maps ingested without standard errors) would yield meaningless
  # group standard errors. Refuse rather than silently mislead.
  synthetic_var <- isTRUE(if (inherits(x, "gds")) x$metadata$synthetic_var else NULL) ||
    isTRUE(tryCatch(metadata(plan)$synthetic_var, error = function(e) NULL)) ||
    isTRUE(plan$source$probe$metadata$synthetic_var)
  if (synthetic_var) {
    red <- reducer
    # For a registered reducer, its declared inputs are authoritative: an
    # unweighted reducer such as ols:voxelwise does not consume `var` and must
    # not be blocked (the default weights = "1/var" is ignored by it). Only fall
    # back to the weight scheme for unknown/legacy reducers.
    needs_var <- if (!is.null(red)) {
      "var" %in% (red$requires %||% character()) ||
        (identical(red$name, "perm:onesample") && identical(weights, "1/var"))
    } else {
      weights %in% c("1/var", "n_eff")
    }
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
  if (!is.null(data)) {
    plan <- .attach_reducer_data(plan, data)
  }
  formula_spec <- NULL
  if (!is.null(formula)) {
    f <- .formula_from_spec(formula)
    if (!is_observation_level && !is.null(reducer$model_contract) &&
        !isTRUE(reducer$model_contract$uses_X)) {
      stop(
        "Reducer '", .normalize_reducer_name(method),
        "' does not consume a design matrix; remove `formula`.",
        call. = FALSE
      )
    }
    formula_spec <- paste(deparse(f), collapse = " ")
  }
  add_op(
    plan,
    op_reduce(
      method,
      weights,
      by = by,
      options = options,
      formula = formula_spec
    )
  )
}

.attach_reducer_data <- function(plan, data) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  current_subjects <- as.character(.plan_subjects_for_model_matrix(plan))
  if (!length(current_subjects)) {
    stop("Cannot attach reducer data without a known subject axis.", call. = FALSE)
  }
  if ("subject" %in% names(data)) {
    ids <- as.character(data$subject)
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("`data$subject` must contain unique non-missing subject IDs.", call. = FALSE)
    }
    rownames(data) <- ids
  } else if (is.null(rownames(data)) || !all(current_subjects %in% rownames(data))) {
    if (nrow(data) != length(current_subjects)) {
      stop(
        "`data` must have matching rownames, a `subject` column, or one row per current subject.",
        call. = FALSE
      )
    }
    rownames(data) <- current_subjects
  }
  plan$meta$col_data <- .align_col_data_for_subjects(
    data,
    current_subjects,
    warn_extra = FALSE,
    context = "reduce"
  )
  plan
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
