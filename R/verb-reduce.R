#' Reduce across subjects (meta-analysis)
#'
#' @param x A [`gds_plan`][as_plan()], [`gds_source`], or realised [`gds`]
#' @param method Reduction method. Built-ins include `"fixed"`, `"random"`,
#'   `"stouffer"`, and `"fisher"`. Registry-backed reducers include
#'   meta-analytic regression reducers such as `"meta:fe_reg"` and the
#'   restricted repeated-measures Gaussian LMM reducers `"lmm:ri"`,
#'   `"lmm:ri_slope1"`, `"lmm:ri_knownvar"`, and
#'   `"lmm:ri_slope1_knownvar"`, the parametric two-sample Welch test
#'   `"ttest:welch"`, plus permutation reducers
#'   `"perm:onesample"` and `"perm:twosample"`.
#' @param weights Weighting scheme (`"1/var"`, `"n_eff"`, `"equal"`, `"custom"`).
#'   Registered reducers apply their own weighting, so only the default
#'   `"1/var"` is accepted for every reducer. `"equal"` is accepted for
#'   unweighted and evidence-combining reducers (e.g. `"ols:voxelwise"`,
#'   `"stouffer"`); other schemes raise an error unless the reducer lists them in
#'   `options_schema$weights`. Unknown `method` names also raise an error.
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
#'   [`gds_plan`][as_plan()] you can pipe directly into `reduce()`.
#' - Start from an in-memory result with a realised [`gds`] returned by
#'   [compute()] or created directly with [new_gds()]. `reduce()` will convert
#'   it with [as_plan()] for you.
#' - You can also pass a low-level [`gds_source`] created by [gds_source()], but
#'   most users do not need this because [gds()] creates the source binding
#'   automatically.
#'
#' If you want that conversion to be explicit, use [as_plan()]
#' before calling `reduce()`.
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
#' - `"perm:onesample"` performs an unweighted one-sample sign-flip t test
#'   for each sample and contrast. It does not accept subject weights: omit
#'   `weights` or pass `weights = "equal"`; an explicit `"1/var"`, `"n_eff"` or
#'   `"custom"` is refused rather than silently ignored (#21).
#' - `"perm:twosample"` performs an unweighted two-sample label-permutation
#'   t test. It can infer the tested two-level group column from `formula`
#'   and `col_data`, or use `options$group` directly.
#' - Common options include `n_perm`, `seed`, `alternative`, and for
#'   `"perm:twosample"` `variance = "welch"` or `"pooled"`.
#' - Outputs include `t_g`, parametric `p_g`, permutation `p_perm`, and
#'   max-statistic family-wise `p_fwer`. The maximum is taken in the direction
#'   of `alternative` (max t for `"greater"`, max -t for `"less"`, max |t| for
#'   `"two.sided"`), so `p_fwer >= p_perm` always holds (#22).
#'
#' Spatial (TFCE and cluster) inference with `"perm:onesample"`:
#' - On a voxel space ([space_voxel()]), `options$spatial = "tfce"`,
#'   `"cluster"` or `c("tfce", "cluster")` adds whole-brain family-wise
#'   inference from the same subject-level sign flips. Every null
#'   distribution (max-T, max-TFCE, largest cluster mass and extent) comes from
#'   one sign matrix in one pass, so the corrected maps are mutually
#'   consistent. Each sign flip multiplies a whole subject map by +1 or -1,
#'   preserving its spatial correlation (unlike voxel-wise sign flipping).
#' - The spatial statistic is the z-equivalent of each voxel's t
#'   (`options$spatial_stat = "z"`, the default; same one-tail probability
#'   with df = n - 1, so voxels with different subject counts are on one
#'   scale) or the t statistic itself (`"t"`, as in FSL randomise).
#' - TFCE: `tfce_H` (default 2), `tfce_E` (0.5) and `tfce_dh` (0.1); the
#'   integral is evaluated at heights `dh, 2 dh, ...` up to the map maximum,
#'   exactly as in [neurothresh::tfce_transform()].
#' - Clusters: the cluster-forming threshold is `options$cluster_threshold` on
#'   the scale of the spatial statistic, or `options$cluster_p` (default
#'   0.001), a one-tailed p converted with that scale (`qnorm` for z, `qt` with
#'   df = n_subjects - 1 for t). Clusters are tested by excess mass, the sum of
#'   `stat - threshold` over the cluster (`cluster_stat = "mass"`, default),
#'   or by extent (`"extent"`).
#' - `options$connectivity` is 6 (default), 18 or 26 (26 is FSL's choice).
#' - Tails: `"greater"` uses positive clusters/TFCE, `"less"` negative ones,
#'   `"two.sided"` both; the two-sided null takes each permutation's maximum
#'   over both tails.
#' - Missing data: subjects with non-finite values are dropped voxel by voxel.
#'   A voxel with fewer than `options$min_subjects` (default 2) finite
#'   subjects, or zero variance, has no statistic: it is excluded from the
#'   mask for TFCE and clustering (it never exceeds a threshold and does not
#'   connect its neighbours), and its spatial assays are `NA`.
#' - p-values are `(1 + #{null >= observed}) / (1 + n_null)`; the observed
#'   (all +1) sign vector is never part of the null. With generated signs,
#'   `n_perm` null flips are drawn; a user-supplied `options$signs` whose first
#'   row is all +1 has that row excluded.
#' - `options$threads` sets the OpenMP team for this call (default: the
#'   package setting, see `options(fmrigds.threads)` / `FMRIGDS_THREADS` /
#'   `OMP_NUM_THREADS`). Results do not depend on the thread count.
#' - Extra assays: `z_g` (z-equivalent of `t_g`), `n_eff` (finite subjects),
#'   and with `"tfce"` `tfce` (TFCE score in the tested tail, signed like
#'   `t_g`) and `p_tfce_fwer`; with `"cluster"` `cluster_id` (0 outside
#'   clusters; clusters numbered by decreasing mass) and `p_cluster_fwer`
#'   (the cluster's FWER p, 1 outside clusters). All export with
#'   [write_nifti_assays()].
#' - The attachment `metadata(x)$attachments[["reduce/perm:onesample/contrast=<c>"]]`
#'   holds the null maxima (`null$max_t`, `null$tfce`, `null$cluster_mass`,
#'   `null$cluster_extent`), the settings, and a cluster table with columns
#'   `cluster_id`, `sign`, `size`, `mass`, `p_fwer_mass`, `p_fwer_extent`,
#'   `peak_stat`, `peak_t`, `peak_sample`, grid indices `peak_i/j/k` and
#'   world coordinates `peak_x/y/z` (from the space affine).
#'
#' Welch two-sample t test (`"ttest:welch"`):
#' - Unweighted, per sample and contrast, identical to
#'   `stats::t.test(y[g == 2], y[g == 1], var.equal = FALSE)` where `g` is the
#'   two-level group: the effect is second level minus first level. The group
#'   is the last (or `options$term`) column of the design built from `formula`
#'   (e.g. `~ group`), or `options$group` directly.
#' - Non-finite subjects are dropped per sample; samples with fewer than
#'   `options$min_group` (default and minimum 2) finite subjects in either
#'   group are `NA`. `options$alternative` sets the tail of `p_g`.
#' @section Variance floor:
#' The inverse-variance reducers `"meta:fe"`/`"fixed"`, `"meta:re"`/`"random"`,
#' `"meta:fe_reg"` and `"meta:re_reg"` share one validity rule and one floor,
#' in both the R and compiled kernels, which front ends may rely on:
#' - a subject enters a sample only if its effect and variance are finite and
#'   the variance is strictly positive (zero, negative, `NA` and infinite
#'   variances exclude that subject for that sample);
#' - a positive variance below `eps = 1e-12` is floored to `eps` before
#'   weighting, i.e. the weight is `1 / max(var, eps)` and never exceeds
#'   `1e12` (random-effects weights are `1 / (max(var, eps) + tau2)`).
#'
#' `options = list(eps = ...)` overrides the floor. Rescale effects whose
#' sampling variances are legitimately below `1e-12`, since the floor would
#' otherwise equalise their weights.
#'
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
#'   `p_g`, `p_perm`, `p_fwer`; with `options$spatial`, `"perm:onesample"`
#'   adds `z_g`, `n_eff`, `tfce`, `p_tfce_fwer`, `cluster_id` and
#'   `p_cluster_fwer` (see above).
#' - `"ttest:welch"`: `beta_g` (mean difference), `var_g`, `se_g`, `t_g`,
#'   Welch--Satterthwaite `df`, `p_g`, the normal-equivalent `z_g` (same
#'   p-value as `t_g`), and per-group finite counts `n0`, `n1`; `z`/`p` mirror
#'   `z_g`/`p_g`.
#'
#' @return Updated plan
#' @seealso [gds()], [compute()], [as_plan()], [new_gds()], [gds_source()]
#' @export
reduce <- function(x,
                   method = c("fixed", "random", "stouffer", "fisher"),
                   weights = c("1/var", "n_eff", "equal", "custom"),
                   by = "contrast",
                   options = list(),
                   formula = NULL,
                   data = NULL) {
  weights_explicit <- !missing(weights)
  # Allow registry-backed methods in addition to legacy names
  method <- as.character(method)[1]
  # Support alias: ivw -> 1/var
  w_in <- as.character(weights)[1]
  if (identical(w_in, "ivw")) w_in <- "1/var"
  weights <- match.arg(w_in, c("1/var", "n_eff", "equal", "custom"))
  if (!is.list(options)) {
    stop("`options` must be a list", call. = FALSE)
  }
  reducer <- get_reducer(.normalize_reducer_name(method))
  if (is.null(reducer)) {
    stop(.unknown_reducer_message(method), call. = FALSE)
  }
  .check_reducer_weights(reducer, method, weights, explicit = weights_explicit)
  # [[ ]], not $: `options$weights` partially matches an option such as
  # `weights_custom`, which would leave the requested scheme unforwarded.
  if (is.character(reducer$options_schema[["weights"]]) && is.null(options[["weights"]])) {
    options[["weights"]] <- weights
  }
  plan <- as_plan(x)

  # Guard: variance-weighted reducers on a synthetic unit-variance placeholder
  # (beta/stat maps ingested without standard errors) would yield meaningless
  # group standard errors. Refuse rather than silently mislead.
  synthetic_var <- isTRUE(if (inherits(x, "gds")) x$metadata$synthetic_var else NULL) ||
    isTRUE(tryCatch(metadata(plan)$synthetic_var, error = function(e) NULL)) ||
    isTRUE(plan$source$probe$metadata$synthetic_var)
  if (synthetic_var) {
    red <- reducer
    # The registered reducer's declared inputs are authoritative: an
    # unweighted reducer such as ols:voxelwise does not consume `var` and must
    # not be blocked (the default weights = "1/var" is ignored by it).
    needs_var <- "var" %in% (red$requires %||% character())
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

.reducer_aliases <- c(
  fixed = "meta:fe",
  random = "meta:re",
  stouffer = "combine:stouffer",
  fisher = "combine:fisher"
)

.unknown_reducer_message <- function(method) {
  sprintf(
    "Unknown reducer '%s'. Available reducers: %s (aliases: %s). Register custom reducers with register_reducer().",
    method,
    paste(list_reducers(), collapse = ", "),
    paste(sprintf("%s -> %s", names(.reducer_aliases), .reducer_aliases), collapse = ", ")
  )
}

# Weight schemes other than the default "1/var" are only honoured when the
# reducer declares them. Registered reducers never receive the `weights`
# argument directly, so silently accepting "n_eff"/"custom" (or "equal" for a
# variance-weighted reducer) would run a different analysis than requested.
#
# A reducer declares support either through its model_contract weight_mode
# ("unweighted"/"evidence" reducers are inherently equal-weighted, so
# weights = "equal" is accepted) or by listing the supported schemes in an
# `options_schema$weights` entry, in which case the scheme is forwarded as
# `options$weights`.
.check_reducer_weights <- function(reducer, method, weights, explicit = FALSE) {
  mode <- reducer$model_contract$weight_mode %||% NA_character_
  declared <- reducer$options_schema[["weights"]] %||% NULL
  if (identical(weights, "1/var")) {
    # The default "1/var" is tolerated when omitted; asking for it explicitly
    # from an unweighted reducer would silently run an unweighted analysis.
    if (isTRUE(explicit) && identical(mode, "unweighted") &&
        !(is.character(declared) && "1/var" %in% declared)) {
      stop(sprintf(
        paste0(
          "weights = \"1/var\" is not supported by reducer '%s', which is unweighted ",
          "(every subject counts equally). Omit `weights` or use weights = \"equal\"; ",
          "for inverse-variance weighting use a meta-analytic reducer such as \"random\"."
        ),
        reducer$name
      ), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  if (is.character(declared) && weights %in% declared) return(invisible(TRUE))
  if (identical(weights, "equal") && mode %in% c("unweighted", "evidence")) {
    return(invisible(TRUE))
  }
  stop(sprintf(
    paste0(
      "weights = \"%s\" is not supported by reducer '%s' (weight mode: %s); ",
      "registered reducers apply their own weighting and would silently ignore it. ",
      "Use the default weights = \"1/var\"%s."
    ),
    weights, reducer$name, if (is.na(mode)) "unspecified" else mode,
    if (identical(weights, "equal")) ", or an unweighted reducer such as \"ols:voxelwise\" or \"stouffer\"" else ""
  ), call. = FALSE)
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
#' @noRd
reduce_eager <- function(x, ...) {
  compute(reduce(x, ...))
}
