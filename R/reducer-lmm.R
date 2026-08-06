#' Restricted repeated-measures Gaussian LMM reducers
#'
#' `fmrigds` includes four specialized repeated-measures Gaussian mixed-model
#' reducers that are optimized for the case where every sample shares the same
#' observation layout and the same fixed/random design.
#'
#' - `method = "lmm:ri"` fits a random-intercept model.
#' - `method = "lmm:ri_slope1"` fits a random intercept plus one within-subject
#'   slope supplied through `options$slope`.
#' - `method = "lmm:ri_knownvar"` and `"lmm:ri_slope1_knownvar"` fit the
#'   corresponding variance-aware models using the input `var` assay as known
#'   diagonal sampling variance.
#'
#' These reducers operate on the full contrast axis jointly. Repeated-measure
#' metadata therefore needs to be attached as contrast-level metadata, either
#' with [with_contrast_data()] or during tabular ingestion via
#' `gds(..., contrast_data_cols = ...)`.
#'
#' Supported contract:
#' - one grouping factor only
#' - Gaussian outcomes only
#' - shared observation layout and shared fixed/random design across samples
#' - unweighted reducers support `theta_mode = "pooled"` for one shared
#'   random-effect-to-residual covariance shape across samples or
#'   `theta_mode = "voxelwise"` for sample-specific parameters
#' - known-variance reducers currently support voxelwise covariance estimation
#'   only
#' - no crossed random effects and no general `lmer()` random-effects grammar
#'
#' @section Input variance contract:
#' The unweighted `lmm:ri` and `lmm:ri_slope1` reducers require `beta` only. An
#' input `var` or `se` assay is ignored by their likelihood. Their `sigma2` and
#' `vc_*` outputs are second-level variance components estimated from the
#' supplied effects.
#'
#' The `*_knownvar` reducers require aligned `beta` and strictly positive `var`
#' assays. For each sample, they fit marginal covariance
#' `diag(var) + vc_resid * I + Z G Z'`. Thus, `var` is treated as known sampling
#' variance, `vc_resid` (`sigma2`) is additional observation-level heterogeneity,
#' and `G` contains the subject random-effect covariance. Only diagonal sampling
#' covariance is currently accepted; correlated first-level estimates require a
#' future full sampling-covariance input rather than pretending their marginal
#' variances are independent.
#'
#' Fixed-effect standard errors use the fitted covariance as a plug-in estimate,
#' and `p_coef:*` uses residual degrees of freedom (`n - p`). Satterthwaite and
#' Kenward--Roger corrections are not applied.
#'
#' @section Output assays:
#' Both reducers expand fixed effects per model term (`<term>` = design-matrix
#' column names of `formula`): `coef:<term>`, `se_coef:<term>`, `t_coef:<term>`,
#' `p_coef:<term>`. In addition:
#' - `lmm:ri` provides `sigma2`, `vc_intercept`, `vc_resid`, `df_res`, `logLik`,
#'   `converged` (1/0), and `lambda` (random-intercept-to-residual variance
#'   ratio).
#' - `lmm:ri_slope1` provides `sigma2`, `vc_intercept`, `vc_slope`,
#'   `vc_cov_intercept_slope`, `vc_resid`, `df_res`, `logLik`, `converged`,
#'   `lambda_intercept`, `lambda_slope`, `lambda_cov_intercept_slope`, and
#'   `corr_intercept_slope`.
#' - `lmm:ri_knownvar` provides `sigma2`, `vc_intercept`, `vc_resid`, and
#'   `sampling_var_mean`; `lmm:ri_slope1_knownvar` additionally provides
#'   `vc_slope`, `vc_cov_intercept_slope`, and `corr_intercept_slope`. The
#'   relative `lambda_*` assays are intentionally omitted because known sampling
#'   variance breaks the single residual-scale parameterization.
#' Pass these names verbatim (including the colon and any parentheses) to
#' [write_nifti_assays()] / [write_out()].
#'
#' @section Usage:
#' ```r
#' rm_df <- data.frame(
#'   sample = rep(c("ROI_1", "ROI_2"), each = 12),
#'   subject = rep(paste0("sub-", sprintf("%02d", 1:4)), each = 3, times = 2),
#'   contrast = rep(c("pre", "mid", "post"), times = 8),
#'   time = rep(c(-1, 0, 1), times = 8),
#'   beta = c(
#'     0.40, 0.80, 1.20, 0.55, 0.95, 1.35, 0.65, 1.05, 1.45, 0.85, 1.25, 1.65,
#'     0.70, 1.20, 1.70, 0.85, 1.35, 1.85, 0.95, 1.45, 1.95, 1.15, 1.65, 2.15
#'   )
#' )
#'
#' g <- gds(rm_df, contrast_data_cols = "time") |>
#'   reduce(
#'     method = "lmm:ri_slope1",
#'     formula = ~ time,
#'     options = list(
#'       slope = "time",
#'       covariance = "diag",
#'       fit = "REML",
#'       theta_mode = "voxelwise"
#'     )
#'   ) |>
#'   compute()
#'
#' names(assays(g))
#' assay(g, "coef:time")
#'
#' rm_df_known <- transform(rm_df, var = 0.04)
#' g_known <- gds(rm_df_known, contrast_data_cols = "time") |>
#'   reduce(
#'     method = "lmm:ri_slope1_knownvar",
#'     formula = ~ time,
#'     options = list(slope = "time", covariance = "diag", fit = "REML")
#'   ) |>
#'   compute()
#' ```
#'
#' @return This is a documentation page (no return value); these reducers are
#'   invoked via [reduce()].
#' @name reducer-lmm
#' @keywords models mixed-models reduce
NULL
