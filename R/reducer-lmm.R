#' Restricted repeated-measures Gaussian LMM reducers
#'
#' `fmrigds` includes two specialized repeated-measures Gaussian mixed-model
#' reducers that are optimized for the case where every sample shares the same
#' observation layout and the same fixed/random design.
#'
#' - `method = "lmm:ri"` fits a random-intercept model.
#' - `method = "lmm:ri_slope1"` fits a random intercept plus one within-subject
#'   slope supplied through `options$slope`.
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
#' - `theta_mode = "pooled"` for one shared variance structure across samples
#' - `theta_mode = "voxelwise"` for sample-specific variance parameters
#' - no crossed random effects and no general `lmer()` random-effects grammar
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
#'   ),
#'   var = 0.04
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
#' ```
#'
#' @return This is a documentation page (no return value); these reducers are
#'   invoked via [reduce()].
#' @name reducer-lmm
#' @keywords models mixed-models reduce
NULL
