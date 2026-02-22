#' Voxelwise OLS reducer (per-sample GLM across subjects)
#'
#' The reducer identified by `method = "ols:voxelwise"` performs unweighted
#' ordinary least squares across subjects independently for each sample and
#' contrast. It accepts a model formula via `reduce(..., formula = ~ 1 + cov)`
#' (design is built from `with_col_data()` or `gds(..., col_data=...)`).
#'
#' Outputs are expanded into per-term assays to preserve the 3-axis invariant:
#' - `coef:<term>` and `se_coef:<term>` (one assay per model term)
#' - `t_coef:<term>` and `p_coef:<term>` convenience statistics
#' - `sigma2` (residual variance) and `df_res` (residual degrees of freedom)
#'
#' Optionally, a packed triangular covariance for the coefficient vector can be
#' attached per contrast using `options = list(return_cov = "tri")`. Retrieve
#' it from a realised GDS via [coef_cov_tri()].
#'
#' @section Usage:
#' ```r
#' # Attach subject-level covariates; build design via formula
#' plan <- gds("group.csv", col_data = subj_cov) |
#'   reduce(method = "ols:voxelwise", formula = ~ 1 + age + sex,
#'          options = list(return_cov = "tri"))
#' g <- compute(plan)
#'
#' # Per-term assays
#' assay(g, "coef:(Intercept)"); assay(g, "coef:age"); assay(g, "coef:sex")
#'
#' # Packed covariance for first contrast
#' cov_info <- coef_cov_tri(g, contrast = contrasts(g)[1])
#' cov_info$terms   # term order
#' cov_info$cov_tri # L x samples packed upper triangle
#'
#' # Stack all coefficients into [samples x terms x contrasts]
#' coef_arr <- coef_array(g)
#' ```
#'
#' @name reducer-ols-voxelwise
#' @keywords models regression reduce
NULL

