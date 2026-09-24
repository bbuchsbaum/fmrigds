#' Define a mask policy
#'
#' @param scope Scope of mask computation. `"group"` checks coverage across
#'   subjects separately for each sample and contrast. `"subject"` checks
#'   coverage across contrasts for each subject and requires every subject to
#'   meet the rule.
#' @param rule Mask rule ("intersection", "union", "threshold", "custom")
#' @param threshold Minimum fraction of finite observations for the threshold
#'   rule, between 0 and 1 (inclusive). With `zero_is_missing = TRUE`, only
#'   finite nonzero observations count. The denominator is the number of
#'   subjects present at this point in the plan for `scope = "group"`.
#' @param custom Custom function returning one non-missing logical value per
#'   sample. It receives the assays after any zero-to-missing conversion.
#' @param zero_is_missing Logical; default `FALSE` preserves zeros as valid
#'   data, including binary 0/1 observations. Set to `TRUE` only when exact
#'   zeros denote background or unavailable observations. Zeros in `beta`
#'   (or the first assay when `beta` is absent) become `NA` in all assays at
#'   the same sample/subject/contrast positions before coverage is computed.
#'   Nonzero values, however small, are preserved.
#'
#' @details
#' Apply [mask()] before a group reducer or [examine_group()]. Masking removes
#' samples that fail the rule in every contrast. For group masks, a sample
#' retained for some contrasts has all observations set to `NA` in failing
#' contrasts, so those contrasts do not enter downstream tests or FDR families.
#' Coverage is based on the effect assay, not variance validity; reducers still
#' enforce their own variance, minimum-observation, and design requirements.
#' Reducers requiring complete data may exclude partially observed samples.
#' A threshold of `1 / 3` requires at least 7 nonzero observations among 20
#' subjects when `zero_is_missing = TRUE`; the other observations are missing
#' in the fit, not zero-valued contributions. Use subject coverage masks when
#' available to distinguish genuine zero effects from missing measurements.
#'
#' @examples
#' background_policy <- MaskPolicy(
#'   rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE
#' )
#' # Select subjects first, then mask, then fit:
#' # fit <- subject_gds |> mask(background_policy) |> one_sample() |> compute()
#'
#' @return A mask policy object
#' @name MaskPolicy
#' @export
MaskPolicy <- function(scope = c("group", "subject"),
                       rule = c("intersection", "union", "threshold", "custom"),
                       threshold = 0.95,
                       custom = NULL,
                       zero_is_missing = FALSE) {
  scope <- match.arg(scope)
  rule <- match.arg(rule)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold > 1) {
    stop("`threshold` must be a finite number between 0 and 1", call. = FALSE)
  }
  if (!is.logical(zero_is_missing) || length(zero_is_missing) != 1L ||
      is.na(zero_is_missing)) {
    stop("`zero_is_missing` must be TRUE or FALSE", call. = FALSE)
  }
  structure(
    list(scope = scope, rule = rule, threshold = threshold, custom = custom,
         zero_is_missing = zero_is_missing),
    class = "gds_mask_policy"
  )
}

#' Apply a mask policy lazily
#'
#' @param x Plan, source, or realised GDS
#' @param policy Mask policy created by [MaskPolicy()]
#'
#' @return Updated plan
#' @export
mask <- function(x, policy = MaskPolicy()) {
  if (!inherits(policy, "gds_mask_policy")) {
    stop("`policy` must be a gds_mask_policy", call. = FALSE)
  }
  plan <- as_plan(x)
  add_op(plan, op_mask_policy(policy))
}

#' Eagerly apply mask policy and compute immediately
#'
#' @param x Plan, source, or realized GDS
#' @param ... Arguments passed to mask(), then compute()
#'
#' @return A realized GDS object
#' @export
mask_eager <- function(x, ...) {
  compute(mask(x, ...))
}
