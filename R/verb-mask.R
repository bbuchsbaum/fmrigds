#' Define a mask policy
#'
#' @param scope Scope of mask computation ("group" or "subject")
#' @param rule Mask rule ("intersection", "union", "threshold", "custom")
#' @param threshold Fraction for threshold rule
#' @param custom Custom function returning logical mask
#'
#' @return A mask policy object
#' @name MaskPolicy
#' @export
MaskPolicy <- function(scope = c("group", "subject"),
                       rule = c("intersection", "union", "threshold", "custom"),
                       threshold = 0.95,
                       custom = NULL) {
  scope <- match.arg(scope)
  rule <- match.arg(rule)
  structure(
    list(scope = scope, rule = rule, threshold = threshold, custom = custom),
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
