#' Align subjects into a consensus space
#'
#' @param x A plan, source, or realised GDS
#' @param family A [`MapFamily`] object describing per-subject transforms
#'
#' @return Updated plan
#' @export
align <- function(x, family) {
  plan <- as_plan(x)
  fam <- if (is.character(family)) {
    get_map_family(plan, family)
  } else {
    if (!inherits(family, "gds_map_family")) {
      stop("`family` must be a MapFamily or the name of a registered family", call. = FALSE)
    }
    family
  }

  existing <- plan$meta$map_families %||% list()
  if (is.null(existing[[fam$name]])) {
    plan <- register_map(plan, fam, overwrite = FALSE)
  }

  add_op(plan, op_align_to_group(family = fam, family_name = fam$name))
}

#' Eagerly align subjects and compute immediately
#'
#' @param x Plan, source, or realized GDS
#' @param ... Arguments passed to align(), then compute()
#'
#' @return A realized GDS object
#' @export
align_eager <- function(x, ...) {
  compute(align(x, ...))
}
