#' Attach a custom weight array to a GDS
#'
#' @param g A realised GDS
#' @param name Assay-like name for the weights (e.g., "n_eff" or "w_custom")
#' @param array Numeric 3D array matching `sample × subject × contrast`
#'
#' @return Updated GDS with weights stored in assays
#' @export
#' @examples
#' beta <- array(0, dim = c(2,2,1)); var <- array(1, dim = c(2,2,1))
#' g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("a","b")), c("s1","s2"), "c1")
#' w <- array(1, dim = c(2,2,1))
#' g2 <- attach_weight(g, "w_custom", w)
attach_weight <- function(g, name, array) {
  stopifnot(inherits(g, "gds"))
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  stopifnot(is.array(array), length(dim(array)) == 3L)
  dims <- dim(assays(g)[[1L]])
  if (!identical(dim(array), dims)) stop("Weight array dims must match GDS assays", call. = FALSE)
  g$assays[[name]] <- array
  g
}

#' Convenience helper to use a stored weight array for reduction
#'
#' @param g A realised GDS containing the weight assay
#' @param name Assay name to use as custom weights
#'
#' @return A list(list) suitable to splice into reduce():
#'   list(weights = "custom", options = list(custom_weights = <array>))
#' @export
#' @examples
#' # opts <- use_weight(g2, "w_custom")
use_weight <- function(g, name) {
  stopifnot(inherits(g, "gds"))
  w <- assays(g)[[name]]
  if (is.null(w)) stop("Weight assay not found: ", name, call. = FALSE)
  list(weights = "custom", options = list(custom_weights = w))
}
