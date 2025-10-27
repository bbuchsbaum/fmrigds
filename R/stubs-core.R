# Core noun constructors ---------------------------------------------------

# Lazy verbs ---------------------------------------------------------------

#' @export
align <- function(...) {
  stop("align() plan verb not yet implemented.")
}

#' @export
mask <- function(...) {
  stop("mask() plan verb not yet implemented.")
}

#' @export
reduce <- function(...) {
  stop("reduce() plan verb not yet implemented.")
}

#' @export
subset_eager <- function(...) {
  compute(subset(...))
}

#' @export
derive_eager <- function(...) {
  compute(derive(...))
}
