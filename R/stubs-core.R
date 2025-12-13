# Core noun constructors ---------------------------------------------------

# Lazy verbs ---------------------------------------------------------------

## Intentionally left blank: real reduce() is implemented in R/verb-reduce.R

#' Eagerly subset a GDS and compute immediately
#'
#' @param ... Arguments passed to subset(), then compute()
#'
#' @return A realized GDS object
#' @export
subset_eager <- function(...) compute(subset(...))

#' Eagerly derive statistics and compute immediately
#'
#' @param ... Arguments passed to derive(), then compute()
#'
#' @return A realized GDS object
#' @export
derive_eager <- function(...) compute(derive(...))
