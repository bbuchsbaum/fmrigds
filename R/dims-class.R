# gds_dims: dim triplet with friendly subsetting for tests

#' Construct a GDS dims vector
#' @keywords internal
gds_dims <- function(sample, subject, contrast) {
  v <- c(sample = as.integer(sample), subject = as.integer(subject), contrast = as.integer(contrast))
  class(v) <- c("gds_dims", class(v))
  v
}

#' Subset method strips names to avoid attribute sensitivity in tests
#' @keywords internal
`[.gds_dims` <- function(x, i, ...) {
  out <- NextMethod("[")
  unname(out)
}

