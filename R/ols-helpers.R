#' Stack coefficient assays into an array
#'
#' Returns a samples × terms × contrasts array by stacking per-term coefficient
#' assays (e.g., coef:(Intercept), coef:age, ...) along a new terms dimension.
#'
#' @param gds A realised GDS
#' @param prefix Assay name prefix to stack (default: "coef:")
#'
#' @return A 3D numeric array `samples × terms × contrasts`
#' @export
coef_array <- function(gds, prefix = "coef:") {
  stopifnot(inherits(gds, "gds"))
  all_names <- names(assays(gds))
  term_assays <- all_names[startsWith(all_names, prefix)]
  if (!length(term_assays)) stop("No coefficient assays with prefix '", prefix, "'", call. = FALSE)
  # terms inferred from suffix after prefix
  terms <- sub(paste0("^", prefix), "", term_assays)
  base <- assay(gds, term_assays[[1]])
  dims <- dim(base)
  n_samples <- dims[1]; n_contrasts <- dims[3]
  out <- array(NA_real_, dim = c(n_samples, length(term_assays), n_contrasts),
               dimnames = list(NULL, terms, dimnames(base)[[3]]))
  for (ti in seq_along(term_assays)) {
    a <- assay(gds, term_assays[[ti]])
    out[, ti, ] <- a[, 1, ]
  }
  out
}

#' Retrieve packed covariance triangles for OLS coefficients
#'
#' @param gds A realised GDS
#' @param method Reducer method key (default: "ols:voxelwise")
#' @param contrast Contrast name to retrieve
#'
#' @return A list with fields `type`, `terms`, `pack`, and `cov_tri` (L × samples)
#' @export
coef_cov_tri <- function(gds, method = "ols:voxelwise", contrast) {
  stopifnot(inherits(gds, "gds"))
  att <- gds$metadata$attachments %||% list()
  key <- paste0("reduce/", method, "/contrast=", contrast)
  out <- att[[key]]
  if (is.null(out)) {
    # try index form
    idx <- match(contrast, gds$contrasts)
    if (!is.na(idx)) {
      out <- att[[paste0("reduce/", method, "/contrast=", idx)]]
    }
  }
  if (is.null(out)) stop("No packed covariance found for reduce/", method, " (contrast=", contrast, ")", call. = FALSE)
  out
}
