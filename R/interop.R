#' Coerce common R objects into a GDS
#'
#' Provides a lightweight front door for external packages holding data in
#' memory (lists, arrays, data.frames) to construct a validated GDS without
#' going through on-disk adapters.
#'
#' @param x Object to coerce (list/array/data.frame)
#' @param ... Additional arguments passed to methods
#'
#' @return A `gds` object
#' @export
as_gds <- function(x, ...) UseMethod("as_gds")

#' @rdname as_gds
#' @param space Optional space descriptor for the sample axis
#' @param subjects Optional subject identifiers
#' @param contrasts Optional contrast identifiers
#' @param col_data Optional subject-level covariates (rownames must match subjects)
#' @param row_data Optional sample-level metadata
#' @param metadata Optional metadata list merged into [`gds_metadata()`]
#' @export
as_gds.list <- function(x,
                        space = NULL,
                        subjects = NULL,
                        contrasts = NULL,
                        col_data = NULL,
                        row_data = NULL,
                        metadata = list(),
                        ...) {
  assays <- x
  # Upcast any 2-D arrays to 3-D [sample x subject x 1]
  assays <- lapply(assays, function(a) {
    if (is.array(a) && length(dim(a)) == 2L) {
      dn <- dimnames(a)
      a <- array(a, dim = c(nrow(a), ncol(a), 1L))
      dimnames(a) <- list(dn[[1L]], dn[[2L]], NULL)
    }
    a
  })
  # Validate and infer dimensions
  assays <- .validate_gds_assays(assays)
  dims <- dim(assays[[1L]])
  dn <- dimnames(assays[[1L]])
  subjects <- subjects %||% (dn[[2L]] %||% as.character(seq_len(dims[2L])))
  contrasts <- contrasts %||% (dn[[3L]] %||% paste0("contrast", seq_len(dims[3L])))
  sample_labels <- dn[[1L]] %||% as.character(seq_len(dims[1L]))
  space <- space %||% space_sample_labels(sample_labels)

  new_gds(
    assays = assays,
    space = space,
    subjects = as.character(subjects),
    contrasts = as.character(contrasts),
    col_data = col_data,
    row_data = row_data,
    metadata = metadata
  )
}

#' @rdname as_gds
#' @param assay_name Name to use for the single-array input (default: "beta")
#' @export
as_gds.array <- function(x,
                         space = NULL,
                         subjects = NULL,
                         contrasts = NULL,
                         assay_name = "beta",
                         col_data = NULL,
                         row_data = NULL,
                         metadata = list(),
                         ...) {
  if (length(dim(x)) == 2L) {
    x <- array(x, dim = c(nrow(x), ncol(x), 1L))
  } else if (length(dim(x)) != 3L) {
    stop("Array input must be 2-D or 3-D; canonical is [sample × subject × contrast]", call. = FALSE)
  }
  lst <- setNames(list(x), assay_name)
  as_gds.list(lst,
              space = space,
              subjects = subjects,
              contrasts = contrasts,
              col_data = col_data,
              row_data = row_data,
              metadata = metadata)
}

#' @rdname as_gds
#' @param mapping Optional list mapping column names for axes and assays.
#' Default: list(sample="sample", subject="subject", contrast="contrast",
#' assays=list(beta="beta", var="var"))
#' @param space Optional space for sample axis; defaults to sample labels from data
#' @export
as_gds.data.frame <- function(x,
                              mapping = NULL,
                              space = NULL,
                              ...) {
  mapping <- utils::modifyList(
    list(sample = "sample", subject = "subject", contrast = "contrast",
         assays = list(beta = "beta", var = "var", se = NULL),
         roi = NULL),
    mapping %||% list()
  )
  effect_cols <- mapping$assays
  effect_cols <- effect_cols[!vapply(effect_cols, is.null, logical(1L))]

  # Reuse tabular adapter's probe/read logic on in-memory data
  handle <- list(data = x, path = "<memory>")
  probe <- .tabular_probe(
    handle,
    effect_cols = effect_cols,
    subject_col = mapping$subject,
    sample_col = mapping$sample,
    roi_col = mapping$roi,
    contrast_col = mapping$contrast,
    space = space
  )
  arrs <- .tabular_read(
    handle,
    assays = names(effect_cols),
    effect_cols = effect_cols,
    subject_col = mapping$subject,
    sample_col = mapping$sample,
    roi_col = mapping$roi,
    contrast_col = mapping$contrast
  )
  new_gds(arrs, probe$space, probe$subjects, probe$contrasts)
}

#' Build a design matrix from attached col_data
#'
#' Convenience helper to construct a model matrix using subject-level
#' covariates already attached to a plan via [with_col_data()] or stored in a
#' realised GDS.
#'
#' @param x A `gds` or `gds_plan` (or an object coercible via [as_plan()])
#' @param formula Model formula (or string convertible via [stats::as.formula])
#'
#' @return A numeric model matrix
#' @export
model_matrix <- function(x, formula) {
  f <- tryCatch(if (inherits(formula, "formula")) formula else stats::as.formula(formula), error = function(e) NULL)
  if (is.null(f)) stop("`formula` must be a valid model formula or string", call. = FALSE)

  if (inherits(x, "gds")) {
    subjects <- x$subjects
    cd <- x$col_data
  } else {
    plan <- as_plan(x)
    subjects <- plan$source$probe$subjects %||% plan$meta$subjects
    cd <- plan$meta$col_data %||% NULL
  }
  if (is.null(subjects) || !length(subjects)) {
    stop("Cannot determine subject identifiers from `x`", call. = FALSE)
  }
  if (is.null(cd) || !nrow(cd)) {
    stop("No col_data found; attach with with_col_data() first", call. = FALSE)
  }
  # Align and build model matrix
  aligned <- .align_col_data_for_subjects(cd, subjects, warn_extra = FALSE, context = "model_matrix")
  dat <- data.frame(subject = subjects, aligned, check.names = FALSE)
  stats::model.matrix(f, data = dat)
}
