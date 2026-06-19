# Internal helper to align/validate col_data against subjects
.align_col_data_for_subjects <- function(col_data, subjects, warn_extra = TRUE, context = "plan") {
  if (is.null(col_data) || is.null(subjects)) return(col_data)
  if (!is.data.frame(col_data)) {
    stop("col_data must be a data.frame", call. = FALSE)
  }
  subjects <- as.character(subjects)
  if (!length(subjects)) return(col_data)
  rn <- rownames(col_data)
  if (is.null(rn) || anyNA(rn)) {
    stop("col_data must have rownames matching subjects", call. = FALSE)
  }
  if (any(duplicated(rn))) {
    stop("col_data rownames must be unique", call. = FALSE)
  }
  missing <- setdiff(subjects, rn)
  if (length(missing)) {
    stop(sprintf("col_data is missing subjects: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  extra <- setdiff(rn, subjects)
  if (length(extra) && warn_extra) {
    warning(sprintf("Dropping extra rows in col_data: %s", paste(extra, collapse = ", ")), call. = FALSE)
  }
  col_data[subjects, , drop = FALSE]
}

.align_row_data_for_samples <- function(row_data, sample_labels, n_samples = NULL, warn_extra = TRUE, context = "plan") {
  if (is.null(row_data)) return(row_data)
  if (!is.data.frame(row_data)) {
    stop("row_data must be a data.frame", call. = FALSE)
  }

  if (is.null(sample_labels)) {
    n_samples <- n_samples %||% nrow(row_data)
    if (nrow(row_data) != n_samples) {
      stop("row_data must have one row per sample", call. = FALSE)
    }
    return(row_data)
  }

  sample_labels <- as.character(sample_labels)
  rn <- rownames(row_data)
  if (!is.null(rn) && !anyNA(rn) && !anyDuplicated(rn)) {
    missing <- setdiff(sample_labels, rn)
    if (length(missing)) {
      stop(sprintf("row_data is missing samples: %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
    extra <- setdiff(rn, sample_labels)
    if (length(extra) && warn_extra) {
      warning(sprintf("Dropping extra rows in row_data: %s", paste(extra, collapse = ", ")), call. = FALSE)
    }
    return(row_data[sample_labels, , drop = FALSE])
  }

  if (nrow(row_data) != length(sample_labels)) {
    stop("row_data must have one row per sample", call. = FALSE)
  }

  row_data
}

.align_contrast_data_for_contrasts <- function(contrast_data, contrasts, warn_extra = TRUE, context = "plan") {
  if (is.null(contrast_data) || is.null(contrasts)) return(contrast_data)
  if (!is.data.frame(contrast_data)) {
    stop("contrast_data must be a data.frame", call. = FALSE)
  }
  contrasts <- as.character(contrasts)
  if (!length(contrasts)) return(contrast_data)
  rn <- rownames(contrast_data)
  if (is.null(rn) || anyNA(rn)) {
    stop("contrast_data must have rownames matching contrasts", call. = FALSE)
  }
  if (any(duplicated(rn))) {
    stop("contrast_data rownames must be unique", call. = FALSE)
  }
  missing <- setdiff(contrasts, rn)
  if (length(missing)) {
    stop(sprintf("contrast_data is missing contrasts: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  extra <- setdiff(rn, contrasts)
  if (length(extra) && warn_extra) {
    warning(sprintf("Dropping extra rows in contrast_data: %s", paste(extra, collapse = ", ")), call. = FALSE)
  }
  contrast_data[contrasts, , drop = FALSE]
}

#' Adapter front door for building a plan
#'
#' @param source Source specification (paths, etc.)
#' @param format Optional adapter name or "auto"
#' @param prefer Optional preferred adapter when multiple match
#' @param ... Adapter-specific arguments forwarded to `open()`/`probe()`.
#'   Common metadata arguments include `col_data`, `row_data`, and
#'   `contrast_data`. For stacked repeated-measures tabular data,
#'   `contrast_data_cols` can name one or more columns to collapse into
#'   contrast-level metadata during ingestion.
#'
#' @return A [`gds_plan`]
#' @export
gds <- function(source,
                format = c("auto", ls(.adapter_registry)),
                prefer = NULL,
                ...) {
  format <- match.arg(format)
  adapter_name <- if (identical(format, "auto")) detect_adapter(source, prefer) else format
  adapter <- get_adapter(adapter_name)

  dots <- list(...)
  # Allow callers to pass subject-level covariates for meta-regression via col_data
  col_data <- NULL
  if ("col_data" %in% names(dots)) {
    col_data <- dots$col_data
    dots$col_data <- NULL
  }
  row_data <- NULL
  if ("row_data" %in% names(dots)) {
    row_data <- dots$row_data
    dots$row_data <- NULL
  }
  contrast_data <- NULL
  if ("contrast_data" %in% names(dots)) {
    contrast_data <- dots$contrast_data
    dots$contrast_data <- NULL
  }
  # Capture optional temporal policy for parcellated ingestion
  temporal_policy <- dots$temporal_policy %||% NULL
  contrast_matrix <- dots$contrast_matrix %||% NULL
  contrast_names  <- dots$contrast_names %||% NULL
  dots$temporal_policy <- NULL
  dots$contrast_matrix <- NULL
  dots$contrast_names  <- NULL

  if (identical(adapter_name, "nifti")) {
    source <- .nifti_attach_source_metadata(source, dots)
  }

  handle <- adapter$open(source)
  probe_result <- do.call(adapter$probe, c(list(handle), dots, list(temporal_policy = temporal_policy, contrast_matrix = contrast_matrix, contrast_names = contrast_names)))
  adapter$close(handle)

  adapter_columns <- probe_result$columns
  probe_result$columns <- NULL
  probe_col_data <- probe_result$col_data %||% NULL
  probe_result$col_data <- NULL
  probe_row_data <- probe_result$row_data %||% NULL
  probe_result$row_data <- NULL
  probe_contrast_data <- probe_result$contrast_data %||% NULL
  probe_result$contrast_data <- NULL

  src <- gds_source(adapter_name, source, probe_result)
  plan <- gds_plan(source = src)
  plan$meta$adapter_columns <- adapter_columns
  plan$meta$map_families <- probe_result$maps %||% list()
  plan$meta$subjects <- probe_result$subjects
  # Legacy-friendly alias for tests expecting plan$metadata$dims
  plan$metadata <- list(dims = probe_result$dims)
  subjects <- probe_result$subjects
  plan_col_data <- if (!is.null(col_data)) col_data else probe_col_data
  if (!is.null(plan_col_data) && !is.null(subjects)) {
    plan_col_data <- .align_col_data_for_subjects(plan_col_data, subjects, warn_extra = TRUE)
  }
  if (!is.null(plan_col_data)) plan$meta$col_data <- plan_col_data
  sample_labels <- probe_result$space$labels %||% NULL
  plan_row_data <- if (!is.null(row_data)) row_data else probe_row_data
  if (!is.null(plan_row_data)) {
    plan_row_data <- .align_row_data_for_samples(
      plan_row_data,
      sample_labels = sample_labels,
      n_samples = probe_result$dims[1L],
      warn_extra = TRUE
    )
    plan$meta$row_data <- plan_row_data
  }
  contrasts <- probe_result$contrasts
  plan_contrast_data <- if (!is.null(contrast_data)) contrast_data else probe_contrast_data
  if (!is.null(plan_contrast_data) && !is.null(contrasts)) {
    plan_contrast_data <- .align_contrast_data_for_contrasts(plan_contrast_data, contrasts, warn_extra = TRUE)
    plan$meta$contrast_data <- plan_contrast_data
  }
  if (!is.null(temporal_policy)) plan$meta$temporal_policy <- temporal_policy
  if (!is.null(contrast_matrix)) plan$meta$contrast_matrix <- contrast_matrix
if (!is.null(contrast_names)) plan$meta$contrast_names <- contrast_names
  plan
}

#' Attach subject-level covariates to a plan or GDS
#'
#' @param x Plan, source, or realised GDS
#' @param col_data Data frame keyed by subject identifiers (rownames)
#'
#' @return A [`gds_plan`] (for plans/sources) or a [`gds`] (for realised objects)
#' @export
with_col_data <- function(x, col_data) {
  if (inherits(x, "gds")) {
    aligned <- .align_col_data_for_subjects(col_data, x$subjects, warn_extra = TRUE, context = "gds")
    x$col_data <- aligned
    return(x)
  }
  plan <- as_plan(x)
  subjects <- plan$source$probe$subjects %||% plan$meta$subjects
  aligned <- .align_col_data_for_subjects(col_data, subjects, warn_extra = TRUE)
  plan$meta$col_data <- aligned
  plan
}

#' Attach sample-level metadata to a plan or GDS
#'
#' @param x Plan, source, or realised GDS
#' @param row_data Data frame keyed by sample identifiers (rownames) when available
#'
#' @return A [`gds_plan`] (for plans/sources) or a [`gds`] (for realised objects)
#' @export
with_row_data <- function(x, row_data) {
  if (inherits(x, "gds")) {
    sample_ids <- sample_labels(x)
    aligned <- .align_row_data_for_samples(row_data, sample_labels = sample_ids, n_samples = dim(assays(x)[[1L]])[1L], warn_extra = TRUE, context = "gds")
    x$row_data <- aligned
    return(x)
  }

  plan <- as_plan(x)
  probe <- plan$source$probe
  sample_ids <- probe$space$labels %||% NULL
  n_samples <- probe$dims[1L] %||% NULL
  aligned <- .align_row_data_for_samples(row_data, sample_labels = sample_ids, n_samples = n_samples, warn_extra = TRUE)
  plan$meta$row_data <- aligned
  plan
}

#' Attach contrast-level metadata to a plan or GDS
#'
#' @param x Plan, source, or realised GDS
#' @param contrast_data Data frame keyed by contrast identifiers (rownames)
#'
#' @return A [`gds_plan`] (for plans/sources) or a [`gds`] (for realised objects)
#' @export
with_contrast_data <- function(x, contrast_data) {
  if (inherits(x, "gds")) {
    aligned <- .align_contrast_data_for_contrasts(contrast_data, x$contrasts, warn_extra = TRUE, context = "gds")
    info <- x$metadata$contrast_info %||% list()
    info$data <- aligned
    x$metadata$contrast_info <- info
    return(x)
  }

  plan <- as_plan(x)
  contrasts <- plan$source$probe$contrasts %||% plan$meta$contrasts
  aligned <- .align_contrast_data_for_contrasts(contrast_data, contrasts, warn_extra = TRUE)
  plan$meta$contrast_data <- aligned
  plan
}
