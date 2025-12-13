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

#' Adapter front door for building a plan
#'
#' @param source Source specification (paths, etc.)
#' @param format Optional adapter name or "auto"
#' @param prefer Optional preferred adapter when multiple match
#' @param ... Adapter-specific arguments forwarded to `open()`/`probe()`
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
  # Capture optional temporal policy for parcellated ingestion
  temporal_policy <- dots$temporal_policy %||% NULL
  contrast_matrix <- dots$contrast_matrix %||% NULL
  contrast_names  <- dots$contrast_names %||% NULL
  dots$temporal_policy <- NULL
  dots$contrast_matrix <- NULL
  dots$contrast_names  <- NULL

  handle <- adapter$open(source)
  probe_result <- do.call(adapter$probe, c(list(handle), dots, list(temporal_policy = temporal_policy, contrast_matrix = contrast_matrix, contrast_names = contrast_names)))
  adapter$close(handle)

  adapter_columns <- probe_result$columns
  probe_result$columns <- NULL
  probe_col_data <- probe_result$col_data %||% NULL
  probe_result$col_data <- NULL

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
