#' Tidy long-form export for GDS
#'
#' Public wrapper around the internal tidy frame builder used by CSV/Parquet
#' exporters. Returns a data.frame; if you prefer a tibble, wrap with
#' `tibble::as_tibble()`.
#'
#' @param gds A realised GDS object
#' @param assays Character vector of assay names to include (default: all
#'   assays in the GDS)
#' @param drop_na If TRUE, drop rows with all-NA assays
#' @param include_col_data If TRUE, join columns from `col_data`
#'
#' @return A data.frame in long form with columns `sample`, `subject`,
#'   `contrast`, plus one column per requested assay.
#' @export
gds_to_tibble <- function(gds,
                          assays = NULL,
                          drop_na = FALSE,
                          include_col_data = FALSE) {
  if (is.null(assays)) assays <- names(fmrigds::assays(gds))
  .gds_tidy_frame(gds, stats = assays, drop_na = drop_na, include_col_data = include_col_data)
}

