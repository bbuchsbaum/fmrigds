#' Find NIfTI map files with optional subject/contrast hooks
#'
#' Convenience helper to discover NIfTI map files on disk using a flexible
#' regular-expression pattern. This is intended for quick "scooping up" of
#' first-level maps spread across a directory tree (e.g. `.*/maps/AUC.nii`).
#'
#' The function does not impose any file naming convention. Instead, callers
#' can provide optional `subject_fun` and `contrast_fun` hooks to derive
#' subject and contrast labels from file paths.
#'
#' @param root Root directory to search.
#' @param pattern Regular expression passed to [base::list.files()] to match
#'   candidate NIfTI files (e.g. `"AUC\\\\.nii(\\\\.gz)?$"`).
#' @param recursive Logical; whether to search subdirectories (default `TRUE`).
#' @param subject_fun Optional function applied to the character vector of file
#'   paths, expected to return a character vector of subject identifiers.
#' @param contrast_fun Optional function applied to the character vector of
#'   file paths, expected to return a character vector of contrast labels.
#'
#' @return A data.frame with columns `file`, `subject`, and `contrast`.
#' @export
find_maps <- function(root,
                      pattern = "\\.nii(\\.gz)?$",
                      recursive = TRUE,
                      subject_fun = NULL,
                      contrast_fun = NULL) {
  if (!dir.exists(root)) stop("Root directory does not exist: ", root, call. = FALSE)
  files <- list.files(root, pattern = pattern, recursive = recursive, full.names = TRUE)
  if (!length(files)) {
    return(data.frame(
      file = character(),
      subject = character(),
      contrast = character(),
      stringsAsFactors = FALSE
    ))
  }
  subject <- rep(NA_character_, length(files))
  contrast <- rep(NA_character_, length(files))
  if (!is.null(subject_fun)) {
    subject <- as.character(subject_fun(files))
  }
  if (!is.null(contrast_fun)) {
    contrast <- as.character(contrast_fun(files))
  }
  data.frame(
    file = as.character(files),
    subject = subject,
    contrast = contrast,
    stringsAsFactors = FALSE
  )
}

#' Construct a GDS from NIfTI maps discovered on disk
#'
#' This helper builds on the NIfTI adapter to materialise a realised GDS from
#' a set of NIfTI map files, typically obtained via [find_maps()]. The heavy
#' lifting (reading NIfTI, spatial metadata, masking) is delegated to the
#' existing NIfTI adapter and the `neuroim2` package.
#'
#' Subject and contrast labels can be injected post-hoc using the `subject` and
#' `contrast` columns when available in `maps`. This avoids hard-coding any
#' filename convention while still allowing callers to impose semantics via
#' hooks.
#'
#' This helper is **single-contrast**: every file in `maps$file` is loaded as
#' one map per subject along a single contrast axis. A `contrast` column is used
#' only to *label* that one contrast (all non-missing values must agree); it
#' does not pivot a long subject x contrast table into a grid. For multi-contrast
#' workflows, call `gds_from_nifti_maps()` once per contrast and combine the
#' results, or assemble the source manually with [nifti_source()].
#'
#' Beta-only raw maps are accepted. When no variance or standard-error maps are
#' available, the realised GDS contains a synthetic `var` assay filled with 1 and
#' a warning is emitted. This placeholder variance is useful for satisfying the
#' GDS assay contract, but it is not meaningful subject-level uncertainty. For
#' group analyses over subject-level searchlight or stat-map metrics, prefer an
#' unweighted reducer such as `method = "ols:voxelwise"` instead of
#' fixed/random-effects meta-analysis.
#'
#' @param maps A data.frame (or tibble) with at least a `file` column of
#'   NIfTI paths. Optional `subject` and `contrast` columns can be used to
#'   relabel the realised GDS.
#' @param mask Optional mask passed through to [gds()] when constructing the
#'   NIfTI plan. Can be a NIfTI file path or a `neuroim2::NeuroVol`, matching
#'   the NIfTI adapter contract.
#' @param ... Additional arguments forwarded to [gds()] when building the
#'   NIfTI-backed plan.
#'
#' @return A realised GDS object.
#' @export
#' @examples
#' \dontrun{
#' maps <- data.frame(
#'   file = c("1001_era_partition/maps/naive_diag_mean.nii.gz",
#'            "1002_era_partition/maps/naive_diag_mean.nii.gz"),
#'   subject = c("1001", "1002"),
#'   contrast = "naive_diag_mean"
#' )
#'
#' group_info <- data.frame(
#'   group = c("control", "patient"),
#'   row.names = maps$subject
#' )
#'
#' g <- gds_from_nifti_maps(maps) |>
#'   with_col_data(group_info)
#'
#' fit <- reduce(g, method = "ols:voxelwise", formula = ~ group) |>
#'   compute()
#' }
gds_from_nifti_maps <- function(maps, mask = NULL, ...) { # nocov start
  if (is.null(maps) || !nrow(maps)) {
    stop("`maps` must contain at least one file", call. = FALSE)
  }
  if (is.null(maps$file)) stop("`maps` must have a `file` column", call. = FALSE)
  files <- as.character(maps$file)
  if (!all(file.exists(files))) {
    missing <- files[!file.exists(files)]
    stop("Some files in `maps$file` do not exist; first few missing: ",
         paste(utils::head(missing, 3L), collapse = ", "),
         call. = FALSE)
  }

  src <- nifti_source(beta = files)
  plan <- gds(src, format = "nifti", mask = mask, ...)
  g <- compute(plan)

  # Optional relabelling using hooks output, if provided. Use relabel_subjects()
  # so the stored col_data rownames and assay dimnames stay consistent with the
  # new labels (otherwise downstream reduce()/formula handling fails with
  # "col_data is missing subjects").
  if (!is.null(maps$subject)) {
    subj <- as.character(maps$subject)
    if (length(subj) == length(g$subjects) && !anyDuplicated(subj)) {
      g <- relabel_subjects(g, stats::setNames(subj, g$subjects))
    }
  }
  if (!is.null(maps$contrast)) {
    contr <- as.character(maps$contrast)
    non_missing <- contr[!is.na(contr)]
    # For a one-contrast NIfTI layout, per-file contrast metadata must agree.
    if (length(g$contrasts) == 1L && length(unique(non_missing)) == 1L) {
      g$contrasts <- unique(non_missing)
    }
  }
  g
} # nocov end
