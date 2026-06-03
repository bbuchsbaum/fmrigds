#' Create a GDS from subject-level scalar NIfTI maps
#'
#' `gds_from_scalar_maps()` is a semantic wrapper for beta-only/stat-map
#' workflows: each NIfTI is a subject-level metric map rather than a first-level
#' beta image with meaningful uncertainty. The resulting GDS uses the NIfTI
#' adapter's unit-variance placeholder behavior for beta-only sources.
#'
#' @param files Character vector of NIfTI paths, a data.frame with a `file`
#'   column, or an [image_catalog()].
#' @param subject Optional subject identifiers. Required for character `files`
#'   unless `files` is named; ignored when `files` already has a `subject`
#'   column.
#' @param contrast Contrast label(s). A scalar value is recycled across files.
#'   Ignored when `files` already has a `contrast` column unless explicitly
#'   supplied.
#' @param col_data Optional subject-level covariates. If a `subject` column is
#'   present and row names are absent, it is used as the row key.
#' @param mask Optional mask passed to [gds_from_nifti_maps()].
#' @param ... Additional arguments forwarded to [gds_from_nifti_maps()].
#'
#' @return A realised [`gds`] object.
#' @export
#' @examples
#' \dontrun{
#' g <- gds_from_scalar_maps(
#'   files = c("sub-01/maps/auc.nii.gz", "sub-02/maps/auc.nii.gz"),
#'   subject = c("sub-01", "sub-02"),
#'   contrast = "auc",
#'   col_data = participants
#' )
#' fit <- group_ols(g, ~ group) |> compute()
#' }
gds_from_scalar_maps <- function(files,
                                 subject = NULL,
                                 contrast = "map",
                                 col_data = NULL,
                                 mask = NULL,
                                 ...) {
  maps <- .scalar_maps_frame(
    files,
    subject = subject,
    contrast = contrast,
    subject_supplied = !missing(subject),
    contrast_supplied = !missing(contrast)
  )
  g <- gds_from_nifti_maps(maps, mask = mask, ...)
  col_data <- col_data %||% .scalar_maps_col_data(files)
  if (!is.null(col_data)) {
    g <- with_col_data(g, .normalise_keyed_col_data(col_data, .gds_subjects(g)))
  }
  g
}

#' @rdname gds_from_scalar_maps
#' @export
as_scalar_map_gds <- gds_from_scalar_maps

#' Voxelwise OLS helpers for scalar-map group analyses
#'
#' These helpers wrap `reduce(method = "ols:voxelwise")` for subject-level map
#' workflows. They return a lazy GDS plan; call [compute()] to materialise the
#' fitted maps.
#'
#' @param x A realised [`gds`] or [`gds_plan`].
#' @param formula One-sided model formula for [group_ols()].
#' @param col_data Optional subject-level covariates to attach before fitting.
#' @param options Reducer options passed to [reduce()].
#' @param ... Additional arguments passed to [reduce()].
#'
#' @return A [`gds_plan`] with an OLS reduce node.
#' @export
#' @examples
#' \dontrun{
#' fit <- group_ols(g, ~ group) |> compute()
#' one <- one_sample(g) |> compute()
#' two <- two_sample(g, group = "diagnosis", baseline = "control") |> compute()
#' }
group_ols <- function(x,
                      formula = ~ 1,
                      col_data = NULL,
                      options = list(),
                      ...) {
  if (!is.null(col_data)) {
    x <- with_col_data(x, .normalise_keyed_col_data(col_data, .gds_subjects(as_plan(x))))
  }
  reduce(
    x,
    method = "ols:voxelwise",
    weights = "equal",
    formula = formula,
    options = options,
    ...
  )
}

#' @rdname group_ols
#' @param subset Optional expression evaluated in `col_data(x)` to select
#'   subjects before fitting the one-sample model.
#' @export
one_sample <- function(x,
                       subset = NULL,
                       col_data = NULL,
                       options = list(),
                       ...) {
  if (!is.null(col_data)) {
    x <- with_col_data(x, .normalise_keyed_col_data(col_data, .gds_subjects(as_plan(x))))
  }
  plan <- as_plan(x)
  if (!missing(subset) && !is.null(substitute(subset))) {
    plan <- .subset_plan_by_col_expr(plan, substitute(subset), parent.frame())
  }
  group_ols(plan, formula = ~ 1, options = options, ...)
}

#' @rdname group_ols
#' @param group Name of the grouping column in `col_data(x)`.
#' @param baseline Optional reference level for the grouping factor.
#' @param level Optional non-baseline level to compare against `baseline`. When
#'   omitted, all non-baseline levels are retained in the model.
#' @export
two_sample <- function(x,
                       group = "group",
                       baseline = NULL,
                       level = NULL,
                       col_data = NULL,
                       options = list(),
                       ...) {
  plan <- as_plan(x)
  if (!is.null(col_data)) {
    plan <- with_col_data(plan, .normalise_keyed_col_data(col_data, .gds_subjects(plan)))
  }
  subjects_in_model <- .plan_subjects_for_model_matrix(plan)
  cd <- .gds_col_data(plan)
  if (is.null(cd)) {
    stop("two_sample() requires col_data; attach it with with_col_data() or pass `col_data`.", call. = FALSE)
  }
  cd <- .align_col_data_for_subjects(cd, subjects_in_model, warn_extra = FALSE, context = "two_sample")
  if (!group %in% names(cd)) {
    stop("Grouping column not found in col_data: ", group, call. = FALSE)
  }

  grp <- cd[[group]]
  if (!is.null(level)) {
    keep_levels <- c(baseline, level)
    if (is.null(baseline)) {
      stop("`baseline` is required when `level` is supplied.", call. = FALSE)
    }
    keep <- grp %in% keep_levels
    if (!any(keep)) stop("No subjects match requested two-sample levels.", call. = FALSE)
    cd <- cd[keep, , drop = FALSE]
    grp <- cd[[group]]
    plan <- subset(plan, subject = rownames(cd))
  }

  grp <- factor(grp)
  if (!is.null(baseline)) {
    if (!baseline %in% levels(grp)) {
      stop("Baseline level not present in `", group, "`: ", baseline, call. = FALSE)
    }
    grp <- stats::relevel(grp, ref = baseline)
  }
  cd[[group]] <- grp
  plan <- with_col_data(plan, cd)
  group_ols(plan, formula = stats::reformulate(group), options = options, ...)
}

#' Write image-like GDS assays as NIfTI files
#'
#' Writes each selected assay/subject/contrast image in a voxel-space GDS to a
#' separate NIfTI file and returns a manifest describing written and skipped
#' outputs. This is intended for model outputs such as `coef:*`, `t_coef:*`,
#' and `p_coef:*` assays.
#'
#' @param g A realised [`gds`] with voxel space.
#' @param out_dir Output directory.
#' @param prefix Optional filename prefix.
#' @param assays Character vector of assays to write. Defaults to all assays.
#' @param subjects Optional subject labels to write. Defaults to all subjects.
#' @param contrasts Optional contrast labels to write. Defaults to all contrasts.
#' @param sanitize Logical; if `TRUE`, make filenames filesystem-friendly.
#' @param overwrite Logical; if `FALSE`, existing files are skipped.
#'
#' @return A data.frame manifest with `assay`, `subject`, `contrast`, `path`,
#'   `written`, and `skipped_reason` columns.
#' @export
#' @examples
#' \dontrun{
#' fit <- group_ols(g, ~ group) |> compute()
#' manifest <- write_nifti_assays(fit, "group_maps", prefix = "auc")
#' }
write_nifti_assays <- function(g,
                               out_dir,
                               prefix = NULL,
                               assays = NULL,
                               subjects = NULL,
                               contrasts = NULL,
                               sanitize = TRUE,
                               overwrite = FALSE) {
  stopifnot(inherits(g, "gds"))
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("The 'RNifti' package is required for NIfTI export.", call. = FALSE)
  }
  sp <- space(g)
  if (!inherits(sp, "space_voxel")) {
    stop("write_nifti_assays() requires a GDS with voxel space.", call. = FALSE)
  }

  available_assays <- names(.gds_assays(g))
  assays_to_write <- assays %||% available_assays
  missing_assays <- setdiff(assays_to_write, available_assays)
  present_assays <- intersect(assays_to_write, available_assays)

  subject_labels <- .gds_subjects(g)
  contrast_labels <- .gds_contrasts(g)
  subj_idx <- .match_axis(subjects %||% subject_labels, subject_labels, "subjects")
  con_idx <- .match_axis(contrasts %||% contrast_labels, contrast_labels, "contrasts")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  add_row <- function(assay_name, subject, contrast, path, written, skipped_reason) {
    rows[[length(rows) + 1L]] <<- data.frame(
      assay = assay_name,
      subject = subject,
      contrast = contrast,
      path = path,
      written = written,
      skipped_reason = skipped_reason,
      stringsAsFactors = FALSE
    )
  }

  for (nm in missing_assays) {
    add_row(nm, NA_character_, NA_character_, NA_character_, FALSE, "assay not found")
  }

  for (nm in present_assays) {
    arr <- assay(g, nm)
    if (!.is_image_assay(arr, sp)) {
      add_row(nm, NA_character_, NA_character_, NA_character_, FALSE, "not an image assay")
      next
    }
    for (j in subj_idx) {
      for (k in con_idx) {
        file_name <- .nifti_assay_filename(
          assay = nm,
          subject = subject_labels[j],
          contrast = contrast_labels[k],
          prefix = prefix,
          include_subject = length(subj_idx) > 1L,
          include_contrast = length(con_idx) > 1L,
          sanitize = sanitize
        )
        path <- file.path(out_dir, file_name)
        if (file.exists(path) && !isTRUE(overwrite)) {
          add_row(nm, subject_labels[j], contrast_labels[k], path, FALSE, "file exists")
          next
        }
        vol <- .array_to_nifti_volume(arr[, j, k], sp)
        .write_voxel_nifti(vol, sp$affine, path)
        add_row(nm, subject_labels[j], contrast_labels[k], path, TRUE, NA_character_)
      }
    }
  }

  if (!length(rows)) {
    return(data.frame(
      assay = character(),
      subject = character(),
      contrast = character(),
      path = character(),
      written = logical(),
      skipped_reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

.scalar_maps_frame <- function(files,
                               subject = NULL,
                               contrast = "map",
                               subject_supplied = !is.null(subject),
                               contrast_supplied = !is.null(contrast)) {
  if (inherits(files, "image_catalog")) {
    maps <- files$metadata
    if (!"file" %in% names(maps)) maps$file <- files$files
  } else if (is.data.frame(files)) {
    maps <- files
  } else {
    paths <- as.character(files)
    if (is.null(subject)) {
      subject <- names(files)
      if (is.null(subject) || any(!nzchar(subject))) {
        subject <- tools::file_path_sans_ext(basename(paths))
        subject <- sub("\\.nii$", "", subject)
      }
    }
    maps <- data.frame(file = paths, subject = as.character(subject), stringsAsFactors = FALSE)
  }

  if (!"file" %in% names(maps)) {
    stop("Scalar map inputs must include a `file` column.", call. = FALSE)
  }
  n <- nrow(maps)
  if (!"subject" %in% names(maps) || subject_supplied) {
    if (is.null(subject)) stop("`subject` is required when map metadata has no subject column.", call. = FALSE)
    if (length(subject) != n) stop("Length of `subject` must match number of files.", call. = FALSE)
    maps$subject <- as.character(subject)
  }
  if (!"contrast" %in% names(maps) || contrast_supplied) {
    if (length(contrast) == 1L) contrast <- rep(contrast, n)
    if (length(contrast) != n) stop("Length of `contrast` must be 1 or match number of files.", call. = FALSE)
    maps$contrast <- as.character(contrast)
  }
  maps[, c("file", "subject", "contrast"), drop = FALSE]
}

.scalar_maps_col_data <- function(files) {
  if (!inherits(files, "image_catalog")) return(NULL)
  meta <- files$metadata
  if (!"subject" %in% names(meta)) return(NULL)
  subject_meta <- meta[!duplicated(meta$subject), , drop = FALSE]
  rownames(subject_meta) <- subject_meta$subject
  drop_cols <- c("file", "basename", "relpath", "contrast")
  keep_cols <- setdiff(names(subject_meta), drop_cols)
  subject_meta <- subject_meta[, keep_cols, drop = FALSE]
  if (ncol(subject_meta) > 1L) subject_meta else NULL
}

.normalise_keyed_col_data <- function(col_data, subjects) {
  if (is.null(col_data)) return(NULL)
  if (!is.data.frame(col_data)) stop("col_data must be a data.frame", call. = FALSE)
  if ((is.null(rownames(col_data)) || identical(rownames(col_data), as.character(seq_len(nrow(col_data))))) &&
      "subject" %in% names(col_data)) {
    rownames(col_data) <- as.character(col_data$subject)
  }
  .align_col_data_for_subjects(col_data, subjects, warn_extra = TRUE, context = "scalar_maps")
}

.subset_plan_by_col_expr <- function(plan, expr, env) {
  subjects_in_model <- .plan_subjects_for_model_matrix(plan)
  cd <- .gds_col_data(plan)
  if (is.null(cd)) stop("`subset` requires col_data.", call. = FALSE)
  cd <- .align_col_data_for_subjects(cd, subjects_in_model, warn_extra = FALSE, context = "one_sample")
  keep <- eval(expr, envir = cd, enclos = env)
  if (!is.logical(keep) || length(keep) != nrow(cd)) {
    stop("`subset` must evaluate to one logical value per subject.", call. = FALSE)
  }
  subset(plan, subject = rownames(cd)[keep])
}

.gds_assays <- function(x) assays(x)
.gds_subjects <- function(x) subjects(x)
.gds_contrasts <- function(x) contrasts(x)
.gds_col_data <- function(x) col_data(x)

.is_image_assay <- function(arr, sp) {
  is.array(arr) &&
    length(dim(arr)) == 3L &&
    dim(arr)[1L] == length(sp$mask_idx %||% seq_len(prod(sp$dim)))
}

.array_to_nifti_volume <- function(vec, sp) {
  n_vox <- prod(sp$dim)
  if (identical(sp$storage, "packed") && !is.null(sp$mask_idx)) {
    full <- numeric(n_vox)
    full[sp$mask_idx] <- vec
  } else {
    full <- vec
  }
  array(full, dim = sp$dim)
}

.nifti_assay_filename <- function(assay,
                                  subject,
                                  contrast,
                                  prefix = NULL,
                                  include_subject = FALSE,
                                  include_contrast = FALSE,
                                  sanitize = TRUE) {
  parts <- c(prefix, assay)
  if (include_subject) parts <- c(parts, subject)
  if (include_contrast) parts <- c(parts, contrast)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (sanitize) parts <- vapply(parts, .sanitize_filename_part, character(1))
  paste0(paste(parts, collapse = "_"), ".nii.gz")
}

.sanitize_filename_part <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "map" else x
}

.match_axis <- function(values, labels, axis_name) {
  idx <- match(values, labels)
  if (anyNA(idx)) {
    stop("Unknown ", axis_name, ": ", paste(values[is.na(idx)], collapse = ", "), call. = FALSE)
  }
  idx
}
