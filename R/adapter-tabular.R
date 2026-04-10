#' Register tabular adapter (CSV/TSV/Parquet)
#'
#' Provides a reader for tabular summary data (ROIs/parcels/basis) by mapping
#' columns to assays and axes.
#'
#' @keywords internal
register_tabular_adapter <- function() {
  register_adapter(
    name = "tabular",
    detect = .tabular_detect,
    open = .tabular_open,
    probe = .tabular_probe,
    read = .tabular_read,
    close = .tabular_close
  )
}

.tabular_detect <- function(source) {
  if (is.data.frame(source)) return(0.9)
  if (is.character(source) && length(source) == 1L) {
    ext <- tolower(tools::file_ext(source))
    if (ext %in% c("csv", "tsv", "parquet")) return(0.8)
  }
  FALSE
}

.tabular_open <- function(source, mode = "r", ...) {
  if (!identical(mode, "r")) stop("tabular adapter is read-only", call. = FALSE)
  if (is.data.frame(source)) {
    return(list(path = "<memory>", data = source, ext = "memory"))
  }
  if (!file.exists(source)) stop("File does not exist: ", source, call. = FALSE)

  ext <- tolower(tools::file_ext(source))
  data <- switch(ext,
    csv = data.table::fread(source, data.table = FALSE, ...),
    tsv = data.table::fread(source, sep = "\t", data.table = FALSE, ...),
    parquet = {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("arrow package required to read parquet files", call. = FALSE)
      }
      as.data.frame(arrow::read_parquet(source))
    },
    stop("Unsupported tabular file extension: ", ext, call. = FALSE)
  )

  list(path = source, data = data, ext = ext)
}

.tabular_probe <- function(handle,
                           effect_cols = list(beta = "beta", var = "var"),
                           subject_col = "subject",
                           sample_col = "sample",
                           contrast_col = "contrast",
                           contrast_data_cols = NULL,
                           roi_col = NULL,
                           space = NULL,
                           ...) {
  # Alias: roi_col overrides sample_col if provided
  sample_col <- roi_col %||% sample_col
  df <- handle$data
  # Required axes columns must exist
  axes_needed <- c(subject_col, sample_col, contrast_col)
  missing_axes <- setdiff(axes_needed, names(df))
  if (length(missing_axes)) {
    stop("Missing required columns in tabular data: ", paste(missing_axes, collapse = ", "), call. = FALSE)
  }
  # Effect columns are optional; keep only those present
  present_effects <- effect_cols[vapply(effect_cols, function(nm) nm %in% names(df), logical(1))]
  if (!length(present_effects)) {
    stop("No effect columns present in tabular data.", call. = FALSE)
  }

  subjects <- unique(df[[subject_col]])
  samples <- unique(df[[sample_col]])
  contrasts <- unique(df[[contrast_col]])
  contrast_data <- .tabular_extract_contrast_data(
    df,
    contrast_col = contrast_col,
    contrast_data_cols = contrast_data_cols,
    contrasts = contrasts
  )

  dims <- gds_dims(sample = length(samples), subject = length(subjects), contrast = length(contrasts))
  space <- space %||% space_sample_labels(labels = samples)

  out <- list(
    assays = names(present_effects),
    dims = dims,
    subjects = as.character(subjects),
    contrasts = as.character(contrasts),
    space = space,
    maps = list(),
    metadata = list(
      schema_version = "0.1.0",
      source_file = handle$path
    ),
    contrast_data = contrast_data,
    columns = list(
      effect_cols = present_effects,
      subject_col = subject_col,
      sample_col = sample_col,
      contrast_col = contrast_col,
      contrast_data_cols = contrast_data_cols
    )
  )
  probe_contract(out)
}

.tabular_read <- function(handle,
                          assays,
                          block = NULL,
                          effect_cols = NULL,
                          subject_col = NULL,
                          sample_col = NULL,
                          contrast_col = NULL,
                          roi_col = NULL,
                          ...) {
  df <- handle$data
  # Fallback to probe-provided columns when args are NULL
  if (is.null(effect_cols) || is.null(subject_col) || is.null(sample_col) || is.null(contrast_col)) {
    # We don't have direct access to probe_result here; compute() passes columns explicitly.
    # If missing, try to infer sensible defaults
    effect_cols <- effect_cols %||% list(beta = "beta", var = "var")
    subject_col <- subject_col %||% "subject"
    sample_col <- sample_col %||% "sample"
    contrast_col <- contrast_col %||% "contrast"
  }
  # Alias: roi_col overrides sample column
  sample_col <- roi_col %||% sample_col
  cols <- effect_cols

  samples <- unique(df[[sample_col]])
  subjects <- unique(df[[subject_col]])
  contrasts <- unique(df[[contrast_col]])

  if (!is.null(block)) {
    if (!is.null(block$sample)) {
      bs <- block$sample
      if (is.logical(bs)) {
        samples <- samples[which(bs)]
      } else {
        samples <- samples[as.integer(bs)]
      }
      df <- df[df[[sample_col]] %in% samples, , drop = FALSE]
    }
    if (!is.null(block$subject)) {
      bsj <- block$subject
      if (is.logical(bsj)) {
        subjects <- subjects[which(bsj)]
      } else {
        subjects <- subjects[as.integer(bsj)]
      }
      df <- df[df[[subject_col]] %in% subjects, , drop = FALSE]
    }
    if (!is.null(block$contrast)) {
      bc <- block$contrast
      if (is.logical(bc)) {
        contrasts <- contrasts[which(bc)]
      } else {
        contrasts <- contrasts[as.integer(bc)]
      }
      df <- df[df[[contrast_col]] %in% contrasts, , drop = FALSE]
    }
  }

  dims <- c(length(samples), length(subjects), length(contrasts))
  dimnames <- list(samples, subjects, contrasts)

  # Warn on duplicate keys
  key <- paste(df[[sample_col]], df[[subject_col]], df[[contrast_col]], sep = "\t")
  if (anyDuplicated(key)) {
    warning("Duplicate (sample,subject,contrast) rows found; using first occurrence.", call. = FALSE)
    df <- df[!duplicated(key), , drop = FALSE]
  }

  result <- lapply(assays, function(name) {
    col <- cols[[name]]
    arr <- array(NA_real_, dim = dims, dimnames = dimnames)
    idx <- cbind(
      match(df[[sample_col]], samples),
      match(df[[subject_col]], subjects),
      match(df[[contrast_col]], contrasts)
    )
    arr[idx] <- df[[col]]
    arr
  })
  names(result) <- assays
  result
}

.tabular_close <- function(handle) {
  invisible(NULL)
}

.tabular_extract_contrast_data <- function(df,
                                           contrast_col,
                                           contrast_data_cols = NULL,
                                           contrasts = unique(df[[contrast_col]])) {
  cols <- as.character(contrast_data_cols %||% character())
  if (!length(cols)) return(NULL)

  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop("Missing contrast_data columns in tabular data: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  contrasts_chr <- as.character(contrasts)
  split_idx <- split(seq_len(nrow(df)), factor(as.character(df[[contrast_col]]), levels = contrasts_chr))
  out <- setNames(vector("list", length(cols)), cols)

  for (nm in cols) {
    values <- df[[nm]]
    if (is.factor(values)) values <- as.character(values)
    col_out <- lapply(seq_along(contrasts_chr), function(i) {
      idx <- split_idx[[i]]
      .tabular_constant_by_group(values[idx], column = nm, group = contrasts_chr[[i]])
    })
    if (is.logical(values)) {
      out[[nm]] <- as.logical(unlist(col_out, use.names = FALSE))
    } else if (is.integer(values)) {
      out[[nm]] <- as.integer(unlist(col_out, use.names = FALSE))
    } else if (is.numeric(values)) {
      out[[nm]] <- as.numeric(unlist(col_out, use.names = FALSE))
    } else {
      out[[nm]] <- as.character(unlist(col_out, use.names = FALSE))
    }
  }

  data.frame(out, row.names = contrasts_chr, check.names = FALSE, stringsAsFactors = FALSE)
}

.tabular_constant_by_group <- function(x, column, group) {
  if (!length(x)) return(NA)
  if (all(is.na(x))) return(NA)
  non_missing <- x[!is.na(x)]
  ref <- non_missing[[1L]]
  if (any(is.na(x)) || any(non_missing != ref)) {
    stop(
      "contrast_data column '", column, "' must be constant within contrast '", group, "'",
      call. = FALSE
    )
  }
  ref
}

register_builtin_adapters <- function() {
  if (!"memory" %in% ls(.adapter_registry)) register_memory_adapter()
  if (!"tabular" %in% ls(.adapter_registry)) register_tabular_adapter()
  if (!"nifti" %in% ls(.adapter_registry)) register_nifti_adapter()
  if (!"h5" %in% ls(.adapter_registry)) register_h5_adapter()
  if (!"fmristore" %in% ls(.adapter_registry)) register_fmristore_adapter()
  if (!"nftab" %in% ls(.adapter_registry) &&
      requireNamespace("neurotabs", quietly = TRUE)) {
    register_nftab_adapter()
  }
}
