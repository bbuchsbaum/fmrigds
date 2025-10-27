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
  if (is.character(source) && length(source) == 1L) {
    ext <- tolower(tools::file_ext(source))
    if (ext %in% c("csv", "tsv", "parquet")) return(0.8)
  }
  FALSE
}

.tabular_open <- function(source, mode = "r", ...) {
  if (!identical(mode, "r")) stop("tabular adapter is read-only", call. = FALSE)
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
                           space = NULL,
                           ...) {
  df <- handle$data
  cols_needed <- c(unlist(effect_cols), subject_col, sample_col, contrast_col)
  missing <- setdiff(cols_needed, names(df))
  if (length(missing)) {
    stop("Missing required columns in tabular data: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  subjects <- unique(df[[subject_col]])
  samples <- unique(df[[sample_col]])
  contrasts <- unique(df[[contrast_col]])

  dims <- c(sample = length(samples), subject = length(subjects), contrast = length(contrasts))
  space <- space %||% space_parcels(labels = samples)

  list(
    assays = names(effect_cols),
    dims = dims,
    subjects = as.character(subjects),
    contrasts = as.character(contrasts),
    space = space,
    maps = list(),
    metadata = list(
      schema_version = "0.1.0",
      source_file = handle$path
    ),
    columns = list(
      effect_cols = effect_cols,
      subject_col = subject_col,
      sample_col = sample_col,
      contrast_col = contrast_col
    )
  )
}

.tabular_read <- function(handle,
                          assays,
                          block = NULL,
                          effect_cols = NULL,
                          subject_col = NULL,
                          sample_col = NULL,
                          contrast_col = NULL,
                          ...) {
  df <- handle$data
  cols <- effect_cols

  samples <- unique(df[[sample_col]])
  subjects <- unique(df[[subject_col]])
  contrasts <- unique(df[[contrast_col]])

  if (!is.null(block) && !is.null(block$sample)) {
    idx <- block$sample
    samples <- samples[idx]
    df <- df[df[[sample_col]] %in% samples, , drop = FALSE]
  }

  dims <- c(length(samples), length(subjects), length(contrasts))
  dimnames <- list(samples, subjects, contrasts)

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

register_builtin_adapters <- function() {
  if (!"tabular" %in% ls(.adapter_registry)) register_tabular_adapter()
  if (!"nifti" %in% ls(.adapter_registry)) register_nifti_adapter()
}
