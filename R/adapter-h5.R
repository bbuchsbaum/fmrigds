register_h5_adapter <- function() {
  register_adapter(
    name = "h5",
    detect = .h5_detect,
    open = .h5_open,
    probe = .h5_probe,
    read = .h5_read,
    close = .h5_close
  )
}

.h5_detect <- function(source) {
  if (!is.character(source) || length(source) != 1L) return(FALSE)
  if (!file.exists(source)) return(FALSE)
  ext <- tolower(tools::file_ext(source))
  if (!ext %in% c("h5", "hdf5")) return(FALSE)
  if (!requireNamespace("hdf5r", quietly = TRUE)) return(FALSE)
  h5 <- hdf5r::H5File$new(source, mode = "r")
  on.exit(h5$close())
  if (h5$exists("/gds")) return(1.0)
  FALSE
}

.h5_open <- function(source, mode = "r", ...) {
  list(file = hdf5r::H5File$new(source, mode = mode), path = source)
}

.h5_probe <- function(handle, ...) {
  file_path <- handle$path
  h5 <- handle$file
  subjects <- .h5_read_dataset(h5, "/gds/axes/subjects")
  contrasts <- .h5_read_dataset(h5, "/gds/axes/contrasts")
  # Schema version validation
  schema_full <- tryCatch(as.character(.h5_read_dataset(h5, "/gds/version")), error = function(e) NA_character_)
  schema_num <- "0.1"
  if (is.character(schema_full) && length(schema_full) == 1L && !is.na(schema_full)) {
    if (startsWith(schema_full, "gds-h5/")) {
      schema_num <- sub("^gds-h5/", "", schema_full)
      major <- strsplit(schema_num, ".", fixed = TRUE)[[1]][1]
      if (!identical(major, "0")) {
        stop(sprintf("Unsupported HDF5 schema version '%s' (major=%s). Supported major is 0.x. How to fix: write with a compatible fmrigds version or upgrade reader.", schema_full, major), call. = FALSE)
      }
    }
  }
  h5 <- handle$file
  assays_group <- h5$open("/gds/assays")
  on.exit(assays_group$close(), add = TRUE)
  assay_names <- assays_group$ls()
  if (is.data.frame(assay_names)) {
    assay_names <- assay_names$name
  }
  first_path <- paste0("/gds/assays/", assay_names[[1]])
  dims_raw <- .h5_dataset_dims(h5, first_path)
  space_group <- h5$open("/gds/space")
  on.exit(space_group$close(), add = TRUE)
  type <- as.character(.h5_read_attr(space_group, "type"))
  space <- switch(type,
    voxel = space_voxel(
      dim = .h5_read_dataset(h5, "/gds/space/voxel/dim"),
      affine = matrix(
        .h5_read_dataset(h5, "/gds/space/voxel/affine"),
        nrow = 4,
        byrow = TRUE
      ),
      mask_idx = if (space_group$exists("voxel/mask_idx")) .h5_read_dataset(h5, "/gds/space/voxel/mask_idx") else NULL,
      storage = "packed"
    ),
    parcels = space_parcels(.h5_read_dataset(h5, "/gds/space/parcels/labels")),
    sample_labels = {
      labels <- if (space_group$exists("sample_labels/labels")) {
        .h5_read_dataset(h5, "/gds/space/sample_labels/labels")
      } else {
        subjects
      }
      space_sample_labels(labels)
    },
    stop("Unsupported space type in HDF5: ", type, call. = FALSE)
  )
  maps <- .h5_read_alignments(h5)
  col_data <- .h5_read_subjects_table(h5, subjects)
  out <- list(
    assays = assay_names,
    dims = gds_dims(sample = dims_raw[1], subject = dims_raw[2], contrast = dims_raw[3]),
    subjects = as.character(subjects),
    contrasts = as.character(contrasts),
    space = space,
    maps = maps,
    metadata = list(schema_version = schema_num, source_file = handle$path),
    columns = list(),
    col_data = col_data
  )
  probe_contract(out)
}

.h5_read <- function(handle,
                     assays,
                     block = NULL,
                     ...) {
  h5 <- handle$file
  assays_group <- h5$open("/gds/assays")
  on.exit(assays_group$close())
  first_path <- paste0("/gds/assays/", assays[[1]])
  dims <- .h5_dataset_dims(h5, first_path)
  samples_idx <- if (!is.null(block) && !is.null(block$sample)) block$sample else seq_len(dims[1])
  subject_idx <- if (!is.null(block) && !is.null(block$subject)) block$subject else TRUE
  contrast_idx <- if (!is.null(block) && !is.null(block$contrast)) block$contrast else TRUE
  index <- list(samples_idx, subject_idx, contrast_idx)
  result <- lapply(assays, function(name) {
    path <- paste0("/gds/assays/", name)
    data <- .h5_read_dataset(h5, path, index = index)
    array(data, dim = c(length(samples_idx), dims[2], dims[3]))
  })
  names(result) <- assays
  result
}

.h5_close <- function(handle) {
  handle$file$close()
  invisible(NULL)
}
