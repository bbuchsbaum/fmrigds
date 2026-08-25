register_h5_adapter <- function() {
  register_adapter(
    name = "h5",
    detect = .h5_detect,
    open = .h5_open,
    probe = .h5_probe,
    read = .h5_read,
    close = .h5_close,
    capabilities = list(
      sample_blocks = TRUE,
      subject_blocks = TRUE,
      contrast_blocks = TRUE,
      persistent_handle = TRUE,
      preferred_axis = "sample",
      cheap_revisit = TRUE
    )
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
  if (!is.character(source) || length(source) != 1L) {
    stop(
      "The HDF5 adapter accepts a single native fmrigds .h5/.hdf5 file ",
      "(one file = one GDS), not a vector of per-subject files. To combine ",
      "per-subject maps use nifti_source()/as_gds(), or write a single GDS first.",
      call. = FALSE
    )
  }
  pre_read_stat <- .source_stat(source)
  list(
    file = hdf5r::H5File$new(source, mode = mode),
    path = source,
    pre_read_stat = pre_read_stat
  )
}

.h5_probe <- function(handle, source_identity = "sha256", ...) {
  file_path <- handle$path
  h5 <- handle$file
  source_identity <- .normalize_source_identity_policy(source_identity)
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
    voxel = {
      mask_idx <- if (space_group$exists("voxel/mask_idx")) {
        .h5_read_dataset(h5, "/gds/space/voxel/mask_idx")
      } else NULL
      space_voxel(
        dim = .h5_read_dataset(h5, "/gds/space/voxel/dim"),
        affine = matrix(
          .h5_read_dataset(h5, "/gds/space/voxel/affine"),
          nrow = 4,
          byrow = TRUE
        ),
        mask_idx = mask_idx,
        storage = if (is.null(mask_idx)) "dense" else "packed"
      )
    },
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
  sample_ids <- if (inherits(space, "space_sample_labels") || inherits(space, "space_parcels")) {
    as.character(space$labels)
  } else if (!is.null(space$mask_idx)) {
    as.character(space$mask_idx)
  } else {
    as.character(seq_len(dims_raw[1]))
  }
  row_data <- .h5_read_samples_table(h5, sample_ids)
  contrast_data <- .h5_read_contrasts_table(h5, as.character(contrasts))
  stored_metadata <- tryCatch({
    if (!.h5_safe_exists(h5, "/gds/metadata/json")) {
      list()
    } else {
      json <- as.character(.h5_read_dataset(h5, "/gds/metadata/json"))
      jsonlite::fromJSON(paste(json, collapse = ""), simplifyVector = FALSE)
    }
  }, error = function(e) {
    warning("Could not decode stored GDS metadata: ", conditionMessage(e), call. = FALSE)
    list()
  })
  # Original receipts describe the data embedded in the container; when the
  # container is re-opened, only its new container receipt is live-verified.
  stored_metadata <- .deactivate_source_verification(stored_metadata)
  identity_records <- list(list(
    path = file_path,
    role = "h5-container",
    ordinal = 1L,
    kind = "file",
    expected_stat = handle$pre_read_stat
  ))
  metadata <- .merge_gds_metadata(gds_metadata(), stored_metadata)
  metadata$source_file <- handle$path
  out <- list(
    assays = assay_names,
    dims = gds_dims(sample = dims_raw[1], subject = dims_raw[2], contrast = dims_raw[3]),
    subjects = as.character(subjects),
    contrasts = as.character(contrasts),
    space = space,
    maps = maps,
    metadata = metadata,
    source_identity_records = identity_records,
    source_identity_policy = source_identity,
    columns = list(),
    col_data = col_data,
    row_data = row_data,
    contrast_data = contrast_data
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
  samples_idx <- .normalize_block_index(block$sample %||% NULL, dims[1], "sample")
  subject_idx <- .normalize_block_index(block$subject %||% NULL, dims[2], "subject")
  contrast_idx <- .normalize_block_index(block$contrast %||% NULL, dims[3], "contrast")
  index <- list(samples_idx, subject_idx, contrast_idx)
  subjects <- as.character(.h5_read_dataset(h5, "/gds/axes/subjects"))[subject_idx]
  contrasts <- as.character(.h5_read_dataset(h5, "/gds/axes/contrasts"))[contrast_idx]
  sample_labels <- if (.h5_safe_exists(h5, "/gds/space/parcels/labels")) {
    as.character(.h5_read_dataset(h5, "/gds/space/parcels/labels"))[samples_idx]
  } else if (.h5_safe_exists(h5, "/gds/space/sample_labels/labels")) {
    as.character(.h5_read_dataset(h5, "/gds/space/sample_labels/labels"))[samples_idx]
  } else if (.h5_safe_exists(h5, "/gds/space/voxel/mask_idx")) {
    as.character(.h5_read_dataset(h5, "/gds/space/voxel/mask_idx"))[samples_idx]
  } else {
    as.character(samples_idx)
  }
  result <- lapply(assays, function(name) {
    path <- paste0("/gds/assays/", name)
    data <- .h5_read_dataset(h5, path, index = index)
    array(
      data,
      dim = c(length(samples_idx), length(subject_idx), length(contrast_idx)),
      dimnames = list(sample_labels, subjects, contrasts)
    )
  })
  names(result) <- assays
  result
}

.h5_close <- function(handle) {
  handle$file$close()
  invisible(NULL)
}
