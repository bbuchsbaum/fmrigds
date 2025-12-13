write_gds_h5 <- function(gds,
                         file,
                         options = list()) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package required to write HDF5", call. = FALSE)
  }
  overwrite <- options$overwrite %||% TRUE
  if (!overwrite && file.exists(file)) {
    stop("File exists and overwrite = FALSE: ", file, call. = FALSE)
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  h5 <- hdf5r::H5File$new(file, mode = "w")
  on.exit(h5$close())

  g <- h5$create_group("gds")
  on.exit(g$close(), add = TRUE)
  # Write a small version stamp for validators/readers
  g$create_dataset("version", "gds-h5/0.1")

  axes <- g$create_group("axes")
  on.exit(axes$close(), add = TRUE)
  axes$create_dataset("subjects", gds$subjects)
  axes$create_dataset("contrasts", gds$contrasts)

  if (!is.null(gds$col_data) && nrow(gds$col_data) && ncol(gds$col_data)) {
    cd <- gds$col_data
    if (is.null(rownames(cd))) {
      rownames(cd) <- gds$subjects
    }
    cd <- cd[gds$subjects, , drop = FALSE]
    subjects_tbl <- axes$create_group("subjects_table")
    on.exit(subjects_tbl$close(), add = TRUE)
    schema <- lapply(colnames(cd), function(nm) {
      column <- cd[[nm]]
      prepared <- .prepare_col_data_column(column)
      ds <- subjects_tbl$create_dataset(nm, prepared$data)
      ds$close()
      list(name = nm, type = prepared$type)
    })
    schema_json <- jsonlite::toJSON(schema, auto_unbox = TRUE, null = "null")
    subjects_tbl$create_dataset("schema", schema_json)
  }

  space_group <- g$create_group("space")
  on.exit(space_group$close(), add = TRUE)
  space_group$create_attr("type", gds$space$type)
  if (inherits(gds$space, "space_voxel")) {
    voxel_group <- space_group$create_group("voxel")
    on.exit(voxel_group$close(), add = TRUE)
    voxel_group$create_dataset("dim", gds$space$dim)
    voxel_group$create_dataset("affine", gds$space$affine)
    if (!is.null(gds$space$mask_idx)) {
      voxel_group$create_dataset("mask_idx", gds$space$mask_idx)
    }
  }
  if (inherits(gds$space, "space_parcels")) {
    parcel_group <- space_group$create_group("parcels")
    on.exit(parcel_group$close(), add = TRUE)
    parcel_group$create_dataset("labels", gds$space$labels)
  }
  if (inherits(gds$space, "space_sample_labels")) {
    sl_group <- space_group$create_group("sample_labels")
    on.exit(sl_group$close(), add = TRUE)
    sl_labels <- gds$space$labels %||% seq_len(dim(assays(gds)[[1]])[1])
    sl_group$create_dataset("labels", sl_labels)
  }

  assays_group <- g$create_group("assays")
  on.exit(assays_group$close(), add = TRUE)
  dims <- dim(assays(gds)[[1]])
  chunk <- c(min(1024, dims[1]), 1, 1)
  for (name in names(assays(gds))) {
    ds <- assays_group$create_dataset(name, assays(gds)[[name]], chunk_dims = chunk)
    ds$close()
  }

  metadata_group <- g$create_group("metadata")
  on.exit(metadata_group$close(), add = TRUE)
  meta_json <- jsonlite::toJSON(.serialize_metadata(gds$metadata), auto_unbox = TRUE, null = "null")
  metadata_group$create_dataset("json", meta_json)

  prov_group <- g$create_group("provenance")
  on.exit(prov_group$close(), add = TRUE)
  log_entries <- gds$metadata$provenance$log %||% character()
  prov_group$create_dataset("log", log_entries)

  families <- gds$metadata$map_families %||% list()
  if (length(families)) {
    align_group <- g$create_group("alignments")
    on.exit(align_group$close(), add = TRUE)
    for (nm in names(families)) {
      fam_group <- align_group$create_group(nm)
      serialized <- .serialize_map_family_lines(families[[nm]])
      fam_group$create_dataset("serialized", serialized)
      fam_group$close()
    }
  }
}

.serialize_metadata <- function(meta) {
  if (inherits(meta, "numeric_version")) return(as.character(meta))
  if (is.list(meta)) return(lapply(meta, .serialize_metadata))
  meta
}

.serialize_map_family_lines <- function(family) {
  out <- character()
  con <- textConnection("out", "w", local = TRUE)
  serialize(family, con, ascii = TRUE)
  close(con)
  out
}

.prepare_col_data_column <- function(x) {
  if (is.factor(x)) {
    return(list(data = as.character(x), type = "character"))
  }
  if (is.logical(x)) {
    return(list(data = as.integer(x), type = "logical"))
  }
  if (is.integer(x)) {
    return(list(data = as.integer(x), type = "integer"))
  }
  if (is.numeric(x)) {
    return(list(data = as.numeric(x), type = "numeric"))
  }
  if (is.character(x)) {
    return(list(data = x, type = "character"))
  }
  # Fallback: coerce to character
  list(data = as.character(x), type = "character")
}
