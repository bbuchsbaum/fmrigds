# Internal HDF5 helpers shared by adapters

.h5_safe_exists <- function(h5, path) {
  tryCatch({ obj <- h5$open(path); obj$close(); TRUE }, error = function(e) FALSE)
}

.h5_read_dataset <- function(h5, path, index = NULL) {
  ds <- h5$open(path)
  on.exit(ds$close())
  if (is.null(index)) {
    ds$read(drop = FALSE)
  } else {
    ds$read(args = index, drop = FALSE)
  }
}

.h5_dataset_dims <- function(h5, path) {
  ds <- h5$open(path)
  on.exit(ds$close())
  as.integer(ds$dims)
}

.h5_read_attr <- function(obj, name) {
  attr <- obj$attr_open(name)
  on.exit(attr$close())
  attr$read()
}

.h5_first_dataset_in_group <- function(h5, group_path) {
  # Prefer common dataset names
  for (nm in c("embedding", "data", "matrix", "values")) {
    cand <- paste0(group_path, "/", nm)
    if (.h5_safe_exists(h5, cand)) return(cand)
  }
  grp <- h5$open(group_path)
  on.exit(grp$close(), add = TRUE)
  lsdf <- tryCatch(grp$ls(), error = function(e) NULL)
  if (is.null(lsdf)) return(NA_character_)
  entries <- if (is.data.frame(lsdf)) lsdf$name else lsdf
  if (!length(entries)) return(NA_character_)
  for (nm in entries) {
    cand <- paste0(group_path, "/", nm)
    # Only return datasets (objects that can be opened and read)
    if (.h5_safe_exists(h5, cand)) return(cand)
  }
  NA_character_
}

.deserialize_map_family_lines <- function(lines) {
  con <- textConnection(lines, "r", local = TRUE)
  on.exit(close(con))
  unserialize(con)
}

.h5_read_alignments <- function(h5) {
  if (!.h5_safe_exists(h5, "/gds/alignments")) return(list())
  grp <- h5$open("/gds/alignments")
  on.exit(grp$close())
  entries <- grp$ls()
  names <- if (is.data.frame(entries)) entries$name else entries
  grp$close()
  on.exit(NULL)
  out <- list()
  for (nm in names) {
    ds <- h5$open(paste0("/gds/alignments/", nm, "/serialized"))
    lines <- ds$read()
    ds$close()
    out[[nm]] <- .deserialize_map_family_lines(lines)
  }
  out
}

.h5_read_axis_table <- function(h5, path, ids, what = "metadata") {
  if (!.h5_safe_exists(h5, path)) return(NULL)
  tbl <- h5$open(path)
  on.exit(tbl$close())
  if (!tbl$exists("schema")) return(NULL)
  schema_ds <- tbl$open("schema")
  on.exit(schema_ds$close(), add = TRUE)
  schema <- tryCatch(jsonlite::fromJSON(schema_ds$read(), simplifyDataFrame = TRUE), error = function(e) NULL)
  schema_ds$close()
  if (is.null(schema)) return(NULL)
  if (is.data.frame(schema)) {
    entries <- split(schema, seq_len(nrow(schema)))
  } else {
    entries <- schema
  }
  ids <- as.character(ids)
  out <- data.frame(row.names = ids)
  for (entry in entries) {
    nm <- entry$name %||% entry[["name"]]
    type <- entry$type %||% entry[["type"]]
    if (!length(nm) || !tbl$exists(nm)) next
    ds <- tbl$open(nm)
    val <- ds$read()
    ds$close()
    if (!length(val)) {
      out[[nm]] <- NA
      next
    }
    val <- switch(type,
      logical = as.logical(val),
      integer = as.integer(val),
      numeric = as.numeric(val),
      character = as.character(val),
      as.character(val)
    )
    if (length(val) != length(ids)) {
      warning(sprintf("Ignoring %s column '%s' with unexpected length", what, nm), call. = FALSE)
      next
    }
    out[[nm]] <- val
  }
  if (!ncol(out)) return(NULL)
  out
}

.h5_read_subjects_table <- function(h5, subjects) {
  .h5_read_axis_table(h5, "/gds/axes/subjects_table", subjects, what = "col_data")
}

.h5_read_samples_table <- function(h5, samples) {
  .h5_read_axis_table(h5, "/gds/axes/samples_table", samples, what = "row_data")
}

.h5_read_contrasts_table <- function(h5, contrasts) {
  .h5_read_axis_table(h5, "/gds/axes/contrasts_table", contrasts, what = "contrast_data")
}
