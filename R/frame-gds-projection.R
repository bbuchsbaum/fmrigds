.one_projection_column <- function(value, name, allow_null = FALSE) {
  if (is.null(value) && allow_null) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop("`", name, "` must be one non-empty column name",
         if (allow_null) " or NULL" else "", ".", call. = FALSE)
  }
  value
}

.projection_axis_ids <- function(frame, column, supplied, label) {
  observations <- fmridataset::observations(frame)
  if (is.null(column)) {
    ids <- supplied %||% if (identical(label, "subject")) "meta" else "contrast1"
    if (!is.character(ids) || length(ids) != 1L || is.na(ids) || !nzchar(ids)) {
      stop("A singleton `", label, "s` value is required when `", label,
           " = NULL`.", call. = FALSE)
    }
    return(list(ids = ids, index = rep.int(1L, nrow(observations))))
  }
  if (!column %in% names(observations)) {
    stop("The ", label, " column '", column,
         "' is absent from frame observations.", call. = FALSE)
  }
  values <- observations[[column]]
  if (is.list(values) || !is.null(dim(values)) || anyNA(values)) {
    stop("The ", label, " column must contain complete scalar keys.",
         call. = FALSE)
  }
  values <- as.character(values)
  registry <- fmridataset::entities(frame)
  matching_entities <- fmridataset::entity_names(registry)[vapply(
    registry,
    function(value) identical(fmridataset::entity_key(value), column),
    logical(1)
  )]
  inferred <- if (length(matching_entities) == 1L) {
    fmridataset::entity_ids(registry[[matching_entities]])
  } else if (is.factor(observations[[column]])) {
    levels(observations[[column]])
  } else {
    unique(values)
  }
  inferred <- inferred[inferred %in% values]
  ids <- supplied %||% inferred
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids) || !setequal(ids, unique(values))) {
    stop("`", label, "s` must contain each observed ", label,
         " key exactly once.", call. = FALSE)
  }
  list(ids = ids, index = match(values, ids))
}

.projection_entity_data <- function(frame, column, ids, label) {
  if (is.null(column)) {
    out <- data.frame(ids, stringsAsFactors = FALSE)
    names(out) <- paste0(label, "_id")
    rownames(out) <- ids
    return(out)
  }
  registry <- fmridataset::entities(frame)
  matching <- fmridataset::entity_names(registry)[vapply(
    registry,
    function(value) identical(fmridataset::entity_key(value), column),
    logical(1)
  )]
  if (length(matching) == 1L) {
    value <- fmridataset::entity_data(registry[[matching]])
    value <- as.data.frame(value[match(ids, as.character(value[[column]])), , drop = FALSE])
  } else {
    value <- data.frame(ids, stringsAsFactors = FALSE)
    names(value) <- column
  }
  rownames(value) <- ids
  value
}

.projection_dtype_bytes <- function(dtype) {
  switch(
    dtype,
    logical = 1,
    int8 = 1,
    uint8 = 1,
    int16 = 2,
    uint16 = 2,
    int32 = 4,
    uint32 = 4,
    float32 = 4,
    int64 = 8,
    uint64 = 8,
    float64 = 8,
    8
  )
}

.frame_to_gds_space <- function(frame) {
  spatial <- fmridataset::space(frame)
  ids <- fmridataset::feature_ids(frame)
  if (inherits(spatial, "volume_space")) {
    dense <- identical(spatial$support, seq_len(prod(spatial$dim)))
    return(space_voxel(
      dim = spatial$dim,
      affine = spatial$affine,
      mask_idx = if (dense) NULL else spatial$support,
      storage = if (dense) "dense" else "packed",
      template_id = spatial$template
    ))
  }
  if (inherits(spatial, "parcel_space")) {
    return(space_parcels(
      labels = ids,
      lookup = as.data.frame(fmridataset::features(frame))
    ))
  }
  space_sample_labels(ids)
}

.projection_assay_names <- function(frame_names, assay_map) {
  if (is.null(assay_map)) assay_map <- c(variance = "var", std_error = "se")
  if (!is.character(assay_map) || is.null(names(assay_map)) ||
      anyNA(assay_map) || any(!nzchar(assay_map)) ||
      anyNA(names(assay_map)) || any(!nzchar(names(assay_map)))) {
    stop("`assay_map` must be a named character vector from frame to GDS names.",
         call. = FALSE)
  }
  out <- frame_names
  mapped <- match(frame_names, names(assay_map))
  out[!is.na(mapped)] <- unname(assay_map[mapped[!is.na(mapped)]])
  if (anyDuplicated(out)) {
    stop("`assay_map` creates duplicate GDS assay names.", call. = FALSE)
  }
  setNames(out, frame_names)
}

#' Project a rectangular frame to the legacy GDS cube
#'
#' This compatibility projection is deliberately strict. The selected
#' observation annotations must define one complete, unique Cartesian product
#' of subject and contrast keys. Ragged frames remain frames and are rejected
#' rather than padded with implicit slabs.
#'
#' @param x An `fmri_frame` or synchronized frame view.
#' @param subject,contrast Observation columns defining the rectangular axes.
#'   Either may be `NULL` only for an explicit singleton axis.
#' @param subjects,contrasts Optional explicit axis order or singleton label.
#'   A `NULL` subject column defaults to a singleton `"meta"` axis; a `NULL`
#'   contrast column defaults to `"contrast1"`.
#' @param assay_names Assays to realize into the cube.
#' @param assay_map Named character vector mapping frame assay names to GDS
#'   assay names. By default `variance` becomes `var` and `std_error` becomes
#'   `se`.
#' @param memory_budget Maximum total bytes realized across projected assays.
#' @param metadata Additional GDS metadata.
#' @param ... Reserved for future use.
#' @return A realized legacy `gds` compatibility projection.
#' @export
as_gds.fmri_frame <- function(
  x,
  subject = "subject_id",
  contrast = "contrast_id",
  subjects = NULL,
  contrasts = NULL,
  assay_names = names(fmridataset::assays(x)),
  assay_map = NULL,
  memory_budget = 2 * 1024^3,
  metadata = list(),
  ...
) {
  subject <- .one_projection_column(subject, "subject", allow_null = TRUE)
  contrast <- .one_projection_column(contrast, "contrast", allow_null = TRUE)
  if (!is.null(subject) && identical(subject, contrast)) {
    stop("`subject` and `contrast` must name distinct observation columns.",
         call. = FALSE)
  }
  if (nrow(x) < 1L) {
    stop("GDS projection requires at least one frame observation.",
         call. = FALSE)
  }
  if (!is.numeric(memory_budget) || length(memory_budget) != 1L ||
      is.na(memory_budget) || !is.finite(memory_budget) || memory_budget <= 0) {
    stop("`memory_budget` must be one positive finite byte count.", call. = FALSE)
  }
  available <- names(fmridataset::assays(x))
  if (!is.character(assay_names) || !length(assay_names) || anyNA(assay_names) ||
      any(!nzchar(assay_names)) || anyDuplicated(assay_names)) {
    stop("`assay_names` must contain unique non-empty assay names.", call. = FALSE)
  }
  missing_assays <- setdiff(assay_names, available)
  if (length(missing_assays)) {
    stop("Frame assays not found: ", paste(missing_assays, collapse = ", "), ".",
         call. = FALSE)
  }
  subject_axis <- .projection_axis_ids(x, subject, subjects, "subject")
  contrast_axis <- .projection_axis_ids(x, contrast, contrasts, "contrast")
  pair <- data.frame(
    subject = subject_axis$index,
    contrast = contrast_axis$index
  )
  if (anyDuplicated(pair)) {
    stop("Frame observations must define unique subject-contrast pairs.",
         call. = FALSE)
  }
  expected_n <- length(subject_axis$ids) * length(contrast_axis$ids)
  if (nrow(pair) != expected_n) {
    stop("Frame observations must form a complete Cartesian grid for GDS projection.",
         call. = FALSE)
  }
  grid <- expand.grid(
    subject = seq_along(subject_axis$ids),
    contrast = seq_along(contrast_axis$ids),
    KEEP.OUT.ATTRS = FALSE
  )
  observation_order <- match(
    paste(grid$subject, grid$contrast, sep = "\r"),
    paste(pair$subject, pair$contrast, sep = "\r")
  )
  if (anyNA(observation_order)) {
    stop("Frame observations must form a complete Cartesian grid for GDS projection.",
         call. = FALSE)
  }
  descriptors <- fmridataset::assays(x)[assay_names]
  required_bytes <- sum(vapply(descriptors, function(value) {
    nrow(x) * ncol(x) * .projection_dtype_bytes(value$dtype)
  }, numeric(1)))
  if (required_bytes > memory_budget) {
    stop(
      "GDS projection requires ", format(required_bytes, scientific = FALSE),
      " bytes, above `memory_budget`.", call. = FALSE
    )
  }
  mapped_names <- .projection_assay_names(assay_names, assay_map)
  arrays <- lapply(assay_names, function(name) {
    matrix_value <- fmridataset::collect_assay(
      x, assay = name, memory_budget = memory_budget, force = TRUE
    )[observation_order, , drop = FALSE]
    array(
      t(matrix_value),
      dim = c(ncol(x), length(subject_axis$ids), length(contrast_axis$ids)),
      dimnames = list(
        fmridataset::feature_ids(x), subject_axis$ids, contrast_axis$ids
      )
    )
  })
  names(arrays) <- unname(mapped_names)
  observation_data <- as.data.frame(fmridataset::observations(x))[
    observation_order, , drop = FALSE
  ]
  feature_data <- as.data.frame(fmridataset::features(x))
  rownames(feature_data) <- fmridataset::feature_ids(x)
  projection <- list(
    schema_version = 1L,
    subject_column = subject,
    contrast_column = contrast,
    observation_ids = fmridataset::observation_ids(x)[observation_order],
    observation_data = observation_data,
    feature_ids = fmridataset::feature_ids(x),
    feature_space = fmridataset::space(x),
    frame_assay_names = assay_names,
    gds_assay_names = unname(mapped_names),
    frame_metadata = x$metadata %||% x$base$metadata %||% list(),
    frame_provenance = x$provenance %||% x$base$provenance %||% NULL
  )
  metadata <- utils::modifyList(metadata, list(frame_projection = projection))
  g <- new_gds(
    assays = arrays,
    space = .frame_to_gds_space(x),
    subjects = subject_axis$ids,
    contrasts = contrast_axis$ids,
    col_data = .projection_entity_data(x, subject, subject_axis$ids, "subject"),
    row_data = feature_data,
    metadata = metadata
  )
  contrast_values <- .projection_entity_data(
    x, contrast, contrast_axis$ids, "contrast"
  )
  with_contrast_data(g, contrast_values)
}

.gds_to_frame_space <- function(x, feature_ids, feature_data) {
  projection <- metadata(x)$frame_projection %||% NULL
  if (!is.null(projection$feature_space)) {
    spatial <- projection$feature_space
    if (!identical(fmridataset::feature_ids(spatial), feature_ids)) {
      stop("Stored frame feature IDs do not align with the GDS sample axis.",
           call. = FALSE)
    }
    return(spatial)
  }
  spatial <- space(x)
  if (inherits(spatial, "space_voxel")) {
    support <- spatial$mask_idx %||% seq_len(prod(spatial$dim))
    return(fmridataset::volume_space(
      dim = spatial$dim,
      affine = spatial$affine,
      support = support,
      template = spatial$template_id
    ))
  }
  fmridataset::index_space(
    length(feature_ids),
    ids = feature_ids,
    namespace = paste0("fmrigds-", digest::digest(feature_ids, algo = "xxhash64")),
    data = feature_data
  )
}

#' Convert a legacy rectangular GDS cube to a canonical frame
#'
#' @param x A realized legacy `gds` object.
#' @param subject,contrast Names for the generated observation key columns.
#'   `NULL` recovers the names stored by [as_gds()] when available, otherwise
#'   defaults to `"subject_id"` and `"contrast_id"`.
#' @param assay_map Named character vector mapping legacy GDS assay names to
#'   canonical frame names.
#' @param ... Reserved for future use.
#' @return An `fmri_frame` with one row per subject-contrast cell.
#' @exportS3Method fmridataset::as_fmri_frame gds
as_fmri_frame.gds <- function(
  x,
  subject = NULL,
  contrast = NULL,
  assay_map = c(var = "variance", se = "std_error"),
  ...
) {
  projection <- metadata(x)$frame_projection %||% NULL
  subject <- subject %||% projection$subject_column %||% "subject_id"
  contrast <- contrast %||% projection$contrast_column %||% "contrast_id"
  subject <- .one_projection_column(subject, "subject")
  contrast <- .one_projection_column(contrast, "contrast")
  if (identical(subject, contrast)) {
    stop("`subject` and `contrast` must name distinct observation columns.",
         call. = FALSE)
  }
  grid <- expand.grid(
    subject_value = subjects(x),
    contrast_value = contrasts(x),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  observation_data <- if (!is.null(projection$observation_data)) {
    as.data.frame(projection$observation_data)
  } else {
    data.frame(stringsAsFactors = FALSE)
  }
  if (nrow(observation_data) != nrow(grid)) {
    observation_data <- data.frame(row.names = seq_len(nrow(grid)))
  }
  observation_data[[subject]] <- grid$subject_value
  observation_data[[contrast]] <- grid$contrast_value
  observation_ids <- projection$observation_ids %||%
    paste(grid$subject_value, grid$contrast_value, sep = "::")
  if (anyDuplicated(observation_ids)) {
    observation_ids <- sprintf("gds-observation-%06d", seq_len(nrow(grid)))
  }
  observation_data$.obs_id <- observation_ids
  observation_data <- observation_data[c(
    ".obs_id", subject, contrast,
    setdiff(names(observation_data), c(".obs_id", subject, contrast))
  )]

  gds_names <- names(assays(x))
  frame_names <- gds_names
  if (!is.null(projection$frame_assay_names) &&
      length(projection$frame_assay_names) == length(gds_names) &&
      identical(projection$gds_assay_names, gds_names)) {
    frame_names <- projection$frame_assay_names
  } else {
    mapped <- match(gds_names, names(assay_map))
    frame_names[!is.na(mapped)] <- unname(assay_map[mapped[!is.na(mapped)]])
  }
  if (anyDuplicated(frame_names)) {
    stop("GDS assay mapping creates duplicate frame assay names.", call. = FALSE)
  }
  frame_assays <- lapply(assays(x), function(value) {
    dims <- dim(value)
    matrix(aperm(value, c(2L, 3L, 1L)), nrow = dims[2L] * dims[3L],
           ncol = dims[1L])
  })
  names(frame_assays) <- frame_names

  feature_data <- as.data.frame(row_data(x))
  feature_ids <- projection$feature_ids %||% if (inherits(space(x), "space_voxel")) {
    support <- space(x)$mask_idx %||% seq_len(prod(space(x)$dim))
    paste0("voxel-", support)
  } else {
    sample_labels(x)
  }
  if (length(feature_ids) != nrow(feature_data)) {
    stop("GDS sample metadata does not align with the feature axis.", call. = FALSE)
  }
  feature_data$.feature_id <- feature_ids
  feature_data <- feature_data[c(
    ".feature_id", setdiff(names(feature_data), ".feature_id")
  )]
  spatial <- .gds_to_frame_space(x, feature_ids, feature_data)

  subject_data <- as.data.frame(col_data(x))
  subject_data[[subject]] <- subjects(x)
  subject_data <- subject_data[c(subject, setdiff(names(subject_data), subject))]
  contrast_values <- contrast_data(x)
  if (is.null(contrast_values)) {
    contrast_values <- data.frame(row.names = contrasts(x))
  }
  contrast_values <- as.data.frame(contrast_values)
  contrast_values[[contrast]] <- contrasts(x)
  contrast_values <- contrast_values[c(
    contrast, setdiff(names(contrast_values), contrast)
  )]
  entities <- list(
    subject = fmridataset::entity_frame(
      subject_data, key = subject, entity_type = "subject"
    ),
    contrast = fmridataset::entity_frame(
      contrast_values, key = contrast, entity_type = "contrast"
    )
  )
  provenance <- fmridataset::provenance_graph(fmridataset::provenance_record(
    "fmrigds::as_fmri_frame",
    inputs = list(
      subjects = subjects(x),
      contrasts = contrasts(x),
      assays = gds_names
    ),
    outputs = list(
      observation_ids = observation_ids,
      feature_ids = feature_ids,
      space_digest = fmridataset::space_digest(spatial)
    ),
    software = list(package = "fmrigds", version = .pkg_version())
  ))
  fmridataset::fmri_frame(
    assays = frame_assays,
    observations = observation_data,
    features = fmridataset::feature_axis(feature_data, space = spatial),
    entities = entities,
    relations = list(
      observation_subject = fmridataset::key_relation(
        subject, target = "subject"
      ),
      observation_contrast = fmridataset::key_relation(
        contrast, target = "contrast"
      )
    ),
    active_assay = frame_names[[1L]],
    metadata = projection$frame_metadata %||% list(
      legacy_gds = metadata(x)
    ),
    provenance = projection$frame_provenance %||% provenance
  )
}
