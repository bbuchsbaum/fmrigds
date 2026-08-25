#' Register a map family on a plan or realised GDS
#'
#' @param x A [`gds_plan`] or realised [`gds`]
#' @param family A [`MapFamily`] object
#' @param overwrite Whether to overwrite an existing family with the same name
#'
#' @return Updated object (`gds_plan` or `gds`)
#' @export
register_map <- function(x, family, overwrite = FALSE) {
  if (!inherits(family, "gds_map_family")) {
    stop("`family` must be a MapFamily-created object", call. = FALSE)
  }

  if (inherits(x, "gds")) {
    .validate_map_family_subjects(family, x$subjects)
    maps <- x$metadata$map_families %||% list()
    if (!overwrite && family$name %in% names(maps)) {
      stop("Map family '", family$name, "' already registered; set overwrite = TRUE to replace.", call. = FALSE)
    }
    maps[[family$name]] <- family
    x$metadata$map_families <- maps
    x$metadata <- .ensure_nonfile_source(x$metadata, "memory", "map-registration-input")
    x$metadata <- add_provenance_node(
      x$metadata,
      "register_map",
      params = list(family = family$name, type = family$type),
      inputs = .provenance_current_inputs(x$metadata)
    )
    return(x)
  }

  plan <- as_plan(x)
  subjects <- plan$meta$subjects %||% plan$source$probe$subjects
  if (!is.null(subjects)) {
    .validate_map_family_subjects(family, subjects)
  }

  maps <- plan$meta$map_families %||% list()
  if (!overwrite && family$name %in% names(maps)) {
    stop("Map family '", family$name, "' already registered on plan; set overwrite = TRUE to replace.", call. = FALSE)
  }
  maps[[family$name]] <- family
  plan$meta$map_families <- maps
  plan
}

#' List registered map families
#'
#' @param x A plan or realised GDS
#'
#' @return Character vector of family names
#' @export
list_map_families <- function(x) {
  registry <- .map_registry(x)
  names(registry)
}

#' Retrieve a registered map family by name
#'
#' @param x A plan or realised GDS
#' @param name Map family name
#'
#' @return A [`MapFamily`] object
#' @export
get_map_family <- function(x, name) {
  registry <- .map_registry(x)
  fam <- registry[[name]]
  if (is.null(fam)) {
    stop("No map family named '", name, "' is registered on this object.", call. = FALSE)
  }
  fam
}

.map_registry <- function(x) {
  if (inherits(x, "gds")) {
    return(x$metadata$map_families %||% list())
  }
  if (inherits(x, "gds_plan")) {
    maps <- list()
    probe_maps <- x$source$probe$maps %||% list()
    if (length(probe_maps)) maps[names(probe_maps)] <- probe_maps
    meta_maps <- x$meta$map_families %||% list()
    if (length(meta_maps)) maps[names(meta_maps)] <- meta_maps
    return(maps)
  }
  stop("Object must be a gds or gds_plan", call. = FALSE)
}

.validate_map_family_subjects <- function(family, subjects) {
  fam_subj <- names(family$by_subject)
  if (is.null(fam_subj) || any(!nzchar(fam_subj))) {
    stop("Map family must have subject-specific operators with names", call. = FALSE)
  }
  if (!setequal(fam_subj, subjects)) {
    stop("Subject mismatch between map family and data", call. = FALSE)
  }
  invisible(NULL)
}
