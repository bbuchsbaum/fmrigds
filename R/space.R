#' Create a voxel space descriptor
#'
#' @param dim Integer vector of length 3 (x, y, z)
#' @param affine 4x4 affine matrix
#' @param mask_bitmap Optional logical array matching `dim`
#' @param mask_idx Optional integer vector of active voxel indices (mask order)
#' @param storage Either "dense" or "packed"
#' @param template_id Optional template identifier string
#'
#' @return An object of class `c("space_voxel", "gds_space")`
#' @export
space_voxel <- function(dim,
                        affine,
                        mask_bitmap = NULL,
                        mask_idx = NULL,
                        storage = c("dense", "packed"),
                        template_id = NULL) {
  storage <- match.arg(storage)
  .validate_voxel_inputs(dim, affine, mask_bitmap, mask_idx)

  if (!is.null(mask_bitmap) && is.null(mask_idx)) {
    mask_idx <- which(as.vector(mask_bitmap))
  }

  structure(
    list(
      type = "voxel",
      dim = as.integer(dim),
      affine = affine,
      mask_bitmap = mask_bitmap,
      mask_idx = if (!is.null(mask_idx)) as.integer(mask_idx) else mask_idx,
      storage = storage,
      template_id = template_id
    ),
    class = c("Space", "space_voxel", "gds_space")
  )
}

#' Alias for space_voxel
#'
#' @rdname space_voxel
#' @export
space_voxels <- space_voxel

#' Create a parcels/ROI space descriptor
#'
#' @param labels Character vector of parcel labels
#' @param lookup Optional lookup information (data frame or list)
#' @param membership Optional mapping (list or sparse matrix) from parcels to voxels
#'
#' @return Space object
#' @export
space_parcels <- function(labels,
                          lookup = NULL,
                          membership = NULL) {
  if (!is.character(labels) || !length(labels)) {
    stop("`labels` must be a non-empty character vector", call. = FALSE)
  }

  structure(
    list(
      type = "parcels",
      labels = labels,
      lookup = lookup,
      membership = membership
    ),
    class = c("Space", "space_parcels", "gds_space")
  )
}

#' Create a surface space descriptor
#'
#' @param vertices Matrix (n x 3) of vertex coordinates
#' @param faces Integer matrix (m x 3) of triangle indices (1-based)
#' @param hemi Hemisphere identifier ("L", "R", "LR")
#' @param template_id Optional template identifier
#'
#' @return Space object
#' @export
space_surface <- function(vertices,
                          faces,
                          hemi,
                          template_id = NULL) {
  if (!is.matrix(vertices) || ncol(vertices) != 3L) {
    stop("`vertices` must be a matrix with 3 columns", call. = FALSE)
  }
  if (!is.matrix(faces) || ncol(faces) != 3L) {
    stop("`faces` must be a matrix with 3 columns", call. = FALSE)
  }
  if (!hemi %in% c("L", "R", "LR")) {
    stop("`hemi` must be one of 'L', 'R', 'LR'", call. = FALSE)
  }

  structure(
    list(
      type = "surface",
      vertices = vertices,
      faces = faces,
      hemi = hemi,
      template_id = template_id
    ),
    class = c("Space", "space_surface", "gds_space")
  )
}

#' Create a latent basis space descriptor
#'
#' @param k Number of components
#' @param basis_name Optional name for the basis (e.g., "ICA")
#' @param projector Optional matrix/function for voxel projection
#' @param voxel_space Optional reference voxel space descriptor
#'
#' @return Space object
#' @export
space_basis <- function(k,
                        basis_name = NULL,
                        projector = NULL,
                        voxel_space = NULL) {
  if (!is.numeric(k) || length(k) != 1L || k <= 0) {
    stop("`k` must be a positive scalar", call. = FALSE)
  }

  if (!is.null(projector) && !(is.matrix(projector) || is.function(projector))) {
    stop("`projector` must be a matrix or function", call. = FALSE)
  }

  structure(
    list(
      type = "basis",
      k = as.integer(k),
      basis_name = basis_name,
      projector = projector,
      voxel_space = voxel_space
    ),
    class = c("Space", "space_basis", "gds_space")
  )
}

.validate_voxel_inputs <- function(dim, affine, mask_bitmap, mask_idx) {
  if (!is.numeric(dim) || length(dim) != 3L) {
    stop("`dim` must be numeric of length 3", call. = FALSE)
  }
  if (!is.matrix(affine) || !all(dim(affine) == c(4L, 4L))) {
    stop("`affine` must be a 4x4 matrix", call. = FALSE)
  }
  if (!is.null(mask_bitmap)) {
    if (!is.logical(mask_bitmap)) {
      stop("`mask_bitmap` must be logical", call. = FALSE)
    }
    if (!identical(dim(mask_bitmap), as.integer(dim))) {
      stop("`mask_bitmap` must have dimensions matching `dim`", call. = FALSE)
    }
  }
  if (!is.null(mask_idx)) {
    if (!is.numeric(mask_idx)) {
      stop("`mask_idx` must be numeric", call. = FALSE)
    }
    if (any(mask_idx < 1L) || any(mask_idx > prod(dim))) {
      stop("`mask_idx` contains indices out of range", call. = FALSE)
    }
  }
  invisible(NULL)
}
#' Create a simple label space for tabular samples
#'
#' @param labels Character vector of sample labels
#' @return Space object of type "sample_labels"
#' @export
space_sample_labels <- function(labels) {
  if (!is.character(labels) || !length(labels)) {
    stop("`labels` must be a non-empty character vector", call. = FALSE)
  }
  structure(
    list(
      type = "sample_labels",
      labels = labels
    ),
    class = c("Space", "space_sample_labels", "gds_space")
  )
}
