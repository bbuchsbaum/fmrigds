#' Create a voxel space from a NIfTI file
#'
#' @param path Path to a NIfTI file (.nii or .nii.gz)
#' @param mask Optional path to a mask NIfTI or a neuroim2::NeuroVol
#'
#' @return A `space_voxel` object
#' @export
#' @examples
#' # if (requireNamespace("neuroim2", quietly = TRUE)) {
#' #   sp <- space_from_nifti(system.file("extdata", "tiny.nii", package = "fmrigds"))
#' # }
space_from_nifti <- function(path, mask = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("neuroim2 package is required for space_from_nifti", call. = FALSE)
  }
  if (!file.exists(path)) stop("NIfTI file does not exist: ", path, call. = FALSE)
  hdr <- neuroim2::read_header(path)
  dim_all <- as.integer(hdr@dims)
  if (length(dim_all) < 3L) {
    stop("NIfTI dimensionality must be at least 3", call. = FALSE)
  }
  spatial_dim <- if (length(dim_all) >= 3L) dim_all[1:3] else dim_all

  img <- if (length(dim_all) == 3L) {
    neuroim2::read_vol(path)
  } else if (length(dim_all) == 4L) {
    neuroim2::read_vec(path)
  } else {
    ns <- asNamespace("neuroim2")
    if (exists("read_hyper_vec", envir = ns, mode = "function", inherits = FALSE)) {
      get("read_hyper_vec", envir = ns, inherits = FALSE)(path)
    } else {
      stop(
        "Reading 5D+ NIfTI requires neuroim2::read_hyper_vec(). ",
        "Please update neuroim2.", call. = FALSE
      )
    }
  }
  nspace <- neuroim2::space(img)
  affine <- neuroim2::trans(nspace)

  mask_bitmap <- NULL
  if (!is.null(mask)) {
    mask_vol <- if (is.character(mask)) neuroim2::read_vol(mask) else mask
    if (!inherits(mask_vol, "NeuroVol")) stop("mask must be a path or a NeuroVol", call. = FALSE)
    mask_bitmap <- as.array(mask_vol) > 0
  }
  space_voxel(dim = spatial_dim, affine = affine, mask_bitmap = mask_bitmap, storage = if (is.null(mask_bitmap)) "dense" else "packed")
}

#' Subset a space by sample indices
#'
#' @param space A `gds_space` object
#' @param idx Integer or logical index for the sample axis
#' @param pack Logical; when TRUE and `space` is a dense voxel space, return a
#'   packed voxel space using `idx` as the new mask indices.
#'
#' @return A new space object with samples restricted
#' @export
#' @examples
#' sp <- space_parcels(letters[1:5])
#' space_subset(sp, 2:4)
space_subset <- function(space, idx, pack = FALSE) {
  if (inherits(space, "space_parcels")) {
    if (is.logical(idx)) idx <- which(idx)
    space_parcels(labels = space$labels[idx], lookup = space$lookup, membership = NULL)
  } else if (inherits(space, "space_sample_labels")) {
    if (is.logical(idx)) idx <- which(idx)
    space_sample_labels(labels = space$labels[idx])
  } else if (inherits(space, "space_voxel")) {
    if (is.logical(idx)) idx <- which(idx)
    if (identical(space$storage, "packed") && !is.null(space$mask_idx)) {
      space$mask_idx <- as.integer(space$mask_idx[idx])
      # keep original mask_bitmap and dim; mask_bitmap may not reflect subset
      return(space)
    }
    # dense voxel space: optionally produce a packed subset
    if (isTRUE(pack)) {
      space$mask_idx <- as.integer(idx)
      space$mask_bitmap <- NULL
      space$storage <- "packed"
      return(space)
    }
    warning("space_subset: dense voxel space cannot be reduced; returning original space (set pack=TRUE to convert)", call. = FALSE)
    space
  } else {
    stop("Unsupported space type for subsetting", call. = FALSE)
  }
}

#' Compute a common mask between two spaces
#'
#' @param g1 First GDS or space
#' @param g2 Second GDS or space
#' @param rule Either "intersection" or "union"
#'
#' @return A list with integer indices `idx1` and `idx2` for subsetting each
#' @export
#' @examples
#' m1 <- array(c(TRUE, FALSE, TRUE, TRUE), c(2, 2, 1))
#' m2 <- array(c(TRUE, TRUE, FALSE, TRUE), c(2, 2, 1))
#' sp1 <- space_voxel(c(2,2,1), diag(4), mask_bitmap = m1,
#'   storage = "packed")
#' sp2 <- space_voxel(c(2,2,1), diag(4), mask_bitmap = m2,
#'   storage = "packed")
#' common_mask(sp1, sp2, rule = "intersection")
common_mask <- function(g1, g2, rule = c("intersection", "union")) {
  rule <- match.arg(rule)
  s1 <- if (inherits(g1, "gds_space")) g1 else g1$space
  s2 <- if (inherits(g2, "gds_space")) g2 else g2$space
  if (!(inherits(s1, "space_voxel") && inherits(s2, "space_voxel"))) {
    stop("common_mask requires voxel spaces", call. = FALSE)
  }
  if (!identical(s1$dim, s2$dim)) stop("Voxel dims differ; cannot combine", call. = FALSE)
  m1 <- s1$mask_idx %||% seq_len(prod(s1$dim))
  m2 <- s2$mask_idx %||% seq_len(prod(s2$dim))
  set <- if (identical(rule, "intersection")) intersect(m1, m2) else union(m1, m2)
  # map back to packed sample indices for each
  idx1 <- match(set, m1, nomatch = 0L)
  idx1 <- as.integer(idx1[idx1 > 0L])
  idx2 <- match(set, m2, nomatch = 0L)
  idx2 <- as.integer(idx2[idx2 > 0L])
  list(idx1 = idx1, idx2 = idx2)
}

#' Apply a common mask to two realised GDS objects
#'
#' @param g1 First GDS
#' @param g2 Second GDS
#' @param rule Mask rule ("intersection" or "union")
#'
#' @return A list with subsetted `g1`, `g2`, and indices `idx1`, `idx2`
#' @export
apply_common_mask <- function(g1, g2, rule = c("intersection", "union")) {
  rule <- match.arg(rule)
  stopifnot(inherits(g1, "gds"), inherits(g2, "gds"))
  cm <- common_mask(g1, g2, rule = rule)
  s1 <- cm$idx1; s2 <- cm$idx2
  # subset arrays
  subset_arrays <- function(arrs, idx) lapply(arrs, function(a) a[idx, , , drop = FALSE])
  g1$assays <- subset_arrays(g1$assays, s1)
  g2$assays <- subset_arrays(g2$assays, s2)
  # update spaces
  g1$space <- space_subset(g1$space, s1, pack = TRUE)
  g2$space <- space_subset(g2$space, s2, pack = TRUE)
  list(g1 = g1, g2 = g2, idx1 = s1, idx2 = s2)
}
