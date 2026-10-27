#' Assert compatibility between two spaces
#'
#' @param g1 First GDS or space
#' @param g2 Second GDS or space
#' @param fields Fields to compare exactly. By default, compares type, voxel
#'   dimensions and affine, and available identity metadata. An explicit list
#'   overrides these defaults and can request a weaker metadata-only check.
#'
#' @return TRUE if compatible; otherwise errors
#' @export
#' @examples
#' sp1 <- space_voxel(c(2,2,2), diag(4), storage = "dense")
#' sp2 <- space_voxel(c(2,2,2), diag(4), storage = "dense")
#' assert_compatible_spaces(sp1, sp2)
assert_compatible_spaces <- function(g1, g2, fields = c("type", "dim", "affine", "template_id")) {
  s1 <- if (inherits(g1, "gds_space")) g1 else g1$space
  s2 <- if (inherits(g2, "gds_space")) g2 else g2$space
  if (missing(fields)) {
    fields <- .default_space_compat_fields(s1, s2)
  }
  for (f in fields) {
    v1 <- s1[[f]] %||% NA
    v2 <- s2[[f]] %||% NA
    if (!identical(v1, v2)) {
      stop(sprintf("Space field '%s' mismatch: %s vs %s", f,
                   paste(format(v1), collapse = ", "),
                   paste(format(v2), collapse = ", ")), call. = FALSE)
    }
  }
  TRUE
}

.default_space_compat_fields <- function(s1, s2) {
  fields <- "type"

  if (inherits(s1, "space_voxel") && inherits(s2, "space_voxel")) {
    fields <- c(fields, "dim", "affine")
  }

  optional <- c(
    "template_id",
    "support_id",
    "grid_id",
    "topology_id",
    "mesh_id",
    "hemi",
    "reference_space"
  )

  for (field in optional) {
    if (!is.null(s1[[field]]) || !is.null(s2[[field]])) {
      fields <- c(fields, field)
    }
  }

  unique(fields)
}

#' Harmonise contrast names in a GDS
#'
#' @param g A realised GDS
#' @param map Named character vector mapping old -> new contrast names
#'   (unique, nonempty source names and nonempty destination names). The full
#'   resulting contrast axis must be unique; renaming does not aggregate data.
#'
#' @return Updated GDS with renamed contrasts
#' @export
#' @examples
#' beta <- array(0, dim = c(2,1,2), dimnames = list(NULL, "s1", c("a","b")))
#' var <- beta
#' g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("x","y")), "s1", c("a","b"))
#' harmonise_contrasts(g, c(a = "A", b = "B"))
harmonise_contrasts <- function(g, map) {
  stopifnot(inherits(g, "gds"))
  stopifnot(is.character(map), !is.null(names(map)))
  old <- names(map)
  new <- unname(map)
  if (anyNA(old) || any(!nzchar(old)) || anyDuplicated(old) ||
      anyNA(new) || any(!nzchar(new))) {
    stop("Contrast map must have unique, nonempty source names and nonempty destination names", call. = FALSE)
  }
  proposed <- unname(vapply(g$contrasts, function(x) if (x %in% old) new[match(x, old)] else x, character(1)))
  if (anyNA(proposed) || any(!nzchar(proposed)) || anyDuplicated(proposed)) {
    stop("Resulting contrast IDs must be unique and nonempty", call. = FALSE)
  }
  # Align metadata using the old IDs before assigning the new ones.
  cd <- contrast_data(g)
  if (!is.null(cd)) {
    cd <- .normalise_contrast_data(cd, g$contrasts)
    rownames(cd) <- proposed
    g$metadata$contrast_info$data <- cd
  }
  g$contrasts <- proposed
  # update dimnames on assays if present
  g$assays <- lapply(g$assays, function(a) {
    dn <- dimnames(a)
    if (!is.null(dn) && length(dn) >= 3L) {
      dn[[3]] <- proposed
      dimnames(a) <- dn
    }
    a
  })
  g
}

#' Relabel subjects in a GDS
#'
#' @param g A realised GDS
#' @param mapping Named character vector mapping old -> new subject ids
#'
#' @return Updated GDS with relabeled subjects
#' @export
#' @examples
#' beta <- array(0, dim = c(2,2,1), dimnames = list(NULL, c("u1","u2"), "c1"))
#' var <- beta
#' g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("x","y")), c("u1","u2"), "c1")
#' relabel_subjects(g, c(u1 = "s1", u2 = "s2"))
relabel_subjects <- function(g, mapping) {
  stopifnot(inherits(g, "gds"))
  stopifnot(is.character(mapping), !is.null(names(mapping)))
  old <- names(mapping)
  new <- unname(mapping)
  if (anyDuplicated(new)) stop("New subject ids must be unique", call. = FALSE)
  # update subjects
  g$subjects <- unname(vapply(g$subjects, function(x) if (x %in% old) new[match(x, old)] else x, character(1)))
  # update col_data rownames
  if (!is.null(g$col_data) && nrow(g$col_data)) {
    rn <- rownames(g$col_data)
    rn2 <- unname(vapply(rn, function(x) if (x %in% old) new[match(x, old)] else x, character(1)))
    rownames(g$col_data) <- rn2
    g$col_data <- g$col_data[g$subjects, , drop = FALSE]
  }
  # update assay dimnames if present
  g$assays <- lapply(g$assays, function(a) {
    dn <- dimnames(a)
    if (!is.null(dn) && length(dn) >= 2L) {
      dn[[2]] <- unname(vapply(dn[[2]], function(x) if (x %in% old) new[match(x, old)] else x, character(1)))
      dimnames(a) <- dn
    }
    a
  })
  g
}
