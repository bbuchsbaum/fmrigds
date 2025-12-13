#' Assert compatibility between two spaces
#'
#' @param g1 First GDS or space
#' @param g2 Second GDS or space
#' @param fields Fields to compare (subset of c("type","dim","template_id"))
#'
#' @return TRUE if compatible; otherwise errors
#' @export
#' @examples
#' sp1 <- space_voxel(c(2,2,2), diag(4), storage = "dense")
#' sp2 <- space_voxel(c(2,2,2), diag(4), storage = "dense")
#' assert_compatible_spaces(sp1, sp2)
assert_compatible_spaces <- function(g1, g2, fields = c("type", "dim", "template_id")) {
  s1 <- if (inherits(g1, "gds_space")) g1 else g1$space
  s2 <- if (inherits(g2, "gds_space")) g2 else g2$space
  for (f in fields) {
    v1 <- s1[[f]] %||% NA
    v2 <- s2[[f]] %||% NA
    if (!identical(v1, v2)) stop(sprintf("Space field '%s' mismatch: %s vs %s", f, format(v1), format(v2)), call. = FALSE)
  }
  TRUE
}

#' Harmonise contrast names in a GDS
#'
#' @param g A realised GDS
#' @param map Named character vector mapping old -> new contrast names
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
  # rename contrasts vector
  g$contrasts <- unname(vapply(g$contrasts, function(x) if (x %in% old) new[match(x, old)] else x, character(1)))
  # update dimnames on assays if present
  g$assays <- lapply(g$assays, function(a) {
    dn <- dimnames(a)
    if (!is.null(dn) && length(dn) >= 3L) {
      dn[[3]] <- unname(vapply(dn[[3]], function(x) if (x %in% old) new[match(x, old)] else x, character(1)))
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
