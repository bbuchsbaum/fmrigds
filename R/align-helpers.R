#' Convenience helpers to create and register alignments
#'
#' @param name Family name
#' @param from Source space
#' @param to Target space
#' @param operators_by_subject Named list of per-subject linear operators (matrices)
#' @param uncertainty Uncertainty rule for variance propagation
#'
#' @return A MapFamily
#' @export
#' @examples
#' from <- space_parcels(c("R1","R2"))
#' to   <- space_parcels(c("G1","G2"))
#' ops <- list(s1 = diag(2), s2 = diag(2))
#' fam <- make_linear_family("lin", from, to, ops)
make_linear_family <- function(name,
                               from,
                               to,
                               operators_by_subject,
                               uncertainty = UncertaintyRule("independent")) {
  MapFamily(
    name = name,
    from_space = from,
    to_space = to,
    type = "linear",
    by_subject = operators_by_subject,
    uncertainty = uncertainty
  )
}

#' Create a warp-based alignment family from on-disk paths
#'
#' @param name Family name
#' @param from Source space
#' @param to Target space
#' @param warp_paths Named character vector of file paths per subject
#' @param loader Function(path) -> operator (e.g., matrix) for a single subject
#' @param uncertainty Uncertainty rule
#'
#' @return A MapFamily (type = deform3d)
#' @export
#' @examples
#' loader <- function(path) as.matrix(utils::read.csv(path, header = FALSE))
#' # paths <- c(s1 = "warp_s1.csv", s2 = "warp_s2.csv")
#' # fam <- make_warp_family("warp", from, to, paths, loader)
make_warp_family <- function(name,
                             from,
                             to,
                             warp_paths,
                             loader,
                             uncertainty = UncertaintyRule("independent")) {
  stopifnot(is.character(warp_paths), !is.null(names(warp_paths)))
  ops <- lapply(warp_paths, loader)
  WarpFamily(name, from_space = from, to_space = to, warps_by_subject = ops, uncertainty = uncertainty)
}

#' Register, list, and get alignments (map families)
#'
#' Syntactic sugar around the map-family registry helpers.
#'
#' @name alignment-registry
NULL

#' @rdname alignment-registry
#' @param x A `gds` or `gds_plan` object
#' @param family A `MapFamily` to register
#' @param overwrite Logical; overwrite existing family with the same name
#' @return The updated object with the family registered
#' @export
register_alignment <- function(x, family, overwrite = FALSE) register_map(x, family, overwrite = overwrite)

#' @rdname alignment-registry
#' @return A character vector of registered alignment family names
#' @export
list_alignments <- function(x) list_map_families(x)

#' @rdname alignment-registry
#' @param name Alignment family name
#' @return The `MapFamily` entry or `NULL` if not found
#' @export
get_alignment <- function(x, name) get_map_family(x, name)
