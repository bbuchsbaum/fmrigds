#' Create a linear map between spaces
#'
#' @param from_space Source space descriptor or name
#' @param to_space Target space descriptor or name
#' @param operator Mapping operator (matrix/function)
#' @param by_subject Optional named list of subject-specific operator descriptors/functions
#' @param traits List of map traits (orthogonal, mass_preserving)
#' @param uncertainty [`UncertaintyRule`] controlling variance propagation
#' @param name Optional map name
#'
#' @return Object of class `c("map_linear", "gds_map")`
#' @export
map_linear <- function(from_space,
                       to_space,
                       operator,
                       by_subject = NULL,
                       traits = list(orthogonal = FALSE, mass_preserving = FALSE),
                       uncertainty = UncertaintyRule("independent"),
                       name = NULL) {
  .validate_space_reference(from_space)
  .validate_space_reference(to_space)
  .validate_operator(operator)
  .validate_subject_operators(by_subject)

  structure(
    list(
      name = name %||% .derive_map_name(from_space, to_space),
      from = from_space,
      to = to_space,
      operator = operator,
      by_subject = by_subject,
      traits = traits,
      uncertainty = uncertainty
    ),
    class = c("map_linear", "gds_map")
  )
}

#' Create a subject-aware map family
#'
#' @param name Family name
#' @param from_space Source space
#' @param to_space Target space
#' @param type Map type identifier
#' @param by_subject Named list of per-subject operator descriptors
#' @param traits Map traits
#' @param uncertainty [`UncertaintyRule`]
#'
#' @return Object of class `c("gds_map_family", "gds_map")`
#' @export
MapFamily <- function(name,
                      from_space,
                      to_space,
                      type = c("linear", "orthogonal", "affine3d", "deform3d", "ot"),
                      by_subject,
                      traits = list(orthogonal = FALSE, mass_preserving = FALSE),
                      uncertainty = UncertaintyRule("independent")) {
  type <- match.arg(type)
  .validate_space_reference(from_space)
  .validate_space_reference(to_space)
  .validate_subject_operators(by_subject, require_named = TRUE)

  structure(
    list(
      name = name,
      from = from_space,
      to = to_space,
      type = type,
      by_subject = by_subject,
      traits = traits,
      uncertainty = uncertainty
    ),
    class = c("gds_map_family", "gds_map")
  )
}

#' Orthogonal map family helper
#' @export
OrthogonalFamily <- function(name,
                             from_space,
                             to_space,
                             matrices_by_subject,
                             uncertainty = UncertaintyRule("independent")) {
  MapFamily(
    name = name,
    from_space = from_space,
    to_space = to_space,
    type = "orthogonal",
    by_subject = matrices_by_subject,
    traits = list(orthogonal = TRUE, mass_preserving = FALSE),
    uncertainty = uncertainty
  )
}

#' Optimal transport family helper
#' @export
OTFamily <- function(name,
                     from_space,
                     to_space,
                     plans_by_subject,
                     uncertainty = UncertaintyRule("independent")) {
  MapFamily(
    name = name,
    from_space = from_space,
    to_space = to_space,
    type = "ot",
    by_subject = plans_by_subject,
    traits = list(orthogonal = FALSE, mass_preserving = TRUE),
    uncertainty = uncertainty
  )
}

#' Deformable warp family helper
#' @export
WarpFamily <- function(name,
                       from_space,
                       to_space,
                       warps_by_subject,
                       uncertainty = UncertaintyRule("independent")) {
  MapFamily(
    name = name,
    from_space = from_space,
    to_space = to_space,
    type = "deform3d",
    by_subject = warps_by_subject,
    traits = list(orthogonal = FALSE, mass_preserving = FALSE),
    uncertainty = uncertainty
  )
}

#' Specify how uncertainty is propagated through a map
#'
#' @param mode One of "independent", "cov_provider", "kernel", "none"
#' @param df_rule Degrees of freedom aggregation rule
#' @param cov_provider Optional covariance provider function
#' @param kernel Optional kernel function
#'
#' @return Uncertainty rule descriptor
#' @export
UncertaintyRule <- function(mode = c("independent", "cov_provider", "kernel", "none"),
                            df_rule = c("satterthwaite", "none"),
                            cov_provider = NULL,
                            kernel = NULL) {
  mode <- match.arg(mode)
  df_rule <- match.arg(df_rule)

  if (mode == "cov_provider" && is.null(cov_provider)) {
    stop("`cov_provider` must be supplied when mode = 'cov_provider'", call. = FALSE)
  }
  if (mode == "kernel" && is.null(kernel)) {
    stop("`kernel` must be supplied when mode = 'kernel'", call. = FALSE)
  }

  structure(
    list(
      mode = mode,
      df_rule = df_rule,
      cov_provider = cov_provider,
      kernel = kernel
    ),
    class = "gds_uncertainty_rule"
  )
}

# -------------------------------------------------------------------------
# Helpers ------------------------------------------------------------------

.validate_space_reference <- function(x) {
  if (!inherits(x, "gds_space") && !is.character(x)) {
    stop("Space references must be 'gds_space' objects or character identifiers", call. = FALSE)
  }
  invisible(NULL)
}

.validate_operator <- function(operator) {
  if (is.null(operator)) return(invisible(NULL))
  if (is.matrix(operator) || is.function(operator)) return(invisible(NULL))
  stop("`operator` must be a matrix or function", call. = FALSE)
}

.validate_subject_operators <- function(by_subject, require_named = FALSE) {
  if (is.null(by_subject)) return(invisible(NULL))
  if (!is.list(by_subject) || !length(by_subject)) {
    stop("`by_subject` must be a non-empty named list", call. = FALSE)
  }
  if (is.null(names(by_subject)) || any(!nzchar(names(by_subject)))) {
    stop("`by_subject` must be named by subject id", call. = FALSE)
  }
  invisible(NULL)
}

.derive_map_name <- function(from_space, to_space) {
  from <- if (inherits(from_space, "gds_space")) from_space$type else as.character(from_space)
  to <- if (inherits(to_space, "gds_space")) to_space$type else as.character(to_space)
  paste0(from, "_to_", to)
}

`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}
