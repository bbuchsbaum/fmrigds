.assay_registry <- new.env(parent = emptyenv())

#' Register an assay type
#'
#' @param name Assay name
#' @param role Semantic category (e.g., "location", "variance", "stdev", "z", "t", "F", "df",
#'   "n_eff", "p", "chi2", "evidence", "log_evidence", "posterior")
#' @param units Optional units string/list
#' @param variance_of Optional name of the assay whose variance this represents
#' @param derive_from Optional character vector of prerequisite assays
#'
#' @return Invisibly, the registered assay definition
#' @export
register_assay <- function(name,
                           role,
                           units = NULL,
                           variance_of = NULL,
                           derive_from = NULL) {
  role <- match.arg(
    role,
    c(
      "location",          # effect scale (beta / log-odds / etc.)
      "variance",          # variance of a location statistic
      "stdev",             # standard deviation of a location statistic
      "z", "t", "F",      # classic test-statistics families
      "df", "n_eff",      # degrees of freedom / effective sample size
      "p", "chi2",        # evidence-on-null metrics
      "evidence",         # generic evidence score (e.g., log likelihood ratio)
      "log_evidence",     # log-scale evidence values (e.g., log BF)
      "posterior"         # posterior summaries (e.g., posterior probability)
    )
  )
  .assay_registry[[name]] <- list(
    name = name,
    role = role,
    units = units,
    variance_of = variance_of,
    derive_from = derive_from
  )
  invisible(.assay_registry[[name]])
}

#' Retrieve assay metadata
#' @export
assay_info <- function(name) {
  .assay_registry[[name]]
}

#' Test whether an assay can be linearly mapped
#' @export
can_map_linear <- function(name) {
  info <- assay_info(name)
  if (is.null(info)) return(FALSE)
  info$role %in% c("location", "variance", "stdev")
}

.register_default_assays <- function() {
  register_assay("beta", role = "location", units = "%BOLD")
  register_assay("var", role = "variance", units = "%BOLD^2", variance_of = "beta")
  register_assay("se", role = "stdev", units = "%BOLD", derive_from = "var")
  register_assay("t", role = "t", derive_from = c("beta", "var", "df"))
  register_assay("z", role = "z", derive_from = c("t", "df"))
  register_assay("F", role = "F", derive_from = c("beta", "var", "df1", "df2"))
  register_assay("df", role = "df")
  register_assay("df1", role = "df")
  register_assay("df2", role = "df")
  register_assay("n_eff", role = "n_eff")
  register_assay("p", role = "p")
  register_assay("chi2", role = "chi2")
  register_assay("logBF", role = "log_evidence", units = "log BF")
}
