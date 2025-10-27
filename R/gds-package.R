#' gdsfmri Package
#'
#' Scaffolding for the Group Data Set (GDS) implementation as specified in
#' `TECHNICAL_SPECIFICATION.md` and summarised in `fmrigds_blueprint.md`.
#'
#' The functions exported from this package will implement:
#'   * Constructors for GDS objects and spaces
#'   * Lazy plan verbs (`gds()`, `subset()`, `derive()`, `align()`, `mask()`,
#'     `map_to()`, `reduce()`, `write_out()`, `compute()`)
#'   * Storage adapters (tabular, NIfTI, HDF5, fmristore)
#'   * Statistical derivations and variance propagation helpers
#'
#' The current file provides package-level documentation placeholder so the
#' skeleton can be installed. Concrete implementations will follow the
#' technical specification.
"_PACKAGE"

#' @keywords internal
".gdsfmri"

.onLoad <- function(libname, pkgname) {
  .register_default_assays()
  register_builtin_adapters()
}
