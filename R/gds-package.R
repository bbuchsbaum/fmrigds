#' fmrigds
#'
#' Lazy, format-agnostic group analysis for fMRI.
#'
#' `fmrigds` provides a common Group Data Set (GDS) abstraction for first-level
#' fMRI outputs together with a lazy, plan-based workflow for group analysis.
#' The package centers on a compact public grammar:
#'
#' - `gds()` to open supported sources and create a plan
#' - `subset()`, `derive()`, `align()`, `mask()`, `map_to()`, `reduce()`,
#'   `posthoc()`, and `write_out()` to describe analysis steps
#' - `compute()` to execute the plan and return a realised GDS
#'
#' Core capabilities include:
#'   * Constructors for GDS objects and spaces
#'   * Storage adapters (tabular, NIfTI, HDF5, fmristore)
#'   * Statistical derivations and variance propagation helpers
#'   * Reducer and post-hoc registries for group and meta-analytic workflows
#'   * Provenance-aware export and inspection helpers
#'
#' See the package reference and vignettes for supported workflows and release
#' guidance.
#' @importFrom stats setNames
#' @importFrom utils head
#' @useDynLib fmrigds, .registration = TRUE
"_PACKAGE"

.onLoad <- function(libname, pkgname) { # nocov start
  .register_default_assays()
  register_builtin_adapters()
  if (exists("register_core_reducers", mode = "function")) {
    register_core_reducers()
  }
  if (exists("register_lmm_reducers", mode = "function")) {
    register_lmm_reducers()
  }
  if (exists(".register_builtin_posthoc", mode = "function")) {
    .register_builtin_posthoc()
  }
  if (exists(".set_threads_from_option", mode = "function")) {
    .set_threads_from_option()
  }
} # nocov end
