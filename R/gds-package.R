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

# Internal mutable package state (once-per-session warning flags, etc.).
.gds_state <- new.env(parent = emptyenv())

# Emit the synthetic unit-variance warning. With `once = TRUE` it fires at most
# once per session (used by the lazy NIfTI read path, which re-reads the `var`
# assay on every compute()). The message keeps the legacy substring so existing
# `fixed=` matchers still match, and points users at the unweighted path.
.warn_synthetic_variance <- function(once = FALSE) {
  if (once && isTRUE(.gds_state$warned_unit_variance)) {
    return(invisible(FALSE))
  }
  if (once) .gds_state$warned_unit_variance <- TRUE
  warning(
    "No variance or SE provided; using unit variance as a synthetic placeholder. ",
    "This is only valid for unweighted reducers such as method = \"ols:voxelwise\" ",
    "(e.g. one_sample()/group_ols()); variance-weighted reducers (fixed/random/meta:*) ",
    "will refuse it. Supply real SE (e.g. nifti_source(se = ...)) for meta-analysis.",
    call. = FALSE
  )
  invisible(TRUE)
}

# Test/inspection helper: reset the once-per-session flag.
.reset_synthetic_variance_warning <- function() {
  .gds_state$warned_unit_variance <- NULL
  invisible(NULL)
}

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
  if (exists(".register_neuropublish_method", mode = "function")) {
    .register_neuropublish_method()
  }
} # nocov end
