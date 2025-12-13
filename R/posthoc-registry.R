# Post-hoc registry -----------------------------------------------------------

.gds_posthoc <- new.env(parent = emptyenv())

#' Register a post-hoc method
#'
#' Register or replace a post-hoc handler identified by `name`.
#'
#' The handler `fun` is called as `fun(arrays, opts)` where `arrays` is a
#' named list of assays shaped `sample × subject × contrast`. It should return
#' either (a) a named list of new/updated assays to merge (e.g., `list(q = ...)`),
#' or (b) the full `arrays` list with modifications. Implementations must
#' preserve the input array dimensions.
#'
#' @param name String identifier, e.g., "fdr:bh"
#' @param fun Function(arrays, opts) -> list of arrays (delta or full arrays)
#' @param requires Character vector of required input assays, e.g., c("p")
#' @param provides Character vector of output assay names, e.g., c("q")
#' @param overwrite Logical; if FALSE and a handler with `name` exists, error
#' @export
register_posthoc <- function(name, fun, requires, provides, overwrite = TRUE) {
  stopifnot(is.character(name), length(name) == 1L, is.function(fun))
  if (!isTRUE(overwrite) && !is.null(.gds_posthoc[[name]])) {
    stop(sprintf("Post-hoc method '%s' already exists; set overwrite = TRUE to replace.", name), call. = FALSE)
  }
  .gds_posthoc[[name]] <- list(
    name = name,
    fun = fun,
    requires = as.character(requires),
    provides = as.character(provides)
  )
  invisible(name)
}

#' Get post-hoc method by name
#'
#' @param name Method identifier (e.g., "fdr:bh")
#' @return A list with fields `name`, `fun`, `requires`, `provides`, or `NULL` if not found
#' @export
get_posthoc <- function(name) {
  .gds_posthoc[[name]]
}

#' List registered post-hoc methods
#' @export
list_posthoc <- function() {
  sort(ls(.gds_posthoc))
}

#' Unregister a post-hoc method
#'
#' @param name Method identifier to remove
#' @return Invisibly, TRUE if removed, FALSE if not found
#' @export
unregister_posthoc <- function(name) {
  if (is.null(.gds_posthoc[[name]])) return(invisible(FALSE))
  rm(list = name, envir = .gds_posthoc)
  invisible(TRUE)
}

# Built-ins: FDR BH/BY --------------------------------------------------------

.posthoc_fdr_template <- function(method) {
  force(method)
  function(arrays, opts) {
    if (is.null(arrays$p)) {
      # try to derive p from z or t
      if (!is.null(arrays$z)) {
        arrays$p <- derive_p(list(z = arrays$z))
      } else if (!is.null(arrays$t) && !is.null(arrays$df)) {
        arrays$p <- derive_p(list(t = arrays$t, df = arrays$df))
      } else {
        stop("post-hoc FDR requires p values (or z/t with df)", call. = FALSE)
      }
    }
    p <- arrays$p
    dims <- dim(p)
    out <- array(NA_real_, dim = dims)
    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        pv <- p[, j, k]
        adj <- stats::p.adjust(pv, method = method)
        out[, j, k] <- adj
      }
    }
    list(q = out)
  }
}

.register_builtin_posthoc <- function() {
  register_posthoc("fdr:bh", .posthoc_fdr_template("BH"), requires = c("p"), provides = c("q"))
  register_posthoc("fdr:by", .posthoc_fdr_template("BY"), requires = c("p"), provides = c("q"))
  # Register spatial FDR if available
  if (exists(".posthoc_spatial_fdr", mode = "function")) {
    register_posthoc("fdr:spatial", .posthoc_spatial_fdr(), requires = c("p"), provides = c("q"))
  }
}

# initialize defaults when package loads
NULL
