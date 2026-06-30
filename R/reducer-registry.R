# Reducer registry -----------------------------------------------------------

.gds_reducers <- new.env(parent = emptyenv())

#' Register a reducer kernel
#'
#' @param name String identifier, e.g., "meta:fe"
#' @param fun Function(beta, var, X, z, p, df, df1, df2, opts) -> named list
#' @param requires Character vector of required inputs, e.g., c("beta","var")
#' @param provides Character vector of outputs to be written
#' @param options_schema Optional schema for options
#' @param input_shape Reducer execution mode: `"contrastwise"` (default) or
#'   `"joint_contrast"` for reducers that consume the full contrast axis jointly
#' @return Invisibly, the registered reducer `name`.
#' @export
register_reducer <- function(name,
                             fun,
                             requires,
                             provides,
                             options_schema = list(),
                             input_shape = c("contrastwise", "joint_contrast")) {
  stopifnot(is.character(name), length(name) == 1L, is.function(fun))
  input_shape <- match.arg(input_shape)
  .gds_reducers[[name]] <- list(
    name = name,
    fun = fun,
    requires = as.character(requires),
    provides = as.character(provides),
    options_schema = options_schema,
    input_shape = input_shape
  )
  invisible(name)
}

#' Get reducer by name
#'
#' @param name Reducer name
#'
#' @return Reducer list or NULL if not found
#' @export
get_reducer <- function(name) {
  .gds_reducers[[name]]
}

#' List registered reducers
#' @return A sorted character vector of registered reducer names.
#' @seealso [list_assays()] to preview a reducer's output assays.
#' @export
list_reducers <- function() {
  sort(ls(.gds_reducers))
}

# Internal: map legacy method names to registry ids
.normalize_reducer_name <- function(name) {
  switch(name,
    fixed = "meta:fe",
    random = "meta:re",
    stouffer = "combine:stouffer",
    fisher = "combine:fisher",
    name
  )
}
