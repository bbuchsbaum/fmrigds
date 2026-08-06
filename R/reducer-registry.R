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
#' @param model_contract Optional declarative description of the reducer's
#'   design, weighting, missingness, variance, estimand, and deletion semantics.
#'   Custom reducers may omit this, in which case model-conditioned examination
#'   diagnostics are reported as unavailable.
#' @param diagnostics Optional diagnostic implementation declaration. This is a
#'   list containing `fun`, `capabilities`, and `modes`; `fun` may be `NULL`
#'   while a capability is declared but not yet implemented.
#' @return Invisibly, the registered reducer `name`.
#' @export
register_reducer <- function(name,
                             fun,
                             requires,
                             provides,
                             options_schema = list(),
                             input_shape = c("contrastwise", "joint_contrast"),
                             model_contract = NULL,
                             diagnostics = NULL) {
  stopifnot(is.character(name), length(name) == 1L, is.function(fun))
  input_shape <- match.arg(input_shape)
  model_contract <- .normalize_reducer_model_contract(model_contract)
  diagnostics <- .normalize_reducer_diagnostics(diagnostics)
  .gds_reducers[[name]] <- list(
    name = name,
    fun = fun,
    requires = as.character(requires),
    provides = as.character(provides),
    options_schema = options_schema,
    input_shape = input_shape,
    model_contract = model_contract,
    diagnostics = diagnostics
  )
  invisible(name)
}

.normalize_reducer_model_contract <- function(contract = NULL) {
  if (is.null(contract)) return(NULL)
  if (!is.list(contract)) {
    stop("Reducer model_contract must be NULL or a list.", call. = FALSE)
  }
  defaults <- list(
    uses_X = FALSE,
    estimands = "none",
    weight_mode = "model_specific",
    missingness = "model_specific",
    synthetic_variance = "allow_effect_only",
    deletion = "unsupported"
  )
  unknown <- setdiff(names(contract), names(defaults))
  if (length(unknown)) {
    stop("Unknown model_contract fields: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- utils::modifyList(defaults, contract)
  if (!is.logical(out$uses_X) || length(out$uses_X) != 1L || is.na(out$uses_X)) {
    stop("model_contract uses_X must be TRUE or FALSE.", call. = FALSE)
  }
  .contract_choice <- function(value, choices, field) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% choices) {
      stop(
        "model_contract ", field, " must be one of: ",
        paste(choices, collapse = ", "), ".",
        call. = FALSE
      )
    }
    value
  }
  out$estimands <- .contract_choice(
    out$estimands,
    c("none", "intercept", "linear"),
    "estimands"
  )
  out$weight_mode <- .contract_choice(
    out$weight_mode,
    c("inverse_variance", "unweighted", "evidence", "model_specific"),
    "weight_mode"
  )
  out$missingness <- .contract_choice(
    out$missingness,
    c("samplewise", "complete_case", "model_specific"),
    "missingness"
  )
  out$synthetic_variance <- .contract_choice(
    out$synthetic_variance,
    c("forbid", "allow_effect_only", "allow"),
    "synthetic_variance"
  )
  out$deletion <- .contract_choice(
    out$deletion,
    c("closed_form", "hat_matrix", "tau2_fixed_full", "model_specific", "unsupported"),
    "deletion"
  )
  out
}

.normalize_reducer_diagnostics <- function(diagnostics = NULL) {
  if (is.null(diagnostics)) return(NULL)
  if (!is.list(diagnostics)) {
    stop("Reducer diagnostics must be NULL or a list.", call. = FALSE)
  }
  defaults <- list(fun = NULL, capabilities = character(), modes = character())
  unknown <- setdiff(names(diagnostics), names(defaults))
  if (length(unknown)) {
    stop("Unknown reducer diagnostics fields: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- utils::modifyList(defaults, diagnostics)
  if (!is.null(out$fun) && !is.function(out$fun)) {
    stop("diagnostics fun must be NULL or a function.", call. = FALSE)
  }
  for (field in c("capabilities", "modes")) {
    value <- out[[field]]
    if (!is.character(value) || anyNA(value) || any(!nzchar(value))) {
      stop("diagnostics ", field, " must be a character vector.", call. = FALSE)
    }
    out[[field]] <- unique(value)
  }
  out
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
