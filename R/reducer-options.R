# Reducer options validation utilities (internal)

#' Validate reducer options against schema
#' @keywords internal
validate_reducer_options <- function(schema, opts) {
  opts <- opts %||% list()
  if (is.null(schema) || !length(schema)) return(opts)

  for (nm in names(schema)) {
    allowed <- schema[[nm]]
    # enum-like schema: character vector of allowed values
    if (is.character(allowed)) {
      if (is.null(opts[[nm]])) {
        # default to first allowed value
        opts[[nm]] <- allowed[[1L]]
      } else if (!opts[[nm]] %in% allowed) {
        stop(sprintf("Invalid option '%s' value '%s'. Allowed: %s. How to fix: choose one of the allowed values.",
                     nm, as.character(opts[[nm]]), paste(allowed, collapse = ", ")), call. = FALSE)
      }
    }
  }
  opts
}

