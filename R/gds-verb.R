#' Adapter front door for building a plan
#'
#' @param source Source specification (paths, etc.)
#' @param format Optional adapter name or "auto"
#' @param prefer Optional preferred adapter when multiple match
#' @param ... Adapter-specific arguments forwarded to `open()`/`probe()`
#'
#' @return A [`gds_plan`]
#' @export
gds <- function(source,
                format = c("auto", ls(.adapter_registry)),
                prefer = NULL,
                ...) {
  format <- match.arg(format)
  adapter_name <- if (identical(format, "auto")) detect_adapter(source, prefer) else format
  adapter <- get_adapter(adapter_name)

  handle <- adapter$open(source, ...)
  probe_result <- adapter$probe(handle, ...)
  adapter$close(handle)

  adapter_columns <- probe_result$columns
  probe_result$columns <- NULL

  src <- gds_source(adapter_name, source, probe_result)
  plan <- gds_plan(source = src)
  plan$meta$adapter_columns <- adapter_columns
  plan
}
