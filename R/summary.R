#' Explain objects and plans
#'
#' Human-readable summary of a realised GDS or a lazy `gds_plan`.
#'
#' @param x A `gds` or `gds_plan`
#' @return Invisibly returns `x` (and prints a summary)
#' @export
explain <- function(x) UseMethod("explain")

.space_brief <- function(space) {
  if (inherits(space, "space_voxel")) {
    paste0("voxel dim=", paste(space$dim, collapse = "x"), 
           ", storage=", space$storage %||% "dense",
           if (!is.null(space$template_id)) paste0(", template=", space$template_id) else "")
  } else if (inherits(space, "space_parcels")) {
    n <- length(space$labels %||% list())
    paste0("parcels n=", n)
  } else if (inherits(space, "space_sample_labels")) {
    paste0("sample_labels n=", length(space$labels))
  } else if (inherits(space, "gds_space")) {
    paste0(space$type)
  } else {
    "unknown_space"
  }
}

#' @export
explain.gds <- function(x) {
  dims <- dim(assays(x)[[1]])
  cat("GDS\n")
  cat("  dims: ", paste(dims, collapse = " x "), " [sample x subject x contrast]", "\n", sep = "")
  cat("  space: ", .space_brief(x$space), "\n", sep = "")
  cat("  subjects: ", length(x$subjects), if (length(x$subjects)) paste0(" (", paste(head(x$subjects, 3), collapse = ", "), if (length(x$subjects) > 3) ", ..." else "", ")") else "", "\n", sep = "")
  cat("  contrasts: ", length(x$contrasts), if (length(x$contrasts)) paste0(" (", paste(head(x$contrasts, 3), collapse = ", "), if (length(x$contrasts) > 3) ", ..." else "", ")") else "", "\n", sep = "")
  ai <- lapply(names(assays(x)), assay_info)
  roles <- vapply(ai, function(a) a$role %||% "?", character(1))
  cat("  assays: ", paste(paste0(names(assays(x)), "[", roles, "]"), collapse = ", "), "\n", sep = "")
  fams <- names(x$metadata$map_families %||% list())
  if (length(fams)) cat("  map families: ", paste(fams, collapse = ", "), "\n", sep = "")
  if (!is.null(x$metadata$design_mats)) {
    for (i in seq_along(x$metadata$design_mats)) {
      dm <- x$metadata$design_mats[[i]]
      cols <- dm$columns %||% character()
      cat("  design[", i, "]: ", paste(cols, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}

#' @export
explain.gds_plan <- function(x) {
  probe <- x$source$probe
  cat("GDS plan\n")
  cat("  adapter: ", x$source$adapter, "\n", sep = "")
  if (!is.null(probe)) {
    dims <- probe$dims
    cat("  source dims: ", paste(as.integer(dims), collapse = " x "), " [sample x subject x contrast]", "\n", sep = "")
    cat("  space: ", .space_brief(probe$space), "\n", sep = "")
    cat("  assays: ", paste(probe$assays, collapse = ", "), "\n", sep = "")
  }
  cat("  nodes: ", length(x$nodes), "\n", sep = "")
  if (length(x$nodes)) {
    for (i in seq_along(x$nodes)) {
      n <- x$nodes[[i]]
      summary <- switch(n$op,
        reduce = paste0("reduce:", n$method),
        map = paste0("map combine=", n$combine %||% "none"),
        mask_policy = paste0("mask:", n$policy$rule),
        subset_axis = "subset",
        align_to_group = paste0("align:", n$family_name %||% n$family$name %||% "(anonymous)"),
        write = paste0("write:", n$format),
        n$op
      )
      cat("    ", i, ". ", summary, "\n", sep = "")
    }
  }
  invisible(x)
}

#' Validate a GDS or plan
#'
#' @param x A `gds`, `gds_plan`, or `image_catalog`
#' @param ... Additional arguments passed to methods
#' @return TRUE if valid; otherwise errors with a descriptive message
#' @export
validate <- function(x, ...) UseMethod("validate")

#' @export
validate.gds <- function(x, ...) {
  assays <- assays(x)
  .validate_gds_assays(assays)
  .validate_gds_dims(assays, x$subjects, x$contrasts)
  .validate_gds_space(x$space, dim(assays[[1L]])[1L])
  if (!is.null(x$col_data)) {
    rn <- rownames(x$col_data)
    if (is.null(rn) || !all(x$subjects %in% rn)) {
      stop("col_data must have rownames matching subjects", call. = FALSE)
    }
  }
  TRUE
}

#' @export
validate.gds_plan <- function(x, ...) {
  # Basic adapter + probe availability
  if (is.null(get_adapter(x$source$adapter))) stop("Unknown adapter: ", x$source$adapter, call. = FALSE)
  if (is.null(x$source$probe)) stop("Plan is missing probe metadata", call. = FALSE)
  # Check reduce nodes reference registered reducers
  for (n in x$nodes) {
    if (identical(n$op, "reduce")) {
      name <- n$method
      red <- get_reducer(.normalize_reducer_name(name))
      if (is.null(red)) stop("Unknown reducer: ", name, call. = FALSE)
    }
  }
  TRUE
}

#' Explain a plan's operations in a tidy table
#'
#' @param plan A `gds_plan`
#' @return A data.frame with one row per node (printed and returned invisibly)
#' @export
explain_plan <- function(plan) {
  plan <- as_plan(plan)
  rows <- lapply(seq_along(plan$nodes), function(i) {
    n <- plan$nodes[[i]]
    data.frame(
      index = i,
      op = n$op,
      summary = switch(n$op,
        reduce = paste0("method=", n$method),
        map = paste0("combine=", n$combine %||% "none"),
        mask_policy = paste0("rule=", n$policy$rule),
        subset_axis = "subset",
        align_to_group = paste0("family=", n$family_name %||% n$family$name %||% "(anonymous)"),
        write = paste0("format=", n$format),
        ""
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame(index = integer(), op = character(), summary = character())
  print(out)
  invisible(out)
}

#' Preview a small block through the plan
#'
#' Executes the plan with a small sample block to quickly inspect shapes or
#' values. If `assays` is provided, returns raw arrays from the adapter stage
#' (pre-operations); otherwise returns a realised GDS of the small block after
#' applying plan nodes.
#'
#' @param plan A `gds_plan` or object coercible via [as_plan()]
#' @param n Number of samples to preview (from start)
#' @param assays Optional assays to request as raw arrays
#'
#' @return Either a list of arrays (when `assays` provided) or a small `gds`
#' @export
preview <- function(plan, n = 3, assays = NULL) {
  plan <- as_plan(plan)
  n_total <- as.integer(plan$source$probe$dims[["sample"]] %||% dim(assays(compute(plan))[[1]])[1])
  n <- min(as.integer(n), n_total)
  blk <- list(sample = seq_len(n))
  if (!is.null(assays)) {
    return(compute(plan, block = blk, assays = assays))
  }
  compute(plan, block = blk)
}

