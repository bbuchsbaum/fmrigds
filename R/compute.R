#' Materialise a plan into a realised GDS
#'
#' @param x Plan or source
#' @param sink Storage sink ("memory" only at present)
#' @param ... Reserved for future use
#'
#' @return A realised [`gds`] object
#' @export
compute <- function(x, sink = c("memory"), ...) {
  sink <- match.arg(sink)
  plan <- as_plan(x)
  plan <- optimize_plan(plan)

  adapter <- get_adapter(plan$source$adapter)
  handle <- adapter$open(plan$source$source)
  on.exit({
    if (!is.null(adapter$close)) adapter$close(handle)
  }, add = TRUE)

  assays <- plan$source$probe$assays
  columns <- plan$meta$adapter_columns %||% list(effect_cols = NULL, subject_col = NULL, sample_col = NULL, contrast_col = NULL)

  arrays <- adapter$read(
    handle,
    assays = assays,
    block = NULL,
    effect_cols = columns$effect_cols,
    subject_col = columns$subject_col,
    sample_col = columns$sample_col,
    contrast_col = columns$contrast_col,
    mask_idx = plan$source$probe$mask_idx,
    spatial_dim = plan$source$probe$spatial_dim,
    n_contrasts = plan$source$probe$n_contrasts
  )

  node_result <- .apply_plan_nodes(arrays, plan, plan$source$probe$space)
  arrays <- node_result$arrays

  subjects <- plan$source$probe$subjects
  contrasts <- plan$source$probe$contrasts
  space <- node_result$space
  subset_info <- node_result$subset
  if (!is.null(subset_info$subjects)) subjects <- subjects[subset_info$subjects]
  if (!is.null(subset_info$contrasts)) contrasts <- contrasts[subset_info$contrasts]
  if (!is.null(subset_info$samples)) space <- .subset_space(space, subset_info$samples)

  gds <- new_gds(
    assays = arrays,
    space = space,
    subjects = subjects,
    contrasts = contrasts,
    metadata = plan$source$probe$metadata
  )

  gds$metadata$provenance$digest <- digest_plan(plan)
  gds
}

#' Compute a stable digest for a plan
#' @export
digest_plan <- function(plan) {
  canonical <- list(
    source_hash = plan$source$hash,
    nodes = lapply(plan$nodes, canonicalize_node)
  )
  digest::digest(canonical, algo = "xxhash64")
}

#' Canonicalize an operation node
#' @export
canonicalize_node <- function(node) {
  params <- node[names(node) != "op"]
  params <- params[order(names(params))]
  list(op = node$op, params = params)
}

# -------------------------------------------------------------------------
# Helpers ------------------------------------------------------------------

.apply_plan_nodes <- function(arrays, plan, space) {
  subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
  current_space <- space
  for (node in plan$nodes) {
    op <- node$op
    if (op == "subset_axis") {
      res <- .apply_subset_node(arrays, plan, node)
      arrays <- res$arrays
      subset_info <- res$subset
    } else if (op == "derive") {
      arrays <- execute_derive(arrays, node$what, node$options)
    } else if (op == "map") {
      res <- apply_map_to(node, arrays)
      arrays <- res$arrays
      current_space <- res$space
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else {
      # TODO: implement derive/map/reduce execution in future sprints
    }
  }
  list(arrays = arrays, subset = subset_info, space = current_space)
}

.apply_subset_node <- function(arrays, plan, node) {
  sample <- node$sample
  subject <- node$subject
  contrast <- node$contrast

  dims <- dim(arrays[[1]])
  samples_idx <- if (is.null(sample)) seq_len(dims[1]) else .coerce_index(sample, plan$source$probe$space, dims[1])
  subjects_idx <- if (is.null(subject)) seq_len(dims[2]) else match(subject, plan$source$probe$subjects)
  contrasts_idx <- if (is.null(contrast)) seq_len(dims[3]) else match(contrast, plan$source$probe$contrasts)

  if (any(is.na(subjects_idx))) stop("Unknown subject in subset operation", call. = FALSE)
  if (any(is.na(contrasts_idx))) stop("Unknown contrast in subset operation", call. = FALSE)

  arrays <- lapply(arrays, function(a) a[samples_idx, subjects_idx, contrasts_idx, drop = FALSE])
  list(arrays = arrays, subset = list(samples = samples_idx, subjects = subjects_idx, contrasts = contrasts_idx))
}

.coerce_index <- function(idx, space, n_samples) {
  if (is.numeric(idx)) return(idx)
  if (is.logical(idx)) return(which(idx))
  if (is.character(idx) && inherits(space, "space_parcels")) {
    return(match(idx, space$labels))
  }
  stop("Unsupported index type for samples", call. = FALSE)
}

.subset_space <- function(space, idx) {
  if (inherits(space, "space_parcels")) {
    space$labels <- space$labels[idx]
  }
  space
}

# Placeholder optimizer ----------------------------------------------------

#' Plan optimizer (placeholder)
#' @export
optimize_plan <- function(plan) {
  plan
}
