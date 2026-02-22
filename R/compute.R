#' Materialise a plan into a realised GDS
#'
#' @param x Plan or source
#' @param sink Storage sink ("memory" for in-memory data, "h5" to write a GDS HDF5 store)
#' @param path Output path when `sink = "h5"`
#' @param sink_options Optional list of sink-specific options
#' @param block Optional list specifying block indices to request from the adapter
#' @param assays Optional character vector to return raw arrays instead of GDS object
#' @param ... Reserved for future use
#'
#' @return A realised [`gds`] object
#' @name compute
#' @export
compute <- function(x,
                    sink = c("memory", "h5"),
                    path = NULL,
                    sink_options = list(),
                    block = NULL,
                    assays = NULL,
                    ...) {
  sink <- match.arg(sink)
  plan <- as_plan(x)
  plan <- optimize_plan(plan)

  adapter <- get_adapter(plan$source$adapter)
  handle <- adapter$open(plan$source$source)
  on.exit({
    if (!is.null(adapter$close)) adapter$close(handle)
  }, add = TRUE)

  all_assays <- plan$source$probe$assays
  columns <- plan$meta$adapter_columns %||% list(effect_cols = NULL, subject_col = NULL, sample_col = NULL, contrast_col = NULL)

  # If assays explicitly requested, return raw arrays (adapter-level read)
  if (!is.null(assays)) {
    return(adapter$read(
      handle,
      assays = assays,
      block = block,
      effect_cols = columns$effect_cols,
      subject_col = columns$subject_col,
      sample_col = columns$sample_col,
      contrast_col = columns$contrast_col,
      mask_idx = plan$source$probe$mask_idx,
      spatial_dim = plan$source$probe$spatial_dim,
      n_contrasts = plan$source$probe$n_contrasts,
      temporal_policy = plan$meta$temporal_policy %||% NULL,
      contrast_matrix = plan$meta$contrast_matrix %||% NULL,
      contrast_names = plan$meta$contrast_names %||% NULL
    ))
  }

  arrays <- adapter$read(
    handle,
    assays = all_assays,
    block = block,
    effect_cols = columns$effect_cols,
    subject_col = columns$subject_col,
    sample_col = columns$sample_col,
    contrast_col = columns$contrast_col,
    mask_idx = plan$source$probe$mask_idx,
    spatial_dim = plan$source$probe$spatial_dim,
    n_contrasts = plan$source$probe$n_contrasts,
    temporal_policy = plan$meta$temporal_policy %||% NULL,
    contrast_matrix = plan$meta$contrast_matrix %||% NULL,
    contrast_names = plan$meta$contrast_names %||% NULL
  )

  node_result <- .apply_plan_nodes(
    arrays,
    plan,
    plan$source$probe$space,
    plan$source$probe$subjects,
    col_data = plan$meta$col_data %||% NULL
  )
  arrays <- node_result$arrays

  subjects <- node_result$subjects %||% plan$source$probe$subjects
  contrasts <- plan$source$probe$contrasts
  space <- node_result$space
  subset_info <- node_result$subset
  if (!is.null(subset_info$subjects)) subjects <- subjects[subset_info$subjects]
  if (!is.null(subset_info$contrasts)) contrasts <- contrasts[subset_info$contrasts]
  if (!is.null(subset_info$samples)) space <- .subset_space(space, subset_info$samples)
  # Harmonize col_data with current subjects; drop if mismatched (e.g., after reduce to group-level)
  col_data_final <- plan$meta$col_data %||% NULL
  if (!is.null(col_data_final)) {
    rn <- rownames(col_data_final)
    if (length(subjects) == 1L && identical(subjects, "meta")) {
      col_data_final <- NULL
    } else if (is.null(rn) || !all(subjects %in% rn)) {
      col_data_final <- NULL
    } else {
      col_data_final <- col_data_final[subjects, , drop = FALSE]
    }
  }

  gds <- new_gds(
    assays = arrays,
    space = space,
    subjects = subjects,
    contrasts = contrasts,
    col_data = col_data_final,
    metadata = plan$source$probe$metadata
  )

  if (length(node_result$designs)) {
    existing_designs <- gds$metadata$design_mats %||% list()
    gds$metadata$design_mats <- c(existing_designs, node_result$designs)
  }
  if (length(node_result$attachments)) {
    existing_att <- gds$metadata$attachments %||% list()
    gds$metadata$attachments <- utils::modifyList(existing_att, node_result$attachments)
  }

  registry <- list()
  if (length(plan$source$probe$maps %||% list())) {
    registry[names(plan$source$probe$maps)] <- plan$source$probe$maps
  }
  if (length(plan$meta$map_families %||% list())) {
    registry[names(plan$meta$map_families)] <- plan$meta$map_families
  }
  if (length(registry)) {
    gds$metadata$map_families <- registry
  }

  for (node in plan$nodes) {
    params <- switch(node$op,
      subset_axis = list(sample = node$sample, subject = node$subject, contrast = node$contrast),
      align_to_group = list(
        family = node$family_name %||% node$family$name %||% NA_character_,
        type = node$type %||% node$family$type %||% NA_character_
      ),
      mask_policy = list(scope = node$policy$scope, rule = node$policy$rule, threshold = node$policy$threshold),
      map = list(combine = node$combine, uncertainty = node$uncertainty$mode),
      reduce = list(method = node$method, weights = node$weights, by = node$by),
      posthoc = list(method = node$method),
      write = list(path = node$path, format = node$format),
      list()
    )
    gds$metadata <- add_provenance_node(gds$metadata, node$op, params)
  }
  gds$metadata$provenance$digest <- digest_plan(plan)

  # Re-run plan nodes that depend on col_data (e.g., build X from formula) if needed
  # Build X during reduce execution using gds$col_data if available.

  if (sink == "h5") {
    if (is.null(path) || !nzchar(path)) {
      stop("`path` must be provided when sink = 'h5'", call. = FALSE)
    }
    .write_export(gds, "h5", path, sink_options)
  }

  if (length(node_result$writes)) {
    for (node in node_result$writes) {
      .write_export(gds, node$format, node$path, node$options)
    }
  }

  gds
}

#' Compute a stable digest for a plan
#'
#' @param plan A gds_plan object
#'
#' @return Character digest hash
#' @export
digest_plan <- function(plan) {
  canonical <- list(
    source_hash = plan$source$hash,
    nodes = lapply(plan$nodes, canonicalize_node)
  )
  digest::digest(canonical, algo = "xxhash64")
}

#' Canonicalize an operation node
#'
#' @param node Operation node list
#'
#' @return Canonicalized node list
#' @export
canonicalize_node <- function(node) {
  params <- node[names(node) != "op"]
  params <- params[order(names(params))]
  list(op = node$op, params = params)
}

# -------------------------------------------------------------------------
# Helpers ------------------------------------------------------------------

.apply_plan_nodes <- function(arrays, plan, space, subjects, col_data = NULL) {
  subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
  current_space <- space
  current_subjects <- subjects
  current_contrasts <- plan$source$probe$contrasts
  writes <- list()
  designs <- list()
  attachments <- list()
  registry <- list()
  source_maps <- plan$source$probe$maps %||% list()
  if (length(source_maps)) registry[names(source_maps)] <- source_maps
  meta_maps <- plan$meta$map_families %||% list()
  if (length(meta_maps)) registry[names(meta_maps)] <- meta_maps

  for (node in plan$nodes) {
    op <- node$op
    if (op == "subset_axis") {
      res <- .apply_subset_node(arrays, plan, node)
      arrays <- res$arrays
      subset_info <- res$subset
      if (!is.null(subset_info$subjects)) {
        current_subjects <- current_subjects[subset_info$subjects]
      }
      if (!is.null(subset_info$contrasts)) {
        current_contrasts <- current_contrasts[subset_info$contrasts]
      }
    } else if (op == "derive") {
      arrays <- execute_derive(arrays, node$what, node$options)
    } else if (op == "align_to_group") {
      family <- node$family
      if (is.null(family)) {
        family <- registry[[node$family_name]]
      }
      if (is.null(family)) {
        stop("Align operation requires a registered map family '", node$family_name, "'", call. = FALSE)
      }
      res <- apply_align(family, arrays, current_subjects, current_space)
      arrays <- res$arrays
      current_space <- res$space
      current_subjects <- current_subjects
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else if (op == "mask_policy") {
      res <- apply_mask_policy(node, arrays, current_space)
      arrays <- res$arrays
      current_space <- res$space
      subset_info <- res$subset
    } else if (op == "map") {
      res <- apply_map_to(node, arrays)
      arrays <- res$arrays
      current_space <- res$space
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else if (op == "reduce") {
      res <- apply_reduce(node, arrays, node$weights, current_subjects, col_data)
      arrays <- res$arrays
      current_subjects <- res$subjects %||% "meta"
       if (!is.null(res$design_info)) {
        designs <- c(designs, list(res$design_info))
      }
      if (!is.null(res$attachments) && length(res$attachments)) {
        attachments <- c(attachments, res$attachments)
      }
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else if (op == "posthoc") {
      res <- apply_posthoc(
        node,
        arrays,
        context = list(space = current_space, subjects = current_subjects, contrasts = current_contrasts)
      )
      arrays <- res$arrays
      if (!is.null(res$attachments) && length(res$attachments)) {
        attachments <- c(attachments, res$attachments)
      }
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else if (op == "write") {
      writes <- c(writes, list(node))
    } else {
      # TODO: implement derive/map/reduce execution in future sprints
    }
  }
  list(
    arrays = arrays,
    subset = subset_info,
    space = current_space,
    subjects = current_subjects,
    writes = writes,
    designs = designs,
    attachments = attachments
  )
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

## optimize_plan is implemented in R/plan-optimizer.R
