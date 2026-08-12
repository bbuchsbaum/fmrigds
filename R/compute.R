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
  if (inherits(x, "fmri_group_plan")) {
    if (!identical(match.arg(sink), "memory") || !is.null(path) ||
        length(sink_options) || !is.null(block) || !is.null(assays)) {
      stop(
        "Frame group plans currently support only compute(plan) to an in-memory result frame.",
        call. = FALSE
      )
    }
    return(.compute_group_plan(x))
  }
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

  read_assays <- .initial_assays_for_plan(plan, all_assays)
  arrays <- adapter$read(
    handle,
    assays = read_assays,
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
    col_data = col_data(plan),
    row_data = row_data(plan),
    contrast_data = contrast_data(plan)
  )
  arrays <- node_result$arrays

  subjects <- node_result$subjects %||% plan$source$probe$subjects
  contrasts <- node_result$contrasts %||% plan$source$probe$contrasts
  space <- node_result$space
  subset_info <- node_result$subset
  row_data_final <- if ("row_data" %in% names(node_result)) {
    node_result$row_data
  } else {
    row_data(plan)
  }
  if (!("row_data" %in% names(node_result)) && !is.null(row_data_final) && !is.null(subset_info$samples)) {
    row_data_final <- row_data_final[subset_info$samples, , drop = FALSE]
  }
  # Harmonize col_data with current subjects; drop if mismatched (e.g., after reduce to group-level)
  col_data_final <- col_data(plan)
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
    row_data = row_data_final,
    metadata = plan$source$probe$metadata
  )

  contrast_data_final <- if ("contrast_data" %in% names(node_result)) {
    node_result$contrast_data
  } else {
    contrast_data(plan)
  }
  if (!is.null(contrast_data_final)) {
    if (!is.null(contrasts) && nrow(contrast_data_final) == length(contrasts)) {
      rownames(contrast_data_final) <- contrasts
    }
    info <- gds$metadata$contrast_info %||% list()
    info$data <- contrast_data_final
    gds$metadata$contrast_info <- info
  }

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
  if (inherits(plan, "fmri_group_plan")) {
    design_digest <- digest::digest(
      list(
        observation_ids = plan$design$observation_ids,
        model_matrix = multidesign::model_matrix(plan$design),
        term_data = multidesign::term_data(plan$design),
        grouping_data = multidesign::grouping_data(plan$design)
      ),
      algo = "xxhash64"
    )
    return(digest::digest(
      list(
        schema_version = plan$schema_version,
        method = plan$method,
        reducer_digest = plan$reducer_digest,
        observation_ids = plan$observation_ids,
        feature_ids = plan$feature_ids,
        space_digest = plan$space_digest,
        estimate = plan$estimate,
        variance = plan$variance,
        estimate_source_fingerprint = plan$estimate_source_fingerprint,
        variance_source_fingerprint = plan$variance_source_fingerprint,
        design = design_digest,
        block_size = plan$block_size,
        memory_budget = plan$memory_budget,
        options = plan$options
      ),
      algo = "xxhash64"
    ))
  }
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
  params <- node[!names(node) %in% c("op", "node_id")]
  params <- params[order(names(params))]
  list(op = node$op, params = params)
}

# -------------------------------------------------------------------------
# Helpers ------------------------------------------------------------------

.initial_assays_for_plan <- function(plan, all_assays) {
  if (!"beta" %in% all_assays) {
    return(all_assays)
  }
  nodes <- plan$nodes %||% list()
  if (!length(nodes)) {
    return(all_assays)
  }

  for (node in nodes) {
    if (identical(node$op, "subset_axis")) {
      next
    }
    if (identical(node$op, "reduce") &&
        identical(.normalize_reducer_name(node$method), "ols:voxelwise")) {
      return("beta")
    }
    return(all_assays)
  }

  all_assays
}

.apply_plan_nodes <- function(arrays, plan, space, subjects, col_data = NULL, row_data = NULL, contrast_data = NULL) {
  subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
  current_space <- space
  current_subjects <- subjects
  current_contrasts <- plan$source$probe$contrasts
  current_row_data <- row_data
  current_contrast_data <- contrast_data
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
      res <- .apply_subset_node(
        arrays,
        node,
        current_space,
        current_subjects,
        current_contrasts
      )
      arrays <- res$arrays
      subset_info <- res$subset
      if (!is.null(subset_info$samples) && !is.null(current_row_data)) {
        current_row_data <- current_row_data[subset_info$samples, , drop = FALSE]
      }
      if (!is.null(subset_info$samples)) {
        current_space <- .subset_space(current_space, subset_info$samples)
      }
      if (!is.null(subset_info$subjects)) {
        current_subjects <- current_subjects[subset_info$subjects]
      }
      if (!is.null(subset_info$contrasts)) {
        current_contrasts <- current_contrasts[subset_info$contrasts]
        if (!is.null(current_contrast_data)) {
          current_contrast_data <- current_contrast_data[current_contrasts, , drop = FALSE]
        }
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
      current_row_data <- NULL
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else if (op == "mask_policy") {
      res <- apply_mask_policy(node, arrays, current_space)
      arrays <- res$arrays
      current_space <- res$space
      subset_info <- res$subset
      if (!is.null(subset_info$samples) && !is.null(current_row_data)) {
        current_row_data <- current_row_data[subset_info$samples, , drop = FALSE]
      }
    } else if (op == "map") {
      res <- apply_map_to(node, arrays)
      arrays <- res$arrays
      current_space <- res$space
      current_row_data <- NULL
      subset_info <- list(samples = NULL, subjects = NULL, contrasts = NULL)
    } else if (op == "reduce") {
      res <- apply_reduce(node, arrays, node$weights, current_subjects, col_data, current_contrast_data, current_contrasts)
      arrays <- res$arrays
      current_subjects <- res$subjects %||% "meta"
      current_contrasts <- res$contrasts %||% current_contrasts
      current_contrast_data <- res$contrast_data %||% current_contrast_data
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
        context = list(
          space = current_space,
          subjects = current_subjects,
          contrasts = current_contrasts,
          row_data = current_row_data
        )
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
    contrasts = current_contrasts,
    row_data = current_row_data,
    contrast_data = current_contrast_data,
    writes = writes,
    designs = designs,
    attachments = attachments
  )
}

.apply_subset_node <- function(arrays,
                               node,
                               current_space,
                               current_subjects,
                               current_contrasts) {
  sample <- node$sample
  subject <- node$subject
  contrast <- node$contrast

  dims <- dim(arrays[[1]])
  samples_idx <- if (is.null(sample)) seq_len(dims[1]) else .coerce_index(sample, current_space, dims[1])
  subjects_idx <- if (is.null(subject)) {
    seq_len(dims[2])
  } else {
    .coerce_named_index(subject, current_subjects, axis = "subject")
  }
  contrasts_idx <- if (is.null(contrast)) {
    seq_len(dims[3])
  } else {
    .coerce_named_index(contrast, current_contrasts, axis = "contrast")
  }

  if (any(is.na(subjects_idx))) stop("Unknown subject in subset operation", call. = FALSE)
  if (any(is.na(contrasts_idx))) stop("Unknown contrast in subset operation", call. = FALSE)

  arrays <- lapply(arrays, function(a) a[samples_idx, subjects_idx, contrasts_idx, drop = FALSE])
  list(arrays = arrays, subset = list(samples = samples_idx, subjects = subjects_idx, contrasts = contrasts_idx))
}

.coerce_index <- function(idx, space, n_samples) {
  if (is.numeric(idx)) return(.validate_positive_index(idx, n_samples, "sample"))
  if (is.logical(idx)) {
    if (length(idx) != n_samples) {
      stop("Logical sample index must match the current sample axis.", call. = FALSE)
    }
    return(which(idx))
  }
  if (is.character(idx) && (inherits(space, "space_parcels") || inherits(space, "space_sample_labels"))) {
    out <- match(idx, as.character(space$labels))
    if (anyNA(out)) stop("Unknown sample in subset operation", call. = FALSE)
    return(out)
  }
  stop("Unsupported index type for samples", call. = FALSE)
}

.coerce_named_index <- function(idx, values, axis) {
  if (is.numeric(idx)) return(.validate_positive_index(idx, length(values), axis))
  if (is.logical(idx)) {
    if (length(idx) != length(values)) {
      stop("Logical ", axis, " index must match the current axis.", call. = FALSE)
    }
    return(which(idx))
  }
  if (is.character(idx)) {
    out <- match(idx, values)
    if (anyNA(out)) stop("Unknown ", axis, " in subset operation", call. = FALSE)
    return(out)
  }
  stop("Unsupported index type for ", axis, call. = FALSE)
}

.validate_positive_index <- function(idx, n, axis) {
  if (anyNA(idx) || any(!is.finite(idx)) || any(idx != as.integer(idx)) ||
      any(idx < 1L) || any(idx > n)) {
    stop("Invalid positional ", axis, " index for current axis.", call. = FALSE)
  }
  as.integer(idx)
}

.subset_space <- function(space, idx) {
  if (is.null(idx)) {
    return(space)
  }
  if (exists("space_subset", mode = "function")) {
    return(space_subset(space, idx))
  }
  if (inherits(space, "space_parcels") || inherits(space, "space_sample_labels")) {
    space$labels <- space$labels[idx]
  }
  space
}

## optimize_plan is implemented in R/plan-optimizer.R
