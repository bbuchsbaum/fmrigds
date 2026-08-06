# Conservative compiler for terminal group examinations -------------------

compile_examination_plan <- function(x,
                                     method = NULL,
                                     formula = NULL,
                                     options = list()) {
  plan <- .ensure_plan_node_ids(as_plan(x))
  validate.gds_plan(plan)

  allowed <- c(
    "subset_axis", "derive", "align_to_group", "mask_policy", "map",
    "reduce", "posthoc", "write"
  )
  operations <- vapply(plan$nodes, function(node) node$op, character(1))
  unknown <- setdiff(unique(operations), allowed)
  if (length(unknown)) {
    stop("Unsupported examination plan operation: ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  discarded_writes <- plan$nodes[operations == "write"]
  analysis_nodes <- plan$nodes[operations != "write"]
  analysis_ops <- vapply(analysis_nodes, function(node) node$op, character(1))
  reducer_positions <- which(analysis_ops == "reduce")

  if (length(reducer_positions) > 1L) {
    stop("Group examination v0.1 requires exactly one reducer; found ",
         length(reducer_positions), ".", call. = FALSE)
  }
  if (length(reducer_positions) == 1L &&
      (!is.null(method) || !is.null(formula) || length(options))) {
    stop(
      "Group examination inherits method, formula, and options from the existing reducer; ",
      "model overrides are allowed only when the plan has no reducer.",
      call. = FALSE
    )
  }
  if (!length(reducer_positions)) {
    if (is.null(method)) {
      stop("Group examination requires exactly one reducer or an explicit method.", call. = FALSE)
    }
    reducer <- op_reduce(
      method = method,
      weights = NULL,
      by = "contrast",
      options = options,
      formula = formula
    )
    reducer$node_id <- .make_plan_node_id(plan, reducer, length(plan$nodes) + 1L)
    prefix <- analysis_nodes
    conclusion_tail <- list()
  } else {
    reducer_position <- reducer_positions[[1L]]
    reducer <- analysis_nodes[[reducer_position]]
    prefix <- if (reducer_position > 1L) {
      analysis_nodes[seq_len(reducer_position - 1L)]
    } else {
      list()
    }
    conclusion_tail <- if (reducer_position < length(analysis_nodes)) {
      analysis_nodes[seq.int(reducer_position + 1L, length(analysis_nodes))]
    } else {
      list()
    }
  }

  prefix_ops <- vapply(prefix, function(node) node$op, character(1))
  if ("posthoc" %in% prefix_ops) {
    stop("A post-hoc node cannot appear before the examination reducer.", call. = FALSE)
  }
  tail_ops <- vapply(conclusion_tail, function(node) node$op, character(1))
  invalid_tail <- setdiff(unique(tail_ops), "posthoc")
  if (length(invalid_tail)) {
    stop(
      "Only post-hoc conclusion nodes may appear after the reducer; found ",
      paste(invalid_tail, collapse = ", "), ".",
      call. = FALSE
    )
  }

  axis_resolution <- .resolve_examination_axes(plan, prefix)
  classes <- lapply(prefix, .examination_node_class)
  caps <- .adapter_capabilities(plan$source$adapter)
  strategy <- .choose_examination_scan_strategy(
    classes,
    caps,
    axis_resolution,
    plan$source$probe$dims
  )

  model_spec <- list(
    method = reducer$method,
    formula = reducer$formula %||% NULL,
    options = reducer$options %||% list(),
    weights = reducer$weights %||% NULL,
    by = reducer$by %||% NULL,
    node_id = reducer$node_id
  )
  fingerprint <- .source_fingerprint(plan$source)

  structure(
    list(
      plan = plan,
      prefix = prefix,
      execution_prefix = if (identical(strategy$strategy, "direct")) {
        axis_resolution$execution_prefix
      } else {
        prefix
      },
      reducer = reducer,
      conclusion_tail = conclusion_tail,
      discarded_writes = discarded_writes,
      model_spec = model_spec,
      scan_strategy = strategy$strategy,
      scan_reasons = strategy$reasons,
      prefix_classes = classes,
      adapter_capabilities = caps,
      axis_selection = axis_resolution$selection,
      node_ids = vapply(plan$nodes, function(node) node$node_id, character(1)),
      source_plan_digest = digest_plan(plan),
      origin_source_hash = plan$source$hash,
      source_fingerprint = fingerprint,
      source_fingerprint_digest = fingerprint$digest
    ),
    class = "gds_examination_plan"
  )
}

.examination_node_class <- function(node) {
  switch(
    node$op,
    subset_axis = list(class = "block_local", support_query = TRUE, reason = "axis selection"),
    derive = list(class = "block_local", support_query = TRUE, reason = "elementwise derivation"),
    mask_policy = list(
      class = "materializing",
      support_query = FALSE,
      reason = "data-derived sample indexing requires staged finalization"
    ),
    map = list(
      class = if (is.matrix(node$map) || inherits(node$map, "Matrix")) {
        "support_aware_linear"
      } else {
        "materializing"
      },
      support_query = FALSE,
      reason = "map has no target-to-source support query"
    ),
    align_to_group = list(
      class = "support_aware_linear",
      support_query = FALSE,
      reason = "alignment has no target-to-source support query"
    ),
    list(
      class = "materializing",
      support_query = FALSE,
      reason = paste("unsupported prefix", node$op)
    )
  )
}

.choose_examination_scan_strategy <- function(classes, caps, axes, source_dims) {
  reasons <- character()
  class_names <- vapply(classes, function(x) x$class, character(1))

  if (isTRUE(caps$requires_staging)) {
    reasons <- c(reasons, "adapter requires staging")
  }
  if (isTRUE(axes$unsafe)) {
    reasons <- c(reasons, axes$reason)
  }
  materializing <- which(class_names == "materializing")
  if (length(materializing)) {
    reasons <- c(
      reasons,
      vapply(classes[materializing], function(x) x$reason, character(1))
    )
  }
  support_aware <- which(class_names == "support_aware_linear")
  unsupported_support <- support_aware[
    !vapply(classes[support_aware], function(x) x$support_query, logical(1))
  ]
  if (length(unsupported_support)) {
    reasons <- c(
      reasons,
      vapply(classes[unsupported_support], function(x) x$reason, character(1))
    )
  }
  if (!isTRUE(caps$sample_blocks) || !isTRUE(caps$persistent_handle)) {
    reasons <- c(reasons, "adapter lacks persistent sample-block scanning")
  }
  if (length(axes$selection$subject) != source_dims[["subject"]] &&
      !isTRUE(caps$subject_blocks)) {
    reasons <- c(reasons, "adapter cannot select subject blocks")
  }
  if (length(axes$selection$contrast) != source_dims[["contrast"]] &&
      !isTRUE(caps$contrast_blocks)) {
    reasons <- c(reasons, "adapter cannot select contrast blocks")
  }

  if (length(reasons)) {
    return(list(strategy = "stage", reasons = unique(reasons)))
  }
  if (length(support_aware)) {
    return(list(strategy = "target_block", reasons = character()))
  }
  list(strategy = "direct", reasons = character())
}

.resolve_examination_axes <- function(plan, prefix) {
  dims <- plan$source$probe$dims
  selection <- list(
    sample = seq_len(dims[["sample"]]),
    subject = seq_len(dims[["subject"]]),
    contrast = seq_len(dims[["contrast"]])
  )
  labels <- list(
    sample = .plan_sample_labels(plan),
    subject = as.character(plan$source$probe$subjects),
    contrast = as.character(plan$source$probe$contrasts)
  )
  execution_prefix <- list()
  axis_barrier <- FALSE

  for (node in prefix) {
    if (identical(node$op, "subset_axis")) {
      if (axis_barrier) {
        return(list(
          selection = selection,
          execution_prefix = prefix,
          unsafe = TRUE,
          reason = "subset after an axis/data-dependent prefix node requires staging"
        ))
      }
      for (axis in c("sample", "subject", "contrast")) {
        selector <- node[[axis]] %||% NULL
        if (is.null(selector)) next
        resolved <- .resolve_current_axis_selector(
          selector,
          selection[[axis]],
          labels[[axis]],
          axis
        )
        selection[[axis]] <- resolved$source
        labels[[axis]] <- resolved$labels
      }
    } else {
      execution_prefix <- c(execution_prefix, list(node))
      if (node$op %in% c("mask_policy", "map", "align_to_group")) {
        axis_barrier <- TRUE
      }
    }
  }

  list(
    selection = selection,
    execution_prefix = execution_prefix,
    unsafe = FALSE,
    reason = character()
  )
}

.resolve_current_axis_selector <- function(selector, source, labels, axis) {
  n <- length(source)
  if (is.character(selector)) {
    if (is.null(labels)) {
      stop("Character ", axis, " selection requires labels.", call. = FALSE)
    }
    position <- match(selector, labels)
    if (anyNA(position)) {
      stop("Unknown ", axis, " in examination subset.", call. = FALSE)
    }
  } else if (is.logical(selector)) {
    if (length(selector) != n) {
      stop("Logical ", axis, " subset must match the current axis.", call. = FALSE)
    }
    position <- which(selector)
  } else {
    position <- .validate_positive_index(selector, n, axis)
  }
  list(
    source = as.integer(source[position]),
    labels = if (is.null(labels)) NULL else labels[position]
  )
}

.plan_sample_labels <- function(plan) {
  space <- plan$source$probe$space
  if (inherits(space, "space_parcels") || inherits(space, "space_sample_labels")) {
    return(as.character(space$labels))
  }
  if (!is.null(space$mask_idx)) return(as.character(space$mask_idx))
  as.character(seq_len(plan$source$probe$dims[["sample"]]))
}
