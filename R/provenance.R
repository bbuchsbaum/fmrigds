# Portable provenance ------------------------------------------------------

.provenance_schema_version <- "1.0.0"
.provenance_digest_schema_version <- "1.0"

.next_portable_occurrence <- function(kind) {
  counter <- as.integer(.gds_state$portable_occurrence_counter %||% 0L) + 1L
  .gds_state$portable_occurrence_counter <- counter
  list(
    kind = as.character(kind),
    timestamp = .iso_utc(Sys.time()),
    process = as.integer(Sys.getpid()),
    counter = counter
  )
}

.iso_utc <- function(x = Sys.time()) {
  format(as.POSIXct(x, tz = "UTC"), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

.canonical_portable_value <- function(x, path = "value") {
  if (is.null(x)) return(NULL)
  if (inherits(x, "POSIXt")) return(.iso_utc(x))
  if (inherits(x, "numeric_version")) return(as.character(x))
  if (inherits(x, "formula")) return(paste(deparse(x), collapse = " "))
  if (is.factor(x)) x <- as.character(x)
  if (is.environment(x) || is.function(x) || is.symbol(x) || is.language(x) || isS4(x)) {
    stop("Non-portable value at `", path, "` (", paste(class(x), collapse = "/"), ").", call. = FALSE)
  }
  if (is.data.frame(x) || is.matrix(x) || (is.array(x) && length(dim(x)) > 0L)) {
    stop("Matrices, arrays, and data frames are not portable at `", path, "`.", call. = FALSE)
  }
  if (is.raw(x)) return(as.integer(x))
  if (is.atomic(x)) {
    if (anyNA(x) || (is.numeric(x) && any(!is.finite(x)))) {
      stop("Missing or non-finite scalar at `", path, "` is not portable.", call. = FALSE)
    }
    if (!is.null(names(x)) && length(x)) {
      if (any(!nzchar(names(x))) || anyDuplicated(names(x))) {
        stop("Named values at `", path, "` need unique non-empty names.", call. = FALSE)
      }
      idx <- order(enc2utf8(names(x)), method = "radix")
      out <- lapply(idx, function(i) .canonical_portable_value(x[[i]], paste0(path, ".", names(x)[i])))
      names(out) <- names(x)[idx]
      return(out)
    }
    return(unname(x))
  }
  if (is.list(x)) {
    nm <- names(x)
    if (!is.null(nm)) {
      if (any(!nzchar(nm)) || anyDuplicated(nm)) {
        stop("Object at `", path, "` needs unique non-empty names.", call. = FALSE)
      }
      idx <- order(enc2utf8(nm), method = "radix")
      out <- lapply(idx, function(i) .canonical_portable_value(x[[i]], paste0(path, ".", nm[i])))
      names(out) <- nm[idx]
      return(out)
    }
    return(lapply(seq_along(x), function(i) .canonical_portable_value(x[[i]], paste0(path, "[", i, "]"))))
  }
  stop("Unsupported portable value at `", path, "`.", call. = FALSE)
}

.canonical_portable_json <- function(x) {
  canonical <- .canonical_portable_value(x)
  enc2utf8(as.character(jsonlite::toJSON(
    canonical,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    force = TRUE
  )))
}

.portable_sha256 <- function(x) {
  bytes <- charToRaw(.canonical_portable_json(x))
  list(
    algorithm = "sha256",
    value = tolower(digest::digest(bytes, algo = "sha256", serialize = FALSE)),
    schemaVersion = .provenance_digest_schema_version
  )
}

.operation_semantic_id <- function(op_name) {
  slug <- tolower(gsub("[^a-z0-9]+", "-", as.character(op_name)))
  slug <- gsub("(^-+|-+$)", "", slug)
  if (!nzchar(slug)) slug <- "operation"
  paste0("org.bbuchsbaum.fmrigds.operation/", slug)
}

.portable_operation_params <- function(params) {
  if (is.null(params)) return(list())
  if (!is.list(params) || (length(params) && (is.null(names(params)) || any(!nzchar(names(params)))))) {
    stop("Provenance `params` must be a named list.", call. = FALSE)
  }
  keep <- !vapply(params, function(x) {
    is.null(x) || (length(x) == 1L && is.atomic(x) && (is.na(x) || (is.numeric(x) && !is.finite(x))))
  }, logical(1L))
  params <- params[keep]
  .canonical_portable_value(params, "params")
}

.activity_descriptor <- function(op_name, params, inputs, software) {
  list(
    schemaVersion = .provenance_schema_version,
    semanticId = .operation_semantic_id(op_name),
    params = .portable_operation_params(params),
    inputs = as.character(inputs),
    software = .canonical_portable_value(software, "software")
  )
}

.new_activity_record <- function(op_name,
                                 params,
                                 inputs,
                                 id,
                                 timestamp = Sys.time(),
                                 software = list(
                                   package = "fmrigds",
                                   version = .pkg_version(),
                                   RVersion = as.character(getRversion())
                                 )) {
  descriptor <- .activity_descriptor(op_name, params, inputs, software)
  digest <- .portable_sha256(descriptor)
  structure(
    list(
      id = as.character(id),
      semanticId = descriptor$semanticId,
      op = as.character(op_name),
      params = descriptor$params,
      inputs = descriptor$inputs,
      timestamp = .iso_utc(timestamp),
      software = descriptor$software,
      digest = digest,
      # Compatibility view for callers that inspected the former hash field.
      hash = paste0("sha256:", digest$value)
    ),
    class = "gds_provenance_node"
  )
}

.empty_provenance <- function() {
  list(
    schemaVersion = .provenance_schema_version,
    entities = list(),
    graph = list(),
    log = character(),
    digest = NULL,
    graphReceipt = list(status = "unavailable", reason = "no provenance records"),
    status = "unavailable",
    availability = list(graph = "unavailable", sourceIdentity = "unavailable"),
    roots = character(),
    heads = character()
  )
}

.is_portable_activity <- function(node) {
  is.list(node) && !is.null(node$id) && !is.null(node$semanticId) &&
    !is.null(node$inputs) && !is.null(node$digest$value)
}

.provenance_graph_problems <- function(provenance, check_declared = TRUE) {
  entities <- provenance$entities %||% list()
  activities <- provenance$graph %||% list()
  problems <- character()
  entity_ids <- if (length(entities)) vapply(entities, function(x) x$id %||% "", character(1L)) else character()
  activity_ids <- if (length(activities)) vapply(activities, function(x) x$id %||% "", character(1L)) else character()
  all_ids <- c(entity_ids, activity_ids)

  if (any(!nzchar(all_ids))) problems <- c(problems, "every entity and activity must have a non-empty id")
  if (anyDuplicated(all_ids)) problems <- c(problems, "entity and activity ids must be globally unique")
  if (length(activities) && !all(vapply(activities, .is_portable_activity, logical(1L)))) {
    return(unique(c(problems, "legacy activities without portable ids/digests cannot form a complete graph")))
  }

  for (node in activities) {
    inputs <- as.character(node$inputs %||% character())
    if (anyDuplicated(inputs)) {
      problems <- c(problems, paste0("activity '", node$id, "' repeats an input"))
    }
    if (node$id %in% inputs) {
      problems <- c(problems, paste0("activity '", node$id, "' has a self-edge"))
    }
    missing <- setdiff(inputs, all_ids)
    if (length(missing)) {
      problems <- c(problems, paste0(
        "activity '", node$id, "' has missing input endpoint(s): ",
        paste(missing, collapse = ", ")
      ))
    }
  }

  # Kahn's algorithm over activity-to-activity edges.
  if (length(activity_ids) && !anyDuplicated(activity_ids)) {
    indegree <- setNames(integer(length(activity_ids)), activity_ids)
    children <- setNames(vector("list", length(activity_ids)), activity_ids)
    for (node in activities) {
      parents <- intersect(as.character(node$inputs %||% character()), activity_ids)
      indegree[[node$id]] <- length(parents)
      for (parent in parents) children[[parent]] <- c(children[[parent]], node$id)
    }
    queue <- names(indegree)[indegree == 0L]
    visited <- 0L
    while (length(queue)) {
      id <- queue[[1L]]
      queue <- queue[-1L]
      visited <- visited + 1L
      for (child in children[[id]]) {
        indegree[[child]] <- indegree[[child]] - 1L
        if (indegree[[child]] == 0L) queue <- c(queue, child)
      }
    }
    if (visited != length(activity_ids)) problems <- c(problems, "activity graph contains a cycle")

    if (isTRUE(check_declared)) {
      expected_roots <- sort(names(indegree)[vapply(activities, function(node) {
        !any(as.character(node$inputs %||% character()) %in% activity_ids)
      }, logical(1L))])
      consumed <- unique(unlist(lapply(activities, function(node) {
        intersect(as.character(node$inputs %||% character()), activity_ids)
      }), use.names = FALSE))
      expected_heads <- sort(setdiff(activity_ids, consumed))
      if (!identical(sort(as.character(provenance$roots %||% character())), expected_roots)) {
        problems <- c(problems, "declared provenance roots do not match graph degree")
      }
      if (!identical(sort(as.character(provenance$heads %||% character())), expected_heads)) {
        problems <- c(problems, "declared provenance heads do not match graph degree")
      }
    }
  }
  unique(problems)
}

#' Validate an fmrigds provenance graph
#'
#' Checks globally unique source/activity IDs, resolvable input endpoints,
#' duplicate/self edges, acyclicity, and declared root/head degree receipts.
#' Legacy list-only provenance is readable but is reported as incomplete.
#'
#' @param x A provenance list, a GDS metadata list, or a realised `gds`.
#' @param error Whether to stop on the first combined validation report.
#' @return Invisibly `TRUE` when valid; otherwise the character problems when
#'   `error = FALSE`.
#' @export
validate_provenance_graph <- function(x, error = TRUE) {
  raw_provenance <- if (inherits(x, "gds")) {
    x$metadata$provenance
  } else if (is.list(x) && !is.null(x$provenance)) {
    x$provenance
  } else {
    x
  }
  declared_roots <- raw_provenance$roots %||% NULL
  declared_heads <- raw_provenance$heads %||% NULL
  provenance <- .normalize_provenance(raw_provenance)
  # Normalization refreshes derived receipts for ordinary consumers. A
  # validator must instead check declarations as supplied, not silently repair
  # a corrupt roots/heads receipt before comparing graph degree.
  if (!is.null(declared_roots)) provenance$roots <- as.character(declared_roots)
  if (!is.null(declared_heads)) provenance$heads <- as.character(declared_heads)
  problems <- .provenance_graph_problems(provenance, check_declared = TRUE)
  if (length(problems)) {
    if (isTRUE(error)) stop(paste(problems, collapse = "; "), call. = FALSE)
    return(problems)
  }
  invisible(TRUE)
}

.graph_receipt_descriptor <- function(provenance) {
  entities <- .portable_source_entities(provenance$entities %||% list())
  if (length(entities)) {
    entities <- entities[order(vapply(entities, `[[`, character(1L), "id"), method = "radix")]
  }
  activities <- lapply(provenance$graph %||% list(), function(node) {
    list(
      id = node$id,
      semanticId = node$semanticId,
      params = node$params,
      inputs = as.character(node$inputs),
      software = node$software,
      digest = node$digest
    )
  })
  if (length(activities)) {
    activities <- activities[order(vapply(activities, `[[`, character(1L), "id"), method = "radix")]
  }
  edges <- unlist(lapply(activities, function(node) {
    lapply(node$inputs, function(input) list(from = input, to = node$id))
  }), recursive = FALSE)
  if (length(edges)) {
    keys <- vapply(edges, function(edge) paste(edge$from, edge$to, sep = "\r"), character(1L))
    edges <- edges[order(keys, method = "radix")]
  }
  list(
    schemaVersion = .provenance_schema_version,
    entities = entities,
    activities = activities,
    edges = edges
  )
}

.refresh_provenance <- function(provenance) {
  provenance$schemaVersion <- provenance$schemaVersion %||% .provenance_schema_version
  provenance$entities <- provenance$entities %||% list()
  provenance$graph <- provenance$graph %||% list()
  provenance$log <- as.character(provenance$log %||% character())

  if (length(provenance$graph) && !all(vapply(provenance$graph, .is_portable_activity, logical(1L)))) {
    provenance$status <- "legacy-incomplete"
    provenance$availability <- list(graph = "legacy-incomplete", sourceIdentity = "unavailable")
    provenance$roots <- character()
    provenance$heads <- character()
    provenance$graphReceipt <- list(
      status = "unavailable",
      reason = "legacy provenance has no explicit activity ids/edges"
    )
    return(provenance)
  }

  activity_ids <- if (length(provenance$graph)) {
    vapply(provenance$graph, `[[`, character(1L), "id")
  } else character()
  provenance$roots <- if (length(activity_ids)) {
    activity_ids[vapply(provenance$graph, function(node) {
      !any(as.character(node$inputs %||% character()) %in% activity_ids)
    }, logical(1L))]
  } else character()
  consumed <- unique(unlist(lapply(provenance$graph, function(node) {
    intersect(as.character(node$inputs %||% character()), activity_ids)
  }), use.names = FALSE))
  provenance$heads <- setdiff(activity_ids, consumed)

  problems <- .provenance_graph_problems(provenance, check_declared = TRUE)
  if (length(problems)) {
    provenance$status <- "invalid"
    provenance$availability <- list(graph = "invalid", sourceIdentity = "unknown")
    provenance$graphReceipt <- list(status = "unavailable", reason = paste(problems, collapse = "; "))
    return(provenance)
  }

  statuses <- if (length(provenance$entities)) {
    vapply(provenance$entities, function(x) x$identityStatus %||% "unavailable", character(1L))
  } else character()
  source_status <- if (!length(statuses)) {
    "unavailable"
  } else if (all(statuses == "verified")) {
    "verified"
  } else if (any(statuses == "changed-during-read")) {
    "changed-during-read"
  } else {
    "incomplete"
  }
  graph_status <- if (length(provenance$graph) || length(provenance$entities)) "available" else "unavailable"
  provenance$status <- if (identical(graph_status, "unavailable")) {
    "unavailable"
  } else if (identical(source_status, "verified")) {
    "complete"
  } else {
    "incomplete"
  }
  provenance$availability <- list(graph = graph_status, sourceIdentity = source_status)
  if (identical(graph_status, "available")) {
    provenance$graphReceipt <- c(list(status = "verified"), .portable_sha256(.graph_receipt_descriptor(provenance)))
  } else {
    provenance$graphReceipt <- list(status = "unavailable", reason = "no provenance records")
  }
  provenance
}

.normalize_provenance <- function(provenance = NULL) {
  if (is.null(provenance) || !length(provenance)) return(.empty_provenance())
  if (!is.list(provenance)) stop("Provenance metadata must be a list.", call. = FALSE)
  out <- .empty_provenance()
  for (nm in names(provenance)) out[[nm]] <- provenance[[nm]]
  # Old metadata used only graph/log/digest. Keep it readable, but never infer
  # ids or edges from list order.
  .refresh_provenance(out)
}

.next_activity_id <- function(provenance, descriptor) {
  ordinal <- length(provenance$graph %||% list()) + 1L
  seed <- .portable_sha256(list(
    ordinal = ordinal,
    descriptor = descriptor,
    occurrence = .next_portable_occurrence("activity")
  ))$value
  sprintf("activity-%04d-%s", ordinal, substr(seed, 1L, 12L))
}

.provenance_current_inputs <- function(metadata) {
  provenance <- .normalize_provenance(metadata$provenance %||% NULL)
  inputs <- as.character(provenance$heads %||% character())
  if (length(provenance$entities)) {
    entity_ids <- vapply(provenance$entities, `[[`, character(1L), "id")
    consumed_entities <- unique(unlist(lapply(provenance$graph, function(node) {
      intersect(as.character(node$inputs %||% character()), entity_ids)
    }), use.names = FALSE))
    inputs <- c(inputs, setdiff(entity_ids, consumed_entities))
  }
  unique(inputs)
}

.append_source_entities <- function(metadata, entities) {
  provenance <- .normalize_provenance(metadata$provenance %||% NULL)
  if (!length(entities)) {
    metadata$provenance <- provenance
    return(metadata)
  }
  existing <- provenance$entities %||% list()
  existing_ids <- if (length(existing)) vapply(existing, `[[`, character(1L), "id") else character()
  for (entity in entities) {
    id <- entity$id
    if (id %in% existing_ids) {
      old <- existing[[match(id, existing_ids)]]
      if (!identical(.canonical_portable_json(.portable_source_entity(old)),
                     .canonical_portable_json(.portable_source_entity(entity)))) {
        suffix <- substr(.portable_sha256(.portable_source_entity(entity))$value, 1L, 8L)
        entity$id <- paste0(substr(id, 1L, 110L), "-", suffix)
      } else {
        next
      }
    }
    existing <- c(existing, list(entity))
    existing_ids <- c(existing_ids, entity$id)
  }
  provenance$entities <- existing
  metadata$provenance <- .refresh_provenance(provenance)
  metadata
}

.ensure_nonfile_source <- function(metadata, kind = "memory", role = "object") {
  provenance <- .normalize_provenance(metadata$provenance %||% NULL)
  if (!length(provenance$entities) && !length(provenance$graph)) {
    metadata$provenance <- provenance
    metadata <- .append_source_entities(metadata, list(.source_nonfile_entity(kind, role)))
  } else {
    metadata$provenance <- provenance
  }
  metadata
}

.portable_activity_equal <- function(x, y) {
  descriptor <- function(node) list(
    id = node$id,
    semanticId = node$semanticId,
    params = node$params,
    inputs = node$inputs,
    software = node$software,
    digest = node$digest
  )
  identical(.canonical_portable_json(descriptor(x)), .canonical_portable_json(descriptor(y)))
}

.rewrite_provenance_ids <- function(provenance, mapping) {
  if (!length(mapping)) return(provenance)
  for (i in seq_along(provenance$entities)) {
    old <- provenance$entities[[i]]$id
    if (old %in% names(mapping)) provenance$entities[[i]]$id <- unname(mapping[[old]])
  }
  for (i in seq_along(provenance$graph)) {
    node <- provenance$graph[[i]]
    if (node$id %in% names(mapping)) node$id <- unname(mapping[[node$id]])
    node$inputs <- vapply(as.character(node$inputs), function(id) mapping[[id]] %||% id, character(1L))
    descriptor <- .activity_descriptor(node$op, node$params, node$inputs, node$software)
    node$digest <- .portable_sha256(descriptor)
    node$hash <- paste0("sha256:", node$digest$value)
    provenance$graph[[i]] <- node
  }
  .refresh_provenance(provenance)
}

.prepare_merge_lineage <- function(metadata, prefix) {
  provenance <- .normalize_provenance(metadata$provenance %||% NULL)
  if (identical(provenance$status, "legacy-incomplete")) {
    legacy <- provenance
    provenance <- .empty_provenance()
    entity <- .source_nonfile_entity(
      "legacy-lineage",
      role = paste0(prefix, "-lineage"),
      reason = "legacy lineage has no addressable endpoints"
    )
    provenance$entities <- list(entity)
    provenance$legacy <- legacy
    provenance <- .refresh_provenance(provenance)
  } else if (!length(provenance$entities) && !length(provenance$graph)) {
    provenance$entities <- list(.source_nonfile_entity(
      "memory",
      role = paste0(prefix, "-lineage"),
      reason = "in-memory lineage has no byte receipt"
    ))
    provenance <- .refresh_provenance(provenance)
  }
  provenance
}

.merge_provenance_lineages <- function(metadata_a, metadata_b) {
  a <- .prepare_merge_lineage(metadata_a, "a")
  b <- .prepare_merge_lineage(metadata_b, "b")
  a_ids <- c(
    vapply(a$entities %||% list(), `[[`, character(1L), "id"),
    vapply(a$graph %||% list(), `[[`, character(1L), "id")
  )
  b_records <- c(b$entities %||% list(), b$graph %||% list())
  b_ids <- if (length(b_records)) vapply(b_records, `[[`, character(1L), "id") else character()
  conflicts <- intersect(a_ids, b_ids)
  mapping <- character()
  for (id in conflicts) {
    a_record <- c(a$entities, a$graph)[[match(id, a_ids)]]
    b_record <- b_records[[match(id, b_ids)]]
    same <- if (.is_portable_activity(a_record) && .is_portable_activity(b_record)) {
      .portable_activity_equal(a_record, b_record)
    } else if (!.is_portable_activity(a_record) && !.is_portable_activity(b_record)) {
      identical(
        .canonical_portable_json(.portable_source_entity(a_record)),
        .canonical_portable_json(.portable_source_entity(b_record))
      )
    } else FALSE
    if (!same) {
      candidate <- paste0("b-", substr(id, 1L, 120L))
      while (candidate %in% c(a_ids, b_ids, unname(mapping))) candidate <- paste0(candidate, "x")
      mapping[id] <- candidate
    }
  }
  b <- .rewrite_provenance_ids(b, mapping)

  out <- .empty_provenance()
  out$entities <- a$entities
  out$graph <- a$graph
  out$log <- a$log
  out$digest <- NULL
  out <- .append_source_entities(list(provenance = out), b$entities)$provenance
  existing_activity_ids <- if (length(out$graph)) vapply(out$graph, `[[`, character(1L), "id") else character()
  for (node in b$graph) {
    if (node$id %in% existing_activity_ids) {
      old <- out$graph[[match(node$id, existing_activity_ids)]]
      if (.portable_activity_equal(old, node)) next
      stop("Cannot merge provenance activities with conflicting id '", node$id, "'.", call. = FALSE)
    }
    out$graph <- c(out$graph, list(node))
    out$log <- c(out$log, b$log[match(node$id, vapply(b$graph, `[[`, character(1L), "id"))] %||% character())
    existing_activity_ids <- c(existing_activity_ids, node$id)
  }
  out <- .refresh_provenance(out)
  list(
    provenance = out,
    parents_a = .provenance_current_inputs(list(provenance = a)),
    parents_b = .provenance_current_inputs(list(provenance = b))
  )
}

.merge_gds_metadata <- function(base, extra) {
  if (is.null(extra)) return(base)
  if (!is.list(extra)) stop("GDS metadata must be a list.", call. = FALSE)
  base_provenance <- .normalize_provenance(base$provenance %||% NULL)
  extra_provenance <- extra$provenance %||% NULL
  extra$provenance <- NULL
  out <- utils::modifyList(base, extra)
  if (is.null(extra_provenance)) {
    out$provenance <- base_provenance
  } else if (!length(base_provenance$entities) && !length(base_provenance$graph)) {
    out$provenance <- .normalize_provenance(extra_provenance)
  } else {
    merged <- .merge_provenance_lineages(
      list(provenance = base_provenance),
      list(provenance = extra_provenance)
    )
    out$provenance <- merged$provenance
  }
  out
}

.portable_scalar_options <- function(options) {
  if (is.null(options) || !length(options)) return(list())
  out <- list()
  for (nm in names(options)) {
    value <- options[[nm]]
    if (is.atomic(value) && length(value) == 1L && !is.na(value) &&
        (!is.numeric(value) || is.finite(value))) {
      out[[nm]] <- unname(value)
    }
  }
  out
}

.portable_design_receipt <- function(X, columns = colnames(X)) {
  if (is.null(X)) return(list(status = "not-applicable", columns = character()))
  X <- as.matrix(X)
  columns <- as.character(columns %||% paste0("X", seq_len(ncol(X))))
  digest <- .portable_sha256(list(
    layout = "row-major",
    dimensions = as.integer(dim(X)),
    columns = columns,
    values = as.numeric(t(unname(X)))
  ))
  list(status = "verified", digest = digest, columns = columns)
}

.portable_plan_node <- function(node) {
  op <- node$op
  id <- node$node_id %||% NULL
  params <- switch(op,
    subset_axis = list(
      sampleSelection = if (is.null(node$sample)) NULL else .portable_sha256(list(node$sample)),
      subjectSelection = if (is.null(node$subject)) NULL else .portable_sha256(list(node$subject)),
      contrastSelection = if (is.null(node$contrast)) NULL else .portable_sha256(list(node$contrast))
    ),
    derive = list(what = as.character(node$what), options = .portable_scalar_options(node$options)),
    align_to_group = list(family = node$family_name %||% node$family$name %||% NULL),
    mask_policy = list(
      scope = node$policy$scope %||% NULL,
      rule = node$policy$rule %||% NULL,
      threshold = node$policy$threshold %||% NULL
    ),
    map = {
      if (is.matrix(node$map)) {
        list(
          mapDigest = .portable_sha256(list(
            dimensions = as.integer(dim(node$map)),
            values = as.numeric(t(unname(node$map)))
          )),
          combine = node$combine %||% NULL,
          uncertainty = node$uncertainty$mode %||% NULL
        )
      } else {
        stop("function-valued map nodes have no portable plan descriptor", call. = FALSE)
      }
    },
    reduce = list(
      method = .normalize_reducer_name(node$method),
      weights = node$weights,
      by = node$by,
      formula = node$formula %||% NULL,
      options = .portable_scalar_options(node$options)
    ),
    posthoc = list(method = node$method, options = .portable_scalar_options(node$options)),
    write = list(format = node$format, options = .portable_scalar_options(node$options)),
    stop("operation '", op, "' has no portable plan descriptor", call. = FALSE)
  )
  list(id = id, op = op, params = params)
}

.portable_plan_receipt <- function(plan) {
  tryCatch({
    descriptor <- list(
      schemaVersion = "1.0",
      adapter = plan$source$adapter,
      sourceEntities = .portable_source_entities(plan$source$probe$metadata$provenance$entities %||% list()),
      nodes = lapply(plan$nodes %||% list(), .portable_plan_node)
    )
    c(list(status = "verified"), .portable_sha256(descriptor))
  }, error = function(e) {
    list(status = "unavailable", reason = conditionMessage(e))
  })
}

.custom_weight_receipt <- function(node) {
  if (!identical(node$weights, "custom")) {
    return(list(mode = node$weights %||% "model-specific", status = "verified"))
  }
  values <- node$options$custom_weights %||% NULL
  if (is.null(values)) return(list(mode = "custom", status = "unavailable", reason = "custom weights missing"))
  numeric_values <- as.numeric(values)
  finite <- numeric_values[is.finite(numeric_values)]
  portable_values <- lapply(numeric_values, function(value) {
    if (is.na(value)) return(list(status = "missing"))
    if (is.infinite(value)) {
      return(list(status = if (value > 0) "positive-infinity" else "negative-infinity"))
    }
    list(status = "finite", value = unname(value))
  })
  list(
    mode = "custom",
    status = "verified",
    digest = .portable_sha256(list(
      dimensions = as.integer(dim(values) %||% length(values)),
      values = portable_values
    )),
    summary = list(
      length = length(numeric_values),
      finite = length(finite),
      min = if (length(finite)) min(finite) else NULL,
      max = if (length(finite)) max(finite) else NULL,
      mean = if (length(finite)) mean(finite) else NULL
    )
  )
}

.build_reduction_receipt <- function(node,
                                     result,
                                     input_subjects,
                                     input_contrasts,
                                     output_contrasts,
                                     output_assays,
                                     plan) {
  reducer_id <- .normalize_reducer_name(node$method)
  reducer <- get_reducer(reducer_id)
  design <- result$design_info$portable %||% list(status = "not-applicable", columns = character())
  terms <- as.character(design$columns %||% character())
  list(
    schemaVersion = "1.0.0",
    status = if (identical(.portable_plan_receipt(plan)$status, "verified")) "portable" else "incomplete",
    reducerId = reducer_id,
    planNodeId = node$node_id,
    formula = node$formula %||% NULL,
    options = .portable_scalar_options(node$options %||% list()),
    weight = .custom_weight_receipt(node),
    modelContract = reducer$model_contract %||% list(status = "unavailable"),
    reducerVersion = as.character(.pkg_version()),
    design = design,
    nInputSubjects = as.integer(length(input_subjects)),
    inputContrasts = as.character(input_contrasts),
    outputContrasts = as.character(output_contrasts),
    terms = terms,
    outputAssays = as.character(output_assays),
    planDigest = .portable_plan_receipt(plan),
    software = list(
      package = "fmrigds",
      version = as.character(.pkg_version()),
      RVersion = as.character(getRversion())
    )
  )
}
