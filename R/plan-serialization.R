# Plan serialization ------------------------------------------------------

#' Save a plan to JSON
#'
#' @param plan Plan or object coercible via [as_plan()]
#' @param file Path to JSON file
#' @name save_plan
#' @export
save_plan <- function(plan, file) {
  plan <- as_plan(plan)
  data <- list(
    source = list(
      adapter = plan$source$adapter,
      source = plan$source$source,
      probe = .serialize_probe(plan$source$probe)
    ),
    meta = .serialize_plan_meta(plan$meta),
    nodes = lapply(plan$nodes, .serialize_node),
    digest = digest_plan(plan)
  )
  jsonlite::write_json(data, file, auto_unbox = TRUE, pretty = TRUE)
  invisible(file)
}

#' Load a plan from JSON
#'
#' @param file JSON file produced by [save_plan]
#' @return A `gds_plan`
#' @export
load_plan <- function(file) {
  data <- jsonlite::read_json(file, simplifyVector = FALSE)
  src <- gds_source(
    data$source$adapter,
    data$source$source,
    probe_result = .deserialize_probe(data$source$probe %||% NULL)
  )
  plan <- gds_plan(src, meta = .deserialize_plan_meta(data$meta %||% list()))
  plan$nodes <- lapply(data$nodes, .deserialize_node)
  if (!is.null(src$probe$dims)) {
    plan$metadata <- list(dims = src$probe$dims)
  }
  plan
}

.serialize_node <- function(node) {
  op <- node$op
  if (op == "subset_axis") {
    out <- list(op = op)
    if (!.is_empty_json_field(node$sample %||% NULL)) out$sample <- node$sample
    if (!.is_empty_json_field(node$subject %||% NULL)) out$subject <- node$subject
    if (!.is_empty_json_field(node$contrast %||% NULL)) out$contrast <- node$contrast
    return(out)
  }
  if (op == "derive") {
    return(node)
  }
  if (op == "map") {
    if (is.matrix(node$map)) {
      node$map <- list(matrix = as.matrix(node$map))
      node$uncertainty <- .serialize_uncertainty(node$uncertainty)
      node$target_space <- .serialize_space(node$target_space)
      return(node)
    }
  }
  if (op == "align_to_group") {
    fam_name <- node$family_name %||% node$family$name %||% NA_character_
    fam_type <- node$family$type %||% node$type %||% NA_character_
    return(list(op = op, family_name = fam_name, type = fam_type))
  }
  if (op == "mask_policy") {
    return(list(op = op, scope = node$policy$scope, rule = node$policy$rule, threshold = node$policy$threshold))
  }
  if (op == "reduce") {
    return(list(
      op = op,
      method = node$method,
      weights = node$weights,
      by = node$by,
      options = node$options %||% list(),
      formula = node$formula %||% NULL
    ))
  }
  if (op == "posthoc") {
    return(list(op = op, method = node$method, options = node$options %||% list()))
  }
  if (op == "write") {
    return(list(op = op, path = node$path, format = node$format, options = node$options %||% list()))
  }
  node
}

.deserialize_node <- function(node) {
  op <- node$op
  if (op == "subset_axis") {
    args <- node[names(node) != "op"]
    args <- lapply(args, function(x) if (.is_empty_json_field(x)) NULL else x)
    return(do.call(op_subset_axis, args))
  }
  if (op == "derive") return(do.call(op_derive, node[names(node) != "op"]))
  if (op == "map" && !is.null(node$map$matrix)) {
    mat <- as.matrix(node$map$matrix)
    uncertainty <- .deserialize_uncertainty(node$uncertainty)
    target_space <- .deserialize_space(node$target_space)
    return(op_map(target_space, mat, uncertainty, node$combine))
  }
  if (op == "reduce") return(op_reduce(node$method, node$weights, node$by, node$options %||% list(), formula = node$formula %||% NULL))
  if (op == "mask_policy") return(op_mask_policy(MaskPolicy(scope = node$scope %||% "group", rule = node$rule %||% "intersection", threshold = node$threshold %||% 0.95)))
  if (op == "align_to_group") return(op_align_to_group(family = NULL, family_name = node$family_name))
  if (op == "posthoc") return(list(op = "posthoc", method = node$method, options = node$options %||% list()))
  if (op == "write") return(op_write(path = node$path, format = node$format, options = node$options %||% list()))
  stop("Cannot deserialize node of type ", op, call. = FALSE)
}

.serialize_uncertainty <- function(x) {
  if (!inherits(x, "gds_uncertainty_rule")) return(x)
  if (!is.null(x$cov_provider) || !is.null(x$kernel)) {
    stop("save_plan() cannot serialize uncertainty rules with function-valued providers", call. = FALSE)
  }
  list(mode = x$mode, df_rule = x$df_rule)
}

.deserialize_uncertainty <- function(x) {
  if (is.null(x)) return(UncertaintyRule("independent"))
  UncertaintyRule(
    mode = x$mode %||% "independent",
    df_rule = x$df_rule %||% "satterthwaite"
  )
}

.serialize_space <- function(x) {
  if (is.null(x) || is.character(x)) return(x)
  if (inherits(x, "space_sample_labels")) {
    return(list(kind = "space_sample_labels", labels = x$labels))
  }
  if (inherits(x, "space_parcels")) {
    return(list(kind = "space_parcels", labels = x$labels, lookup = x$lookup %||% NULL))
  }
  if (inherits(x, "space_voxel")) {
    return(list(
      kind = "space_voxel",
      dim = as.integer(x$dim),
      affine = unclass(x$affine),
      mask_idx = x$mask_idx %||% NULL,
      storage = x$storage %||% "dense",
      template_id = x$template_id %||% NULL
    ))
  }
  stop("save_plan() cannot serialize target_space of class ", paste(class(x), collapse = "/"), call. = FALSE)
}

.deserialize_space <- function(x) {
  if (is.null(x) || is.character(x)) return(x)
  kind <- x$kind %||% NULL
  if (identical(kind, "space_sample_labels")) {
    return(space_sample_labels(as.character(unlist(x$labels, use.names = FALSE))))
  }
  if (identical(kind, "space_parcels")) {
    return(space_parcels(as.character(unlist(x$labels, use.names = FALSE)), lookup = x$lookup %||% NULL))
  }
  if (identical(kind, "space_voxel")) {
    affine <- .json_to_matrix(x$affine, ncol = 4L)
    return(space_voxel(
      dim = as.integer(unlist(x$dim, use.names = FALSE)),
      affine = affine,
      mask_idx = if (is.null(x$mask_idx)) NULL else as.integer(unlist(x$mask_idx, use.names = FALSE)),
      storage = x$storage %||% "dense",
      template_id = x$template_id %||% NULL
    ))
  }
  stop("Cannot deserialize target_space of kind ", kind %||% "<missing>", call. = FALSE)
}

.json_to_matrix <- function(x, ncol) {
  if (is.matrix(x)) return(x)
  if (is.atomic(x)) return(matrix(as.numeric(x), ncol = ncol, byrow = TRUE))
  rows <- lapply(x, function(row) as.numeric(unlist(row, use.names = FALSE)))
  do.call(rbind, rows)
}

.is_empty_json_field <- function(x) {
  is.null(x) || (is.list(x) && !length(x))
}

.serialize_data_frame <- function(x) {
  if (is.null(x)) return(NULL)
  stopifnot(is.data.frame(x))
  cols <- lapply(names(x), function(nm) {
    prepared <- .prepare_col_data_column(x[[nm]])
    list(name = nm, type = prepared$type, data = unname(prepared$data))
  })
  list(
    nrow = nrow(x),
    row_names = rownames(x),
    columns = cols
  )
}

.deserialize_data_frame <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  nrow_x <- as.integer(unlist(x$nrow %||% 0L, use.names = FALSE))
  row_names <- x$row_names %||% NULL
  if (!is.null(row_names)) {
    row_names <- as.character(unlist(row_names, use.names = FALSE))
  }

  cols <- x$columns %||% list()
  if (!length(cols)) {
    out <- data.frame(row.names = if (length(row_names)) row_names else seq_len(nrow_x))
    return(out)
  }

  values <- setNames(vector("list", length(cols)), vapply(cols, function(col) col$name %||% col[["name"]], character(1L)))
  for (col in cols) {
    nm <- col$name %||% col[["name"]]
    type <- col$type %||% col[["type"]] %||% "character"
    data <- col$data %||% col[["data"]] %||% list()
    vec <- switch(type,
      logical = as.logical(unlist(data, use.names = FALSE)),
      integer = as.integer(unlist(data, use.names = FALSE)),
      numeric = as.numeric(unlist(data, use.names = FALSE)),
      character = as.character(unlist(data, use.names = FALSE)),
      as.character(unlist(data, use.names = FALSE))
    )
    values[[nm]] <- vec
  }

  out <- as.data.frame(values, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(row_names)) {
    rownames(out) <- row_names
  }
  out
}

.serialize_map_families <- function(maps) {
  if (is.null(maps) || !length(maps)) return(list())
  lapply(maps, .serialize_map_family_lines)
}

.deserialize_map_families <- function(maps) {
  if (is.null(maps) || !length(maps)) return(list())
  out <- lapply(maps, function(lines) .deserialize_map_family_lines(as.character(unlist(lines, use.names = FALSE))))
  names(out) <- names(maps)
  out
}

.serialize_probe <- function(probe) {
  if (is.null(probe) || !length(probe)) return(NULL)
  list(
    assays = probe$assays %||% NULL,
    dims = unclass(probe$dims %||% NULL),
    subjects = probe$subjects %||% NULL,
    contrasts = probe$contrasts %||% NULL,
    space = .serialize_space(probe$space %||% NULL),
    maps = .serialize_map_families(probe$maps %||% list()),
    metadata = .serialize_metadata(probe$metadata %||% list()),
    columns = probe$columns %||% list(),
    col_data = .serialize_data_frame(probe$col_data %||% NULL),
    row_data = .serialize_data_frame(probe$row_data %||% NULL),
    contrast_data = .serialize_data_frame(probe$contrast_data %||% NULL)
  )
}

.deserialize_probe <- function(probe) {
  if (is.null(probe) || !length(probe)) return(NULL)
  dims <- probe$dims %||% NULL
  dims <- if (is.null(dims)) NULL else {
    vals <- as.integer(unlist(dims, use.names = TRUE))
    if (!length(vals)) {
      NULL
    } else {
    sample_n <- if ("sample" %in% names(vals)) vals[["sample"]] else vals[[1L]]
    subject_n <- if ("subject" %in% names(vals)) vals[["subject"]] else vals[[2L]]
    contrast_n <- if ("contrast" %in% names(vals)) vals[["contrast"]] else vals[[3L]]
    gds_dims(
      sample = sample_n,
      subject = subject_n,
      contrast = contrast_n
    )
    }
  }
  list(
    assays = if (is.null(probe$assays)) NULL else as.character(unlist(probe$assays, use.names = FALSE)),
    dims = dims,
    subjects = if (is.null(probe$subjects)) NULL else as.character(unlist(probe$subjects, use.names = FALSE)),
    contrasts = if (is.null(probe$contrasts)) NULL else as.character(unlist(probe$contrasts, use.names = FALSE)),
    space = .deserialize_space(probe$space %||% NULL),
    maps = .deserialize_map_families(probe$maps %||% list()),
    metadata = probe$metadata %||% list(),
    columns = probe$columns %||% list(),
    col_data = .deserialize_data_frame(probe$col_data %||% NULL),
    row_data = .deserialize_data_frame(probe$row_data %||% NULL),
    contrast_data = .deserialize_data_frame(probe$contrast_data %||% NULL)
  )
}

.serialize_plan_meta <- function(meta) {
  if (is.null(meta)) return(list())
  out <- meta
  if ("col_data" %in% names(out)) out$col_data <- .serialize_data_frame(out$col_data)
  if ("row_data" %in% names(out)) out$row_data <- .serialize_data_frame(out$row_data)
  if ("contrast_data" %in% names(out)) out$contrast_data <- .serialize_data_frame(out$contrast_data)
  out
}

.deserialize_plan_meta <- function(meta) {
  if (is.null(meta)) return(list())
  out <- meta
  if ("col_data" %in% names(out)) out$col_data <- .deserialize_data_frame(out$col_data)
  if ("row_data" %in% names(out)) out$row_data <- .deserialize_data_frame(out$row_data)
  if ("contrast_data" %in% names(out)) out$contrast_data <- .deserialize_data_frame(out$contrast_data)
  out
}
