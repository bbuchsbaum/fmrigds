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
    source = list(adapter = plan$source$adapter, source = plan$source$source),
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
  src <- gds_source(data$source$adapter, data$source$source)
  plan <- gds_plan(src)
  plan$nodes <- lapply(data$nodes, .deserialize_node)
  plan
}

.serialize_node <- function(node) {
  op <- node$op
  if (op %in% c("subset_axis", "derive")) {
    return(node)
  }
  if (op == "map") {
    if (is.matrix(node$map)) {
      node$map <- list(matrix = as.matrix(node$map))
      node$uncertainty <- list(mode = node$uncertainty$mode)
      node$target_space <- NULL
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
    return(list(op = op, method = node$method, weights = node$weights, by = node$by))
  }
  if (op == "posthoc") {
    return(list(op = op, method = node$method, options = node$options %||% list()))
  }
  node
}

.deserialize_node <- function(node) {
  op <- node$op
  if (op == "subset_axis") return(do.call(op_subset_axis, node[names(node) != "op"]))
  if (op == "derive") return(do.call(op_derive, node[names(node) != "op"]))
  if (op == "map" && !is.null(node$map$matrix)) {
    mat <- as.matrix(node$map$matrix)
    uncertainty <- UncertaintyRule(node$uncertainty$mode %||% "independent")
    return(op_map(NULL, mat, uncertainty, node$combine))
  }
  if (op == "reduce") return(op_reduce(node$method, node$weights, node$by, list()))
  if (op == "mask_policy") return(op_mask_policy(MaskPolicy(scope = node$scope %||% "group", rule = node$rule %||% "intersection", threshold = node$threshold %||% 0.95)))
  if (op == "align_to_group") return(op_align_to_group(family = NULL, family_name = node$family_name))
  if (op == "posthoc") return(list(op = "posthoc", method = node$method, options = node$options %||% list()))
  stop("Cannot deserialize node of type ", op, call. = FALSE)
}
