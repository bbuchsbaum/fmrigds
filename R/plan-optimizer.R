.merge_subset <- function(a, b) {
  merge_field <- function(old, new) {
    if (is.null(old)) return(new)
    if (is.null(new)) return(old)
    intersect(old, new)
  }
  list(
    op = "subset_axis",
    sample = merge_field(a$sample, b$sample),
    subject = merge_field(a$subject, b$subject),
    contrast = merge_field(a$contrast, b$contrast)
  )
}

.combine_subsets <- function(nodes) {
  combined <- NULL
  out <- list()
  for (node in nodes) {
    if (node$op == "subset_axis") {
      combined <- if (is.null(combined)) node else .merge_subset(combined, node)
    } else {
      if (!is.null(combined)) {
        out <- c(out, list(combined))
        combined <- NULL
      }
      out <- c(out, list(node))
    }
  }
  if (!is.null(combined)) out <- c(list(combined), out)
  out
}

.coalesce_derives <- function(nodes) {
  out <- list()
  pending <- NULL
  for (node in nodes) {
    if (node$op == "derive") {
      if (is.null(pending)) {
        pending <- node
      } else {
        pending$what <- unique(c(pending$what, node$what))
        pending$options <- utils::modifyList(pending$options, node$options)
      }
    } else {
      if (!is.null(pending)) {
        out <- c(out, list(pending))
        pending <- NULL
      }
      out <- c(out, list(node))
    }
  }
  if (!is.null(pending)) out <- c(out, list(pending))
  out
}

.reorder_nodes <- function(nodes) {
  priority <- c(
    subset_axis = 1,
    align_to_group = 2,
    mask_policy = 3,
    derive = 4,
    map = 5,
    reduce = 6,
    posthoc = 7,
    write = 8
  )
  order_idx <- order(vapply(nodes, function(n) priority[[n$op]] %||% 50, numeric(1)), seq_along(nodes))
  nodes[order_idx]
}

.fuse_masks <- function(nodes) {
  out <- list()
  for (node in nodes) {
    if (node$op == "mask_policy" && length(out) && out[[length(out)]]$op == "mask_policy") {
      existing <- out[[length(out)]]$policy
      if (!is.list(existing) || inherits(existing[[1]], "gds_mask_policy")) existing <- list(existing)
      out[[length(out)]]$policy <- c(existing, list(node$policy))
    } else {
      out <- c(out, list(node))
    }
  }
  out
}

#' Plan optimizer applying rewrite rules
#' @keywords internal
optimize_plan <- function(plan) {
  nodes <- plan$nodes
  if (!length(nodes)) return(plan)
  nodes <- .combine_subsets(nodes)
  nodes <- .coalesce_derives(nodes)
  nodes <- .fuse_masks(nodes)
  nodes <- .reorder_nodes(nodes)
  plan$nodes <- nodes
  plan
}
