# Coverage tests for verb-align, verb-map-to, verb-mask, verb-reduce, verb-write,
# align-exec, map-exec, mask-exec

# ===========================================================================
# Helpers
# ===========================================================================
.plan_gds <- function(ns = 3, nsubj = 2, nc = 2) {
  beta <- array(seq_len(ns * nsubj * nc), dim = c(ns, nsubj, nc))
  var  <- array(0.1, dim = c(ns, nsubj, nc))
  new_gds(
    assays    = list(beta = beta, var = var),
    space     = space_sample_labels(paste0("s", seq_len(ns))),
    subjects  = paste0("sub", seq_len(nsubj)),
    contrasts = paste0("c", seq_len(nc))
  ) |> as_plan()
}

# ===========================================================================
# verb-reduce.R: reduce()
# ===========================================================================

test_that("reduce() adds a reduce node to the plan", {
  p <- .plan_gds()
  p2 <- reduce(p, "fixed")
  expect_s3_class(p2, "gds_plan")
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last_node$op, "reduce")
})

test_that("reduce() accepts method and weights", {
  p <- .plan_gds()
  p2 <- reduce(p, method = "random", weights = "1/var")
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last_node$method, "random")
  expect_equal(last_node$weights, "1/var")
})

test_that("reduce_eager() computes immediately", {
  p <- .plan_gds()
  result <- reduce_eager(p, "fixed")
  expect_s3_class(result, "gds")
  expect_equal(dim(result$assays$beta)[2], 1L)
})

# ===========================================================================
# verb-align.R: align()
# ===========================================================================

test_that("align() adds an align node to the plan", {
  p <- .plan_gds()
  from_sp <- space_sample_labels(c("s1", "s2", "s3"))
  to_sp <- space_sample_labels(c("g1", "g2", "g3"))
  M <- diag(3)
  fam <- MapFamily("test", from_sp, to_sp, by_subject = list(sub1 = M, sub2 = M))
  p2 <- align(p, fam)
  expect_s3_class(p2, "gds_plan")
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last_node$op, "align_to_group")
})

test_that("align_eager() computes immediately", {
  p <- .plan_gds()
  from_sp <- space_sample_labels(c("s1", "s2", "s3"))
  to_sp <- space_sample_labels(c("g1", "g2", "g3"))
  M <- diag(3)
  fam <- MapFamily("test", from_sp, to_sp, by_subject = list(sub1 = M, sub2 = M))
  result <- align_eager(p, fam)
  expect_s3_class(result, "gds")
})

# ===========================================================================
# verb-map-to.R: map_to()
# ===========================================================================

test_that("map_to() adds a map node to the plan", {
  p <- .plan_gds()
  target <- space_sample_labels(c("t1", "t2"))
  M <- matrix(c(1, 0, 0, 1, 0, 1), nrow = 2, ncol = 3)
  mp <- map_linear(space_sample_labels(c("s1", "s2", "s3")), target, M)
  p2 <- map_to(p, target, mp)
  expect_s3_class(p2, "gds_plan")
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last_node$op, "map")
})

test_that("map_to_eager() computes immediately", {
  p <- .plan_gds()
  target <- space_sample_labels(c("t1", "t2"))
  M <- matrix(c(1, 0, 0, 1, 0, 1), nrow = 2, ncol = 3)
  mp <- map_linear(space_sample_labels(c("s1", "s2", "s3")), target, M)
  result <- map_to_eager(p, target, mp)
  expect_s3_class(result, "gds")
  expect_equal(dim(result$assays$beta)[1], 2L)
})

# ===========================================================================
# verb-mask.R: mask()
# ===========================================================================

test_that("mask() adds a mask node to the plan", {
  p <- .plan_gds()
  pol <- MaskPolicy(scope = "group", rule = "intersection")
  p2 <- mask(p, pol)
  expect_s3_class(p2, "gds_plan")
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last_node$op, "mask_policy")
})

test_that("mask_eager() computes immediately", {
  beta <- array(c(1, NA, 3, 4, NA, 6,
                  1, 2, 3, 4, 5, 6), dim = c(3, 2, 2))
  var <- array(0.1, dim = c(3, 2, 2))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = c("s1", "s2"), contrasts = c("c1", "c2")
  )
  p <- as_plan(g)
  pol <- MaskPolicy(scope = "group", rule = "intersection")
  result <- mask_eager(p, pol)
  expect_s3_class(result, "gds")
  # Intersection removes samples with any NA
  expect_true(dim(result$assays$beta)[1] <= 3)
})

# ===========================================================================
# verb-write.R: write_out()
# ===========================================================================

test_that("write_out() adds a write node to the plan", {
  p <- .plan_gds()
  p2 <- write_out(p, path = "/tmp/test.csv", format = "csv")
  expect_s3_class(p2, "gds_plan")
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last_node$op, "write")
})

# ===========================================================================
# align-exec.R: apply_align
# ===========================================================================

test_that("apply_align applies per-subject transforms", {
  beta <- array(c(1, 2, 3, 4, 5, 6), dim = c(3, 2, 1))
  var  <- array(0.1, dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var)

  M <- diag(3)  # identity
  from_sp <- space_sample_labels(c("a", "b", "c"))
  to_sp <- space_sample_labels(c("a", "b", "c"))
  fam <- MapFamily("id", from_sp, to_sp, by_subject = list(sub1 = M, sub2 = M))

  result <- apply_align(fam, arrays, c("sub1", "sub2"), from_sp)
  expect_equal(dim(result$arrays$beta), c(3, 2, 1))
  # Identity transform should preserve values
  expect_equal(result$arrays$beta[, 1, 1], c(1, 2, 3))
})

# ===========================================================================
# map-exec.R: apply_map_to
# ===========================================================================

test_that("apply_map_to applies linear map", {
  beta <- array(c(1, 2, 3, 4, 5, 6), dim = c(3, 2, 1))
  var  <- array(0.1, dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var)

  # Map 3 samples -> 2 samples with a simple averaging matrix
  M <- matrix(c(1, 0, 0, 1, 0, 0), nrow = 2, ncol = 3)
  target <- space_sample_labels(c("t1", "t2"))
  mp <- map_linear(space_sample_labels(c("a", "b", "c")), target, M)
  unc <- UncertaintyRule("independent")

  node <- list(op = "map", target_space = target, map = mp, uncertainty = unc)
  result <- apply_map_to(node, arrays)
  expect_equal(dim(result$arrays$beta)[1], 2L)
})

# ===========================================================================
# mask-exec.R: apply_mask_policy
# ===========================================================================

test_that("apply_mask_policy with group intersection", {
  beta <- array(c(1, NA, 3, 4, 5, 6), dim = c(3, 2, 1))
  var  <- array(0.1, dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var)
  space <- space_sample_labels(c("a", "b", "c"))

  pol <- MaskPolicy(scope = "group", rule = "intersection")
  node <- list(op = "mask_policy", policy = pol)
  result <- apply_mask_policy(node, arrays, space)

  # Sample "b" has NA in subject 1, so intersection removes it
  expect_true(dim(result$arrays$beta)[1] < 3)
})

test_that("apply_mask_policy with group union", {
  beta <- array(c(1, NA, 3, NA, 5, 6), dim = c(3, 2, 1))
  var  <- array(0.1, dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var)
  space <- space_sample_labels(c("a", "b", "c"))

  pol <- MaskPolicy(scope = "group", rule = "union")
  node <- list(op = "mask_policy", policy = pol)
  result <- apply_mask_policy(node, arrays, space)

  # Union keeps samples that are finite in ANY subject
  expect_true(dim(result$arrays$beta)[1] >= 2)
})

# ===========================================================================
# MaskPolicy constructor
# ===========================================================================

test_that("MaskPolicy creates valid policy", {
  pol <- MaskPolicy(scope = "group", rule = "intersection")
  expect_equal(pol$scope, "group")
  expect_equal(pol$rule, "intersection")
})

test_that("MaskPolicy accepts threshold rule with value", {
  pol <- MaskPolicy(scope = "group", rule = "threshold", threshold = 0.5)
  expect_equal(pol$threshold, 0.5)
})
