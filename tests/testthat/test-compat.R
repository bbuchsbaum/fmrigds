test_that("assert_compatible_spaces checks fields", {
  sp1 <- space_voxel(dim = c(2,2,2), affine = diag(4), storage = "dense")
  sp2 <- space_voxel(dim = c(2,2,2), affine = diag(4), storage = "dense")
  expect_true(assert_compatible_spaces(sp1, sp2))
  sp3 <- space_voxel(dim = c(2,2,3), affine = diag(4), storage = "dense")
  expect_error(assert_compatible_spaces(sp1, sp3), "dim")
})

test_that("harmonise_contrasts renames contrasts and dimnames", {
  beta <- array(runif(2*2*2), dim = c(2,2,2), dimnames = list(NULL, c("s1","s2"), c("a","b")))
  var <- array(runif(2*2*2), dim = c(2,2,2), dimnames = dimnames(beta))
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("x","y")), c("s1","s2"), c("a","b"))
  g2 <- harmonise_contrasts(g, c(a = "A", b = "B"))
  expect_equal(contrasts(g2), c("A","B"))
  expect_equal(dimnames(assay(g2, "beta"))[[3]], c("A","B"))
})

test_that("relabel_subjects updates subjects, col_data, and dimnames", {
  beta <- array(runif(2*2*1), dim = c(2,2,1), dimnames = list(NULL, c("u1","u2"), "c1"))
  var <- array(runif(2*2*1), dim = c(2,2,1), dimnames = dimnames(beta))
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("x","y")), c("u1","u2"), c("c1"), col_data = data.frame(age=c(1,2), row.names = c("u1","u2")))
  g2 <- relabel_subjects(g, c(u1 = "s1", u2 = "s2"))
  expect_equal(subjects(g2), c("s1","s2"))
  expect_equal(rownames(col_data(g2)), c("s1","s2"))
  expect_equal(dimnames(assay(g2, "beta"))[[2]], c("s1","s2"))
})

test_that("attach_weight and use_weight integrate with reduce", {
  set.seed(1)
  beta <- array(rnorm(3*3*1), dim = c(3,3,1))
  var <- array(runif(3*3*1, min=0.05, max=0.2), dim = c(3,3,1))
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(paste0("i",1:3)), paste0("s",1:3), "c1")
  # Custom weights favor subject 1 strongly
  w <- array(1, dim = dim(beta)); w[,1,] <- 100
  g <- attach_weight(g, "w_custom", w)
  opts <- use_weight(g, "w_custom")
  plan <- reduce(as_plan(gds_source("tabular", list(path="dummy"), list(space=g$space, subjects=g$subjects, contrasts=g$contrasts, assays = names(assays(g)), dims = dim(beta)))), method = "fixed", weights = opts$weights, options = opts$options)
  # Execute by bypassing adapter read using compute(assays=...) pattern
  # Instead, call reducer directly to ensure no I/O
  res <- fmrigds:::apply_reduce(list(method = "meta:fe", options = opts$options), arrays = assays(g), weights = opts$weights, subjects = subjects(g))
  expect_true("beta" %in% names(res$arrays))
})

