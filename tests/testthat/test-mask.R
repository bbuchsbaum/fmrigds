test_that("mask policy removes samples below threshold", {
  beta <- array(c(1, NA, 2, 3, 4, 5), c(3, 2, 1))
  var <- array(1, dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var)
  space <- space_voxel(c(3, 1, 1), diag(4), mask_idx = 1:3, storage = "packed")

  policy <- MaskPolicy(scope = "group", rule = "threshold", threshold = 0.75)
  res <- apply_mask_policy(op_mask_policy(policy), arrays, space)

  expect_equal(dim(res$arrays$beta)[1], 2)
  expect_equal(res$space$mask_idx, c(1, 3))
})
