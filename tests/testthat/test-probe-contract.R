test_that("memory adapter probe respects contract", {
  set.seed(1)
  beta <- array(rnorm(4*3*2), dim = c(4,3,2))
  var  <- array(rexp(4*3*2),  dim = c(4,3,2))
  arrs <- list(beta = beta, var = var)
  # Force memory adapter path
  pl <- gds(arrs, format = "memory")
  pr <- pl$source$probe
  expect_true(is.list(pr))
  expect_true(is.character(pr$assays))
  expect_true(inherits(pr$dims, "gds_dims"))
  expect_equal(length(pr$subjects), pr$dims[2])
  expect_equal(length(pr$contrasts), pr$dims[3])
  expect_true(inherits(pr$space, "gds_space"))
})

