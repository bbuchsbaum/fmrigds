test_that("fixed-effects reduction matches weighted mean", {
  beta <- array(c(1, 2, 1.5, 2.5), c(2, 2, 1))
  var <- array(c(0.25, 0.5, 0.5, 1), c(2, 2, 1))
  df <- array(30, c(2, 2, 1))
  arrays <- list(beta = beta, var = var, df = df)
  node <- op_reduce(method = "fixed", weights = "1/var", by = "contrast", options = list())

  res <- fmrigds:::apply_reduce(node, arrays, NULL, NULL)
  arrays_out <- res$arrays
  expected <- c(
    sum(c(4, 2) * c(1, 1.5)) / sum(c(4, 2)),
    sum(c(2, 1) * c(2, 2.5)) / sum(c(2, 1))
  )
  expect_equal(arrays_out$beta[, 1, 1], expected)
})

test_that("random-effects reduction computes tau2", {
  beta <- array(c(1, 4, 1, 5), c(2, 2, 1))
  var <- array(c(1, 1, 1, 1), c(2, 2, 1))
  df <- array(20, c(2, 2, 1))
  arrays <- list(beta = beta, var = var, df = df)
  node <- op_reduce(method = "random", weights = "1/var", by = "contrast", options = list())

  res <- fmrigds:::apply_reduce(node, arrays, NULL, NULL)
  expect_true(all(res$arrays$var[, 1, 1] > 0))
})

test_that("evidence-only reduction (Stouffer)", {
  z <- array(c(1.96, -1.96), c(1, 2, 1))
  arrays <- list(z = z)
  node <- op_reduce(method = "stouffer", weights = "equal", by = "contrast", options = list())
  res <- fmrigds:::apply_reduce(node, arrays, NULL, NULL)
  expect_equal(res$arrays$z[1, 1, 1], 0, tolerance = 1e-8)
})

test_that("meta FE alternative tail handling works", {
  beta <- array(c(1, 2), c(1, 2, 1))
  var  <- array(c(1, 1), c(1, 2, 1))
  arrays <- list(beta = beta, var = var)
  # FE two-sided
  node_ts <- op_reduce(method = "meta:fe", weights = "1/var", by = "contrast", options = list(alternative = "two.sided"))
  res_ts <- fmrigds:::apply_reduce(node_ts, arrays, NULL, NULL)
  # FE greater (one-sided upper)
  node_gt <- op_reduce(method = "meta:fe", weights = "1/var", by = "contrast", options = list(alternative = "greater"))
  res_gt <- fmrigds:::apply_reduce(node_gt, arrays, NULL, NULL)
  # FE less (one-sided lower)
  node_lt <- op_reduce(method = "meta:fe", weights = "1/var", by = "contrast", options = list(alternative = "less"))
  res_lt <- fmrigds:::apply_reduce(node_lt, arrays, NULL, NULL)

  zval <- as.numeric(res_ts$arrays$z_g[1, 1, 1])
  expect_true(is.finite(zval))
  expect_lt(as.numeric(res_gt$arrays$p_g[1, 1, 1]), as.numeric(res_ts$arrays$p_g[1, 1, 1]))
  expect_gt(as.numeric(res_lt$arrays$p_g[1, 1, 1]), as.numeric(res_ts$arrays$p_g[1, 1, 1]))
})

test_that("Lancaster infers df weights from df", {
  # p-values for two subjects and one sample
  p <- array(c(0.05, 0.10), c(1, 2, 1))
  df <- array(c(10, 20), c(1, 2, 1))
  arrays <- list(p = p, df = df)
  node <- op_reduce(method = "combine:lancaster", weights = "equal", by = "contrast", options = list())
  res <- fmrigds:::apply_reduce(node, arrays, NULL, NULL)
  # df should be sum of 2*wi across subjects (wi ~ df_i)
  expect_true(all(res$arrays$df[, 1, 1] >= 2 * sum(round(c(10, 20)))))
  expect_true(is.finite(res$arrays$p[, 1, 1]))
})
