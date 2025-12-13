test_that("core_meta_fe_kernel() performs fixed-effects meta-analysis", {
  # Create test data: 4 subjects × 3 contrasts
  set.seed(123)
  beta <- matrix(rnorm(12, mean = 1), nrow = 4, ncol = 3)
  var <- matrix(rexp(12, rate = 10), nrow = 4, ncol = 3)

  result <- core_meta_fe_kernel(beta, var)

  # Should return meta-analysis results
  expect_true("beta_g" %in% names(result))
  expect_true("var_g" %in% names(result))
  expect_true("se_g" %in% names(result))
  expect_true("z_g" %in% names(result))
  expect_true("p_g" %in% names(result))
  expect_true("Q" %in% names(result))
  expect_true("I2" %in% names(result))

  # Check dimensions
  expect_equal(length(result$beta_g), 3)
  expect_equal(length(result$var_g), 3)
  expect_equal(length(result$Q), 3)
})

test_that("core_meta_re_dl_kernel() performs random-effects meta-analysis", {
  # Create test data with heterogeneity
  set.seed(124)
  beta <- matrix(rnorm(12, mean = c(1, 2, 0.5)), nrow = 4, ncol = 3, byrow = TRUE)
  var <- matrix(rexp(12, rate = 10), nrow = 4, ncol = 3)

  result <- core_meta_re_dl_kernel(beta, var)

  # Should return RE results including tau2
  expect_true("beta_g" %in% names(result))
  expect_true("tau2" %in% names(result))
  expect_true("Q" %in% names(result))
  expect_true("I2" %in% names(result))

  # tau2 should be non-negative
  expect_true(all(result$tau2 >= 0))
})

test_that("core_stouffer_kernel() combines z-scores", {
  # Create test z-scores
  set.seed(125)
  z <- matrix(rnorm(12), nrow = 4, ncol = 3)

  result <- core_stouffer_kernel(z = z)

  # Should return combined z and p-values
  expect_true("z_g" %in% names(result))
  expect_true("p_g" %in% names(result))
  expect_equal(length(result$z_g), 3)
})

test_that("core_stouffer_kernel() handles weighted combination", {
  set.seed(126)
  z <- matrix(rnorm(12), nrow = 4, ncol = 3)
  weights <- c(1, 2, 1, 2)  # Different weights per subject

  result <- core_stouffer_kernel(z = z, opts = list(weights = weights))

  expect_true("z_g" %in% names(result))
  expect_true("p_g" %in% names(result))
})

test_that("core_fisher_kernel() combines p-values", {
  # Create test p-values
  set.seed(127)
  p <- matrix(runif(12), nrow = 4, ncol = 3)

  result <- core_fisher_kernel(p = p)

  # Should return chi-square statistic and p-value
  expect_true("chi2" %in% names(result))
  expect_true("df" %in% names(result))
  expect_true("p_g" %in% names(result))
  expect_equal(length(result$chi2), 3)
})

test_that("core_lancaster_kernel() combines p-values with weights", {
  # Create test p-values and df weights
  set.seed(128)
  p <- matrix(runif(12), nrow = 4, ncol = 3)
  dfw <- c(10, 15, 20, 12)  # df per subject

  result <- core_lancaster_kernel(p = p, dfw = dfw)

  # Should return chi-square and p-value
  expect_true("chi2" %in% names(result))
  expect_true("df" %in% names(result))
  expect_true("p_g" %in% names(result))
})

test_that("meta-analysis kernels handle NA values", {
  # Test with missing data
  beta <- matrix(c(1, NA, 2, 3), nrow = 4, ncol = 1)
  var <- matrix(c(0.1, NA, 0.2, 0.15), nrow = 4, ncol = 1)

  result_fe <- core_meta_fe_kernel(beta, var)
  result_re <- core_meta_re_dl_kernel(beta, var)

  # Should complete without error and return valid results
  expect_true(is.finite(result_fe$beta_g[1]))
  expect_true(is.finite(result_re$beta_g[1]))
})

test_that("meta-analysis kernels respect eps parameter", {
  # Create data with very small variances
  beta <- matrix(1:4, nrow = 4, ncol = 1)
  var <- matrix(rep(1e-20, 4), nrow = 4, ncol = 1)

  # Should not crash with small eps
  result <- core_meta_fe_kernel(beta, var, opts = list(eps = 1e-12))
  expect_true(is.finite(result$beta_g[1]))
})

test_that("stouffer handles missing z-scores", {
  z <- matrix(c(1.5, NA, 2.0, -1.0), nrow = 4, ncol = 1)

  result <- core_stouffer_kernel(z = z)
  expect_true(is.finite(result$z_g[1]))
})

test_that("fisher handles edge case p-values", {
  # Test with p-values near 0 and 1
  p <- matrix(c(0.001, 0.999, 0.5, 0.1), nrow = 4, ncol = 1)

  result <- core_fisher_kernel(p = p)
  expect_true(is.finite(result$chi2[1]))
  expect_true(result$p_g[1] >= 0 && result$p_g[1] <= 1)
})
