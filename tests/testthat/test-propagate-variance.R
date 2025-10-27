test_that("propagate_variance_independent matches analytical result", {
  M <- matrix(c(0.6, 0.8, 0.2, 0.4), nrow = 2, byrow = TRUE)
  beta <- array(c(1, 2), c(2, 1, 1))
  var <- array(c(0.5, 1.0), c(2, 1, 1))

  res <- gdsfmri:::propagate_variance_independent(M, beta, var)

  expected_beta <- M %*% beta[, 1, 1]
  expected_var <- (M^2) %*% var[, 1, 1]

  expect_equal(res$beta[, 1, 1], as.vector(expected_beta))
  expect_equal(res$var[, 1, 1], as.vector(expected_var))
})

test_that("propagate_variance_covariance uses cov_provider", {
  M <- matrix(c(0.6, 0.8), nrow = 1)
  beta <- array(c(1, 2), c(2, 1, 1))
  var <- array(c(0.4, 0.7), c(2, 1, 1))
  Sigma <- matrix(c(0.4, 0.1, 0.1, 0.7), nrow = 2)

  cov_provider <- function(idx) Sigma[idx, idx, drop = FALSE]
  res <- gdsfmri:::propagate_variance_covariance(M, beta, var, cov_provider)

  expect_equal(res$beta[1, 1, 1], sum(M[1, ] * beta[, 1, 1]))
  expect_equal(res$var[1, 1, 1], as.numeric(M[1, ] %*% Sigma %*% M[1, ]), tolerance = 1e-8)
})
