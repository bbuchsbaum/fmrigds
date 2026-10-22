test_that("ols:voxelwise matches lm() per sample with complete and incomplete samples mixed", {
  set.seed(11)
  N <- 12L
  X <- cbind(1, rnorm(N), rep(0:1, N / 2))
  beta <- matrix(rnorm(N * 6), N, 6)
  beta[c(2, 5), 3] <- NA
  beta[1:9, 5] <- NA
  res <- suppressWarnings(ols_voxelwise_cpp(beta, X, return_cov_tri = TRUE))
  for (b in seq_len(ncol(beta))) {
    ok <- is.finite(beta[, b])
    if (sum(ok) < 4L) {
      expect_true(all(is.na(res$coef[, b])))
      expect_equal(res$n_obs[b], sum(ok))
      next
    }
    fit <- stats::lm(beta[ok, b] ~ X[ok, ] - 1)
    V <- stats::vcov(fit)
    expect_equal(res$coef[, b], unname(stats::coef(fit)), tolerance = 1e-10)
    expect_equal(res$se_coef[, b], unname(sqrt(diag(V))), tolerance = 1e-10)
    expect_equal(res$df_res[b], fit$df.residual)
    expect_equal(res$n_obs[b], sum(ok))
    expect_equal(res$cov_tri[, b], unname(V[upper.tri(V, diag = TRUE)][order(
      row(V)[upper.tri(V, diag = TRUE)], col(V)[upper.tri(V, diag = TRUE)])]),
      tolerance = 1e-10)
  }
})
