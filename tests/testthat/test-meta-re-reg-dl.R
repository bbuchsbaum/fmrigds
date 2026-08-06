.dl_meta_reg_reference <- function(y, v, X, eps = 1e-12) {
  ok <- is.finite(y) & is.finite(v) & v > 0 &
    rowSums(!is.finite(X)) == 0L
  y <- y[ok]
  v <- pmax(v[ok], eps)
  X <- X[ok, , drop = FALSE]
  if (nrow(X) < ncol(X) + 1L || qr(X)$rank < ncol(X)) return(NULL)

  W <- diag(1 / v)
  A <- solve(crossprod(X, W %*% X))
  beta_fe <- A %*% crossprod(X, W %*% y)
  residual <- y - drop(X %*% beta_fe)
  Q <- drop(crossprod(residual, W %*% residual))
  P <- W - W %*% X %*% A %*% t(X) %*% W
  df_res <- nrow(X) - ncol(X)
  tau2 <- max(0, (Q - df_res) / sum(diag(P)))

  Wstar <- diag(1 / (v + tau2))
  Astar <- solve(crossprod(X, Wstar %*% X))
  list(
    coef = drop(Astar %*% crossprod(X, Wstar %*% y)),
    se_coef = sqrt(diag(Astar)),
    tau2 = tau2,
    Q = Q,
    df_res = df_res
  )
}

.expect_meta_reg_equal <- function(observed, expected, column = 1L, tolerance = 1e-10) {
  expect_equal(unname(observed$coef[, column]), unname(expected$coef), tolerance = tolerance)
  expect_equal(unname(observed$se_coef[, column]), unname(expected$se_coef), tolerance = tolerance)
  expect_equal(observed$tau2[column], expected$tau2, tolerance = tolerance)
  expect_equal(observed$Q[column], expected$Q, tolerance = tolerance)
  expect_equal(observed$df_res[column], expected$df_res, tolerance = tolerance)
}

test_that("DL meta-regression uses tr(P), including both weight factors", {
  y <- c(0.1, 0.4, 1.2, 1.0, 2.2, 2.8)
  v <- c(0.02, 0.08, 0.15, 0.25, 0.5, 0.9)
  X <- model.matrix(~ c(0, 0, 0, 1, 1, 1))
  expected <- .dl_meta_reg_reference(y, v, X)

  observed_r <- fmrigds:::.meta_re_reg_dl_r(
    matrix(y, ncol = 1L), matrix(v, ncol = 1L), X
  )
  .expect_meta_reg_equal(observed_r, expected)
  expect_equal(observed_r$tau2, 0.209463037073, tolerance = 1e-12)

  observed_cpp <- fmrigds:::meta_re_reg_dl_cpp(
    matrix(y, ncol = 1L), matrix(v, ncol = 1L), X
  )
  .expect_meta_reg_equal(observed_cpp, expected)
})

test_that("DL meta-regression agrees with metafor", {
  skip_if_not_installed("metafor")

  y <- c(0.1, 0.4, 1.2, 1.0, 2.2, 2.8)
  v <- c(0.02, 0.08, 0.15, 0.25, 0.5, 0.9)
  group <- c(0, 0, 0, 1, 1, 1)
  X <- model.matrix(~ group)
  fit <- metafor::rma.uni(yi = y, vi = v, mods = ~ group, method = "DL")

  observed <- fmrigds:::.meta_re_reg_dl_r(
    matrix(y, ncol = 1L), matrix(v, ncol = 1L), X
  )
  expect_equal(unname(drop(observed$coef)), unname(drop(fit$beta)), tolerance = 1e-10)
  expect_equal(unname(drop(observed$se_coef)), unname(fit$se), tolerance = 1e-10)
  expect_equal(observed$tau2, fit$tau2, tolerance = 1e-10)
  expect_equal(observed$Q, fit$QE, tolerance = 1e-10)
  expect_equal(observed$df_res, fit$k - fit$p)
})

test_that("R and C++ meta-regression share validity and missingness rules", {
  X <- cbind(`(Intercept)` = 1, group = c(0, 0, 0, 1, 1, 1, 1))
  beta <- cbind(
    complete = c(0.1, 0.4, 1.2, 1.0, 2.2, 2.8, 2.4),
    missing_effect = c(0.1, NA, 1.2, 1.0, 2.2, 2.8, 2.4),
    invalid_variance = c(0.1, 0.4, 1.2, 1.0, 2.2, 2.8, 2.4),
    insufficient = c(0.1, NA, NA, 1.0, NA, NA, NA)
  )
  var <- cbind(
    complete = c(0.02, 0.08, 0.15, 0.25, 0.5, 0.9, 0.7),
    missing_effect = c(0.02, 0.08, 0.15, 0.25, 0.5, 0.9, 0.7),
    invalid_variance = c(0.02, 0, -1, 0.25, 0.5, 0.9, 0.7),
    insufficient = rep(0.2, 7)
  )

  observed_r <- fmrigds:::.meta_re_reg_dl_r(beta, var, X)
  observed_cpp <- fmrigds:::meta_re_reg_dl_cpp(beta, var, X)
  expect_equal(observed_cpp, observed_r, tolerance = 1e-10)
  expect_true(all(is.na(observed_r$coef[, 4L])))
  expect_true(all(is.na(observed_r$se_coef[, 4L])))
  expect_true(is.na(observed_r$tau2[4L]))
  expect_true(is.na(observed_r$Q[4L]))
  expect_true(is.na(observed_r$df_res[4L]))

  X_missing <- X
  X_missing[2L, 2L] <- NA_real_
  expect_equal(
    fmrigds:::meta_re_reg_dl_cpp(beta, var, X_missing),
    fmrigds:::.meta_re_reg_dl_r(beta, var, X_missing),
    tolerance = 1e-10
  )
})

test_that("rank loss and minimum-subject boundaries are nonestimable", {
  y <- matrix(c(1, 2, 3, 4), ncol = 1L)
  v <- matrix(rep(0.2, 4), ncol = 1L)
  singular_X <- cbind(1, c(0, 1, 2, 3), c(0, 2, 4, 6))

  for (observed in list(
    fmrigds:::.meta_re_reg_dl_r(y, v, singular_X),
    fmrigds:::meta_re_reg_dl_cpp(y, v, singular_X)
  )) {
    expect_true(all(is.na(observed$coef)))
    expect_true(all(is.na(observed$se_coef)))
    expect_true(is.na(observed$tau2))
    expect_true(is.na(observed$Q))
    expect_true(is.na(observed$df_res))
  }

  X <- cbind(1, c(0, 0, 1, 1))
  accepted <- fmrigds:::.meta_re_reg_dl_r(y, v, X, min_subj = 3L)
  rejected <- fmrigds:::.meta_re_reg_dl_r(y, v, X, min_subj = 5L)
  expect_true(all(is.finite(unlist(accepted))))
  expect_true(all(is.na(unlist(rejected))))
})

test_that("case-deletion reference fits preserve validity and expose rank loss", {
  y <- c(0.2, 0.4, 0.7, 1.5, 1.7, 2.0)
  v <- c(0.1, 0.2, 0.15, 0.1, 0.25, 0.2)
  X <- cbind(`(Intercept)` = 1, group = c(0, 0, 0, 1, 1, 1))

  for (i in seq_along(y)) {
    keep <- seq_along(y) != i
    expected <- .dl_meta_reg_reference(y[keep], v[keep], X[keep, , drop = FALSE])
    observed <- fmrigds:::.meta_re_reg_dl_r(
      matrix(y[keep], ncol = 1L),
      matrix(v[keep], ncol = 1L),
      X[keep, , drop = FALSE]
    )
    .expect_meta_reg_equal(observed, expected)
  }

  sparse_X <- cbind(`(Intercept)` = 1, group = c(0, 0, 0, 1))
  keep <- 1:3
  lost <- fmrigds:::.meta_re_reg_dl_r(
    matrix(y[keep], ncol = 1L),
    matrix(v[keep], ncol = 1L),
    sparse_X[keep, , drop = FALSE]
  )
  expect_true(all(is.na(unlist(lost))))
})

test_that("meta-regression validates array and design dimensions", {
  beta <- matrix(1:6, nrow = 3)
  var <- matrix(1, nrow = 3, ncol = 2)
  X <- cbind(1, c(0, 1, 1))

  expect_error(
    fmrigds:::.meta_re_reg_dl_r(beta, var[, 1L, drop = FALSE], X),
    "identical dimensions"
  )
  expect_error(
    fmrigds:::core_meta_re_reg_dl_kernel(beta, var, X[-1L, , drop = FALSE]),
    "one row per subject"
  )
})
