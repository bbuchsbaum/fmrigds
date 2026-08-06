test_that("random-effects screening holds full tau2 fixed explicitly", {
  beta <- matrix(c(0.2, 0.8, 1.1, 0.4, 1.5, 0.7), 6, 1)
  var <- matrix(c(0.1, 0.2, 0.15, 0.3, 0.12, 0.25), 6, 1)
  fit <- fmrigds:::core_meta_re_dl_kernel(beta, var)
  diag <- fmrigds:::.diagnose_meta_re_block(
    fit, beta, var, NULL,
    matrix(1, 1, 1, dimnames = list("pooled_effect", "pooled_effect")),
    list(), examination_control()$tolerance
  )
  tau2 <- fit$tau2[1]
  w <- 1 / (var[, 1] + tau2)
  full_mu <- weighted.mean(beta[, 1], w)
  full_z <- full_mu * sqrt(sum(w))
  for (i in seq_len(nrow(beta))) {
    wm <- w[-i]
    mu_minus <- weighted.mean(beta[-i, 1], wm)
    z_minus <- mu_minus * sqrt(sum(wm))
    expect_equal(diag$expected[i, 1], mu_minus, tolerance = 1e-12)
    expect_equal(
      diag$delta_stat[i, 1, 1],
      unname(full_z - z_minus),
      tolerance = 1e-12
    )
  }
  expect_identical(diag$mode, "tau2_fixed_full")
})

test_that("random-effects meta-regression fixed-tau deletion matches refit weights", {
  X <- cbind("(Intercept)" = 1, group = c(0, 0, 0, 1, 1, 1, 1, 1))
  beta <- matrix(c(0.1, 0.3, -0.1, 1.0, 1.4, 0.8, 1.3, 0.9), 8, 1)
  var <- matrix(seq(0.1, 0.35, length.out = 8), 8, 1)
  reducer <- get_reducer("meta:re_reg")
  fit <- reducer$fun(beta, var, X, NULL, NULL, NULL, NULL, NULL, list())
  estimands <- matrix(c(0, 1), 1, dimnames = list("group", colnames(X)))
  diag <- fmrigds:::.diagnose_meta_re_reg_block(
    fit, beta, var, X, estimands, list(), examination_control()$tolerance
  )
  tau2 <- fit$tau2[1]
  w <- 1 / (var[, 1] + tau2)
  A <- solve(crossprod(X * sqrt(w)))
  theta <- drop(A %*% crossprod(X, w * beta[, 1]))
  full_z <- theta[2] / sqrt(A[2, 2])
  for (i in seq_len(nrow(beta))) {
    Xm <- X[-i, , drop = FALSE]
    wm <- w[-i]
    Am <- solve(crossprod(Xm * sqrt(wm)))
    theta_minus <- drop(Am %*% crossprod(Xm, wm * beta[-i, 1]))
    z_minus <- theta_minus[2] / sqrt(Am[2, 2])
    expect_equal(
      diag$delta_stat[i, 1, 1],
      unname(full_z - z_minus),
      tolerance = 1e-10
    )
  }
  expect_identical(diag$mode, "tau2_fixed_full")
})

test_that("selected exact random-effects meta-regression re-estimates tau2", {
  X <- cbind("(Intercept)" = 1, group = c(0, 0, 0, 0, 1, 1, 1, 1, 1))
  beta <- cbind(
    c(-0.2, 0.1, 0.3, -0.1, 1.0, 1.4, 0.7, 1.3, 0.9),
    c(0.1, -0.1, 0.2, 0, 0.8, 1.2, 0.6, 1.0, 1.4)
  )
  var <- matrix(seq(0.08, 0.4, length.out = 18), 9, 2)
  reducer <- get_reducer("meta:re_reg")
  full <- reducer$fun(beta, var, X, NULL, NULL, NULL, NULL, NULL, list())
  estimands <- matrix(c(0, 1), 1, dimnames = list("group", colnames(X)))
  exact <- fmrigds:::.exact_random_deletion_block(
    beta, var, X, estimands, full,
    selected_index = 2L,
    reducer = reducer,
    options = list(),
    tolerance = examination_control()$tolerance
  )
  deleted <- reducer$fun(
    beta[-2, , drop = FALSE],
    var[-2, , drop = FALSE],
    X[-2, , drop = FALSE],
    NULL, NULL, NULL, NULL, NULL, list()
  )
  deleted_z <- deleted$coef[2, ] / deleted$se_coef[2, ]
  full_z <- full$coef[2, ] / full$se_coef[2, ]
  expect_equal(
    drop(exact$delta_stat[1, 1, ]),
    full_z - deleted_z,
    tolerance = 1e-10
  )
  expect_equal(drop(exact$tau2_deleted[1, ]), deleted$tau2, tolerance = 1e-12)
})
