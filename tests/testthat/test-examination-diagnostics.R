.diagnostic_tolerance <- function() {
  examination_control()$tolerance
}

test_that("fixed-effect prediction and deletion statistics match brute force", {
  beta <- matrix(c(
    1.1, 0.8, -0.4,
    0.9, 1.2, -0.1,
    1.4, 0.7, 0.2,
    1.0, 1.0, -0.3,
    0.7, 0.9, 0.1
  ), nrow = 5, byrow = TRUE)
  var <- matrix(c(0.2, 0.3, 0.5, 0.25, 0.4), 5, 3)
  fit <- fmrigds:::core_meta_fe_kernel(beta, var)
  diag <- fmrigds:::.diagnose_meta_fe_block(
    fit, beta, var, NULL,
    matrix(1, 1, 1, dimnames = list("pooled_effect", "pooled_effect")),
    list(), .diagnostic_tolerance()
  )

  for (b in seq_len(ncol(beta))) {
    for (i in seq_len(nrow(beta))) {
      deleted <- fmrigds:::core_meta_fe_kernel(
        beta[-i, b, drop = FALSE],
        var[-i, b, drop = FALSE]
      )
      expected_resid <- (beta[i, b] - deleted$beta_g) /
        sqrt(var[i, b] + deleted$var_g)
      expect_equal(diag$expected[i, b], deleted$beta_g, tolerance = 1e-12)
      expect_equal(diag$predictive_resid[i, b], expected_resid, tolerance = 1e-12)
      expect_equal(
        diag$delta_effect[i, 1, b],
        fit$beta_g[b] - deleted$beta_g,
        tolerance = 1e-12
      )
      expect_equal(
        diag$delta_stat[i, 1, b],
        fit$z_g[b] - deleted$z_g,
        tolerance = 1e-12
      )
    }
  }
})

test_that("fixed-effect meta-regression deletion matches brute force", {
  X <- cbind("(Intercept)" = 1, age = c(-2, -1, 0, 1, 2, 3, 4))
  beta <- cbind(
    0.5 + 0.3 * X[, "age"] + c(-0.1, 0.2, 0, -0.2, 0.1, 0.15, -0.05),
    -0.2 + 0.1 * X[, "age"] + c(0.1, -0.1, 0.2, 0, -0.2, 0.1, -0.05)
  )
  var <- matrix(seq(0.15, 0.45, length.out = 14), 7, 2)
  reducer <- get_reducer("meta:fe_reg")
  fit <- reducer$fun(beta, var, X, NULL, NULL, NULL, NULL, NULL, list())
  estimands <- rbind(
    age = c(0, 1),
    at_two = c(1, 2)
  )
  colnames(estimands) <- colnames(X)
  diag <- fmrigds:::.diagnose_meta_fe_reg_block(
    fit, beta, var, X, estimands, list(), .diagnostic_tolerance()
  )

  for (b in seq_len(ncol(beta))) {
    full_effect <- drop(estimands %*% fit$coef[, b])
    full_se <- sqrt(rowSums((estimands %*%
      solve(crossprod(X * sqrt(1 / var[, b])))) * estimands))
    full_z <- full_effect / full_se
    for (i in seq_len(nrow(beta))) {
      deleted <- reducer$fun(
        beta[-i, b, drop = FALSE],
        var[-i, b, drop = FALSE],
        X[-i, , drop = FALSE],
        NULL, NULL, NULL, NULL, NULL, list()
      )
      deleted_effect <- drop(estimands %*% deleted$coef[, 1])
      Wm <- 1 / var[-i, b]
      Am <- solve(crossprod(X[-i, , drop = FALSE] * sqrt(Wm)))
      deleted_se <- sqrt(rowSums((estimands %*% Am) * estimands))
      deleted_z <- deleted_effect / deleted_se
      predicted <- sum(X[i, ] * deleted$coef[, 1])
      prediction_var <- var[i, b] + drop(X[i, , drop = FALSE] %*% Am %*% X[i, ])
      expect_equal(diag$expected[i, b], predicted, tolerance = 1e-10)
      expect_equal(
        diag$predictive_resid[i, b],
        (beta[i, b] - predicted) / sqrt(prediction_var),
        tolerance = 1e-10
      )
      expect_equal(
        unname(diag$delta_effect[i, , b]),
        unname(full_effect - deleted_effect),
        tolerance = 1e-10
      )
      expect_equal(
        unname(diag$delta_stat[i, , b]),
        unname(full_z - deleted_z),
        tolerance = 1e-10
      )
    }
  }
})

test_that("OLS externally studentized residual and deleted t match brute force", {
  X <- cbind("(Intercept)" = 1, group = c(0, 0, 0, 1, 1, 1, 1, 1))
  beta <- cbind(
    c(0.2, -0.1, 0.1, 1.2, 0.8, 1.1, 1.0, 1.3),
    c(-0.2, 0.1, 0, 0.4, 0.7, 0.6, 0.5, 0.8)
  )
  fit <- fmrigds:::ols_voxelwise_cpp(beta, X)
  estimands <- matrix(c(0, 1), 1, dimnames = list("group", colnames(X)))
  diag <- fmrigds:::.diagnose_ols_block(
    fit, beta, NULL, X, estimands, list(), .diagnostic_tolerance()
  )
  full_t <- drop(estimands %*% fit$coef) /
    sqrt(drop(estimands %*% solve(crossprod(X)) %*% estimands[1, ]) * fit$sigma2)

  for (b in seq_len(ncol(beta))) {
    for (i in seq_len(nrow(beta))) {
      Xm <- X[-i, , drop = FALSE]
      ym <- beta[-i, b, drop = FALSE]
      deleted <- fmrigds:::ols_voxelwise_cpp(ym, Xm)
      Am <- solve(crossprod(Xm))
      effect_minus <- drop(estimands %*% deleted$coef[, 1])
      se_minus <- sqrt(drop(estimands %*% Am %*% estimands[1, ]) * deleted$sigma2)
      t_minus <- effect_minus / se_minus
      predicted <- sum(X[i, ] * deleted$coef[, 1])
      prediction_se <- sqrt(deleted$sigma2 *
        (1 + drop(X[i, , drop = FALSE] %*% Am %*% X[i, ])))
      expect_equal(diag$expected[i, b], predicted, tolerance = 1e-10)
      expect_equal(
        diag$predictive_resid[i, b],
        (beta[i, b] - predicted) / prediction_se,
        tolerance = 1e-10
      )
      expect_equal(diag$deleted_stat[i, 1, b], unname(t_minus), tolerance = 1e-10)
      expect_equal(
        diag$delta_stat[i, 1, b],
        unname(full_t[b] - t_minus),
        tolerance = 1e-10
      )
    }
  }
})

test_that("deleted rank loss is nonestimable rather than clamped", {
  X <- cbind("(Intercept)" = 1, singleton = c(1, 0, 0, 0, 0))
  beta <- matrix(c(3, 1, 1.1, 0.9, 1.2), 5, 1)
  fit <- fmrigds:::ols_voxelwise_cpp(beta, X)
  estimands <- matrix(c(0, 1), 1, dimnames = list("singleton", colnames(X)))
  diag <- fmrigds:::.diagnose_ols_block(
    fit, beta, NULL, X, estimands, list(), .diagnostic_tolerance()
  )
  expect_false(diag$influence_eligible[1, 1, 1])
  expect_true(is.na(diag$delta_stat[1, 1, 1]))
  expect_true(is.na(diag$predictive_resid[1, 1]))
})

test_that("samplewise missingness uses the same validity pattern as brute force", {
  X <- cbind("(Intercept)" = 1, age = seq(-2, 3, length.out = 7))
  beta <- cbind(
    1 + 0.2 * X[, "age"] + c(0.1, -0.1, 0, 0.2, -0.2, 0.1, 0),
    -0.5 + 0.1 * X[, "age"] + c(0, 0.1, -0.1, 0.2, 0, -0.2, 0.1)
  )
  var <- matrix(0.25, 7, 2)
  beta[2, 1] <- NA_real_
  var[3, 2] <- 0
  reducer <- get_reducer("meta:fe_reg")
  fit <- reducer$fun(beta, var, X, NULL, NULL, NULL, NULL, NULL, list())
  estimands <- matrix(c(0, 1), 1, dimnames = list("age", colnames(X)))
  diag <- fmrigds:::.diagnose_meta_fe_reg_block(
    fit, beta, var, X, estimands, list(), .diagnostic_tolerance()
  )
  expect_false(diag$surprise_eligible[2, 1])
  expect_false(diag$surprise_eligible[3, 2])

  for (b in seq_len(ncol(beta))) {
    valid <- is.finite(beta[, b]) & is.finite(var[, b]) & var[, b] > 0
    full_A <- solve(crossprod(X[valid, , drop = FALSE] * sqrt(1 / var[valid, b])))
    full_z <- drop(estimands %*% fit$coef[, b]) /
      sqrt(drop(estimands %*% full_A %*% estimands[1, ]))
    for (i in which(valid)) {
      keep <- valid
      keep[i] <- FALSE
      deleted <- reducer$fun(
        beta[keep, b, drop = FALSE],
        var[keep, b, drop = FALSE],
        X[keep, , drop = FALSE],
        NULL, NULL, NULL, NULL, NULL, list()
      )
      Am <- solve(crossprod(X[keep, , drop = FALSE] * sqrt(1 / var[keep, b])))
      deleted_z <- drop(estimands %*% deleted$coef[, 1]) /
        sqrt(drop(estimands %*% Am %*% estimands[1, ]))
      expect_equal(
        diag$delta_stat[i, 1, b],
        unname(full_z - deleted_z),
        tolerance = 1e-10
      )
    }
  }
})
