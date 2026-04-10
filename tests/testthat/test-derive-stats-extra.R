# Extra tests for derive-stats.R — covering all derivation paths

# ===========================================================================
# derive_var
# ===========================================================================

test_that("derive_var computes se^2", {
  se <- array(c(0.5, 1.0, 2.0, 0.3), dim = c(2, 1, 2))
  result <- derive_var(list(se = se))
  expect_equal(result, se^2)
})

test_that("derive_var errors without se", {
  expect_error(derive_var(list(beta = array(1, dim = c(1,1,1)))), "standard error")
})

# ===========================================================================
# derive_se
# ===========================================================================

test_that("derive_se computes sqrt(var)", {
  var <- array(c(0.25, 1.0, 4.0, 0.09), dim = c(2, 1, 2))
  result <- derive_se(list(var = var))
  expect_equal(result, sqrt(var))
})

test_that("derive_se errors without var", {
  expect_error(derive_se(list(beta = array(1, dim = c(1,1,1)))), "variance")
})

# ===========================================================================
# derive_t
# ===========================================================================

test_that("derive_t computes beta / sqrt(var)", {
  beta <- array(c(2, 4), dim = c(2, 1, 1))
  var  <- array(c(1, 4), dim = c(2, 1, 1))
  result <- derive_t(list(beta = beta, var = var))
  expect_equal(result[1,1,1], 2.0)
  expect_equal(result[2,1,1], 2.0)
})

test_that("derive_t errors without beta and var", {
  expect_error(derive_t(list(beta = array(1, dim = c(1,1,1)))), "require")
})

# ===========================================================================
# derive_z
# ===========================================================================

test_that("derive_z from t and df", {
  t_val <- array(c(3.0, -2.0), dim = c(2, 1, 1))
  df <- array(c(20, 20), dim = c(2, 1, 1))
  result <- derive_z(list(t = t_val, df = df))
  expect_equal(dim(result), c(2, 1, 1))
  # z should preserve sign
  expect_true(result[1,1,1] > 0)
  expect_true(result[2,1,1] < 0)
})

test_that("derive_z from p and beta", {
  p <- array(c(0.01, 0.05), dim = c(2, 1, 1))
  beta <- array(c(1, -1), dim = c(2, 1, 1))
  result <- derive_z(list(p = p, beta = beta))
  expect_equal(dim(result), c(2, 1, 1))
  expect_true(result[1,1,1] > 0)
  expect_true(result[2,1,1] < 0)
})

test_that("derive_z from F, df1, df2, beta", {
  F_val <- array(c(9.0), dim = c(1, 1, 1))
  df1 <- array(c(1), dim = c(1, 1, 1))
  df2 <- array(c(30), dim = c(1, 1, 1))
  beta <- array(c(1), dim = c(1, 1, 1))
  result <- derive_z(list(F = F_val, df1 = df1, df2 = df2, beta = beta))
  expect_equal(dim(result), c(1, 1, 1))
  expect_true(result[1,1,1] > 0)
})

test_that("derive_z errors without enough inputs", {
  expect_error(derive_z(list(beta = array(1, dim = c(1,1,1)))), "require")
})

test_that("derive_z from p errors without beta for sign", {
  p <- array(0.05, dim = c(1,1,1))
  expect_error(derive_z(list(p = p)), "without effect sign")
})

test_that("derive_z from F errors without beta for sign", {
  expect_error(
    derive_z(list(F = array(4, dim = c(1,1,1)),
                  df1 = array(1, dim = c(1,1,1)),
                  df2 = array(20, dim = c(1,1,1)))),
    "without effect sign"
  )
})

# ===========================================================================
# derive_p
# ===========================================================================

test_that("derive_p from t and df", {
  t_val <- array(c(0, 2.0), dim = c(2, 1, 1))
  df <- array(c(100, 100), dim = c(2, 1, 1))
  result <- derive_p(list(t = t_val, df = df))
  expect_equal(dim(result), c(2, 1, 1))
  expect_equal(result[1,1,1], 1.0, tolerance = 0.01)
  expect_true(result[2,1,1] < 0.05)
})

test_that("derive_p from z", {
  z <- array(c(0, 2.0, -3.0), dim = c(3, 1, 1))
  result <- derive_p(list(z = z))
  expect_equal(dim(result), c(3, 1, 1))
  expect_equal(result[1,1,1], 1.0, tolerance = 0.01)
  expect_true(result[2,1,1] < 0.05)
  expect_true(result[3,1,1] < 0.01)
})

test_that("derive_p from F, df1, df2", {
  F_val <- array(c(0.1, 10.0), dim = c(2, 1, 1))
  df1 <- array(c(1, 1), dim = c(2, 1, 1))
  df2 <- array(c(30, 30), dim = c(2, 1, 1))
  result <- derive_p(list(F = F_val, df1 = df1, df2 = df2))
  expect_true(result[1,1,1] > 0.5)
  expect_true(result[2,1,1] < 0.01)
})

test_that("derive_p from chi2 and df", {
  chi2 <- array(c(1.0, 20.0), dim = c(2, 1, 1))
  df <- array(c(1, 1), dim = c(2, 1, 1))
  result <- derive_p(list(chi2 = chi2, df = df))
  expect_true(result[1,1,1] > 0.1)
  expect_true(result[2,1,1] < 0.001)
})

test_that("derive_p errors without enough inputs", {
  expect_error(derive_p(list(beta = array(1, dim = c(1,1,1)))), "missing required")
})

# ===========================================================================
# execute_derive
# ===========================================================================

test_that("execute_derive dispatches correctly", {
  beta <- array(c(2, 4), dim = c(2, 1, 1))
  var  <- array(c(1, 4), dim = c(2, 1, 1))
  result <- execute_derive(list(beta = beta, var = var), c("se", "t"))
  expect_true("se" %in% names(result))
  expect_true("t" %in% names(result))
  expect_equal(result$se, sqrt(var))
})

test_that("execute_derive skips existing assays by default", {
  beta <- array(1, dim = c(1, 1, 1))
  var  <- array(0.25, dim = c(1, 1, 1))
  se_orig <- array(999, dim = c(1, 1, 1))
  result <- execute_derive(list(beta = beta, var = var, se = se_orig), "se")
  expect_equal(result$se[1,1,1], 999)
})

test_that("execute_derive overwrites with overwrite=TRUE", {
  var  <- array(0.25, dim = c(1, 1, 1))
  se_orig <- array(999, dim = c(1, 1, 1))
  result <- execute_derive(list(var = var, se = se_orig), "se",
                           options = list(overwrite = TRUE))
  expect_equal(result$se[1,1,1], 0.5)
})

test_that("execute_derive errors on unknown target", {
  expect_error(execute_derive(list(), "unknown_stat"), "Unknown derivation")
})
