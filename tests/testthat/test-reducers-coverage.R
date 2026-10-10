# Coverage tests for reduce-exec.R and reducers-core.R

# ===========================================================================
# Helper: small test arrays [sample x subject x contrast]
# ===========================================================================
.mk <- function(vals, ns = 2, nsubj = 3, nc = 1) {
  array(vals, dim = c(ns, nsubj, nc))
}

# ===========================================================================
# .p_from_z
# ===========================================================================

test_that(".p_from_z two-sided", {
  expect_equal(.p_from_z(0), 1.0, tolerance = 0.01)
  expect_true(.p_from_z(2.0) < 0.05)
  expect_equal(.p_from_z(2.0), .p_from_z(-2.0))
})

test_that(".p_from_z one-sided", {
  expect_true(.p_from_z(2.0, "greater") < 0.025)
  expect_true(.p_from_z(-2.0, "less") < 0.025)
  expect_true(.p_from_z(-2.0, "greater") > 0.97)
})

# ===========================================================================
# .colwise_fe
# ===========================================================================

test_that(".colwise_fe produces correct fixed-effects estimates", {
  # 3 subjects, 4 samples
  beta <- matrix(c(1, 2, 3, 2, 3, 4, 1.5, 2.5, 3.5, 1, 1, 1), nrow = 3)
  var  <- matrix(0.1, nrow = 3, ncol = 4)
  result <- .colwise_fe(beta, var)

  expect_true(all(c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2") %in% names(result)))
  expect_equal(length(result$beta_g), 4)
  # With equal variance, FE mean is simply the column mean
  expect_equal(result$beta_g[1], mean(beta[, 1]), tolerance = 1e-6)
})

test_that(".colwise_fe respects min_subj", {
  beta <- matrix(c(1, 2), nrow = 2, ncol = 1)
  var  <- matrix(0.1, nrow = 2, ncol = 1)
  result <- .colwise_fe(beta, var, min_subj = 3L)
  expect_true(is.na(result$beta_g[1]))
})

# ===========================================================================
# .colwise_tau2_dl
# ===========================================================================

test_that(".colwise_tau2_dl returns non-negative tau2", {
  beta <- matrix(c(1, 1.1, 0.9), nrow = 3, ncol = 2)
  var  <- matrix(0.5, nrow = 3, ncol = 2)
  result <- .colwise_tau2_dl(beta, var)
  expect_equal(length(result), 2)
  expect_true(all(result >= 0))
})

# ===========================================================================
# core_meta_fe_kernel
# ===========================================================================

test_that("core_meta_fe_kernel runs FE meta-analysis", {
  beta <- matrix(c(1, 2, 3), nrow = 3, ncol = 5)
  var  <- matrix(0.1, nrow = 3, ncol = 5)
  result <- core_meta_fe_kernel(beta, var)
  expect_equal(length(result$beta_g), 5)
  expect_true(all(!is.na(result$beta_g)))
})

# ===========================================================================
# core_meta_re_dl_kernel
# ===========================================================================

test_that("core_meta_re_dl_kernel runs RE meta-analysis", {
  beta <- matrix(c(1, 5, 3, 2, 6, 4), nrow = 3, ncol = 2)
  var  <- matrix(c(0.1, 0.2, 0.15, 0.1, 0.3, 0.1), nrow = 3, ncol = 2)
  result <- core_meta_re_dl_kernel(beta, var)
  expect_true("tau2" %in% names(result))
  expect_equal(length(result$beta_g), 2)
})

# ===========================================================================
# .stouffer_fallback
# ===========================================================================

test_that(".stouffer_fallback combines z-scores", {
  z <- matrix(c(2.0, 1.5, 1.8), nrow = 3, ncol = 2)
  result <- .stouffer_fallback(z)
  expect_equal(length(result$z_g), 2)
  expect_true(all(result$z_g > 0))
  expect_true(all(result$p_g < 0.05))
})

test_that(".stouffer_fallback with weights", {
  z <- matrix(c(2.0, 1.0, 2.0, 1.0), nrow = 2, ncol = 2)
  result <- .stouffer_fallback(z, weights = c(2, 1))
  expect_equal(length(result$z_g), 2)
})

test_that(".stouffer_fallback with min_subj filtering", {
  z <- matrix(c(2.0), nrow = 1, ncol = 2)
  result <- .stouffer_fallback(z, min_subj = 2L)
  expect_true(all(is.na(result$z_g)))
})

# ===========================================================================
# core_stouffer_kernel
# ===========================================================================

test_that("core_stouffer_kernel dispatches correctly", {
  z <- matrix(c(2, 1.5, 1.8, 2.1, 1.6, 1.9), nrow = 3, ncol = 2)
  result <- core_stouffer_kernel(z = z)
  expect_true("z_g" %in% names(result))
})

# ===========================================================================
# core_fisher_kernel
# ===========================================================================

test_that("core_fisher_kernel combines p-values", {
  p <- matrix(c(0.01, 0.05, 0.1, 0.001, 0.02, 0.04), nrow = 3, ncol = 2)
  result <- core_fisher_kernel(p = p)
  expect_true(all(c("chi2", "df", "p_g") %in% names(result)))
  expect_equal(length(result$p_g), 2)
  # With small p-values, combined should be small
  expect_true(result$p_g[2] < 0.01)
})

test_that("core_fisher_kernel with min_subjects", {
  p <- matrix(0.05, nrow = 1, ncol = 2)
  result <- core_fisher_kernel(p = p, opts = list(min_subjects = 2L))
  expect_true(all(is.na(result$p_g)))
})

# ===========================================================================
# .lancaster_fallback
# ===========================================================================

test_that(".lancaster_fallback combines with df weights", {
  p <- matrix(c(0.01, 0.05, 0.02, 0.04), nrow = 2, ncol = 2)
  dfw <- c(10L, 15L)
  result <- .lancaster_fallback(p, dfw)
  expect_true(all(c("chi2", "df", "p_g") %in% names(result)))
  expect_equal(length(result$p_g), 2)
})

test_that(".lancaster_fallback errors on mismatched dfw length", {
  p <- matrix(0.05, nrow = 2, ncol = 1)
  expect_error(.lancaster_fallback(p, dfw = 1L), "dfw length")
})

# ===========================================================================
# core_lancaster_kernel
# ===========================================================================

test_that("core_lancaster_kernel dispatches correctly", {
  p <- matrix(c(0.01, 0.05), nrow = 2, ncol = 1)
  result <- core_lancaster_kernel(p = p, dfw = c(5L, 10L))
  expect_true("p_g" %in% names(result))
})

# ===========================================================================
# .ensure_required_arrays
# ===========================================================================

test_that(".ensure_required_arrays derives var from se", {
  se <- array(0.5, dim = c(2, 1, 1))
  result <- .ensure_required_arrays(list(se = se), "var")
  expect_true("var" %in% names(result))
  expect_equal(result$var[1,1,1], 0.25)
})

test_that(".ensure_required_arrays derives se from var", {
  var <- array(0.25, dim = c(2, 1, 1))
  result <- .ensure_required_arrays(list(var = var), "se")
  expect_true("se" %in% names(result))
  expect_equal(result$se[1,1,1], 0.5)
})

# ===========================================================================
# register_core_reducers
# ===========================================================================

test_that("register_core_reducers populates reducer registry", {
  register_core_reducers()
  expect_false(is.null(get_reducer("meta:fe")))
  expect_false(is.null(get_reducer("meta:re")))
  expect_false(is.null(get_reducer("combine:stouffer")))
  expect_false(is.null(get_reducer("combine:fisher")))
  expect_false(is.null(get_reducer("combine:lancaster")))
  expect_false(is.null(get_reducer("meta:fe_reg")))
  expect_false(is.null(get_reducer("meta:re_reg")))
  expect_false(is.null(get_reducer("ols:voxelwise")))
})

test_that(".normalize_reducer_name maps aliases", {
  expect_equal(.normalize_reducer_name("fixed"), "meta:fe")
  expect_equal(.normalize_reducer_name("random"), "meta:re")
  expect_equal(.normalize_reducer_name("stouffer"), "combine:stouffer")
  expect_equal(.normalize_reducer_name("fisher"), "combine:fisher")
  # Unknown names pass through
  expect_equal(.normalize_reducer_name("lancaster"), "lancaster")
  expect_equal(.normalize_reducer_name("custom_thing"), "custom_thing")
})

# ===========================================================================
# ols_voxelwise_cpp (R fallback)
# ===========================================================================

test_that("ols_voxelwise_cpp R fallback produces correct OLS", {
  set.seed(42)
  N <- 10; B <- 5
  X <- cbind(1, rnorm(N))
  beta_true <- c(2, 0.5)
  # Y is subjects x samples: each column is one "sample" (voxel)
  Y <- matrix(NA_real_, N, B)
  for (b in seq_len(B)) Y[, b] <- X %*% beta_true + rnorm(N, sd = 0.1)
  result <- ols_voxelwise_cpp(Y, X)

  expect_equal(nrow(result$coef), 2)
  expect_equal(ncol(result$coef), B)
  # Intercept should be close to 2
  expect_equal(mean(result$coef[1, ]), 2.0, tolerance = 0.2)
  expect_equal(length(result$df_res), B)
  expect_true(all(result$df_res == N - 2))
})

test_that("ols_voxelwise_cpp with cov_tri", {
  N <- 5; B <- 3
  X <- cbind(1, 1:N)
  Y <- matrix(rnorm(N * B), N, B)
  result <- ols_voxelwise_cpp(Y, X, return_cov_tri = TRUE)
  expect_true("cov_tri" %in% names(result))
  # 2 params -> 3 upper triangle elements
  expect_equal(nrow(result$cov_tri), 3)
  expect_equal(ncol(result$cov_tri), B)
})

# ===========================================================================
# apply_reduce end-to-end via registered reducers
# ===========================================================================

test_that("apply_reduce uses registered meta:fe reducer", {
  register_core_reducers()
  beta <- .mk(c(1, 2, 1.5, 2.5, 1.2, 2.2))
  var  <- .mk(rep(0.1, 6))
  node <- list(method = "meta:fe", weights = "1/var", options = list(), formula = NULL)
  result <- apply_reduce(node, list(beta = beta, var = var), "1/var", c("s1", "s2", "s3"))
  expect_equal(dim(result$arrays$beta)[2L], 1L)
})

test_that("apply_reduce uses registered combine:stouffer reducer", {
  register_core_reducers()
  z <- .mk(c(2, 1, 2.5, 1.5, 1.8, 1.2))
  node <- list(method = "combine:stouffer", weights = "equal", options = list(), formula = NULL)
  result <- apply_reduce(node, list(z = z), "equal", c("s1", "s2", "s3"))
  expect_true("z" %in% names(result$arrays))
})
