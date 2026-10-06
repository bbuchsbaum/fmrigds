# Regression tests for numerical review fixes ---------------------------------

# 1. Permutation FWER honours one-sided alternatives --------------------------

test_that("perm:onesample FWER uses the requested tail", {
  set.seed(11)
  n_subject <- 12L
  n_sample <- 5L
  beta <- matrix(stats::rnorm(n_subject * n_sample, mean = -2, sd = 0.5), n_subject, n_sample)

  greater <- fmrigds:::core_perm_onesample_kernel(
    beta,
    opts = list(n_perm = 199L, seed = 3L, alternative = "greater")
  )
  # A strongly negative effect is no evidence for "greater".
  expect_true(all(greater$p_perm > 0.9))
  expect_true(all(greater$p_fwer > 0.9))

  less <- fmrigds:::core_perm_onesample_kernel(
    beta,
    opts = list(n_perm = 199L, seed = 3L, alternative = "less")
  )
  expect_true(all(less$p_fwer < 0.05))
  expect_true(all(less$p_fwer >= less$p_perm))
})

test_that("perm:twosample FWER uses the requested tail", {
  set.seed(12)
  group <- rep(c(0L, 1L), each = 8L)
  beta <- matrix(stats::rnorm(16L * 4L), 16L, 4L)
  beta[group == 1L, ] <- beta[group == 1L, ] - 3
  res <- fmrigds:::core_perm_twosample_kernel(
    beta,
    opts = list(group = group, n_perm = 199L, seed = 5L, alternative = "greater")
  )
  expect_true(all(res$t_g < 0))
  expect_true(all(res$p_fwer > 0.9))
})

# 2. Two-sample observed statistic and one p-value convention -----------------

test_that("perm:twosample observed statistic does not depend on include_observed", {
  set.seed(21)
  group <- rep(c(0L, 1L), c(7L, 6L))
  beta <- matrix(stats::rnorm(13L * 3L), 13L, 3L)
  beta[group == 1L, ] <- beta[group == 1L, ] + 1

  with_obs <- fmrigds:::core_perm_twosample_kernel(
    beta,
    opts = list(group = group, n_perm = 99L, seed = 8L, include_observed = TRUE)
  )
  without_obs <- fmrigds:::core_perm_twosample_kernel(
    beta,
    opts = list(group = group, n_perm = 99L, seed = 8L, include_observed = FALSE)
  )
  expected_t <- apply(beta, 2L, function(x) {
    unname(stats::t.test(x[group == 1L], x[group == 0L])$statistic)
  })
  expected_beta <- colMeans(beta[group == 1L, ]) - colMeans(beta[group == 0L, ])

  expect_equal(with_obs$t_g, expected_t, tolerance = 1e-10)
  expect_equal(without_obs$t_g, expected_t, tolerance = 1e-10)
  expect_equal(with_obs$beta_g, unname(expected_beta), tolerance = 1e-10)
  expect_equal(without_obs$beta_g, unname(expected_beta), tolerance = 1e-10)
})

test_that("perm p-values use (1 + #null >= obs) / (1 + n_null) without double counting", {
  set.seed(22)
  n_subject <- 8L
  beta <- matrix(stats::rnorm(n_subject * 2L, mean = 0.3), n_subject, 2L)
  signs <- matrix(sample(c(-1L, 1L), 20L * n_subject, replace = TRUE), 20L, n_subject)
  signs[1L, ] <- 1L # identity relabelling supplied in row 1
  res <- fmrigds:::core_perm_onesample_kernel(beta, opts = list(signs = signs))

  t1 <- function(y) mean(y) / (stats::sd(y) / sqrt(length(y)))
  for (b in seq_len(ncol(beta))) {
    tobs <- t1(beta[, b])
    tnull <- apply(signs[-1L, , drop = FALSE], 1L, function(s) t1(s * beta[, b]))
    expected <- (1 + sum(abs(tnull) >= abs(tobs))) / (1 + length(tnull))
    expect_equal(res$p_perm[[b]], expected)
  }

  # Same convention for the two-sample kernel.
  group <- rep(c(0L, 1L), each = 4L)
  gmat <- t(replicate(15L, sample(group)))
  gmat[1L, ] <- group
  res2 <- fmrigds:::core_perm_twosample_kernel(beta, opts = list(group = group, group_mat = gmat))
  t2 <- function(y, g) unname(stats::t.test(y[g == 1L], y[g == 0L])$statistic)
  tobs <- t2(beta[, 1L], group)
  tnull <- apply(gmat[-1L, , drop = FALSE], 1L, function(g) t2(beta[, 1L], g))
  expect_equal(res2$p_perm[[1L]], (1 + sum(abs(tnull) >= abs(tobs) - 1e-12)) / (1 + length(tnull)))
})

# 3. OLS never reports exactly-determined fits --------------------------------

test_that("ols:voxelwise returns NA when only p finite subjects remain", {
  set.seed(31)
  n_subject <- 6L
  X <- cbind(1, age = stats::rnorm(n_subject))
  beta <- matrix(stats::rnorm(n_subject * 2L), n_subject, 2L)
  beta[3:6, 2L] <- NA_real_ # only 2 finite subjects for p = 2

  res <- suppressWarnings(fmrigds:::ols_voxelwise_cpp(beta, X, min_obs = 2L))
  expect_true(all(is.finite(res$coef[, 1L])))
  expect_equal(res$df_res[[1L]], n_subject - 2)
  expect_true(all(is.na(res$coef[, 2L])))
  expect_true(all(is.na(res$se_coef[, 2L])))
  expect_true(is.na(res$df_res[[2L]]))

  reducer <- fmrigds:::get_reducer("ols:voxelwise")
  out <- suppressWarnings(reducer$fun(
    beta = beta, var = NULL, X = X, z = NULL, p = NULL, df = NULL,
    df1 = NULL, df2 = NULL, opts = list(min_subjects = 2L)
  ))
  expect_true(all(is.na(out$p_coef[, 2L])))
})

# 4. z from t / p keeps precision in the far tail -----------------------------

test_that("derive_z from t and df is finite in the far tail", {
  arrays <- list(t = array(c(12, -12, 0.5), c(3, 1, 1)), df = array(100, c(3, 1, 1)))
  z <- fmrigds:::derive_z(arrays)
  expected <- -stats::qnorm(stats::pt(-12, 100, log.p = TRUE), log.p = TRUE)
  expect_true(all(is.finite(z)))
  expect_equal(z[1, 1, 1], expected, tolerance = 1e-10)
  expect_equal(z[2, 1, 1], -expected, tolerance = 1e-10)
  expect_gt(z[1, 1, 1], 9)
  expect_lt(z[1, 1, 1], 10)
  expect_equal(z[3, 1, 1], stats::qnorm(stats::pt(0.5, 100)), tolerance = 1e-10)

  zp <- fmrigds:::derive_z(list(p = array(1e-20, c(1, 1, 1)), beta = array(-1, c(1, 1, 1))))
  expect_equal(zp[1, 1, 1], -stats::qnorm(5e-21, lower.tail = FALSE))
})

test_that("combine:stouffer keeps far-tail subjects derived from t/df", {
  t_arr <- array(c(12, 1, 1), c(1, 3, 1))
  df_arr <- array(100, c(1, 3, 1))
  g <- new_gds(
    list(t = t_arr, df = df_arr),
    space_sample_labels("a"), c("s1", "s2", "s3"), "c1"
  )
  out <- as_plan(g) |> reduce("combine:stouffer") |> compute()
  z_ind <- -stats::qnorm(stats::pt(-abs(c(12, 1, 1)), 100, log.p = TRUE), log.p = TRUE)
  expect_equal(as.numeric(assay(out, "z_g")[1, 1, 1]), sum(z_ind) / sqrt(3), tolerance = 1e-8)
})

# 5. Lancaster with tiny p-values ----------------------------------------------

test_that("Lancaster combiner stays finite for tiny p-values", {
  p <- matrix(c(1e-17, 0.5), nrow = 2L, ncol = 1L)
  res <- fmrigds:::core_lancaster_kernel(p = p, dfw = c(1L, 1L))
  expected <- stats::qchisq(1e-17, 2, lower.tail = FALSE) + stats::qchisq(0.5, 2, lower.tail = FALSE)
  expect_true(is.finite(res$chi2[[1L]]))
  expect_equal(res$chi2[[1L]], expected, tolerance = 1e-10)

  fb <- fmrigds:::.lancaster_fallback(p, dfw = c(1L, 1L))
  expect_equal(fb$chi2[[1L]], expected, tolerance = 1e-10)
})

# 6. Covariance propagation uses the per-subject variances ---------------------

test_that("propagate_variance_covariance scales correlation by each subject's variance", {
  M <- matrix(c(0.5, 0.5), nrow = 1)
  beta <- array(c(1, 2, 3, 4), c(2, 2, 1))
  var <- array(c(1, 1, 4, 9), c(2, 2, 1)) # subject 2 has larger variances
  R <- matrix(c(1, 0.3, 0.3, 1), 2)
  res <- fmrigds:::propagate_variance_covariance(M, beta, var, function(idx) R[idx, idx, drop = FALSE])

  expected <- vapply(1:2, function(j) {
    D <- diag(sqrt(var[, j, 1]))
    as.numeric(M %*% D %*% R %*% D %*% t(M))
  }, numeric(1))
  expect_equal(res$var[1, , 1], expected)
  expect_false(isTRUE(all.equal(res$var[1, 1, 1], res$var[1, 2, 1])))
  expect_equal(res$beta[1, , 1], c(1.5, 3.5))
})

# 7. OpenMP thread defaults ------------------------------------------------------

test_that("thread count defaults to at most 2 and honours option and limit", {
  old_opt <- options(fmrigds.threads = NULL)
  old_env <- Sys.getenv(c("FMRIGDS_THREADS", "OMP_THREAD_LIMIT"), unset = NA)
  on.exit({
    options(old_opt)
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(old_env[nm]))
    }
    fmrigds:::.set_threads_from_option()
  }, add = TRUE)
  Sys.unsetenv(c("FMRIGDS_THREADS", "OMP_THREAD_LIMIT"))
  expect_equal(fmrigds:::.resolve_thread_count(), 2L)

  Sys.setenv(FMRIGDS_THREADS = "4")
  expect_equal(fmrigds:::.resolve_thread_count(), 4L)

  options(fmrigds.threads = 3L)
  expect_equal(fmrigds:::.resolve_thread_count(), 3L)
  expect_no_error(fmrigds:::.set_threads_from_option())

  Sys.setenv(OMP_THREAD_LIMIT = "1")
  expect_equal(fmrigds:::.resolve_thread_count(), 1L)
  expect_no_error(fmrigds:::.set_threads_from_option())
})

# 8. Smaller robustness fixes ---------------------------------------------------

test_that("lmm:ri lambda optimizer errors clearly when every candidate fails", {
  Y <- matrix(stats::rnorm(8), 8, 1)
  X <- cbind(1, rep(1, 8)) # singular fixed-effects design
  expect_error(
    suppressWarnings(fmrigds:::.optimize_lmm_ri_theta(Y, X, n_repeat = 2L, fit = "REML")),
    "failed to optimize the random-intercept variance ratio"
  )
})

test_that("lmm:ri reports non-convergence when lambda is pinned at the search bound", {
  set.seed(82)
  subjects <- sprintf("s%02d", 1:10)
  contrasts <- c("c1", "c2", "c3")
  u <- stats::rnorm(10, sd = 5)
  beta <- array(NA_real_, c(2, 10, 3))
  beta[1, , ] <- u + matrix(stats::rnorm(30, sd = 1e-9), 10, 3) # no within-subject noise
  beta[2, , ] <- u + matrix(stats::rnorm(30, sd = 1), 10, 3)
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))), space_sample_labels(c("a", "b")), subjects, contrasts)
  out <- reduce(as_plan(g), method = "lmm:ri", formula = ~ 1,
                options = list(theta_mode = "voxelwise")) |> compute()
  expect_equal(as.numeric(assay(out, "converged")[, 1, 1]), c(0, 1))
})

test_that("lmm:ri_knownvar marks samples with non-positive variance as NA", {
  set.seed(81)
  subjects <- paste0("s", 1:10)
  contrasts <- c("c0", "c1")
  beta <- array(stats::rnorm(2 * 10 * 2), c(2, 10, 2))
  var <- array(0.05, dim(beta))
  var[1, 1, 1] <- 0
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("a", "b")), subjects, contrasts)
  g <- with_contrast_data(g, data.frame(time = c(0, 1), row.names = contrasts))
  out <- NULL
  expect_warning(
    out <- reduce(as_plan(g), method = "lmm:ri_knownvar", formula = ~ time) |> compute(),
    "non-positive sampling variances"
  )
  coef_time <- assay(out, "coef:time")[, 1, 1]
  expect_true(is.na(coef_time[[1L]]))
  expect_true(is.finite(coef_time[[2L]]))
})

test_that("meta:fe_reg honours min_subjects", {
  beta <- matrix(c(1, 2, 3, NA), 4, 1)
  var <- matrix(1, 4, 1)
  X <- matrix(1, 4, 1)
  reducer <- fmrigds:::get_reducer("meta:fe_reg")
  call_fe_reg <- function(opts) reducer$fun(
    beta = beta, var = var, X = X, z = NULL, p = NULL, df = NULL,
    df1 = NULL, df2 = NULL, opts = opts
  )
  expect_true(is.finite(call_fe_reg(list(min_subjects = 2L))$coef[1, 1]))
  expect_true(is.na(call_fe_reg(list(min_subjects = 4L))$coef[1, 1]))
})

test_that("validate_reducer_options accepts vector values and rejects invalid ones", {
  vo <- fmrigds:::validate_reducer_options
  schema <- list(alternative = c("two.sided", "less", "greater"))
  expect_equal(vo(schema, list(alternative = c("less", "greater")))$alternative, c("less", "greater"))
  expect_error(vo(schema, list(alternative = c("less", "bad"))), "Invalid option 'alternative'")
})

test_that("meta_fe_reg_cpp Q skips subjects with non-finite design rows", {
  beta <- matrix(c(1, 2, 3, 4), 4, 1)
  var <- matrix(1, 4, 1)
  X <- cbind(1, c(0, 1, 2, NA))
  res <- fmrigds:::meta_fe_reg_cpp(beta, var, X, min_subj = 2L)
  expect_true(is.finite(res$Q[[1L]]))
  expect_equal(res$df_res[[1L]], 1)
})

# 9. LMM between/within degrees of freedom -------------------------------------

test_that("lmm:ri uses between-subject df for between-subject terms", {
  set.seed(91)
  n_subject <- 20L
  n_rep <- 4L
  subjects <- sprintf("s%02d", seq_len(n_subject))
  contrasts <- paste0("c", seq_len(n_rep))
  grp <- rep(c("A", "B"), each = n_subject / 2)
  cond <- c(0, 1, 2, 3)
  u <- stats::rnorm(n_subject, sd = 0.7)
  beta <- array(NA_real_, c(2, n_subject, n_rep))
  for (s in seq_len(n_subject)) {
    mu <- 0.5 * (grp[s] == "B") + 0.3 * cond + u[s]
    beta[1, s, ] <- mu + stats::rnorm(n_rep, sd = 0.3)
    beta[2, s, ] <- mu + stats::rnorm(n_rep, sd = 0.3)
  }
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_sample_labels(c("a", "b")), subjects, contrasts,
    col_data = data.frame(group = grp, row.names = subjects)
  )
  g <- with_contrast_data(g, data.frame(condition = cond, row.names = contrasts))
  out <- reduce(as_plan(g), method = "lmm:ri", formula = ~ group * condition) |> compute()

  df_group <- as.numeric(assay(out, "df_coef:groupB")[, 1, 1])
  df_int <- as.numeric(assay(out, "df_coef:(Intercept)")[, 1, 1])
  df_cond <- as.numeric(assay(out, "df_coef:condition")[, 1, 1])
  df_inter <- as.numeric(assay(out, "df_coef:groupB:condition")[, 1, 1])
  expect_equal(df_group, c(18, 18))
  expect_equal(df_int, c(18, 18))
  expect_equal(df_cond, c(76, 76))
  expect_equal(df_inter, c(76, 76))
  expect_equal(as.numeric(assay(out, "df_res")[, 1, 1]), c(76, 76))

  t_group <- as.numeric(assay(out, "t_coef:groupB")[, 1, 1])
  expect_equal(
    as.numeric(assay(out, "p_coef:groupB")[, 1, 1]),
    2 * stats::pt(-abs(t_group), df = 18)
  )
})

test_that("map_to Fisher combine stays finite for far-tail z", {
  target <- space_sample_labels("t1")
  from_sp <- space_sample_labels(c("a", "b"))
  mp <- map_linear(from_sp, target, matrix(c(0.5, 0.5), nrow = 1))
  z <- array(c(40, 39), dim = c(2, 1, 1))
  out <- apply_map_to(
    list(op = "map", target_space = target, map = mp,
         uncertainty = UncertaintyRule("none"), combine = "fisher"),
    list(z = z)
  )$arrays
  log_p <- log(2) + pnorm(-abs(c(40, 39)), log.p = TRUE)
  expect_true(is.finite(out$chi2[1, 1, 1]))
  expect_equal(out$chi2[1, 1, 1], -2 * sum(log_p))
  expect_true(is.finite(out$z[1, 1, 1]))
})
