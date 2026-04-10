make_lmm_slope_fixture <- function() {
  subjects <- paste0("sub-", sprintf("%02d", 1:10))
  contrasts <- c("t0", "t1", "t2")
  time <- c(-1, 0, 1)
  b0 <- c(-0.30, -0.18, -0.12, -0.05, 0.00, 0.08, 0.13, 0.18, 0.24, 0.31)
  b1 <- c(-0.15, -0.10, -0.06, -0.02, 0.01, 0.05, 0.08, 0.11, 0.15, 0.19)
  eps1 <- matrix(c(
    0.01, -0.02, 0.00,
    -0.01, 0.02, 0.01,
    0.02, 0.00, -0.01,
    0.00, -0.01, 0.02,
    -0.02, 0.01, 0.00,
    0.01, 0.00, -0.02,
    0.00, 0.02, -0.01,
    -0.01, 0.01, 0.02,
    0.02, -0.01, 0.01,
    0.00, 0.01, -0.01
  ), ncol = 3, byrow = TRUE)
  eps2 <- matrix(c(
    -0.01, 0.00, 0.02,
    0.02, -0.01, 0.01,
    0.00, 0.01, -0.02,
    -0.02, 0.02, 0.00,
    0.01, -0.02, 0.01,
    0.00, 0.01, -0.01,
    -0.01, 0.00, 0.02,
    0.02, -0.02, 0.00,
    0.01, 0.01, -0.02,
    -0.02, 0.00, 0.01
  ), ncol = 3, byrow = TRUE)

  beta <- array(NA_real_, dim = c(2, length(subjects), length(contrasts)))
  beta[1, , ] <- outer(1.1 + b0, rep(1, length(time))) + outer(0.65 + b1, time) + eps1
  beta[2, , ] <- outer(-0.3 + 0.8 * b0, rep(1, length(time))) + outer(0.95 + 1.3 * b1, time) + eps2
  var <- array(0.04, dim = dim(beta))

  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("ROI_1", "ROI_2")),
    subjects = subjects,
    contrasts = contrasts
  )
  g <- with_contrast_data(
    g,
    data.frame(
      time = time,
      row.names = contrasts,
      stringsAsFactors = FALSE
    )
  )

  list(gds = g, beta = beta, subjects = subjects, contrasts = contrasts, time = time)
}

permute_contrast_order <- function(g, perm) {
  g$assays <- lapply(g$assays, function(x) x[, , perm, drop = FALSE])
  g$contrasts <- g$contrasts[perm]
  g
}

test_that("lmm:ri_slope1 with diag covariance returns slope random-effects outputs", {
  fixture <- make_lmm_slope_fixture()

  gout <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  expect_equal(subjects(gout), "meta")
  expect_equal(contrasts(gout), "model")
  expect_true(all(c(
    "coef:(Intercept)",
    "coef:time",
    "vc_intercept",
    "vc_slope",
    "vc_cov_intercept_slope",
    "corr_intercept_slope",
    "lambda_intercept",
    "lambda_slope"
  ) %in% names(assays(gout))))

  expect_true(all(assay(gout, "coef:time")[, 1, 1] > 0.4))
  expect_true(all(assay(gout, "vc_slope")[, 1, 1] > 0))
  expect_equal(as.numeric(assay(gout, "corr_intercept_slope")[, 1, 1]), c(0, 0), tolerance = 1e-8)
})

test_that("lmm:ri_slope1 with full covariance emits non-zero random-effect correlation", {
  fixture <- make_lmm_slope_fixture()

  gout <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "full", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  corr <- assay(gout, "corr_intercept_slope")[, 1, 1]
  expect_true(all(is.finite(corr)))
  expect_true(any(abs(corr) > 0.05))
})

test_that("lmm:ri_slope1 requires a numeric slope column in contrast_data", {
  fixture <- make_lmm_slope_fixture()
  g_bad <- with_contrast_data(
    fixture$gds,
    data.frame(
      time = c("pre", "mid", "post"),
      row.names = fixture$contrasts,
      stringsAsFactors = FALSE
    )
  )

  plan_missing <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(covariance = "diag")
  )
  expect_error(compute(plan_missing), "options\\$slope")

  plan_bad <- reduce(
    as_plan(g_bad),
    method = "lmm:ri_slope1",
    formula = ~ 1,
    options = list(slope = "time", covariance = "diag")
  )
  expect_error(compute(plan_bad), "must be numeric")
})

test_that("lmm:ri_slope1 is numerically close to lme4 for the supported full model subset", {
  skip_if_not_installed("lme4")

  fixture <- make_lmm_slope_fixture()
  beta1 <- fixture$beta[1, , ]
  long_df <- data.frame(
    y = as.numeric(t(beta1)),
    subject = factor(rep(fixture$subjects, each = length(fixture$contrasts)), levels = fixture$subjects),
    time = rep(fixture$time, times = length(fixture$subjects))
  )

  ref <- lme4::lmer(y ~ time + (1 + time | subject), data = long_df, REML = TRUE)
  gout <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "full", fit = "REML", theta_mode = "pooled", center_slope = FALSE)
  ) |> compute()

  fixef_ref <- lme4::fixef(ref)
  expect_equal(as.numeric(assay(gout, "coef:(Intercept)")[1, 1, 1]), unname(fixef_ref[["(Intercept)"]]), tolerance = 1e-4)
  expect_equal(as.numeric(assay(gout, "coef:time")[1, 1, 1]), unname(fixef_ref[["time"]]), tolerance = 1e-4)
})

test_that("lmm:ri_slope1 works from a raw long table with extracted contrast_data", {
  fixture <- make_lmm_slope_fixture()
  dims <- dim(fixture$beta)
  df <- expand.grid(
    sample = c("ROI_1", "ROI_2"),
    subject = fixture$subjects,
    visit = fixture$contrasts,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df$time <- fixture$time[match(df$visit, fixture$contrasts)]

  sample_idx <- match(df$sample, c("ROI_1", "ROI_2"))
  subject_idx <- match(df$subject, fixture$subjects)
  contrast_idx <- match(df$visit, fixture$contrasts)
  df$beta <- fixture$beta[cbind(sample_idx, subject_idx, contrast_idx)]
  df$var <- 0.04

  gout <- gds(df, contrast_col = "visit", contrast_data_cols = "time") |>
    reduce(
      method = "lmm:ri_slope1",
      formula = ~ time,
      options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "pooled")
    ) |>
    compute()

  expect_equal(subjects(gout), "meta")
  expect_equal(contrasts(gout), "model")
  expect_true(all(assay(gout, "coef:time")[, 1, 1] > 0.4))
})

test_that("lmm:ri_slope1 supports voxelwise theta mode with sample-specific slope variance", {
  fixture <- make_lmm_slope_fixture()

  gout <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "voxelwise")
  ) |> compute()

  lambda_slope <- assay(gout, "lambda_slope")[, 1, 1]
  expect_true(all(is.finite(lambda_slope)))
  expect_gt(abs(diff(lambda_slope)), 1e-3)
})

test_that("lmm:ri_slope1 voxelwise theta collapses to pooled theta for identical responses", {
  fixture <- make_lmm_slope_fixture()
  g <- fixture$gds
  beta <- assay(g, "beta")
  beta[2, , ] <- beta[1, , ]
  g$assays$beta <- beta

  pooled <- reduce(
    as_plan(g),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  voxelwise <- reduce(
    as_plan(g),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "voxelwise")
  ) |> compute()

  expect_equal(assay(voxelwise, "lambda_intercept"), assay(pooled, "lambda_intercept"), tolerance = 1e-6)
  expect_equal(assay(voxelwise, "lambda_slope"), assay(pooled, "lambda_slope"), tolerance = 1e-6)
  expect_equal(assay(voxelwise, "coef:(Intercept)"), assay(pooled, "coef:(Intercept)"), tolerance = 1e-6)
  expect_equal(assay(voxelwise, "coef:time"), assay(pooled, "coef:time"), tolerance = 1e-6)
  expect_equal(assay(voxelwise, "vc_slope"), assay(pooled, "vc_slope"), tolerance = 1e-6)
})

test_that("lmm:ri_slope1 is invariant to contrast permutation when contrast_data are aligned by name", {
  fixture <- make_lmm_slope_fixture()
  perm <- c(3, 1, 2)

  g_perm <- permute_contrast_order(fixture$gds, perm)
  g_perm <- with_contrast_data(g_perm, contrast_data(fixture$gds))

  base <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "full", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  shuffled <- reduce(
    as_plan(g_perm),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "full", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  for (nm in names(assays(base))) {
    expect_equal(assay(shuffled, nm), assay(base, nm), tolerance = 1e-5, info = nm)
  }
})

test_that("lmm:ri_slope1 obeys the response scaling law", {
  fixture <- make_lmm_slope_fixture()
  scale_factor <- -2.5
  g_scaled <- fixture$gds
  g_scaled$assays$beta <- scale_factor * assay(g_scaled, "beta")

  base <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  scaled <- reduce(
    as_plan(g_scaled),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "pooled")
  ) |> compute()

  expect_equal(assay(scaled, "coef:(Intercept)"), scale_factor * assay(base, "coef:(Intercept)"), tolerance = 1e-6)
  expect_equal(assay(scaled, "coef:time"), scale_factor * assay(base, "coef:time"), tolerance = 1e-6)
  expect_equal(assay(scaled, "se_coef:(Intercept)"), abs(scale_factor) * assay(base, "se_coef:(Intercept)"), tolerance = 1e-6)
  expect_equal(assay(scaled, "se_coef:time"), abs(scale_factor) * assay(base, "se_coef:time"), tolerance = 1e-6)
  expect_equal(assay(scaled, "t_coef:(Intercept)"), sign(scale_factor) * assay(base, "t_coef:(Intercept)"), tolerance = 1e-6)
  expect_equal(assay(scaled, "t_coef:time"), sign(scale_factor) * assay(base, "t_coef:time"), tolerance = 1e-6)
  expect_lt(
    max(abs(as.numeric(assay(scaled, "p_coef:(Intercept)")) - as.numeric(assay(base, "p_coef:(Intercept)")))),
    1e-8
  )
  expect_lt(
    max(abs(as.numeric(assay(scaled, "p_coef:time")) - as.numeric(assay(base, "p_coef:time")))),
    1e-8
  )
  expect_equal(assay(scaled, "sigma2"), scale_factor^2 * assay(base, "sigma2"), tolerance = 1e-6)
  expect_equal(assay(scaled, "vc_intercept"), scale_factor^2 * assay(base, "vc_intercept"), tolerance = 1e-6)
  expect_equal(assay(scaled, "vc_slope"), scale_factor^2 * assay(base, "vc_slope"), tolerance = 1e-6)
  expect_equal(assay(scaled, "vc_cov_intercept_slope"), scale_factor^2 * assay(base, "vc_cov_intercept_slope"), tolerance = 1e-6)
  expect_equal(assay(scaled, "vc_resid"), scale_factor^2 * assay(base, "vc_resid"), tolerance = 1e-6)
  expect_equal(assay(scaled, "lambda_intercept"), assay(base, "lambda_intercept"), tolerance = 1e-4)
  expect_equal(assay(scaled, "lambda_slope"), assay(base, "lambda_slope"), tolerance = 1e-4)
  expect_equal(assay(scaled, "lambda_cov_intercept_slope"), assay(base, "lambda_cov_intercept_slope"), tolerance = 1e-4)
  expect_equal(assay(scaled, "corr_intercept_slope"), assay(base, "corr_intercept_slope"), tolerance = 1e-6)
})

test_that("lmm:ri_slope1 stays finite at the near-zero random-slope boundary", {
  fixture <- make_lmm_slope_fixture()
  g <- fixture$gds
  beta <- assay(g, "beta")
  time <- fixture$time
  subj_intercept <- seq(-0.25, 0.25, length.out = length(fixture$subjects))
  eps_boundary <- matrix(c(
    -0.002, 0.001, 0.002,
    0.001, -0.001, 0.000,
    -0.001, 0.000, 0.001,
    0.002, -0.001, -0.001,
    0.000, 0.001, -0.001,
    -0.001, 0.002, -0.001,
    0.001, -0.002, 0.001,
    0.000, 0.001, -0.001,
    -0.001, 0.000, 0.002,
    0.001, -0.001, 0.000
  ), ncol = length(time), byrow = TRUE)
  beta[1, , ] <- outer(1.0 + subj_intercept, rep(1, length(time))) +
    outer(rep(0.7, length(fixture$subjects)), time) +
    eps_boundary
  g$assays$beta <- beta

  gout <- reduce(
    as_plan(g),
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(slope = "time", covariance = "diag", fit = "REML", theta_mode = "voxelwise")
  ) |> compute()

  lambda_slope <- assay(gout, "lambda_slope")[, 1, 1]
  converged <- assay(gout, "converged")[, 1, 1]
  expect_true(all(is.finite(lambda_slope)))
  expect_true(all(converged == 1))
  expect_lt(lambda_slope[1], 0.1)
  expect_lt(lambda_slope[1], lambda_slope[2] / 5)
})
