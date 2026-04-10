make_lmm_fixture <- function() {
  subjects <- paste0("sub-", sprintf("%02d", 1:8))
  contrasts <- c("baseline", "task")
  u <- c(-0.5, -0.2, 0.0, 0.1, 0.25, 0.35, -0.15, 0.2)
  eps1 <- matrix(c(
    -0.10, 0.02,
    -0.05, 0.01,
    0.03, -0.02,
    0.04, 0.00,
    -0.02, 0.03,
    0.01, -0.01,
    0.02, 0.04,
    -0.03, -0.02
  ), ncol = 2, byrow = TRUE)
  eps2 <- matrix(c(
    0.02, -0.01,
    -0.04, 0.03,
    0.01, -0.02,
    -0.03, 0.02,
    0.00, 0.01,
    0.03, -0.03,
    -0.01, 0.02,
    0.02, 0.00
  ), ncol = 2, byrow = TRUE)

  beta <- array(NA_real_, dim = c(2, length(subjects), length(contrasts)))
  beta[1, , ] <- cbind(1.2 + u, 1.2 + 0.7 + u) + eps1
  beta[2, , ] <- cbind(-0.4 + u, -0.4 + 1.1 + u) + eps2
  var <- array(0.05, dim = dim(beta))

  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("ROI_1", "ROI_2")),
    subjects = subjects,
    contrasts = contrasts,
    col_data = data.frame(
      cohort = rep(c("A", "B"), each = 4),
      row.names = subjects,
      stringsAsFactors = FALSE
    )
  )
  g <- with_contrast_data(
    g,
    data.frame(
      condition = c(0, 1),
      row.names = contrasts,
      stringsAsFactors = FALSE
    )
  )

  list(gds = g, beta = beta, subjects = subjects, contrasts = contrasts)
}

permute_subject_order <- function(g, perm) {
  g$assays <- lapply(g$assays, function(x) x[, perm, , drop = FALSE])
  g$subjects <- g$subjects[perm]
  g
}

test_that("lmm:ri returns joint repeated-measures estimates with shared theta", {
  fixture <- make_lmm_fixture()

  gout <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri",
    formula = ~ condition,
    options = list(fit = "REML", theta_mode = "pooled")
  ) |> compute()

  expect_equal(subjects(gout), "meta")
  expect_equal(contrasts(gout), "model")
  expect_equal(rownames(contrast_data(gout)), "model")
  expect_equal(contrast_data(gout)$label, "model")

  assay_names <- names(assays(gout))
  expect_true(all(c(
    "coef:(Intercept)",
    "coef:condition",
    "se_coef:(Intercept)",
    "se_coef:condition",
    "t_coef:(Intercept)",
    "t_coef:condition",
    "p_coef:(Intercept)",
    "p_coef:condition",
    "lambda"
  ) %in% assay_names))

  condition_beta <- assay(gout, "coef:condition")[, 1, 1]
  expect_true(all(is.finite(condition_beta)))
  expect_true(all(condition_beta > 0.5))

  lambda <- assay(gout, "lambda")[, 1, 1]
  expect_equal(lambda[1], lambda[2], tolerance = 1e-10)
})

test_that("lmm:ri rejects lmer-style random-effect syntax", {
  fixture <- make_lmm_fixture()
  plan <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri",
    formula = ~ condition + (1 | subject)
  )

  expect_error(compute(plan), "do not support lmer-style random-effects syntax")
})

test_that("lmm:ri is numerically close to lme4 on the supported random-intercept subset", {
  skip_if_not_installed("lme4")

  fixture <- make_lmm_fixture()
  beta1 <- fixture$beta[1, , ]
  long_df <- data.frame(
    y = as.numeric(t(beta1)),
    subject = factor(rep(fixture$subjects, each = length(fixture$contrasts)), levels = fixture$subjects),
    condition = rep(c(0, 1), times = length(fixture$subjects))
  )

  ref <- lme4::lmer(y ~ condition + (1 | subject), data = long_df, REML = TRUE)
  gout <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri",
    formula = ~ condition,
    options = list(fit = "REML", theta_mode = "pooled")
  ) |> compute()

  fixef_ref <- lme4::fixef(ref)
  expect_equal(as.numeric(assay(gout, "coef:(Intercept)")[1, 1, 1]), unname(fixef_ref[["(Intercept)"]]), tolerance = 1e-4)
  expect_equal(as.numeric(assay(gout, "coef:condition")[1, 1, 1]), unname(fixef_ref[["condition"]]), tolerance = 1e-4)
})

test_that("lmm:ri supports voxelwise theta mode with sample-specific lambdas", {
  fixture <- make_lmm_fixture()
  g <- fixture$gds
  beta <- assay(g, "beta")
  beta[2, , ] <- cbind(-0.4 + c(-0.03, -0.02, 0.00, 0.01, 0.01, 0.02, -0.01, 0.02),
                       -0.4 + 0.7 + c(-0.03, -0.02, 0.00, 0.01, 0.01, 0.02, -0.01, 0.02))
  g$assays$beta <- beta

  gout <- reduce(
    as_plan(g),
    method = "lmm:ri",
    formula = ~ condition,
    options = list(fit = "REML", theta_mode = "voxelwise")
  ) |> compute()

  lambda <- assay(gout, "lambda")[, 1, 1]
  expect_true(all(is.finite(lambda)))
  expect_gt(abs(diff(lambda)), 1e-3)
})

test_that("lmm:ri voxelwise theta collapses to pooled theta for identical responses", {
  fixture <- make_lmm_fixture()
  g <- fixture$gds
  beta <- assay(g, "beta")
  beta[2, , ] <- beta[1, , ]
  g$assays$beta <- beta

  pooled <- reduce(
    as_plan(g),
    method = "lmm:ri",
    formula = ~ condition,
    options = list(fit = "REML", theta_mode = "pooled")
  ) |> compute()

  voxelwise <- reduce(
    as_plan(g),
    method = "lmm:ri",
    formula = ~ condition,
    options = list(fit = "REML", theta_mode = "voxelwise")
  ) |> compute()

  expect_equal(assay(voxelwise, "lambda"), assay(pooled, "lambda"), tolerance = 1e-4)
  expect_equal(assay(voxelwise, "coef:(Intercept)"), assay(pooled, "coef:(Intercept)"), tolerance = 1e-6)
  expect_equal(assay(voxelwise, "coef:condition"), assay(pooled, "coef:condition"), tolerance = 1e-6)
  expect_equal(assay(voxelwise, "vc_intercept"), assay(pooled, "vc_intercept"), tolerance = 1e-4)
  expect_lt(
    max(abs(as.numeric(assay(voxelwise, "sigma2")) - as.numeric(assay(pooled, "sigma2")))),
    1e-6
  )
})

test_that("lmm:ri is invariant to subject permutation and col_data row order", {
  fixture <- make_lmm_fixture()
  perm <- c(6, 2, 8, 1, 5, 3, 7, 4)

  g_perm <- permute_subject_order(fixture$gds, perm)
  g_perm <- with_col_data(
    g_perm,
    fixture$gds$col_data[rev(fixture$subjects), , drop = FALSE]
  )

  base <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri",
    formula = ~ condition + cohort,
    options = list(fit = "REML", theta_mode = "pooled")
  ) |> compute()

  shuffled <- reduce(
    as_plan(g_perm),
    method = "lmm:ri",
    formula = ~ condition + cohort,
    options = list(fit = "REML", theta_mode = "pooled")
  ) |> compute()

  for (nm in names(assays(base))) {
    expect_equal(assay(shuffled, nm), assay(base, nm), tolerance = 1e-8, info = nm)
  }
})
