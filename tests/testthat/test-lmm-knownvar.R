make_lmm_knownvar_fixture <- function(kind = c("ri", "ri_slope1"), n_subject = 36L) {
  kind <- match.arg(kind)
  set.seed(if (identical(kind, "ri")) 741L else 913L)
  subjects <- paste0("sub-", sprintf("%03d", seq_len(n_subject)))
  contrasts <- c("t0", "t1", "t2")
  time <- c(-1, 0, 1)
  n_obs <- n_subject * length(time)
  z0 <- rnorm(n_subject)
  b0 <- 0.32 * z0
  b1 <- if (identical(kind, "ri_slope1")) {
    0.16 * (0.35 * z0 + sqrt(1 - 0.35^2) * rnorm(n_subject))
  } else {
    rep(0, n_subject)
  }
  sampling_var <- matrix(
    runif(n_obs, min = 0.008, max = 0.055),
    nrow = n_subject,
    ncol = length(time)
  )
  residual_var <- if (identical(kind, "ri")) 0.045 else 0.035
  noise <- matrix(
    rnorm(n_obs, sd = sqrt(as.numeric(t(sampling_var)) + residual_var)),
    nrow = n_subject,
    ncol = length(time),
    byrow = TRUE
  )
  response <- outer(1.1 + b0, rep(1, length(time))) +
    outer(0.55 + b1, time) + noise

  beta <- array(response, dim = c(1L, n_subject, length(time)))
  var <- array(sampling_var, dim = dim(beta))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels("ROI_1"),
    subjects = subjects,
    contrasts = contrasts
  )
  g <- with_contrast_data(
    g,
    data.frame(time = time, row.names = contrasts, stringsAsFactors = FALSE)
  )
  list(
    gds = g,
    subjects = subjects,
    contrasts = contrasts,
    time = time,
    y = as.numeric(t(response)),
    vi = as.numeric(t(sampling_var)),
    residual_var = residual_var
  )
}

knownvar_long_data <- function(fixture) {
  data.frame(
    yi = fixture$y,
    vi = fixture$vi,
    time = rep(fixture$time, times = length(fixture$subjects)),
    subject = factor(
      rep(fixture$subjects, each = length(fixture$time)),
      levels = fixture$subjects
    ),
    obs = factor(seq_along(fixture$y))
  )
}

test_that("known-variance low-rank kernel matches dense GLS algebra", {
  n_subject <- 4L
  time <- c(-1, 0, 1)
  n_repeat <- length(time)
  X <- cbind(Intercept = 1, time = rep(time, times = n_subject))
  W <- cbind(Intercept = 1, time = time)
  G <- matrix(c(0.12, 0.018, 0.018, 0.035), nrow = 2L)
  L <- t(chol(G))
  U <- W %*% L
  sampling_var <- seq(0.015, 0.07, length.out = nrow(X))
  residual_var <- 0.025
  y <- c(0.1, 0.8, 1.3, -0.2, 0.4, 1.0, 0.3, 0.7, 1.5, -0.1, 0.5, 0.9)

  Z <- kronecker(diag(n_subject), W)
  G_big <- kronecker(diag(n_subject), G)
  V <- diag(sampling_var + residual_var) + Z %*% G_big %*% t(Z)
  Vinv <- solve(V)
  A <- crossprod(X, Vinv %*% X)
  Ainv <- solve(A)
  beta_ref <- Ainv %*% crossprod(X, Vinv %*% y)
  se_ref <- sqrt(diag(Ainv))
  resid_ref <- as.numeric(y - X %*% beta_ref)
  objective_ref <- as.numeric(
    determinant(V, logarithm = TRUE)$modulus +
      determinant(A, logarithm = TRUE)$modulus +
      crossprod(resid_ref, Vinv %*% resid_ref)
  )

  fit <- lmm_knownvar_fit_cpp(
    y,
    sampling_var,
    X,
    U_block = U,
    residual_var = residual_var,
    fit = "REML"
  )
  objective <- lmm_knownvar_objective_cpp(
    y,
    sampling_var,
    X,
    U_block = U,
    residual_var = residual_var,
    fit = "REML"
  )
  objective_ml <- lmm_knownvar_objective_cpp(
    y,
    sampling_var,
    X,
    U_block = U,
    residual_var = residual_var,
    fit = "ML"
  )
  objective_ml_ref <- as.numeric(
    determinant(V, logarithm = TRUE)$modulus +
      crossprod(resid_ref, Vinv %*% resid_ref)
  )

  expect_equal(as.numeric(fit$coef), as.numeric(beta_ref), tolerance = 1e-10)
  expect_equal(as.numeric(fit$se_coef), unname(se_ref), tolerance = 1e-10)
  expect_equal(objective, objective_ref, tolerance = 1e-10)
  expect_equal(objective_ml, objective_ml_ref, tolerance = 1e-10)
})

test_that("known-variance LMM reducers declare and enforce their input contract", {
  ri <- get_reducer("lmm:ri_knownvar")
  slope <- get_reducer("lmm:ri_slope1_knownvar")
  expect_equal(ri$requires, c("beta", "var"))
  expect_equal(slope$requires, c("beta", "var"))
  expect_equal(ri$model_contract$weight_mode, "model_specific")
  expect_equal(ri$model_contract$synthetic_variance, "forbid")

  fixture <- make_lmm_knownvar_fixture("ri")
  expect_error(
    reduce(
      as_plan(fixture$gds),
      method = "lmm:ri_knownvar",
      formula = ~ time,
      options = list(theta_mode = "pooled")
    ) |> compute(),
    "Invalid option 'theta_mode'"
  )
  expect_error(
    reduce(
      as_plan(fixture$gds),
      method = "lmm:ri_slope1_knownvar",
      formula = ~ time
    ) |> compute(),
    "options\\$slope.*lmm:ri_slope1_knownvar"
  )

  bad <- fixture$gds
  bad$assays$var[1, 1, 1] <- 0
  expect_error(
    reduce(as_plan(bad), method = "lmm:ri_knownvar", formula = ~ time) |> compute(),
    "strictly positive sampling variances"
  )

  arrays <- list(
    beta = assay(fixture$gds, "beta"),
    var = assay(fixture$gds, "var")
  )
  attr(arrays$var, "synthetic_unit_variance") <- TRUE
  node <- op_reduce(
    method = "lmm:ri_knownvar",
    weights = "1/var",
    by = "contrast",
    formula = ~ time
  )
  expect_error(
    apply_reduce(
      node,
      arrays,
      weights = "1/var",
      subjects = fixture$subjects,
      contrasts = fixture$contrasts,
      contrast_data = data.frame(time = fixture$time, row.names = fixture$contrasts)
    ),
    "synthetic"
  )
})

test_that("lmm:ri_knownvar agrees with metafor rma.mv", {
  skip_if_not_installed("metafor")
  fixture <- make_lmm_knownvar_fixture("ri")
  dat <- knownvar_long_data(fixture)
  reference <- metafor::rma.mv(
    yi,
    V = vi,
    mods = ~ time,
    random = list(~ 1 | subject, ~ 1 | obs),
    data = dat,
    method = "REML",
    control = list(sigma2.init = c(0.08, 0.03))
  )
  result <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_knownvar",
    formula = ~ time,
    options = list(fit = "REML")
  ) |> compute()

  expect_equal(
    c(assay(result, "coef:(Intercept)"), assay(result, "coef:time")),
    unname(stats::coef(reference)),
    tolerance = 2e-4
  )
  expect_equal(
    c(assay(result, "se_coef:(Intercept)"), assay(result, "se_coef:time")),
    reference$se,
    tolerance = 3e-4
  )
  expect_equal(
    as.numeric(assay(result, "vc_intercept")),
    reference$sigma2[[1L]],
    tolerance = 8e-4
  )
  expect_equal(
    as.numeric(assay(result, "vc_resid")),
    reference$sigma2[[2L]],
    tolerance = 8e-4
  )
})

test_that("known-variance reducers derive var from an aligned se assay", {
  fixture <- make_lmm_knownvar_fixture("ri", n_subject = 24L)
  g_se <- new_gds(
    assays = list(
      beta = assay(fixture$gds, "beta"),
      se = sqrt(assay(fixture$gds, "var"))
    ),
    space = fixture$gds$space,
    subjects = fixture$subjects,
    contrasts = fixture$contrasts
  )
  g_se <- with_contrast_data(
    g_se,
    data.frame(time = fixture$time, row.names = fixture$contrasts)
  )
  fit_one <- function(g) {
    reduce(as_plan(g), method = "lmm:ri_knownvar", formula = ~ time) |> compute()
  }
  from_var <- fit_one(fixture$gds)
  from_se <- fit_one(g_se)
  for (name in names(assays(from_var))) {
    expect_equal(assay(from_se, name), assay(from_var, name), tolerance = 1e-6, info = name)
  }
})

test_that("lmm:ri_slope1_knownvar agrees with metafor general random slopes", {
  skip_if_not_installed("metafor")
  fixture <- make_lmm_knownvar_fixture("ri_slope1")
  dat <- knownvar_long_data(fixture)
  reference <- metafor::rma.mv(
    yi,
    V = vi,
    mods = ~ time,
    random = list(~ time | subject, ~ 1 | obs),
    struct = "GEN",
    data = dat,
    method = "REML",
    control = list(
      sigma2.init = 0.03,
      tau2.init = c(0.08, 0.02),
      rho.init = 0.2
    )
  )
  result <- reduce(
    as_plan(fixture$gds),
    method = "lmm:ri_slope1_knownvar",
    formula = ~ time,
    options = list(slope = "time", covariance = "full", fit = "REML")
  ) |> compute()

  expect_equal(
    c(assay(result, "coef:(Intercept)"), assay(result, "coef:time")),
    unname(stats::coef(reference)),
    tolerance = 5e-4
  )
  expect_equal(
    c(assay(result, "se_coef:(Intercept)"), assay(result, "se_coef:time")),
    reference$se,
    tolerance = 8e-4
  )
  expect_equal(as.numeric(assay(result, "vc_intercept")), reference$tau2[[1L]], tolerance = 2e-3)
  expect_equal(as.numeric(assay(result, "vc_slope")), reference$tau2[[2L]], tolerance = 2e-3)
  expect_equal(as.numeric(assay(result, "corr_intercept_slope")), reference$rho[[1L]], tolerance = 1e-2)
  expect_equal(as.numeric(assay(result, "vc_resid")), reference$sigma2[[1L]], tolerance = 2e-3)
})

test_that("known-variance random-slope fits obey the joint effect-variance scaling law", {
  fixture <- make_lmm_knownvar_fixture("ri_slope1", n_subject = 24L)
  scale_factor <- -3.25
  scaled <- fixture$gds
  scaled$assays$beta <- scale_factor * assay(scaled, "beta")
  scaled$assays$var <- scale_factor^2 * assay(scaled, "var")

  fit_one <- function(g) {
    reduce(
      as_plan(g),
      method = "lmm:ri_slope1_knownvar",
      formula = ~ time,
      options = list(slope = "time", covariance = "full", fit = "REML")
    ) |> compute()
  }
  base <- fit_one(fixture$gds)
  transformed <- fit_one(scaled)

  for (term in c("(Intercept)", "time")) {
    expect_equal(
      assay(transformed, paste0("coef:", term)),
      scale_factor * assay(base, paste0("coef:", term)),
      tolerance = 2e-4
    )
    expect_equal(
      assay(transformed, paste0("se_coef:", term)),
      abs(scale_factor) * assay(base, paste0("se_coef:", term)),
      tolerance = 2e-4
    )
  }
  for (name in c("sigma2", "vc_intercept", "vc_slope", "vc_cov_intercept_slope", "vc_resid")) {
    expect_equal(assay(transformed, name), scale_factor^2 * assay(base, name), tolerance = 3e-3)
  }
  expect_equal(
    assay(transformed, "corr_intercept_slope"),
    assay(base, "corr_intercept_slope"),
    tolerance = 5e-4
  )
})

test_that("known-variance GLS responds to heterogeneous sampling precision", {
  fixture <- make_lmm_knownvar_fixture("ri", n_subject = 20L)
  precise <- fixture$gds
  imprecise <- fixture$gds
  precise$assays$beta[1, 1, 3] <- precise$assays$beta[1, 1, 3] + 2
  imprecise$assays$beta <- precise$assays$beta
  precise$assays$var[1, 1, 3] <- 1e-4
  imprecise$assays$var[1, 1, 3] <- 100

  fit_time <- function(g) {
    as.numeric(assay(
      reduce(as_plan(g), method = "lmm:ri_knownvar", formula = ~ time) |> compute(),
      "coef:time"
    ))
  }
  ordinary <- fit_time(fixture$gds)
  precise_value <- fit_time(precise)
  imprecise_value <- fit_time(imprecise)
  expect_gt(abs(precise_value - ordinary), abs(imprecise_value - ordinary))
})

test_that("non-finite known variances invalidate only their sample", {
  fixture <- make_lmm_knownvar_fixture("ri", n_subject = 20L)
  beta <- array(NA_real_, dim = c(2L, 20L, 3L))
  sampling_var <- array(NA_real_, dim = dim(beta))
  beta[1, , ] <- beta[2, , ] <- assay(fixture$gds, "beta")[1, , ]
  sampling_var[1, , ] <- sampling_var[2, , ] <- assay(fixture$gds, "var")[1, , ]
  sampling_var[2, 1, 1] <- NA_real_
  g <- new_gds(
    assays = list(beta = beta, var = sampling_var),
    space = space_sample_labels(c("ROI_1", "ROI_2")),
    subjects = fixture$subjects,
    contrasts = fixture$contrasts
  )
  g <- with_contrast_data(
    g,
    data.frame(time = fixture$time, row.names = fixture$contrasts)
  )

  result <- reduce(as_plan(g), method = "lmm:ri_knownvar", formula = ~ time) |> compute()
  expect_true(is.finite(as.numeric(assay(result, "coef:time")[1, 1, 1])))
  expect_true(is.na(as.numeric(assay(result, "coef:time")[2, 1, 1])))
  expect_equal(as.numeric(assay(result, "converged")[2, 1, 1]), 0)
})
