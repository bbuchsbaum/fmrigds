.huber_reference_fit <- function(y, X, variance, tuning_constant, start) {
  rho <- function(value) {
    magnitude <- abs(value)
    ifelse(
      magnitude <= tuning_constant,
      0.5 * value^2,
      tuning_constant * magnitude - 0.5 * tuning_constant^2
    )
  }
  objective <- function(coefficients) {
    standardized <- (y - drop(X %*% coefficients)) / sqrt(variance)
    sum(rho(standardized))
  }
  gradient <- function(coefficients) {
    standardized <- (y - drop(X %*% coefficients)) / sqrt(variance)
    psi <- pmax(-tuning_constant, pmin(tuning_constant, standardized))
    -drop(crossprod(X, psi / sqrt(variance)))
  }
  stats::optim(
    start,
    objective,
    gr = gradient,
    method = "BFGS",
    control = list(reltol = 1e-13, maxit = 10000L)
  )
}

.robust_test_settings <- function(...) {
  examination_control(robust = list(...))$robust
}

test_that("robust examination control is explicit and default-off compatible", {
  expect_identical(examination_control(), examination_control(robust = NULL))
  expect_null(examination_control()$robust)

  settings <- examination_control(robust = list())$robust
  expect_identical(settings$method, "huber_ivw")
  expect_equal(settings$tuning_constant, 1.345)
  expect_identical(settings$max_iterations, 100L)
  expect_equal(settings$convergence_tolerance, 1e-8)
  expect_lt(settings$leverage_limit, 1)

  expect_error(
    examination_control(robust = list(method = "bisquare")),
    "method"
  )
  expect_error(
    examination_control(robust = list(tuning_constant = 0)),
    "tuning_constant"
  )
  expect_error(
    examination_control(robust = list(max_iterations = 0L)),
    "max_iterations"
  )
  expect_error(
    examination_control(robust = list(convergence_tolerance = Inf)),
    "convergence_tolerance"
  )
  expect_error(
    examination_control(robust = list(leverage_limit = 1)),
    "leverage_limit"
  )
  expect_error(
    examination_control(robust = list(unknown = TRUE)),
    "Unknown robust settings"
  )
})

test_that("Huber IVW regression agrees with a direct objective oracle", {
  y <- c(-0.2, 0, 0.1, 0.2, 8)
  X <- cbind(intercept = 1, trend = seq_along(y))
  variance <- c(0.8, 1.1, 0.9, 1.2, 1)
  settings <- .robust_test_settings(
    max_iterations = 500L,
    convergence_tolerance = 1e-10
  )
  start <- drop(solve(crossprod(X / sqrt(variance))) %*%
    crossprod(X, y / variance))
  reference <- .huber_reference_fit(
    y, X, variance, settings$tuning_constant, start
  )
  fit <- fmrigds:::.huber_ivw_fit(
    y, X, variance, settings, rank_tolerance = 1e-10
  )

  expect_identical(fit$status, "available")
  expect_identical(reference$convergence, 0L)
  expect_equal(fit$objective, reference$value, tolerance = 1e-8)
  expect_equal(fit$coefficients, unname(reference$par), tolerance = 1e-6)
  expect_gt(fit$downweighted_n, 0L)
  expect_true(all(fit$robust_factor[is.finite(fit$robust_factor)] > 0))
  expect_true(all(fit$robust_factor[is.finite(fit$robust_factor)] <= 1))
})

test_that("Huber IVW reduces to ordinary IVW when every score is quadratic", {
  X <- cbind(intercept = 1, trend = c(-2, -1, 0, 1, 2, 3))
  variance <- c(0.7, 1.2, 0.9, 1.1, 0.8, 1.3)
  y <- drop(X %*% c(0.4, -0.15)) + c(-0.08, 0.06, -0.03, 0.04, -0.05, 0.07)
  ordinary <- stats::lm.wfit(X, y, w = 1 / variance)$coefficients
  fit <- fmrigds:::.huber_ivw_fit(
    y,
    X,
    variance,
    .robust_test_settings(convergence_tolerance = 1e-12),
    1e-10
  )

  expect_identical(fit$status, "available")
  expect_equal(fit$coefficients, unname(ordinary), tolerance = 1e-11)
  expect_identical(fit$downweighted_n, 0L)
  expect_equal(fit$robust_factor, rep(1, nrow(X)))
})

test_that("Huber IVW obeys ordering and scale equivariance", {
  y <- c(-0.3, 0.1, 0.2, 0.4, 6, 0.5)
  X <- cbind(1, c(-2, -1, 0, 1, 2, 3))
  variance <- c(0.7, 1.2, 0.9, 1.1, 0.8, 1.3)
  settings <- .robust_test_settings(
    max_iterations = 500L,
    convergence_tolerance = 1e-10
  )
  fit <- fmrigds:::.huber_ivw_fit(y, X, variance, settings, 1e-10)
  permutation <- c(4, 1, 6, 3, 5, 2)
  reordered <- fmrigds:::.huber_ivw_fit(
    y[permutation], X[permutation, , drop = FALSE], variance[permutation],
    settings, 1e-10
  )
  scaled <- fmrigds:::.huber_ivw_fit(
    7 * y, X, 49 * variance, settings, 1e-10
  )

  expect_identical(fit$status, "available")
  expect_identical(reordered$status, "available")
  expect_identical(scaled$status, "available")
  expect_equal(reordered$coefficients, fit$coefficients, tolerance = 1e-9)
  expect_equal(
    reordered$robust_factor[order(permutation)],
    fit$robust_factor,
    tolerance = 1e-9
  )
  expect_equal(scaled$coefficients, 7 * fit$coefficients, tolerance = 1e-8)
  expect_equal(scaled$robust_factor, fit$robust_factor, tolerance = 1e-9)
  expect_equal(scaled$objective, fit$objective, tolerance = 1e-9)
})

test_that("Huber sensitivity is bounded under isolated and masking contamination", {
  settings <- .robust_test_settings(
    max_iterations = 500L,
    convergence_tolerance = 1e-10
  )
  X <- matrix(1, 10L, 1L)
  variance <- rep(1, 10L)
  moderate <- c(rep(0, 9), 10)
  extreme <- c(rep(0, 9), 1000)
  moderate_fit <- fmrigds:::.huber_ivw_fit(
    moderate, X, variance, settings, 1e-10
  )
  extreme_fit <- fmrigds:::.huber_ivw_fit(
    extreme, X, variance, settings, 1e-10
  )
  masked <- c(rep(0, 8), 8, 8)
  masked_fit <- fmrigds:::.huber_ivw_fit(
    masked, X, variance, settings, 1e-10
  )

  expect_identical(moderate_fit$status, "available")
  expect_identical(extreme_fit$status, "available")
  expect_identical(masked_fit$status, "available")
  expect_lt(abs(extreme_fit$coefficients - moderate_fit$coefficients), 1e-6)
  expect_lt(abs(masked_fit$coefficients), mean(masked) / 2)
  expect_equal(masked_fit$downweighted_n, 2L)
})

test_that("Huber IVW fails closed on numerical boundary states", {
  settings <- .robust_test_settings()
  insufficient <- fmrigds:::.huber_ivw_fit(
    c(1, 2), cbind(1, c(0, 1)), c(1, 1), settings
  )
  rank_deficient <- fmrigds:::.huber_ivw_fit(
    1:5, cbind(1, rep(1, 5)), rep(1, 5), settings
  )
  leverage <- fmrigds:::.huber_ivw_fit(
    c(0, 0, 0, 1), cbind(1, c(0, 0, 0, 100)), rep(1, 4),
    .robust_test_settings(leverage_limit = 0.8)
  )
  nonconverged <- fmrigds:::.huber_ivw_fit(
    c(-0.2, 0, 0.1, 0.2, 8), cbind(1, 1:5), rep(1, 5),
    .robust_test_settings(max_iterations = 1L, convergence_tolerance = 1e-14)
  )

  expect_identical(insufficient$status, "insufficient_samples")
  expect_identical(rank_deficient$status, "rank_deficient")
  expect_identical(leverage$status, "leverage_limit")
  expect_identical(nonconverged$status, "nonconverged")
  expect_null(nonconverged$coefficients)
})

test_that("robust examination is additive and does not alter review semantics", {
  plan <- reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  ordinary <- examine_group(
    plan,
    control = examination_control(block_size = 7L)
  )
  robust <- examine_group(
    plan,
    control = examination_control(block_size = 7L, robust = list())
  )

  expect_false("robust_sensitivity" %in% names(ordinary))
  expect_s3_class(robust$robust_sensitivity, "gds_robust_sensitivity")
  expect_identical(robust$subject_data, ordinary$subject_data)
  expect_identical(robust$contrast_data, ordinary$contrast_data)
  expect_identical(robust$estimand_data, ordinary$estimand_data)
  expect_identical(robust$availability, ordinary$availability)
  expect_identical(
    robust$subject_data$review_status,
    ordinary$subject_data$review_status
  )
  expect_match(robust$robust_sensitivity$interpretation, "does not replace")
  expect_match(robust$robust_sensitivity$interpretation, "outlier probability")
  expect_s3_class(robust$robust_sensitivity$maps, "gds")
  expect_true(all(robust$robust_sensitivity$availability$status == "available"))

  subject_data <- robust$robust_sensitivity$subject_data
  reversed <- subject_data[subject_data$subject == "s10", , drop = FALSE]
  ordinary_subjects <- subject_data[subject_data$subject != "s10", , drop = FALSE]
  expect_lt(
    reversed$mean_downweight_factor,
    min(ordinary_subjects$mean_downweight_factor)
  )
  map_info <- metadata(robust$robust_sensitivity$maps)$examination
  expect_true("status_code" %in% map_info$categorical_assays)
  expect_identical(
    map_info$status_lookup$status,
    names(fmrigds:::.robust_feature_status)
  )

  ordinary_summary <- capture.output(print(summary(ordinary)))
  robust_summary <- capture.output(print(summary(robust)))
  expect_false(any(grepl("robust sensitivity", ordinary_summary, fixed = TRUE)))
  expect_true(any(grepl("robust sensitivity: huber_ivw", robust_summary, fixed = TRUE)))
})

test_that("robust sensitivity is invariant to sample block partitions", {
  plan <- reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  small <- examine_group(
    plan,
    control = examination_control(block_size = 3L, robust = list())
  )$robust_sensitivity
  large <- examine_group(
    plan,
    control = examination_control(block_size = 17L, robust = list())
  )$robust_sensitivity

  expect_equal(small$contrast_data, large$contrast_data, tolerance = 1e-12)
  expect_equal(small$estimand_data, large$estimand_data, tolerance = 1e-12)
  expect_equal(small$subject_data, large$subject_data, tolerance = 1e-12)
  for (name in names(assays(small$maps))) {
    expect_equal(
      assay(small$maps, name),
      assay(large$maps, name),
      tolerance = 1e-12,
      info = name
    )
  }
})

test_that("random-effects sensitivity freezes full-data tau2", {
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:re"),
    control = examination_control(robust = list())
  )
  sensitivity <- exam$robust_sensitivity
  expect_identical(sensitivity$mode, "huber_ivw_tau2_fixed_full")
  expect_identical(
    sensitivity$provenance$tau2_contract,
    "full-data tau2 held fixed"
  )
  expect_true(all(sensitivity$availability$status == "available"))
})

test_that("partial robust availability exposes feature failure counts", {
  g <- .group_examination_fixture()
  g$assays$var[1L, , ] <- NA_real_
  sensitivity <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = examination_control(robust = list())
  )$robust_sensitivity

  expect_identical(sensitivity$contrast_data$status, "partial")
  expect_identical(sensitivity$contrast_data$feature_n, 40L)
  expect_identical(sensitivity$contrast_data$available_n, 39L)
  expect_identical(sensitivity$contrast_data$insufficient_samples_n, 1L)
  expect_identical(sensitivity$availability$status, "partial")
  expect_match(sensitivity$availability$reason, "39 of 40")
  status <- assay(sensitivity$maps, "status_code")[, 1L, 1L]
  expect_identical(
    status[1L],
    unname(fmrigds:::.robust_feature_status[["insufficient_samples"]])
  )
})

test_that("unsupported robust models report availability without fabricating output", {
  exam <- examine_group(
    reduce(
      as_plan(.group_examination_fixture()),
      method = "ols:voxelwise",
      formula = ~ 1
    ),
    estimands = "(Intercept)",
    control = examination_control(robust = list())
  )
  sensitivity <- exam$robust_sensitivity
  expect_true(all(sensitivity$availability$status == "unsupported_reducer"))
  expect_match(sensitivity$availability$reason, "measured-variance meta reducer")
  expect_null(sensitivity$maps)
  expect_true(all(is.na(sensitivity$subject_data$mean_downweight_factor)))
  expect_false(any(exam$subject_data$review_source == "robust", na.rm = TRUE))
})
