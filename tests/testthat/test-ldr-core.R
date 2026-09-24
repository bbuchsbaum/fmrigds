test_that("half-normal marginal evidence agrees with numerical integration", {
  u <- c(-2, 0, 3)
  q <- c(0.5, 2, 4)
  scale <- 1.2
  expected <- vapply(seq_along(u), function(i) {
    integrand <- function(a) {
      exp(a * u[i] - 0.5 * q[i] * a^2) *
        2 * stats::dnorm(a, sd = scale)
    }
    log(stats::integrate(integrand, lower = 0, upper = Inf,
                         rel.tol = 1e-12)$value)
  }, numeric(1L))

  expect_equal(
    .ldr_halfnormal_log_bf(u, q, scale),
    expected,
    tolerance = 1e-10
  )
})

test_that("LDR geometry uses physical millimetres and preserves signed shifts", {
  geometry <- .ldr_patch_geometry(
    center = 13L,
    dims3 = c(5L, 5L, 1L),
    affine = diag(c(2, 3, 4, 1)),
    patch_radius_mm = 4,
    shift_radius_mm = 3
  )

  expect_true(all(geometry$distance_mm <= 4 + 1e-12))
  expect_true(all(geometry$shift_distance_mm <= 3 + 1e-12))
  has_offset <- function(offset) {
    any(rowSums(abs(sweep(geometry$shift_offsets, 2L, offset, "-"))) == 0L)
  }
  expect_true(has_offset(c(1L, 0L, 0L)))
  expect_true(has_offset(c(0L, 1L, 0L)))
  expect_false(has_offset(c(1L, 1L, 0L)))
  expect_equal(sum(.ldr_shift_prior(geometry$shift_distance_mm, 2)), 1,
               tolerance = 1e-12)
  zero_prior <- .ldr_shift_prior(geometry$shift_distance_mm, 0)
  expect_equal(sum(zero_prior > 0), 1L)
  expect_equal(geometry$shift_distance_mm[zero_prior > 0], 0)
})

test_that("cross-fitted LDR distinguishes displacement from aligned signal", {
  shifted <- .ldr_test_fit(.ldr_test_fixture("shifted", seed = 5102L))
  aligned <- .ldr_test_fit(.ldr_test_fixture("aligned", seed = 5103L))

  expect_gt(shifted$activation_t, 3)
  expect_gt(shifted$displacement_t, 3)
  expect_gt(shifted$mni_loss, 2)
  expect_gt(shifted$shift_rms_mm, 0)
  expect_lt(abs(shifted$ordinary_t), shifted$counterfactual_t)

  expect_gt(aligned$activation_t, 3)
  expect_lte(aligned$displacement_score, 0)
  expect_lte(aligned$displacement_t, 0)
  expect_lt(abs(aligned$mni_loss), 1e-8)
})

test_that("signed-amplitude competitor rejects fixed polarity heterogeneity", {
  polarity <- .ldr_test_fit(.ldr_test_fixture("polarity", seed = 5104L))

  expect_gt(polarity$displacement_components[["aligned"]], 0)
  expect_lt(polarity$displacement_components[["heterogeneity"]], 0)
  expect_lt(polarity$displacement_score, 0)
  expect_lt(polarity$displacement_t, 0)
})

test_that("LDR held-out scores are invariant to subject storage order", {
  fixture <- .ldr_test_fixture("shifted", seed = 5105L, n_subject = 24L)
  baseline <- .ldr_test_fit(fixture)
  order <- c(seq(2L, 24L, by = 2L), seq(1L, 23L, by = 2L))
  reordered <- fixture
  reordered$subjects <- fixture$subjects[order]
  reordered$patch$y_a <- fixture$patch$y_a[, order, drop = FALSE]
  reordered$patch$y_b <- fixture$patch$y_b[, order, drop = FALSE]
  reordered$patch$var_a <- fixture$patch$var_a[, order, drop = FALSE]
  reordered$patch$var_b <- fixture$patch$var_b[, order, drop = FALSE]
  reordered$patch$subjects <- fixture$patch$subjects[order]
  permuted <- .ldr_test_fit(reordered)

  expect_equal(permuted$activation_score, baseline$activation_score,
               tolerance = 1e-8)
  expect_equal(permuted$displacement_score, baseline$displacement_score,
               tolerance = 1e-8)
  expect_equal(permuted$mni_loss, baseline$mni_loss, tolerance = 1e-8)
  expect_equal(
    permuted$jitter_log_bf[order(order)],
    baseline$jitter_log_bf,
    tolerance = 1e-8
  )
})

test_that("LDR reference core fails closed on invalid or incomplete patches", {
  expect_error(
    .ldr_patch_geometry(
      center = 1L,
      dims3 = c(5L, 1L, 1L),
      affine = diag(4),
      patch_radius_mm = 2,
      shift_radius_mm = 1
    ),
    "crosses the image boundary"
  )
  expect_error(
    .ldr_patch_geometry(
      center = 3L,
      dims3 = c(5L, 1L, 1L),
      affine = diag(4),
      patch_radius_mm = 2,
      shift_radius_mm = 1,
      mask_idx = c(1L, 2L, 3L, 5L)
    ),
    "outside the active mask"
  )
  expect_error(.ldr_halfnormal_log_bf(1, -1, 1), "q >= 0")
})

test_that("exact LDR ignores missing measurements at zero-weight patch voxels", {
  fixture <- .ldr_test_fixture("shifted", seed = 5107L, n_subject = 18L)
  y_a <- fixture$patch$y_a
  var_a <- fixture$patch$var_a
  y_a[1L, 1L] <- NA_real_
  var_a[1L, 1L] <- NA_real_
  patch <- .ldr_prepare_patch(
    y_a,
    fixture$patch$y_b,
    var_a,
    fixture$patch$var_b,
    fixture$geometry
  )

  fit <- .ldr_crossfit(
    patch,
    fixture$subjects,
    tau_grid_mm = c(1, 2),
    folds = 3L,
    iterations = 6L
  )

  expect_true(is.finite(fit$activation_t))
  expect_true(is.finite(fit$displacement_t))
  expect_true(is.finite(fit$mni_loss))
})

test_that("ROI calibration is reproducible and restores caller RNG state", {
  fixture <- .ldr_test_fixture("shifted", seed = 5106L, n_subject = 12L)
  observed <- .ldr_crossfit(
    fixture$patch,
    fixture$subjects,
    tau_grid_mm = c(1, 2),
    folds = 2L,
    iterations = 4L
  )
  calibrate <- function() {
    .ldr_calibrate(
      fixture$patch,
      subject_ids = fixture$subjects,
      observed = observed,
      tau_grid_mm = c(1, 2),
      folds = 2L,
      iterations = 4L,
      ridge = 1e-4,
      tolerance = 1e-5,
      n_resamples = 2L,
      seed = 91L
    )
  }

  set.seed(8801)
  before <- .Random.seed
  first <- calibrate()
  expect_identical(.Random.seed, before)
  second <- calibrate()

  expect_identical(first$component_p, second$component_p)
  expect_identical(first$null, second$null)
  expect_true(all(first$component_p >= 1 / 3 & first$component_p <= 1))
  expect_equal(
    first$p_ldr,
    max(first$p_activation, first$p_displacement, first$p_loss)
  )
})

test_that("calibrated conjunction rejects aligned and polarity mechanisms", {
  calibrate_mode <- function(mode, seed) {
    fixture <- .ldr_test_fixture(mode, seed = seed, n_subject = 12L)
    observed <- .ldr_crossfit(
      fixture$patch,
      fixture$subjects,
      tau_grid_mm = c(1, 2),
      folds = 2L,
      iterations = 4L
    )
    .ldr_calibrate(
      fixture$patch,
      fixture$subjects,
      observed,
      tau_grid_mm = c(1, 2),
      folds = 2L,
      iterations = 4L,
      ridge = 1e-4,
      tolerance = 1e-5,
      n_resamples = 5L,
      seed = seed + 100L
    )
  }

  aligned <- calibrate_mode("aligned", 71L)
  shifted <- calibrate_mode("shifted", 72L)
  polarity <- calibrate_mode("polarity", 73L)

  expect_equal(aligned$p_displacement, 1)
  expect_equal(aligned$p_ldr, 1)
  expect_equal(polarity$p_displacement, 1)
  expect_equal(polarity$p_ldr, 1)
  expect_equal(shifted$p_ldr, shifted$minimum_p)
})
