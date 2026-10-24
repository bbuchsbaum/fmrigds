.make_experimental_gds <- function(beta, variance, dims3, subjects = NULL,
                                   mask_idx = NULL, use_se = FALSE) {
  if (is.null(subjects)) subjects <- paste0("s", seq_len(dim(beta)[2L]))
  uncertainty <- if (use_se) list(se = sqrt(variance)) else list(var = variance)
  new_gds(
    assays = c(list(beta = beta), uncertainty),
    space = space_voxel(
      dims3, diag(4), mask_idx = mask_idx,
      storage = if (is.null(mask_idx)) "dense" else "packed"
    ),
    subjects = subjects,
    contrasts = dimnames(beta)[[3L]] %||% "task"
  )
}

.shift_synthetic_volume <- function(x, dx = 0L, dy = 0L, dz = 0L) {
  out <- array(0, dim(x))
  dims <- dim(x)
  target <- arrayInd(seq_len(length(x)), dims)
  source <- sweep(target, 2L, c(dx, dy, dz), "+")
  ok <- rowSums(source < 1L | source > rep(dims, each = nrow(source))) == 0L
  source_idx <- source[ok, 1L] +
    (source[ok, 2L] - 1L) * dims[1L] +
    (source[ok, 3L] - 1L) * dims[1L] * dims[2L]
  out[ok] <- x[source_idx]
  out
}

.experimental_pattern <- function(dims3 = c(11L, 11L, 1L)) {
  pattern <- array(0, dims3)
  pattern[6L, 6L, 1L] <- 1
  pattern[5L, 6L, 1L] <- 0.65
  pattern[7L, 6L, 1L] <- 0.45
  pattern[6L, 5L, 1L] <- 0.55
  pattern[6L, 7L, 1L] <- 0.25
  pattern[5L, 7L, 1L] <- 0.35
  pattern
}

.experimental_region <- function(dims3 = c(11L, 11L, 1L)) {
  region <- array(0L, dims3)
  region[3:9, 3:9, 1L] <- 1L
  as.integer(region)
}

.make_shifted_cohort <- function(seed = 4102L, noise_sd = 0.02,
                                 shifts = NULL) {
  set.seed(seed)
  dims3 <- c(11L, 11L, 1L)
  n_subject <- 36L
  if (is.null(shifts)) {
    shifts <- rep(c(-1L, 0L, 1L), length.out = n_subject)
  }
  pattern <- .experimental_pattern(dims3)
  beta <- array(0, c(prod(dims3), n_subject, 1L))
  for (i in seq_len(n_subject)) {
    truth <- .shift_synthetic_volume(pattern, dx = shifts[i])
    beta[, i, 1L] <- as.numeric(truth) +
      stats::rnorm(prod(dims3), sd = noise_sd)
  }
  variance <- array(noise_sd^2, dim(beta))
  list(
    gds = .make_experimental_gds(beta, variance, dims3),
    beta = beta,
    variance = variance,
    dims3 = dims3,
    shifts = shifts,
    region = .experimental_region(dims3)
  )
}

.reference_cancellation_probability <- function(beta, variance, delta,
                                                equivalence, prevalence) {
  w <- 1 / variance
  mean_fe <- sum(w * beta) / sum(w)
  Q <- sum(w * (beta - mean_fe)^2)
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (length(beta) - 1L)) / C)
  wstar <- 1 / (variance + tau2)
  mean_re <- sum(wstar * beta) / sum(wstar)
  se_re <- sqrt(1 / sum(wstar))
  tau <- sqrt(tau2)
  lower <- max(delta + tau * qnorm(prevalence), -equivalence)
  upper <- min(-delta - tau * qnorm(prevalence), equivalence)
  if (tau <= 0 || upper <= lower) return(0)
  pnorm((upper - mean_re) / se_re) - pnorm((lower - mean_re) / se_re)
}

test_that("cancellation probability uses one subject-level map per subject", {
  set.seed(4101)
  n_subject <- 48L
  dims3 <- c(4L, 1L, 1L)
  beta <- array(0, c(prod(dims3), n_subject, 1L))

  opposing <- rep(c(-1, 1), length.out = n_subject)
  coherent <- rep(1, n_subject)
  null <- rnorm(n_subject, sd = 0.1)
  imbalanced <- c(rep(1, 40L), rep(-1, 8L))
  beta[, , 1L] <- rbind(opposing, coherent, null, imbalanced)
  variance <- array(0.01, dim(beta))

  result <- experimental_cancellation(
    .make_experimental_gds(beta, variance, dims3),
    delta = 0.2,
    equivalence = 0.2,
    prevalence = 0.2,
    regions = integer(prod(dims3))
  )

  probability <- assay(result, "cancellation_probability")[, 1L, 1L]
  expect_true(all(is.finite(probability)))
  expect_true(all(probability >= 0 & probability <= 1))
  expect_gt(probability[1L], 0.7)
  expect_lt(probability[2L], 0.05)
  expect_lt(probability[3L], 0.10)
  expect_lt(probability[4L], 0.05)
  expect_true(all(is.na(assay(result, "shift_rescue_fraction"))))

  expected <- .reference_cancellation_probability(
    opposing,
    variance = rep(0.01, n_subject),
    delta = 0.2,
    equivalence = 0.2,
    prevalence = 0.2
  )
  expect_equal(probability[1L], expected, tolerance = 1e-10)
})

test_that("cancellation probability respects sign, scale, and order invariance", {
  n_subject <- 40L
  dims3 <- c(2L, 1L, 1L)
  opposing <- c(rep(-0.9, n_subject / 2L), rep(0.9, n_subject / 2L))
  beta <- array(rep(opposing, each = 2L), c(2L, n_subject, 1L))
  variance <- array(0.02, dim(beta))
  empty_regions <- integer(prod(dims3))
  baseline <- experimental_cancellation(
    .make_experimental_gds(beta, variance, dims3),
    delta = 0.15, equivalence = 0.2, regions = empty_regions
  )

  order <- sample.int(n_subject)
  transformed <- experimental_cancellation(
    .make_experimental_gds(
      -beta[, order, , drop = FALSE],
      variance[, order, , drop = FALSE],
      dims3,
      subjects = paste0("s", order)
    ),
    delta = 0.15, equivalence = 0.2, regions = empty_regions
  )
  scaled <- experimental_cancellation(
    .make_experimental_gds(10 * beta, 100 * variance, dims3),
    delta = 1.5, equivalence = 2, regions = empty_regions
  )

  expect_equal(
    assay(baseline, "cancellation_probability"),
    assay(transformed, "cancellation_probability"),
    tolerance = 1e-10
  )
  expect_equal(
    assay(baseline, "cancellation_probability"),
    assay(scaled, "cancellation_probability"),
    tolerance = 1e-10
  )
})

test_that("automatic support pools nearby locations across subjects", {
  dims3 <- c(5L, 5L, 1L)
  n_subject <- 18L
  probability <- matrix(0, nrow = prod(dims3), ncol = n_subject)
  positions <- as.integer(outer(2:4, (2:4 - 1L) * dims3[1L], "+"))
  for (i in seq_len(n_subject)) {
    probability[positions[(i - 1L) %% length(positions) + 1L], i] <- 1
  }
  spec <- list(dim = dims3, n_full = prod(dims3), mask_idx = seq_len(prod(dims3)))
  support <- fmrigds:::.experimental_local_support(
    probability, spec, radius = 1L, support_probability = 0.95,
    min_subjects = 8L
  )
  center <- 3L + (3L - 1L) * dims3[1L]

  expect_lt(max(rowMeans(probability)), 0.2)
  expect_equal(support[center], 1)
})

test_that("between-subject translations rescue a wandering same-signed activation", {
  cohort <- .make_shifted_cohort()
  result <- experimental_cancellation(
    cohort$gds,
    delta = 0.15,
    equivalence = 0.2,
    prevalence = 0.2,
    shift_radius = 1L,
    patch_radius = 1L,
    regions = cohort$region,
    min_correlation_gain = 0.02
  )

  rescue <- assay(result, "shift_rescue_fraction")[cohort$region == 1L, 1L, 1L]
  expect_true(all(is.finite(rescue)))
  expect_gt(unique(rescue), 0.25)
  expect_lte(unique(rescue), 1)

  # Same-signed spatial dispersion does not require voxelwise sign cancellation.
  probability <- assay(result, "cancellation_probability")[, 1L, 1L]
  expect_lt(max(probability, na.rm = TRUE), 0.1)

  receipt <- metadata(result)$experimental$region_receipts[[1L]][["1"]]
  expect_gte(receipt$n_shifted, length(cohort$shifts) / 2)
  expect_identical(metadata(result)$experimental$shift_scope, "between_subject")

  automatic <- experimental_cancellation(
    cohort$gds,
    delta = 0.15,
    equivalence = 0.2,
    prevalence = 0.2,
    shift_radius = 1L,
    patch_radius = 1L,
    min_correlation_gain = 0.02
  )
  expect_gt(max(assay(automatic, "shift_rescue_fraction"), na.rm = TRUE), 0.25)
})

test_that("spatial rescue is invariant to subject order, sign, and units", {
  cohort <- .make_shifted_cohort(seed = 4104L)
  baseline <- experimental_cancellation(
    cohort$gds, delta = 0.15, regions = cohort$region,
    min_correlation_gain = 0.02
  )
  order <- sample.int(dim(cohort$beta)[2L])
  transformed <- experimental_cancellation(
    .make_experimental_gds(
      -10 * cohort$beta[, order, , drop = FALSE],
      100 * cohort$variance[, order, , drop = FALSE],
      cohort$dims3,
      subjects = paste0("s", order)
    ),
    delta = 1.5,
    equivalence = 1.5,
    regions = cohort$region,
    min_correlation_gain = 0.02
  )

  expect_equal(
    assay(baseline, "shift_rescue_fraction"),
    assay(transformed, "shift_rescue_fraction"),
    tolerance = 1e-10
  )
})

test_that("amplitude and polarity heterogeneity are not spatial rescue", {
  set.seed(4105)
  dims3 <- c(11L, 11L, 1L)
  n_subject <- 36L
  pattern <- .experimental_pattern(dims3)
  region <- .experimental_region(dims3)
  variance <- array(0.02^2, c(prod(dims3), n_subject, 1L))

  amplitude_beta <- array(0, dim(variance))
  amplitude <- seq(0.6, 1.4, length.out = n_subject)
  for (i in seq_len(n_subject)) {
    amplitude_beta[, i, 1L] <- amplitude[i] * as.numeric(pattern) +
      rnorm(prod(dims3), sd = 0.02)
  }
  amplitude_result <- experimental_cancellation(
    .make_experimental_gds(amplitude_beta, variance, dims3),
    delta = 0.15, regions = region, min_correlation_gain = 0.02
  )
  amplitude_rescue <- unique(
    assay(amplitude_result, "shift_rescue_fraction")[region == 1L, 1L, 1L]
  )
  expect_true(is.finite(amplitude_rescue))
  expect_lt(amplitude_rescue, 0.1)

  polarity_beta <- amplitude_beta
  signs <- rep(c(-1, 1), length.out = n_subject)
  for (i in seq_len(n_subject)) {
    polarity_beta[, i, 1L] <- signs[i] * as.numeric(pattern) +
      rnorm(prod(dims3), sd = 0.02)
  }
  polarity_result <- experimental_cancellation(
    .make_experimental_gds(polarity_beta, variance, dims3),
    delta = 0.15, equivalence = 0.2, regions = region,
    min_correlation_gain = 0.02
  )
  center <- 6L + (6L - 1L) * dims3[1L]
  expect_gt(
    assay(polarity_result, "cancellation_probability")[center, 1L, 1L],
    0.5
  )
  polarity_rescue <- unique(
    assay(polarity_result, "shift_rescue_fraction")[region == 1L, 1L, 1L]
  )
  expect_true(is.na(polarity_rescue) || polarity_rescue < 0.1)
})

test_that("unstructured noise is not reported as local spatial rescue", {
  set.seed(4107)
  dims3 <- c(11L, 11L, 1L)
  n_subject <- 36L
  noise_sd <- 0.12
  beta <- array(
    rnorm(prod(dims3) * n_subject, sd = noise_sd),
    c(prod(dims3), n_subject, 1L)
  )
  variance <- array(noise_sd^2, dim(beta))
  region <- .experimental_region(dims3)
  result <- experimental_cancellation(
    .make_experimental_gds(beta, variance, dims3),
    delta = 0.1,
    regions = region
  )
  rescue <- unique(
    assay(result, "shift_rescue_fraction")[region == 1L, 1L, 1L]
  )
  expect_true(is.na(rescue) || rescue < 0.1)

  automatic <- experimental_cancellation(
    .make_experimental_gds(beta, variance, dims3),
    delta = 0.1
  )
  expect_true(all(is.na(assay(automatic, "shift_rescue_fraction"))))
})

test_that("zero search radius cannot manufacture spatial rescue", {
  cohort <- .make_shifted_cohort(seed = 4106L)
  result <- experimental_cancellation(
    cohort$gds,
    delta = 0.15,
    regions = cohort$region,
    shift_radius = 0L,
    min_correlation_gain = 0.02
  )
  rescue <- unique(
    assay(result, "shift_rescue_fraction")[cohort$region == 1L, 1L, 1L]
  )
  expect_equal(rescue, 0)
})

test_that("a known subject displacement is recovered against an external template", {
  dims3 <- c(9L, 9L, 1L)
  pattern <- array(0, dims3)
  pattern[5L, 5L, 1L] <- 1
  pattern[4L, 5L, 1L] <- 0.6
  pattern[5L, 6L, 1L] <- 0.3
  n_subject <- 12L
  beta <- array(rep(as.numeric(pattern), n_subject),
                c(prod(dims3), n_subject, 1L))
  beta[, 1L, 1L] <- as.numeric(.shift_synthetic_volume(pattern, dx = 1L))
  variance <- array(0.001, dim(beta))
  patch_idx <- which(.experimental_region(dims3) == 1L)

  fit <- fmrigds:::.experimental_subject_shift(
    beta[, , 1L], variance[, , 1L], subject = 1L,
    patch_idx = patch_idx, dims3 = dims3,
    offsets = fmrigds:::.experimental_offset_grid(1L, dims3),
    delta = 0.1, translation_penalty = 0,
    mass_retention = c(0.75, 1.25), min_correlation_gain = 0.05,
    min_aligned_correlation = 0.5
  )
  expect_true(fit$shifted)
  expect_identical(fit$shift, c(-1L, 0L, 0L))
  expect_equal(fit$aligned_correlation, 1, tolerance = 1e-12)
})

test_that("experimental cancellation validates scientific and spatial contracts", {
  dims3 <- c(3L, 3L, 1L)
  beta <- array(0, c(prod(dims3), 8L, 1L))
  variance <- array(1, dim(beta))
  x <- .make_experimental_gds(beta, variance, dims3)

  expect_error(experimental_cancellation(x, delta = 0), "delta")
  expect_error(experimental_cancellation(x, delta = 0.1, prevalence = 0.5),
               "prevalence")
  expect_error(experimental_cancellation(
    x, delta = 0.1, support_probability = 0.5
  ), "support_probability")
  expect_error(experimental_cancellation(
    x, delta = 0.1, min_correlation_gain = -0.1
  ), "min_correlation_gain")
  expect_error(experimental_cancellation(
    x, delta = 0.1, mass_retention = c(1, 0.5)
  ), "mass_retention")
  expect_error(experimental_cancellation(
    x, delta = 0.1, regions = 1:3
  ), "regions")

  labels <- new_gds(
    list(beta = beta, var = variance),
    space_sample_labels(paste0("v", seq_len(prod(dims3)))),
    subjects(x),
    "task"
  )
  expect_error(experimental_cancellation(labels, delta = 0.1), "voxel")

  no_variance <- x
  no_variance$assays$var <- NULL
  expect_error(experimental_cancellation(no_variance, delta = 0.1), "var.*se")

  synthetic_variance <- x
  synthetic_variance$metadata$synthetic_var <- TRUE
  expect_error(
    experimental_cancellation(synthetic_variance, delta = 0.1),
    "genuine first-level"
  )

  se_result <- experimental_cancellation(
    .make_experimental_gds(beta, variance, dims3, use_se = TRUE),
    delta = 0.1,
    regions = integer(prod(dims3))
  )
  expect_identical(
    dim(assay(se_result, "cancellation_probability")),
    c(as.integer(prod(dims3)), 1L, 1L)
  )
})

test_that("experimental cancellation preserves packed voxel support", {
  dims3 <- c(5L, 5L, 1L)
  mask_idx <- c(7L, 8L, 9L, 12L, 13L, 14L, 17L, 18L, 19L)
  n_subject <- 24L
  beta <- array(0, c(length(mask_idx), n_subject, 1L))
  opposing <- rep(c(-0.8, 0.8), length.out = n_subject)
  beta[, , 1L] <- matrix(
    rep(opposing, each = length(mask_idx)),
    nrow = length(mask_idx)
  )
  variance <- array(0.01, dim(beta))
  x <- .make_experimental_gds(
    beta, variance, dims3, mask_idx = mask_idx
  )

  result <- experimental_cancellation(
    x, delta = 0.15, equivalence = 0.2,
    regions = integer(length(mask_idx))
  )
  expect_identical(space(result)$mask_idx, mask_idx)
  expect_identical(
    dim(assay(result, "cancellation_probability")),
    c(length(mask_idx), 1L, 1L)
  )
  expect_true(all(is.finite(assay(result, "cancellation_probability"))))
})
