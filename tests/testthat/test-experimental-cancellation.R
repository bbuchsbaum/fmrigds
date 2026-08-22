.make_experimental_split <- function(beta, var, dims3, subjects = NULL) {
  if (is.null(subjects)) subjects <- paste0("s", seq_len(dim(beta)[2L]))
  new_gds(
    assays = list(beta = beta, var = var),
    space = space_voxel(dims3, diag(4), storage = "dense"),
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

test_that("experimental cancellation separates opposing effects from coherence and noise", {
  set.seed(4101)
  n_subject <- 48L
  dims3 <- c(4L, 1L, 1L)
  beta_a <- beta_b <- array(0, c(prod(dims3), n_subject, 1L))

  opposing <- rep(c(-1, 1), length.out = n_subject)
  coherent <- rep(1, n_subject)
  null_a <- rnorm(n_subject, sd = 0.1)
  null_b <- rnorm(n_subject, sd = 0.1)
  imbalanced <- c(rep(1, 40L), rep(-1, 8L))

  beta_a[, , 1L] <- rbind(opposing, coherent, null_a, imbalanced)
  beta_b[, , 1L] <- rbind(opposing, coherent, null_b, imbalanced)
  var_a <- var_b <- array(0.01, dim(beta_a))

  result <- experimental_cancellation(
    .make_experimental_split(beta_a, var_a, dims3),
    .make_experimental_split(beta_b, var_b, dims3),
    delta = 0.2,
    equivalence = 0.2,
    prevalence = 0.2,
    candidate_threshold = 0.5
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
    variance = rep(0.005, n_subject),
    delta = 0.2,
    equivalence = 0.2,
    prevalence = 0.2
  )
  expect_equal(probability[1L], expected, tolerance = 1e-10)
})

test_that("cancellation probability is invariant to subject order and global sign", {
  n_subject <- 40L
  dims3 <- c(2L, 1L, 1L)
  opposing <- c(rep(-0.9, n_subject / 2L), rep(0.9, n_subject / 2L))
  beta <- array(0, c(prod(dims3), n_subject, 1L))
  beta[, , 1L] <- rbind(opposing, opposing)
  var <- array(0.02, dim(beta))
  a <- .make_experimental_split(beta, var, dims3)
  b <- .make_experimental_split(beta, var, dims3)
  baseline <- experimental_cancellation(a, b, delta = 0.15, equivalence = 0.2)

  order <- sample.int(n_subject)
  a_permuted <- a
  b_permuted <- b
  a_permuted$assays <- lapply(a$assays, function(x) x[, order, , drop = FALSE])
  b_permuted$assays <- lapply(b$assays, function(x) x[, order, , drop = FALSE])
  a_permuted$assays$beta <- -a_permuted$assays$beta
  b_permuted$assays$beta <- -b_permuted$assays$beta
  # Variance is not signed.
  b_permuted$assays$var <- b$assays$var[, order, , drop = FALSE]
  a_permuted$subjects <- a$subjects[order]
  b_permuted$subjects <- b$subjects[order]
  signed <- experimental_cancellation(a_permuted, b_permuted,
                                      delta = 0.15, equivalence = 0.2)

  expect_equal(
    assay(baseline, "cancellation_probability"),
    assay(signed, "cancellation_probability"),
    tolerance = 1e-10
  )

  scaled_a <- a
  scaled_b <- b
  scaled_a$assays$beta <- 10 * a$assays$beta
  scaled_b$assays$beta <- 10 * b$assays$beta
  scaled_a$assays$var <- 100 * a$assays$var
  scaled_b$assays$var <- 100 * b$assays$var
  scaled <- experimental_cancellation(
    scaled_a, scaled_b, delta = 1.5, equivalence = 2
  )
  expect_equal(
    assay(baseline, "cancellation_probability"),
    assay(scaled, "cancellation_probability"),
    tolerance = 1e-10
  )
})

test_that("cross-fitted local translations rescue a shifted opponent pattern", {
  set.seed(4102)
  dims3 <- c(9L, 9L, 1L)
  n_subject <- 32L
  base <- array(0, dims3)
  base[5L, 5L, 1L] <- 0.9
  base[4:6, 4:6, 1L] <- -0.9
  base[5L, 5L, 1L] <- 0.9
  shifts <- rep(c(0L, 0L, -1L, 1L), length.out = n_subject)

  beta_a <- beta_b <- array(0, c(prod(dims3), n_subject, 1L))
  for (i in seq_len(n_subject)) {
    truth <- .shift_synthetic_volume(base, dx = shifts[i])
    beta_a[, i, 1L] <- as.numeric(truth) + rnorm(prod(dims3), sd = 0.03)
    beta_b[, i, 1L] <- as.numeric(truth) + rnorm(prod(dims3), sd = 0.03)
  }
  var_a <- var_b <- array(0.03^2, dim(beta_a))
  region <- array(0L, dims3)
  region[3:7, 3:7, 1L] <- 1L

  result <- experimental_cancellation(
    .make_experimental_split(beta_a, var_a, dims3),
    .make_experimental_split(beta_b, var_b, dims3),
    delta = 0.2,
    equivalence = 0.25,
    prevalence = 0.2,
    candidate_threshold = 0.35,
    shift_radius = 1L,
    patch_radius = 1L,
    regions = as.integer(region)
  )

  center <- 5L + (5L - 1L) * dims3[1L]
  probability <- assay(result, "cancellation_probability")[center, 1L, 1L]
  rescue <- assay(result, "shift_rescue_fraction")[center, 1L, 1L]
  expect_gt(probability, 0.5)
  expect_true(is.finite(rescue))
  expect_gt(rescue, 0.5)
  expect_lte(rescue, 1)

  automatic <- experimental_cancellation(
    .make_experimental_split(beta_a, var_a, dims3),
    .make_experimental_split(beta_b, var_b, dims3),
    delta = 0.2,
    equivalence = 0.25,
    prevalence = 0.2,
    candidate_threshold = 0.35,
    shift_radius = 1L,
    patch_radius = 1L
  )
  expect_gt(assay(automatic, "shift_rescue_fraction")[center, 1L, 1L], 0.5)
})

test_that("fixed polarity heterogeneity is not called spatial rescue", {
  set.seed(4103)
  dims3 <- c(9L, 9L, 1L)
  n_subject <- 32L
  base <- array(0, dims3)
  base[4:6, 4:6, 1L] <- -0.8
  base[5L, 5L, 1L] <- 1.4
  signs <- rep(c(-1, 1), length.out = n_subject)

  beta_a <- beta_b <- array(0, c(prod(dims3), n_subject, 1L))
  for (i in seq_len(n_subject)) {
    beta_a[, i, 1L] <- signs[i] * as.numeric(base) + rnorm(prod(dims3), sd = 0.03)
    beta_b[, i, 1L] <- signs[i] * as.numeric(base) + rnorm(prod(dims3), sd = 0.03)
  }
  variance <- array(0.03^2, dim(beta_a))
  region <- array(0L, dims3)
  region[3:7, 3:7, 1L] <- 1L

  result <- experimental_cancellation(
    .make_experimental_split(beta_a, variance, dims3),
    .make_experimental_split(beta_b, variance, dims3),
    delta = 0.2,
    equivalence = 0.25,
    candidate_threshold = 0.35,
    shift_radius = 1L,
    patch_radius = 1L,
    regions = as.integer(region)
  )

  center <- 5L + (5L - 1L) * dims3[1L]
  expect_gt(assay(result, "cancellation_probability")[center, 1L, 1L], 0.5)
  rescue <- assay(result, "shift_rescue_fraction")[center, 1L, 1L]
  expect_true(is.na(rescue) || rescue < 0.2)
})

test_that("experimental cancellation validates scientific and spatial contracts", {
  dims3 <- c(2L, 2L, 1L)
  beta <- array(0, c(prod(dims3), 8L, 1L))
  variance <- array(1, dim(beta))
  a <- .make_experimental_split(beta, variance, dims3)
  b <- .make_experimental_split(beta, variance, dims3)

  expect_error(experimental_cancellation(a, b, delta = 0), "delta")
  expect_error(experimental_cancellation(a, b, delta = 0.1, prevalence = 0.5), "prevalence")

  wrong_subjects <- b
  wrong_subjects$subjects[1L] <- "other"
  expect_error(experimental_cancellation(a, wrong_subjects, delta = 0.1), "subjects")

  labels <- new_gds(
    list(beta = beta, var = variance),
    space_sample_labels(paste0("v", seq_len(prod(dims3)))),
    a$subjects,
    "task"
  )
  expect_error(experimental_cancellation(labels, labels, delta = 0.1), "voxel")
})

test_that("experimental cancellation preserves packed voxel support", {
  dims3 <- c(5L, 5L, 1L)
  mask_idx <- c(7L, 8L, 9L, 12L, 13L, 14L, 17L, 18L, 19L)
  n_subject <- 24L
  beta <- array(0, c(length(mask_idx), n_subject, 1L))
  beta[, , 1L] <- rep(c(-0.8, 0.8), length.out = n_subject)[col(beta[, , 1L])]
  variance <- array(0.01, dim(beta))
  make_packed <- function() {
    new_gds(
      list(beta = beta, var = variance),
      space_voxel(dims3, diag(4), mask_idx = mask_idx, storage = "packed"),
      paste0("s", seq_len(n_subject)),
      "task"
    )
  }

  result <- experimental_cancellation(
    make_packed(), make_packed(), delta = 0.15, equivalence = 0.2,
    candidate_threshold = 1
  )
  expect_identical(space(result)$mask_idx, mask_idx)
  expect_identical(dim(assay(result, "cancellation_probability")),
                   c(length(mask_idx), 1L, 1L))
  expect_true(all(is.finite(assay(result, "cancellation_probability"))))
})
