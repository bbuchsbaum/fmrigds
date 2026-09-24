.ldr_map_gds_fixture <- function(mode = "shifted", seed = 5501L,
                                 packed = FALSE) {
  fixture <- .ldr_test_fixture(
    mode, seed = seed, n_subject = 24L, noise_sd = 0.08
  )
  dims3 <- c(11L, 1L, 1L)
  n_subject <- length(fixture$subjects)
  full_a <- full_b <- matrix(0, prod(dims3), n_subject)
  full_var <- matrix(0.08^2, prod(dims3), n_subject)
  full_a[fixture$geometry$full_idx, ] <- fixture$patch$y_a
  full_b[fixture$geometry$full_idx, ] <- fixture$patch$y_b
  mask_idx <- if (packed) 1:11 else NULL
  sample_idx <- mask_idx %||% seq_len(prod(dims3))
  make_split <- function(beta) {
    beta <- array(
      beta[sample_idx, , drop = FALSE],
      c(length(sample_idx), n_subject, 1L)
    )
    variance <- array(full_var[sample_idx, , drop = FALSE], dim(beta))
    new_gds(
      assays = list(beta = beta, var = variance),
      space = space_voxel(
        dims3,
        diag(c(2, 2, 2, 1)),
        mask_idx = mask_idx,
        storage = if (packed) "packed" else "dense"
      ),
      subjects = fixture$subjects,
      contrasts = "task"
    )
  }
  list(a = make_split(full_a), b = make_split(full_b))
}

test_that("split differences recover a shrunk positive-definite correlation", {
  set.seed(5502L)
  truth <- matrix(c(
    1, 0.55, -0.2,
    0.55, 1, 0.3,
    -0.2, 0.3, 1
  ), 3L, 3L)
  n_subject <- 1200L
  standardized <- t(chol(truth)) %*%
    matrix(rnorm(3L * n_subject), 3L, n_subject)
  y_a <- standardized / 2
  y_b <- -standardized / 2
  variance <- matrix(0.5, 3L, n_subject)
  estimate <- .ldr_noise_correlation(
    y_a, y_b, variance, variance, seq_len(n_subject), shrinkage = 0
  )

  expect_equal(diag(estimate$correlation), rep(1, 3L), tolerance = 1e-12)
  expect_gt(estimate$min_eigenvalue, 0)
  expect_equal(estimate$correlation, truth, tolerance = 0.06)
})

test_that("the physical tangent basis is a least-squares shift derivative", {
  geometry <- .ldr_patch_geometry(
    center = 6L,
    dims3 = c(11L, 1L, 1L),
    affine = diag(c(2, 2, 2, 1)),
    patch_radius_mm = 8,
    shift_radius_mm = 2
  )
  template <- exp(-geometry$coordinates_mm[, 1L]^2 / 8)
  template <- template / sqrt(sum(template^2))
  tangent <- .ldr_tangent_basis(template, geometry)
  shifted <- .ldr_shift_templates(template, geometry$shift_operators)
  nonzero <- geometry$shift_distance_mm > 0
  coordinate <- geometry$shift_world[nonzero, , drop = FALSE] %*%
    tangent$world_basis
  residual <- shifted[, nonzero, drop = FALSE] - template +
    tangent$derivative %*% t(coordinate)

  expect_identical(tangent$rank, 1L)
  expect_lt(max(abs(crossprod(coordinate, t(residual)))), 1e-10)
})

test_that("the tangent fingerprint separates shift from aligned amplitude", {
  shifted <- .ldr_test_fixture(
    "shifted", seed = 5503L, n_subject = 30L
  )
  aligned <- .ldr_test_fixture(
    "aligned", seed = 5503L, n_subject = 30L
  )
  polarity <- .ldr_test_fixture(
    "polarity", seed = 5503L, n_subject = 30L
  )
  fit <- function(x) {
    .ldr_tangent_center(
      x$patch$y_a, x$patch$y_b,
      x$patch$var_a, x$patch$var_b,
      x$geometry, x$subjects,
      folds = 3L,
      covariance_shrinkage = 0.5
    )
  }
  shifted_fit <- fit(shifted)
  aligned_fit <- fit(aligned)
  polarity_fit <- fit(polarity)

  expect_gt(shifted_fit$activation_t, 2)
  expect_gt(shifted_fit$displacement_t, aligned_fit$displacement_t + 4)
  expect_gt(shifted_fit$rho_shift, 0.1)
  expect_lt(aligned_fit$rho_shift, 0.01)
  expect_lt(polarity_fit$rho_shift, 0.01)
  expect_lte(shifted_fit$shift_rms_mm, 2 + 1e-12)
})

test_that("the tangent screen recovers all three physical derivative axes", {
  set.seed(5508L)
  dims3 <- c(9L, 9L, 9L)
  center <- 365L
  geometry <- .ldr_patch_geometry(
    center, dims3, diag(c(2, 2, 2, 1)),
    patch_radius_mm = 4, shift_radius_mm = 2
  )
  feature <- numeric(length(geometry$sample_idx))
  feature[geometry$center] <- 1.2
  feature[geometry$distance_mm == 2] <- -0.35
  n_subject <- 28L
  subject_ids <- paste0("s", seq_len(n_subject))
  shift_index <- rep(
    seq_along(geometry$shift_operators), length.out = n_subject
  )
  simulate <- function(shifted) {
    y_a <- y_b <- matrix(0, length(feature), n_subject)
    for (i in seq_len(n_subject)) {
      truth <- if (shifted) {
        geometry$shift_operators[[shift_index[i]]] %*% feature
      } else feature
      y_a[, i] <- truth + rnorm(length(feature), sd = 0.08)
      y_b[, i] <- truth + rnorm(length(feature), sd = 0.08)
    }
    variance <- matrix(0.08^2, length(feature), n_subject)
    .ldr_tangent_center(
      y_a, y_b, variance, variance, geometry, subject_ids,
      folds = 4L, covariance_shrinkage = 0.5
    )
  }
  shifted <- simulate(TRUE)
  set.seed(5508L)
  aligned <- simulate(FALSE)

  expect_true(all(vapply(
    shifted$folds, `[[`, integer(1L), "tangent_rank"
  ) == 3L))
  expect_gt(shifted$displacement_t, aligned$displacement_t + 6)
  expect_gt(shifted$rho_shift, 0.5)
  expect_lt(aligned$rho_shift, 0.01)
  expect_true(is.na(shifted$shift_rms_mm))
})

test_that("tangent screening is invariant to subject storage order", {
  fixture <- .ldr_test_fixture(
    "shifted", seed = 5504L, n_subject = 24L
  )
  run <- function(x) {
    .ldr_tangent_center(
      x$patch$y_a, x$patch$y_b,
      x$patch$var_a, x$patch$var_b,
      x$geometry, x$subjects,
      folds = 3L,
      covariance_shrinkage = 0.5
    )
  }
  baseline <- run(fixture)
  order <- sample(seq_along(fixture$subjects))
  fixture$subjects <- fixture$subjects[order]
  for (name in c("y_a", "y_b", "var_a", "var_b")) {
    fixture$patch[[name]] <- fixture$patch[[name]][, order, drop = FALSE]
  }
  fixture$patch$subjects <- fixture$patch$subjects[order]
  reordered <- run(fixture)

  expect_equal(reordered$activation_t, baseline$activation_t,
               tolerance = 1e-10)
  expect_equal(reordered$displacement_t, baseline$displacement_t,
               tolerance = 1e-10)
  expect_equal(reordered$mni_loss, baseline$mni_loss, tolerance = 1e-10)
})

test_that("whole-search calibration preserves RNG and max-p invariants", {
  fixture <- .ldr_map_gds_fixture("shifted", seed = 5505L)
  beta_a <- .ldr_matrix_contrast(fixture$a, "beta", 1L)
  beta_b <- .ldr_matrix_contrast(fixture$b, "beta", 1L)
  var_a <- .ldr_matrix_contrast(fixture$a, "var", 1L)
  var_b <- .ldr_matrix_contrast(fixture$b, "var", 1L)
  configurations <- .ldr_search_configurations(c(8, 10), 2)
  units <- .ldr_search_units(
    centers = 5:7,
    configurations = configurations,
    dims3 = c(11L, 1L, 1L),
    affine = diag(c(2, 2, 2, 1))
  )$units
  observed <- .ldr_tangent_scan(
    beta_a, beta_b, var_a, var_b, units,
    subjects(fixture$a), folds = 3L,
    covariance_shrinkage = 0.5
  )
  set.seed(912L)
  before <- .Random.seed
  calibrated <- .ldr_tangent_calibrate(
    beta_a, beta_b, var_a, var_b, units,
    subjects(fixture$a), folds = 3L,
    covariance_shrinkage = 0.5,
    observed = observed,
    n_resamples = 5L,
    seed = 73L
  )

  expect_identical(.Random.seed, before)
  expect_equal(calibrated$ldr_fwer_p, pmax(
    calibrated$activation_fwer_p,
    calibrated$displacement_fwer_p
  ))
  expect_true(all(
    calibrated$activation_fwer_p + 1e-12 >=
      calibrated$pointwise_p$activation,
    na.rm = TRUE
  ))
  expect_true(all(
    calibrated$displacement_fwer_p + 1e-12 >=
      calibrated$pointwise_p$displacement,
    na.rm = TRUE
  ))
  expect_length(calibrated$null_maximum$activation, 5L)
  expect_length(calibrated$null_maximum$displacement, 5L)
})

test_that("corrected map flags shifted rescue but rejects aligned alternatives", {
  run <- function(mode) {
    fixture <- .ldr_map_gds_fixture(mode, seed = 5401L)
    local_displacement_rescue_map(
      fixture$a,
      fixture$b,
      centers = 5:7,
      patch_radius_mm = 8,
      shift_radius_mm = 2,
      folds = 3L,
      min_subjects = 12L,
      covariance_shrinkage = 0.5,
      iterations = 8L,
      n_resamples = 19L,
      seed = 11L,
      alpha = 0.05,
      max_refits = 3L
    )
  }
  shifted <- run("shifted")
  aligned <- run("aligned")
  polarity <- run("polarity")
  null <- run("null")

  expect_gt(assay(shifted, "ldr_flag")[6L, 1L, 1L], 0)
  expect_equal(sum(assay(aligned, "ldr_flag") > 0), 0)
  expect_equal(sum(assay(polarity, "ldr_flag") > 0), 0)
  expect_equal(sum(assay(null, "ldr_flag") > 0), 0)
})

test_that("exact prevalence prevents subgroup activation from becoming rescue", {
  mixed <- .ldr_map_gds_fixture("shifted", seed = 6101L)
  inactive <- .ldr_map_gds_fixture("null", seed = 6102L)
  mixed$a$assays$beta[, 13:24, ] <- inactive$a$assays$beta[, 13:24, ]
  mixed$b$assays$beta[, 13:24, ] <- inactive$b$assays$beta[, 13:24, ]
  result <- local_displacement_rescue_map(
    mixed$a,
    mixed$b,
    centers = 5:7,
    patch_radius_mm = 8,
    shift_radius_mm = 2,
    folds = 3L,
    min_subjects = 12L,
    covariance_shrinkage = 0.5,
    iterations = 8L,
    n_resamples = 19L,
    seed = 61L,
    alpha = 0.05,
    min_prevalence = 0.7,
    max_refits = 3L
  )

  expect_equal(metadata(result)$ldr$admitted_candidates, 1L)
  expect_lt(metadata(result)$ldr$exact_refits[[1L]]$prevalence, 0.7)
  expect_false(metadata(result)$ldr$exact_refits[[1L]]$confirmed)
  expect_equal(sum(assay(result, "ldr_flag") > 0), 0)
})

test_that("a broad aligned feature is not attributed to displacement", {
  fixture <- .ldr_map_gds_fixture("null", seed = 6201L)
  broad <- 0.8 * exp(-((1:11) - 6)^2 / (2 * 2.5^2))
  fixture$a$assays$beta[, , 1L] <-
    fixture$a$assays$beta[, , 1L] + broad
  fixture$b$assays$beta[, , 1L] <-
    fixture$b$assays$beta[, , 1L] + broad
  result <- local_displacement_rescue_map(
    fixture$a,
    fixture$b,
    centers = 5:7,
    patch_radius_mm = 8,
    shift_radius_mm = 2,
    folds = 3L,
    min_subjects = 12L,
    covariance_shrinkage = 0.5,
    iterations = 8L,
    n_resamples = 19L,
    seed = 62L,
    alpha = 0.05,
    max_refits = 3L
  )

  expect_true(all(assay(result, "ordinary_fwer_p")[5:7, 1L, 1L] <= 0.05))
  expect_equal(metadata(result)$ldr$admitted_candidates, 0L)
  expect_equal(sum(assay(result, "ldr_flag") > 0), 0)
})

test_that("spatially correlated null maps remain unflagged", {
  set.seed(6301L)
  n_subject <- 24L
  n_voxel <- 11L
  covariance <- 0.08^2 * 0.65^abs(outer(
    seq_len(n_voxel), seq_len(n_voxel), "-"
  ))
  root <- t(chol(covariance))
  beta_a <- root %*% matrix(rnorm(n_voxel * n_subject), n_voxel)
  beta_b <- root %*% matrix(rnorm(n_voxel * n_subject), n_voxel)
  variance <- matrix(0.08^2, n_voxel, n_subject)
  make_split <- function(beta) {
    new_gds(
      assays = list(
        beta = array(beta, c(n_voxel, n_subject, 1L)),
        var = array(variance, c(n_voxel, n_subject, 1L))
      ),
      space = space_voxel(
        c(11L, 1L, 1L), diag(c(2, 2, 2, 1))
      ),
      subjects = paste0("s", seq_len(n_subject)),
      contrasts = "task"
    )
  }
  result <- local_displacement_rescue_map(
    make_split(beta_a),
    make_split(beta_b),
    centers = 5:7,
    patch_radius_mm = 8,
    shift_radius_mm = 2,
    folds = 3L,
    min_subjects = 12L,
    covariance_shrinkage = 0.25,
    iterations = 4L,
    n_resamples = 19L,
    seed = 63L,
    alpha = 0.05,
    max_refits = 3L
  )

  expect_gte(min(assay(result, "ldr_fwer_p"), na.rm = TRUE), 0.05)
  expect_equal(metadata(result)$ldr$admitted_candidates, 0L)
  expect_equal(sum(assay(result, "ldr_flag") > 0), 0)
})

test_that("the map API records search scope and fails closed", {
  fixture <- .ldr_map_gds_fixture("shifted", seed = 5507L, packed = TRUE)
  result <- local_displacement_rescue_map(
    fixture$a,
    fixture$b,
    centers = NULL,
    patch_radius_mm = 8,
    shift_radius_mm = 2,
    folds = 3L,
    min_subjects = 12L
  )
  expected <- c(
    "ldr_fwer_p", "ldr_flag", "ldr_activation_fwer_p",
    "ldr_displacement_fwer_p", "ordinary_fwer_p",
    "tangent_activation_t", "tangent_displacement_t", "ordinary_t",
    "tangent_mni_loss", "tangent_shift_rms_mm", "rho_shift",
    "mni_loss", "shift_rms_mm"
  )
  expect_identical(names(assays(result)), expected)
  expect_identical(
    metadata(result)$ldr$calibration_scope,
    "all-eligible-active-centers"
  )
  expect_equal(metadata(result)$ldr$requested_center_count, 11L)
  expect_equal(sum(is.finite(assay(result, "tangent_activation_t"))), 3L)
  expect_true(all(assay(result, "ldr_flag") == 0))
  expect_true(all(is.na(assay(result, "ldr_fwer_p"))))
  expect_true(validate_provenance_graph(result))
  provenance <- metadata(result)$provenance
  activity <- provenance$graph[[length(provenance$graph)]]
  expect_identical(activity$op, "local_displacement_rescue_map")
  expect_length(activity$inputs, 2L)

  expect_error(
    local_displacement_rescue_map(
      fixture$a, fixture$b, centers = 1L,
      patch_radius_mm = 8, shift_radius_mm = 2
    ),
    "no center with a complete eligible patch"
  )
  expect_error(
    local_displacement_rescue_map(
      fixture$a, fixture$b, centers = 6L,
      patch_radius_mm = 2, shift_radius_mm = 2
    ),
    "No patch/shift configuration"
  )
  expect_error(
    local_displacement_rescue_map(
      fixture$a, fixture$b, centers = 6L,
      patch_radius_mm = 8, shift_radius_mm = 2,
      min_prevalence = 0
    ),
    "min_prevalence"
  )
})

test_that("a finite exact-refit cap is conservative and receipted", {
  fixture <- .ldr_map_gds_fixture("shifted", seed = 5401L)
  result <- local_displacement_rescue_map(
    fixture$a,
    fixture$b,
    centers = 5:7,
    patch_radius_mm = 8,
    shift_radius_mm = 2,
    folds = 3L,
    min_subjects = 12L,
    covariance_shrinkage = 0.5,
    iterations = 4L,
    n_resamples = 7L,
    seed = 11L,
    alpha = 0.5,
    max_refits = 1L
  )

  expect_gt(metadata(result)$ldr$admitted_candidates, 1L)
  expect_identical(metadata(result)$ldr$completed_refits, 1L)
  expect_true(metadata(result)$ldr$refits_truncated)
  expect_length(metadata(result)$ldr$exact_refits, 1L)
})
