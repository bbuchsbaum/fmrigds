.ldr_fixture_gds <- function(mode = "shifted", seed = 5201L,
                             packed = FALSE) {
  fixture <- .ldr_test_fixture(mode, seed = seed, n_subject = 24L)
  dims3 <- c(11L, 1L, 1L)
  n_subject <- length(fixture$subjects)
  full_a <- full_b <- matrix(0, nrow = prod(dims3), ncol = n_subject)
  full_var <- matrix(0.08^2, nrow = prod(dims3), ncol = n_subject)
  full_a[fixture$geometry$full_idx, ] <- fixture$patch$y_a
  full_b[fixture$geometry$full_idx, ] <- fixture$patch$y_b
  mask_idx <- if (packed) 2:10 else NULL
  sample_idx <- mask_idx %||% seq_len(prod(dims3))
  make_split <- function(beta) {
    beta <- array(beta[sample_idx, , drop = FALSE],
                  c(length(sample_idx), n_subject, 1L))
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

.run_ldr_fixture <- function(x, n_resamples = 0L, seed = 1L) {
  local_displacement_rescue(
    x$a,
    x$b,
    center = 6L,
    patch_radius_mm = 8,
    shift_radius_mm = 2,
    tau_grid_mm = c(1, 2),
    folds = 3L,
    min_subjects = 12L,
    iterations = 8L,
    n_resamples = n_resamples,
    seed = seed
  )
}

test_that("experimental LDR exposes a narrow one-center reference result", {
  result <- .run_ldr_fixture(.ldr_fixture_gds("shifted", seed = 5202L))

  expect_identical(
    names(assays(result)),
    c(
      "ldr_p",
      "mni_loss",
      "shift_rms_mm"
    )
  )
  finite <- vapply(assays(result), function(value) sum(is.finite(value)),
                   integer(1L))
  expect_identical(unname(finite), c(0L, 1L, 1L))
  expect_true(is.na(assay(result, "ldr_p")[6L, 1L, 1L]))
  expect_gt(assay(result, "mni_loss")[6L, 1L, 1L], 0)
  expect_gt(assay(result, "shift_rms_mm")[6L, 1L, 1L], 0)

  receipt <- metadata(result)$ldr
  expect_false(receipt$calibrated)
  expect_gt(receipt$activation_score, 0)
  expect_gt(receipt$displacement_score, 0)
  expect_identical(receipt$scope, "one-center-one-contrast")
  expect_identical(
    receipt$comparators,
    c("no-feature", "aligned-positive", "aligned-signed-amplitude")
  )
  expect_length(receipt$fold_receipts, 3L)
})

test_that("LDR preserves packed voxel coordinates and both provenances", {
  result <- .run_ldr_fixture(
    .ldr_fixture_gds("shifted", seed = 5203L, packed = TRUE)
  )

  expect_identical(space(result)$mask_idx, 2:10)
  expect_true(is.finite(assay(result, "mni_loss")[5L, 1L, 1L]))
  expect_true(validate_provenance_graph(result))
  provenance <- metadata(result)$provenance
  activity <- provenance$graph[[length(provenance$graph)]]
  expect_identical(activity$op, "local_displacement_rescue")
  expect_length(activity$inputs, 2L)
})

test_that("LDR exposes calibrated ROI p only when resampling is requested", {
  result <- .run_ldr_fixture(
    .ldr_fixture_gds("shifted", seed = 5205L),
    n_resamples = 1L,
    seed = 27L
  )

  p_value <- assay(result, "ldr_p")[6L, 1L, 1L]
  expect_true(is.finite(p_value))
  expect_true(p_value >= 0.5 && p_value <= 1)
  expect_true(metadata(result)$ldr$calibrated)
  expect_identical(metadata(result)$ldr$calibration_scope, "prespecified-ROI")
  expect_equal(metadata(result)$ldr$calibration$p_ldr, p_value)
  expect_identical(metadata(result)$ldr$calibration$n_resamples, 1L)
})

test_that("LDR API rejects scientifically invalid scope", {
  x <- .ldr_fixture_gds("shifted", seed = 5204L)

  expect_error(
    local_displacement_rescue(
      x$a, x$b, center = 6L,
      patch_radius_mm = 3, shift_radius_mm = 2
    ),
    "at least twice"
  )
  expect_error(
    local_displacement_rescue(
      x$a, x$b, center = 1L,
      patch_radius_mm = 8, shift_radius_mm = 2
    ),
    "crosses the image boundary"
  )
  expect_error(
    local_displacement_rescue(
      x$a, x$b, center = 6L,
      patch_radius_mm = 8, shift_radius_mm = 2,
      tau_grid_mm = 3
    ),
    "no larger"
  )
  synthetic <- x$a
  synthetic$metadata$synthetic_var <- TRUE
  expect_error(
    local_displacement_rescue(
      synthetic, x$b, center = 6L,
      patch_radius_mm = 8, shift_radius_mm = 2
    ),
    "synthetic variance"
  )
  expect_error(
    local_displacement_rescue(
      x$a, x$b, center = 6L,
      patch_radius_mm = 8, shift_radius_mm = 2,
      n_resamples = -1
    ),
    "nonnegative integer"
  )
})
