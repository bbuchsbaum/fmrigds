# Regression tests for issues reported at github.com/bbuchsbaum/fmrigds
# #1 intercept-only OLS without col_data, #6 affine preservation on NIfTI
# export, #7 per-voxel NaN handling in the ols:voxelwise reducer.

test_that("#1 one_sample()/group_ols(~1) run without col_data", {
  set.seed(1)
  n_samples <- 3L; n_subj <- 8L
  beta <- array(rnorm(n_samples * n_subj, mean = 1.5), dim = c(n_samples, n_subj, 1L))
  g <- new_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space_sample_labels(paste0("v", seq_len(n_samples))),
    paste0("s", seq_len(n_subj)),
    "c1"
  )

  fit <- one_sample(g) |> compute()
  expect_true("coef:(Intercept)" %in% names(assays(fit)))
  # Intercept estimate equals the across-subject mean per sample.
  expect_equal(
    as.numeric(assay(fit, "coef:(Intercept)")[, 1, 1]),
    as.numeric(apply(beta[, , 1], 1, mean)),
    tolerance = 1e-8
  )

  # group_ols(~1) should likewise build the trivial design itself.
  fit2 <- group_ols(g, ~ 1) |> compute()
  expect_true("t_coef:(Intercept)" %in% names(assays(fit2)))
})

test_that("#7 ols:voxelwise applies per-voxel listwise deletion for NaN subjects", {
  set.seed(2)
  n_samples <- 4L; n_subj <- 6L
  beta <- array(rnorm(n_samples * n_subj, mean = 2), dim = c(n_samples, n_subj, 1L))
  beta[1L, 1:5, 1L] <- NaN          # sample 1: only 1 finite subject -> unestimable (NA)
  beta[2L, 1L, 1L]  <- NaN          # sample 2: 5 finite subjects -> estimated from subset
  # sample 3 + 4: all finite
  poison <- beta
  poison[, 6L, ] <- NaN             # subject 6 entirely NaN -> must NOT poison whole map

  g <- new_gds(
    list(beta = poison, var = array(1, dim = dim(poison))),
    space_sample_labels(paste0("v", seq_len(n_samples))),
    paste0("s", seq_len(n_subj)),
    "c1"
  )

  expect_warning(fit <- one_sample(g) |> compute(), "listwise deletion")
  co <- as.numeric(assay(fit, "coef:(Intercept)")[, 1, 1])
  n_obs <- as.numeric(assay(fit, "n_obs")[, 1, 1])

  # An all-NaN subject no longer poisons every voxel.
  expect_true(all(is.finite(co[3:4])))
  # Sample 1 has <2 finite obs (subjects 1-5 NaN, subject 6 NaN) -> NA.
  expect_true(is.na(co[1L]))
  # Sample 2: subjects 1 and 6 NaN -> 4 finite -> estimated.
  expect_true(is.finite(co[2L]))
  expect_equal(n_obs, c(0, 4, 5, 5))
  # Estimated intercept equals the mean of the finite subjects.
  expect_equal(co[2L], mean(poison[2L, , 1L], na.rm = TRUE), tolerance = 1e-8)
})

test_that("#6 write_nifti_assays preserves the spatial affine", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("RNifti")

  aff <- matrix(c(-2, 0, 0, 90,
                   0, 2, 0, -126,
                   0, 0, 2, -72,
                   0, 0, 0, 1), nrow = 4L, byrow = TRUE)
  arr <- array(as.numeric(1:8), dim = c(8L, 1L, 1L))
  g <- new_gds(
    list("t_coef:x" = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2L, 2L, 2L), aff),
    "meta",
    "c1"
  )

  td <- tempfile("nifti-affine-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  manifest <- write_nifti_assays(g, out_dir = td, assays = "t_coef:x")
  path <- manifest$path[manifest$written][1L]
  expect_true(file.exists(path))

  v <- neuroim2::read_vol(path)
  sp <- neuroim2::space(v)
  expect_equal(as.numeric(neuroim2::spacing(sp)), c(2, 2, 2), tolerance = 1e-4)
  expect_equal(as.numeric(neuroim2::origin(sp)), c(90, -126, -72), tolerance = 1e-4)
  expect_equal(unname(neuroim2::trans(sp)), unname(aff), tolerance = 1e-4)
  # scl_slope must be a usable scaling (1), not 0.
  expect_equal(RNifti::niftiHeader(path)$scl_slope, 1)
  # Data values survive the round-trip.
  expect_equal(as.numeric(as.array(v)), as.numeric(1:8), tolerance = 1e-4)
})

test_that("#5 variance-weighted reducers refuse a synthetic unit-variance placeholder", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("RNifti")

  td <- tempfile("syn-var-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  files <- file.path(td, paste0("sub-", 1:3, ".nii"))
  for (i in 1:3) RNifti::writeNifti(array(as.numeric(i), c(2, 2, 2)), files[i])

  suppressWarnings(
    g <- gds_from_scalar_maps(files, subject = paste0("s", 1:3), contrast = "metric")
  )
  expect_true(isTRUE(g$metadata$synthetic_var))

  # Variance-weighted reducers are refused with an actionable error.
  expect_error(reduce(g, method = "random"), "synthetic")
  expect_error(reduce(g, method = "fixed"), "synthetic")
  expect_error(reduce(g, method = "meta:fe"), "synthetic")

  # The unweighted OLS path is allowed.
  expect_s3_class(group_ols(g, ~ 1), "gds_plan")
})

test_that("#6 write_out(format = 'nifti') preserves the spatial affine", {
  skip_if_not_installed("neuroim2")

  aff <- matrix(c(3, 0, 0, 10,
                  0, 3, 0, 20,
                  0, 0, 3, 30,
                  0, 0, 0, 1), nrow = 4L, byrow = TRUE)
  arr <- array(as.numeric(1:8), dim = c(8L, 1L, 1L))
  g <- new_gds(
    list(beta = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2L, 2L, 2L), aff),
    "meta",
    "c1"
  )

  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path), add = TRUE)
  .write_gds_nifti(g, path, options = list(stat = "beta"))

  v <- neuroim2::read_vol(path)
  expect_equal(as.numeric(neuroim2::spacing(neuroim2::space(v))), c(3, 3, 3), tolerance = 1e-4)
  expect_equal(as.numeric(neuroim2::origin(neuroim2::space(v))), c(10, 20, 30), tolerance = 1e-4)
})
