test_that("adapter registry detects tabular files", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("sample,subject,contrast,beta,var", tmp)
  expect_equal(detect_adapter(tmp), "tabular")
})

test_that("nifti adapter produces plan", {
  skip_if_not_installed("RNifti")
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  arr1 <- array(runif(27), dim = c(3, 3, 3))
  arr2 <- array(runif(27), dim = c(3, 3, 3))
  f1 <- file.path(tmp_dir, "sub-01.nii")
  f2 <- file.path(tmp_dir, "sub-02.nii")

  RNifti::writeNifti(arr1, f1)
  RNifti::writeNifti(arr2, f2)

  plan <- gds(c(f1, f2))
  expect_s3_class(plan, "gds_plan")
  expect_error(compute(plan), "If 'beta' is provided", fixed = TRUE)
})

test_that("nifti adapter flattens 5D non-spatial axes into contrasts", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
  ns <- asNamespace("neuroim2")
  if (!exists("read_hyper_vec", envir = ns, mode = "function", inherits = FALSE)) {
    skip("neuroim2::read_hyper_vec() not available")
  }

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  arr1 <- array(seq_len(2 * 2 * 2 * 2 * 3), dim = c(2, 2, 2, 2, 3))
  arr2 <- arr1 + 1000
  f1 <- file.path(tmp_dir, "sub-01_beta.nii")
  f2 <- file.path(tmp_dir, "sub-02_beta.nii")
  RNifti::writeNifti(arr1, f1)
  RNifti::writeNifti(arr2, f2)

  plan <- gds(c(f1, f2), format = "nifti")
  expect_equal(unname(plan$source$probe$dims), c(8L, 2L, 6L))

  beta <- compute(plan, assays = "beta")$beta
  expect_equal(dim(beta), c(8L, 2L, 6L))

  expected1 <- matrix(as.vector(arr1), nrow = 8, ncol = 6)
  expected2 <- matrix(as.vector(arr2), nrow = 8, ncol = 6)
  expect_equal(beta[, 1, ], expected1)
  expect_equal(beta[, 2, ], expected2)
})

test_that("space_from_nifti handles 5D nifti headers", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
  ns <- asNamespace("neuroim2")
  if (!exists("read_hyper_vec", envir = ns, mode = "function", inherits = FALSE)) {
    skip("neuroim2::read_hyper_vec() not available")
  }

  tmp <- tempfile(fileext = ".nii")
  on.exit(unlink(tmp), add = TRUE)
  arr <- array(runif(2 * 3 * 4 * 2 * 2), dim = c(2, 3, 4, 2, 2))
  RNifti::writeNifti(arr, tmp)

  sp <- space_from_nifti(tmp)
  expect_s3_class(sp, "space_voxel")
  expect_equal(as.integer(sp$dim), c(2L, 3L, 4L))
})
