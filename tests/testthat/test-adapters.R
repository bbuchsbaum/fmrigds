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
