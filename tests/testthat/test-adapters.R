test_that("adapter registry detects tabular files", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("sample,subject,contrast,beta,var", tmp)
  expect_equal(detect_adapter(tmp), "tabular")
})

test_that("nifti adapter roundtrip", {
  skip_if_not_installed("neuroim2")
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create test volumes using neuroim2
  arr1 <- array(runif(27), dim = c(3, 3, 3))
  arr2 <- array(runif(27), dim = c(3, 3, 3))
  f1 <- file.path(tmp_dir, "sub-01.nii")
  f2 <- file.path(tmp_dir, "sub-02.nii")

  # Create NeuroSpace for the volumes
  nspace <- neuroim2::NeuroSpace(dim = c(3L, 3L, 3L), spacing = c(1, 1, 1), origin = c(0, 0, 0))

  # Create NeuroVol objects and write
  vol1 <- neuroim2::NeuroVol(arr1, nspace)
  vol2 <- neuroim2::NeuroVol(arr2, nspace)
  neuroim2::write_vol(vol1, f1)
  neuroim2::write_vol(vol2, f2)

  plan <- gds(c(f1, f2))
  g <- compute(plan)

  expect_s3_class(g, "gds")
  expect_equal(dim(assay(g, "beta")), c(27, 2, 1))
})
