test_that("fmristore adapter detects labeled volume layout and reads assays", {
  skip_if_not_installed("hdf5r")

  # Create a minimal labeled volume HDF5 file
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  header <- h5$create_group("header")
  header$create_dataset("dim", c(4L, 2L, 2L, 2L))
  header$close()
  mask <- array(c(1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L), dim = c(2, 2, 2))
  h5$create_dataset("mask", mask)
  labels <- c("cope1", "cope2")
  h5$create_dataset("labels", labels)
  data_grp <- h5$create_group("data")
  # Packed voxels length = sum(mask) = 6
  data_grp$create_dataset("cope1", c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0))
  data_grp$create_dataset("cope1_var", rep(1.0, 6))
  data_grp$create_dataset("cope2", c(7.0, 8.0, 9.0, 10.0, 11.0, 12.0))
  data_grp$create_dataset("cope2_var", rep(2.0, 6))
  data_grp$close(); h5$close()

  # Adapter detection
  name <- detect_adapter(tmp)
  expect_true(name %in% c("fmristore", "h5"))
  plan <- gds(tmp, format = name)
  res <- compute(plan)
  expect_equal(dim(assay(res, "beta")), c(6, 1, 2))
  expect_true(all(!is.na(assay(res, "var"))))
})

test_that("fmristore adapter detects latent layout and reads embeddings", {
  skip_if_not_installed("hdf5r")

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  # basis: [k, V] = [3, 5]
  basis_grp <- h5$create_group("basis")
  basis_grp$create_dataset("basis_matrix", matrix(runif(15), nrow = 3, ncol = 5))
  basis_grp$create_dataset("basis_method", "ICA")
  scans <- h5$create_group("scans")
  scans$create_group("sub-01")$create_dataset("embedding", matrix(c(1,2,3), nrow = 3, ncol = 1))
  scans$create_group("sub-02")$create_dataset("embedding", matrix(c(4,5,6), nrow = 3, ncol = 1))
  scans$close(); h5$close()

  name <- detect_adapter(tmp)
  expect_true(name %in% c("fmristore", "h5"))
  plan <- gds(tmp, format = name)
  res <- compute(plan)
  expect_equal(dim(assay(res, "beta")), c(3, 2, 1))
})
