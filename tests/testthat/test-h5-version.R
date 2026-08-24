test_that("h5 adapter validates schema version", {
  testthat::skip_if_not_installed("hdf5r")
  tf <- tempfile(fileext = ".h5")
  # Some environments may restrict creating HDF5 files; skip if so
  can_create <- try({ h5 <- hdf5r::H5File$new(tf, mode = "w"); h5$close(); TRUE }, silent = TRUE)
  if (inherits(can_create, "try-error") || isFALSE(can_create)) {
    testthat::skip("Unable to create HDF5 files in this environment")
  }
  write_min_h5 <- function(path, version) {
    if (file.exists(path)) {
      try(unlink(path), silent = TRUE)
    }
    h5 <- hdf5r::H5File$new(path, mode = "w")
    on.exit(h5$close_all())
    g <- h5$create_group("gds"); on.exit(g$close(), add = TRUE)
    g$create_dataset("version", version)
    axes <- g$create_group("axes"); on.exit(axes$close(), add = TRUE)
    axes$create_dataset("subjects", c("s1"))
    axes$create_dataset("contrasts", c("c1"))
    sp <- g$create_group("space"); on.exit(sp$close(), add = TRUE)
    sp$create_attr("type", "sample_labels")
    sl <- sp$create_group("sample_labels"); on.exit(sl$close(), add = TRUE)
    sl$create_dataset("labels", c("v1", "v2"))
    assays <- g$create_group("assays"); on.exit(assays$close(), add = TRUE)
    assays$create_dataset("beta", array(0, dim = c(2,1,1)))
  }
  # Unsupported major version should error
  write_min_h5(tf, "gds-h5/1.0")
  expect_error(gds(tf, format = "h5"), "Unsupported HDF5 schema version")
  # Supported 0.x should pass
  write_min_h5(tf, "gds-h5/0.1")
  pl <- gds(tf, format = "h5")
  expect_true(inherits(pl, "gds_plan"))
})
