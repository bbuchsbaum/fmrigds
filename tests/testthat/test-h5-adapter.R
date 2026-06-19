test_that("h5 adapter roundtrip", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  g <- new_gds(
    assays = list(beta = array(1:8, c(4, 2, 1)), var = array(1, c(4, 2, 1))),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  write_gds_h5(g, tmp)

  plan <- gds(tmp)
  expect_s3_class(plan, "gds_plan")
  res <- compute(plan)
  expect_equal(assay(res, "beta")[1, 1, 1], 1)
})

test_that("h5 adapter rejects vector sources with a clear single-file error", {
  skip_if_not_installed("hdf5r")
  expect_error(
    gds(c("a.h5", "b.h5"), format = "h5"),
    "single native fmrigds .h5/.hdf5 file",
    fixed = TRUE
  )
})
