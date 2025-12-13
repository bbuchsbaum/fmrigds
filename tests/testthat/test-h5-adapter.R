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
