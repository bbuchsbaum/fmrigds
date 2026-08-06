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

test_that("h5 adapter honors all block axes and preserves selected labels", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  beta <- array(seq_len(4 * 3 * 2), c(4, 3, 2))
  g <- new_gds(
    assays = list(beta = beta, var = array(1, dim(beta))),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = paste0("s", 1:3),
    contrasts = paste0("c", 1:2)
  )
  write_gds_h5(g, tmp)

  adapter <- get_adapter("h5")
  handle <- adapter$open(tmp)
  on.exit(adapter$close(handle), add = TRUE)
  observed <- adapter$read(
    handle,
    assays = "beta",
    block = list(sample = c(2L, 4L), subject = 2L, contrast = 2L)
  )$beta

  expect_equal(dim(observed), c(2L, 1L, 1L))
  expect_equal(unname(drop(observed)), beta[c(2L, 4L), 2L, 2L])
  expect_equal(dimnames(observed), list(c("roi2", "roi4"), "s2", "c2"))
})
