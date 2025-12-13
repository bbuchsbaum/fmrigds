test_that("compute sink = 'h5' writes a gds store", {
  skip_if_not_installed("hdf5r")

  tmp_src <- tempfile(fileext = ".h5")
  tmp_out <- tempfile(fileext = ".h5")
  on.exit(unlink(c(tmp_src, tmp_out)), add = TRUE)

  g <- new_gds(
    assays = list(beta = array(1:8, c(4, 2, 1)), var = array(1, c(4, 2, 1))),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  write_gds_h5(g, tmp_src)

  plan <- gds(tmp_src)
  res <- compute(plan, sink = "h5", path = tmp_out)
  expect_true(file.exists(tmp_out))

  res2 <- compute(gds(tmp_out))
  expect_equal(assay(res2, "beta"), assay(res, "beta"))

  h5 <- hdf5r::H5File$new(tmp_out, mode = "r")
  on.exit(h5$close(), add = TRUE)
  log_ds <- h5$open("/gds/provenance/log")
  on.exit(log_ds$close(), add = TRUE)
  log_entries <- log_ds$read()
  expect_true(is.character(log_entries))
})

test_that("write_out node materialises an h5 file", {
  skip_if_not_installed("hdf5r")

  tmp_src <- tempfile(fileext = ".h5")
  tmp_out <- tempfile(fileext = ".h5")
  on.exit(unlink(c(tmp_src, tmp_out)), add = TRUE)

  g <- new_gds(
    assays = list(beta = array(rnorm(8), c(4, 2, 1)), var = array(1, c(4, 2, 1))),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  write_gds_h5(g, tmp_src)

  plan <- write_out(gds(tmp_src), path = tmp_out, format = "h5")
  res <- compute(plan)

  expect_true(file.exists(tmp_out))
  saved <- compute(gds(tmp_out))
  expect_equal(assay(saved, "beta"), assay(res, "beta"))
})

test_that("write_out csv export works", {
  skip_if_not_installed("hdf5r")

  tmp_src <- tempfile(fileext = ".h5")
  tmp_out <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp_src, tmp_out)), add = TRUE)

  g <- new_gds(
    assays = list(beta = array(1:8, c(4, 2, 1)), var = array(1, c(4, 2, 1))),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  write_gds_h5(g, tmp_src)

  plan <- write_out(gds(tmp_src), path = tmp_out, format = "csv", options = list(stats = "beta"))
  compute(plan)

  expect_true(file.exists(tmp_out))
  df <- utils::read.csv(tmp_out, stringsAsFactors = FALSE)
  expect_true(all(c("sample", "subject", "contrast", "beta") %in% names(df)))
  expect_equal(nrow(df), 8)
})

test_that("write_out parquet export works when arrow available", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("arrow")

  tmp_src <- tempfile(fileext = ".h5")
  tmp_out <- tempfile(fileext = ".parquet")
  on.exit(unlink(c(tmp_src, tmp_out)), add = TRUE)

  g <- new_gds(
    assays = list(beta = array(1:8, c(4, 2, 1)), var = array(1, c(4, 2, 1))),
    space = space_parcels(paste0("roi", 1:4)),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  write_gds_h5(g, tmp_src)

  plan <- write_out(gds(tmp_src), path = tmp_out, format = "parquet", options = list(stats = "beta"))
  compute(plan)

  expect_true(file.exists(tmp_out))
  tbl <- arrow::read_parquet(tmp_out)
  expect_true("beta" %in% names(tbl))
  expect_equal(nrow(tbl), 8)
})

test_that("write_out nifti export works when RNifti available", {
  skip_if_not_installed("RNifti")

  tmp_out <- tempfile(fileext = ".nii")
  on.exit(unlink(tmp_out), add = TRUE)

  dims <- c(2, 2, 2)
  beta <- array(1:8, c(8, 1, 1))
  var <- array(1, c(8, 1, 1))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_voxel(dim = dims, affine = diag(4), storage = "dense"),
    subjects = "s1",
    contrasts = "c1"
  )

  fmrigds:::.write_gds_nifti(g, tmp_out, options = list(stat = "beta"))

  expect_true(file.exists(tmp_out))
  vol <- RNifti::readNifti(tmp_out)
  expect_equal(dim(vol)[1:3], dims)
  expect_equal(vol[1], 1)
})
