skip_if_not_installed("neuroim2")

.write_t_maps <- function(dir, n_subjects = 4L, dims = c(3L, 3L, 2L), seed = 5) {
  set.seed(seed)
  sp <- neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2))
  vals <- lapply(seq_len(n_subjects), function(i) array(stats::rnorm(prod(dims)), dims))
  files <- file.path(dir, sprintf("sub-%02d_stat-t.nii.gz", seq_len(n_subjects)))
  for (i in seq_len(n_subjects)) neuroim2::write_vol(neuroim2::NeuroVol(vals[[i]], sp), files[i])
  list(files = files, vals = vals)
}

test_that("NIfTI t maps with df give t and df assays", {
  dir <- withr::local_tempdir()
  m <- .write_t_maps(dir)
  df <- c(20, 22, 25, 30)
  g <- compute(gds(list(t = m$files, df = df), format = "nifti"))
  expect_setequal(names(assays(g)), c("t", "df"))
  t_arr <- assay(g, "t")
  df_arr <- assay(g, "df")
  for (i in seq_along(m$files)) {
    expect_equal(as.numeric(t_arr[, i, 1]), as.numeric(m$vals[[i]]), tolerance = 1e-6)
    expect_true(all(df_arr[, i, 1] == df[i]))
  }

  # One df value is used for every subject; nifti_source() builds the same source.
  g1 <- compute(gds(nifti_source(t = m$files, df = 18), format = "nifti"))
  expect_true(all(assay(g1, "df") == 18))
  expect_equal(assay(g1, "t"), t_arr)
})

test_that("NIfTI t maps combine across subjects (vs Stouffer by hand)", {
  dir <- withr::local_tempdir()
  m <- .write_t_maps(dir)
  df <- 20
  src <- gds(list(t = m$files, df = df), format = "nifti")
  fit <- compute(reduce(src, method = "combine:stouffer"))
  # NIfTI stores float32, so the reference uses the t values as read.
  tmat <- assay(compute(src), "t")[, , 1]
  z <- stats::qnorm(stats::pt(tmat, df))
  expect_equal(as.numeric(assay(fit, "z_g")), rowSums(z) / sqrt(ncol(z)), tolerance = 1e-10)
})

test_that("NIfTI t sources are validated", {
  dir <- withr::local_tempdir()
  m <- .write_t_maps(dir)
  expect_error(gds(list(t = m$files), format = "nifti"), "requires `df`")
  expect_error(gds(list(t = m$files, df = c(10, 12)), format = "nifti"), "length 1 or one value per")
  expect_error(gds(list(t = m$files, df = 0), format = "nifti"), "positive")
  expect_error(gds(list(t = m$files, df = 10, beta = m$files), format = "nifti"),
               "cannot be combined")
})
