test_that("adapter registry detects tabular files", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("sample,subject,contrast,beta,var", tmp)
  expect_equal(detect_adapter(tmp), "tabular")
})

test_that("tabular adapter can extract contrast_data from long-table columns", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "sample,subject,visit,time,beta,var",
    "ROI_1,sub-01,baseline,0,0.5,0.04",
    "ROI_1,sub-01,task,1,0.8,0.04",
    "ROI_1,sub-02,baseline,0,0.4,0.05",
    "ROI_1,sub-02,task,1,0.9,0.05"
  ), tmp)

  plan <- gds(tmp, contrast_col = "visit", contrast_data_cols = "time")
  expect_equal(rownames(contrast_data(plan)), c("baseline", "task"))
  expect_equal(as.numeric(contrast_data(plan)$time), c(0, 1))
})

test_that("nifti adapter produces plan", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
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
  fmrigds:::.reset_synthetic_variance_warning()
  expect_warning(
    g <- compute(plan),
    "No variance or SE provided; using unit variance",
    fixed = TRUE
  )
  expect_s3_class(g, "gds")
  expect_true("var" %in% names(assays(g)))
  expect_true(all(assay(g, "var") == 1))
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

test_that("nifti adapter list detection rejects non-NIfTI files", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("a,b\n1,2", tmp)
  expect_false(fmrigds:::.nifti_detect(list(beta = tmp)))
})

test_that("nifti beta and se file sets are aligned by subject key", {
  aligned <- fmrigds:::.nifti_align_file_sets(
    c("sub-02_beta.nii.gz", "sub-01_beta.nii.gz"),
    c("sub-01_se.nii.gz", "sub-02_se.nii.gz")
  )
  expect_equal(basename(aligned$files_beta), c("sub-01_beta.nii.gz", "sub-02_beta.nii.gz"))
  expect_equal(basename(aligned$files_se), c("sub-01_se.nii.gz", "sub-02_se.nii.gz"))
  expect_error(
    fmrigds:::.nifti_align_file_sets("sub-01_beta.nii.gz", c("sub-01_se.nii.gz", "sub-02_se.nii.gz")),
    "same subjects"
  )
})

test_that("nifti explicit subjects bypass filename-key beta/se pairing", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  beta_files <- file.path(tmp_dir, paste0("sub-0", 1:2, "_stat-beta_sm-2.nii"))
  se_files <- file.path(tmp_dir, paste0("sub-0", 1:2, "_stat-se_sm-2.nii"))
  for (i in seq_along(beta_files)) {
    RNifti::writeNifti(array(i, dim = c(2, 2, 2)), beta_files[i])
    RNifti::writeNifti(array(i / 10, dim = c(2, 2, 2)), se_files[i])
  }

  # BIDS-aware keying: beta/se files sharing a sub-<label> entity pair on that
  # entity even though their stat-<token> infixes differ, so no explicit
  # `subject` is required and pairing succeeds.
  g_keyed <- compute(gds(nifti_source(beta = beta_files, se = se_files), format = "nifti"))
  expect_equal(subjects(g_keyed), c("sub-01", "sub-02"))
  expect_equal(as.numeric(assay(g_keyed, "beta")[1, , 1]), c(1, 2))
  expect_equal(as.numeric(assay(g_keyed, "se")[1, , 1]), c(0.1, 0.2))

  src <- nifti_source(
    beta = beta_files,
    se = se_files,
    subject = c("sub-A", "sub-B"),
    contrast = "cue"
  )
  g <- compute(gds(src, format = "nifti"))
  expect_equal(subjects(g), c("sub-A", "sub-B"))
  expect_equal(contrasts(g), "cue")
  expect_equal(as.numeric(assay(g, "beta")[1, , 1]), c(1, 2))
  expect_equal(as.numeric(assay(g, "se")[1, , 1]), c(0.1, 0.2))
})

test_that("nifti gds() subjects argument is stored for compute", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  beta_files <- file.path(tmp_dir, paste0("beta_s0", 1:2, ".nii"))
  se_files <- file.path(tmp_dir, paste0("se_s0", 1:2, ".nii"))
  for (i in seq_along(beta_files)) {
    RNifti::writeNifti(array(10 + i, dim = c(2, 2, 2)), beta_files[i])
    RNifti::writeNifti(array(i, dim = c(2, 2, 2)), se_files[i])
  }

  plan <- gds(
    list(beta = beta_files, se = se_files),
    format = "nifti",
    subjects = c("s1", "s2"),
    contrasts = "effect"
  )
  expect_equal(plan$source$source$subjects, c("s1", "s2"))
  g <- compute(plan)
  expect_equal(subjects(g), c("s1", "s2"))
  expect_equal(contrasts(g), "effect")
  expect_equal(as.numeric(assay(g, "beta")[1, , 1]), c(11, 12))
})
