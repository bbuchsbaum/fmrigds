# Tests for nifti-utils.R: find_maps()

test_that("find_maps discovers files with pattern", {
  td <- tempfile("findmaps-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))

  # Create fake NIfTI files
  dir.create(file.path(td, "sub-01"))
  dir.create(file.path(td, "sub-02"))
  file.create(file.path(td, "sub-01", "stat.nii.gz"))
  file.create(file.path(td, "sub-02", "stat.nii.gz"))
  file.create(file.path(td, "sub-01", "other.txt"))

  result <- find_maps(td, pattern = "\\.nii(\\.gz)?$")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(grepl("nii", result$file)))
  expect_true(all(is.na(result$subject)))
  expect_true(all(is.na(result$contrast)))
})

test_that("find_maps with subject_fun and contrast_fun", {
  td <- tempfile("findmaps-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))

  dir.create(file.path(td, "sub-01"))
  dir.create(file.path(td, "sub-02"))
  file.create(file.path(td, "sub-01", "cope.nii"))
  file.create(file.path(td, "sub-02", "cope.nii"))

  sfun <- function(files) sub(".*(sub-[0-9]+).*", "\\1", files)
  cfun <- function(files) rep("cope1", length(files))

  result <- find_maps(td, subject_fun = sfun, contrast_fun = cfun)
  expect_equal(sort(result$subject), c("sub-01", "sub-02"))
  expect_equal(unique(result$contrast), "cope1")
})

test_that("find_maps returns empty data.frame when no files match", {
  td <- tempfile("findmaps-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))

  result <- find_maps(td, pattern = "\\.nii$")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("find_maps errors on missing root", {
  expect_error(find_maps("/nonexistent/dir"), "does not exist")
})

test_that("gds_from_nifti_maps accepts beta-only maps with unit variance", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")

  td <- tempfile("niftimaps-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  arr1 <- array(1, dim = c(2, 2, 2))
  arr2 <- array(2, dim = c(2, 2, 2))
  f1 <- file.path(td, "sub-01_naive_diag_mean.nii")
  f2 <- file.path(td, "sub-02_naive_diag_mean.nii")
  RNifti::writeNifti(arr1, f1)
  RNifti::writeNifti(arr2, f2)

  maps <- data.frame(
    file = c(f1, f2),
    subject = c("1001", "1002"),
    contrast = "naive_diag_mean",
    stringsAsFactors = FALSE
  )

  fmrigds:::.reset_synthetic_variance_warning()
  expect_warning(
    g <- gds_from_nifti_maps(maps, mask = NULL),
    "No variance or SE provided; using unit variance",
    fixed = TRUE
  )

  expect_s3_class(g, "gds")
  expect_equal(subjects(g), c("1001", "1002"))
  expect_equal(contrasts(g), "naive_diag_mean")
  expect_true(all(assay(g, "var") == 1))
  expect_equal(dim(assay(g, "beta")), c(8L, 2L, 1L))
})

test_that("#2 gds_from_nifti_maps warns when maps$contrast has >1 distinct value", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")

  td <- tempfile("niftimaps-multi-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  f1 <- file.path(td, "sub-01_A.nii")
  f2 <- file.path(td, "sub-02_B.nii")
  RNifti::writeNifti(array(1, dim = c(2, 2, 2)), f1)
  RNifti::writeNifti(array(2, dim = c(2, 2, 2)), f2)

  maps <- data.frame(
    file = c(f1, f2),
    subject = c("1001", "1002"),
    contrast = c("condA", "condB"),
    stringsAsFactors = FALSE
  )

  fmrigds:::.reset_synthetic_variance_warning()
  # Two warnings fire (multi-contrast collapse + beta-only unit variance);
  # capture both and assert the collapse warning is present.
  w <- testthat::capture_warnings(g <- gds_from_nifti_maps(maps, mask = NULL))
  expect_match(w, "single-contrast", all = FALSE)
  # Files are still loaded onto a single contrast axis (not pivoted).
  expect_equal(length(contrasts(g)), 1L)
})

# Tests for BIDS-aware subject keying (adapter-nifti.R)

test_that(".nifti_subject_key keys on the BIDS sub- entity (desc- infix names)", {
  beta <- c("sub-01_task-stroop_space-MNI152NLin2009cAsym_desc-beta_bold.nii.gz",
            "sub-02_task-stroop_space-MNI152NLin2009cAsym_desc-beta_bold.nii.gz")
  se <- sub("desc-beta", "desc-se", beta)
  expect_equal(fmrigds:::.nifti_subject_key(beta), c("sub-01", "sub-02"))
  expect_equal(fmrigds:::.nifti_subject_key(se), c("sub-01", "sub-02"))
  # beta/se now pair (this errored before the BIDS-aware fix)
  al <- fmrigds:::.nifti_align_file_sets(beta, se)
  expect_equal(al$subjects, c("sub-01", "sub-02"))
})

test_that(".nifti_align_file_sets refuses beta/se sharing sub- but differing in other entities", {
  # Same subject, but beta is task-A/run-1 and se is task-B/run-9: pairing on
  # sub- alone would silently mismatch. The non-statistic context guard catches it.
  beta <- c("sub-01_task-A_run-1_desc-beta_bold.nii.gz",
            "sub-02_task-A_run-1_desc-beta_bold.nii.gz")
  se <- c("sub-01_task-B_run-9_desc-se_bold.nii.gz",
          "sub-02_task-B_run-9_desc-se_bold.nii.gz")
  expect_error(
    fmrigds:::.nifti_align_file_sets(beta, se),
    "non-statistic filename entities"
  )
  # Explicit subjects bypass the filename pairing entirely.
  al <- fmrigds:::.nifti_align_file_sets(beta, se, subjects = c("a", "b"))
  expect_equal(al$subjects, c("a", "b"))
})

test_that(".nifti_subject_key falls back to legacy trailing-stat strip", {
  expect_equal(
    fmrigds:::.nifti_subject_key(c("subjA_beta.nii.gz", "subjB_se.nii.gz")),
    c("subjA", "subjB")
  )
  # a BIDS sub- entity with a trailing stat token still keys on the entity
  expect_equal(fmrigds:::.nifti_subject_key("sub-03_beta.nii.gz"), "sub-03")
})
