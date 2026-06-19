# Tests for export-neuroim.R functions (import and export)

# Helper to skip tests when neuroim2 is not available
skip_if_no_neuroim2 <- function() {
  skip_if_not_installed("neuroim2")
}

# =============================================================================
# IMPORT TESTS
# =============================================================================

test_that("as_gds.NeuroVol creates valid GDS from single volume", {
  skip_if_no_neuroim2()

  # Create a simple NeuroVol
  vdim <- c(4, 4, 4)
  arr <- array(rnorm(prod(vdim)), dim = vdim)
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))
  vol <- neuroim2::NeuroVol(arr, nspace)

  gds_obj <- as_gds(vol, subject = "sub01", contrast = "faces")

  expect_s3_class(gds_obj, "gds")
  expect_equal(subjects(gds_obj), "sub01")
  expect_equal(contrasts(gds_obj), "faces")
  expect_true("beta" %in% names(assays(gds_obj)))

  # Space should be voxel
  expect_s3_class(space(gds_obj), "space_voxel")
  expect_equal(space(gds_obj)$dim, vdim)
})

test_that("as_gds.NeuroVol respects mask='none'", {
  skip_if_no_neuroim2()

  vdim <- c(3, 3, 3)
  arr <- array(0, dim = vdim)
  arr[2, 2, 2] <- 1  # Only one non-zero voxel
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))
  vol <- neuroim2::NeuroVol(arr, nspace)

  # With default mask (non-zero), should have 1 sample
  gds_masked <- as_gds(vol, mask = NULL)
  expect_equal(dim(assay(gds_masked, "beta"))[1], 1)  # Only non-zero voxel

  # With mask = "none", should have all voxels
  gds_full <- as_gds(vol, mask = "none")
  expect_equal(dim(assay(gds_full, "beta"))[1], prod(vdim))
})

test_that("as_gds.NeuroVec imports 4D volume along subjects", {
  skip_if_no_neuroim2()

  vdim <- c(4, 4, 4)
  n_subj <- 3
  arr <- array(rnorm(prod(vdim) * n_subj), dim = c(vdim, n_subj))

  # NeuroSpace for 4D: spacing/origin are 3D, dim is 4D
  nspace <- neuroim2::NeuroSpace(
    dim = c(vdim, n_subj),
    spacing = c(1, 1, 1),
    origin = c(0, 0, 0)
  )
  vec4d <- neuroim2::NeuroVec(arr, nspace)

  gds_obj <- as_gds(vec4d,
                    along = "subject",
                    subjects = c("s1", "s2", "s3"),
                    mask = "none")

  expect_equal(subjects(gds_obj), c("s1", "s2", "s3"))
  expect_equal(contrasts(gds_obj), "contrast1")
  expect_equal(dim(assay(gds_obj, "beta")), c(prod(vdim), 3, 1))
})

test_that("as_gds.NeuroVec imports 4D volume along contrasts", {
  skip_if_no_neuroim2()

  vdim <- c(4, 4, 4)
  n_con <- 2
  arr <- array(rnorm(prod(vdim) * n_con), dim = c(vdim, n_con))

  nspace <- neuroim2::NeuroSpace(
    dim = c(vdim, n_con),
    spacing = c(1, 1, 1),
    origin = c(0, 0, 0)
  )
  vec4d <- neuroim2::NeuroVec(arr, nspace)

  gds_obj <- as_gds(vec4d,
                    along = "contrast",
                    contrasts = c("faces", "places"),
                    mask = "none")

  expect_equal(subjects(gds_obj), "subject1")
  expect_equal(contrasts(gds_obj), c("faces", "places"))
  expect_equal(dim(assay(gds_obj, "beta")), c(prod(vdim), 1, 2))
})

test_that("gds_from_neurovols creates GDS from list of volumes", {
  skip_if_no_neuroim2()

  vdim <- c(4, 4, 4)
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))

  # Create named list of beta volumes
  beta_vols <- list(
    `sub-01` = neuroim2::NeuroVol(array(rnorm(prod(vdim)), dim = vdim), nspace),
    `sub-02` = neuroim2::NeuroVol(array(rnorm(prod(vdim)), dim = vdim), nspace)
  )

  # Create named list of variance volumes
  var_vols <- list(
    `sub-01` = neuroim2::NeuroVol(array(abs(rnorm(prod(vdim))) + 0.1, dim = vdim), nspace),
    `sub-02` = neuroim2::NeuroVol(array(abs(rnorm(prod(vdim))) + 0.1, dim = vdim), nspace)
  )

  gds_obj <- gds_from_neurovols(beta_vols, var = var_vols, mask = "none")

  expect_s3_class(gds_obj, "gds")
  expect_equal(subjects(gds_obj), c("sub-01", "sub-02"))
  expect_true("beta" %in% names(assays(gds_obj)))
  expect_true("var" %in% names(assays(gds_obj)))
})

test_that("gds_from_neurovols accepts col_data", {
  skip_if_no_neuroim2()

  vdim <- c(3, 3, 3)
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))

  beta_vols <- list(
    `s1` = neuroim2::NeuroVol(array(rnorm(prod(vdim)), dim = vdim), nspace),
    `s2` = neuroim2::NeuroVol(array(rnorm(prod(vdim)), dim = vdim), nspace)
  )

  covars <- data.frame(
    age = c(25, 30),
    group = c("A", "B"),
    row.names = c("s1", "s2")
  )

  expect_warning(
    gds_obj <- gds_from_neurovols(beta_vols, col_data = covars, mask = "none"),
    "No variance or SE provided; using unit variance",
    fixed = TRUE
  )

  cd <- col_data(gds_obj)
  expect_equal(nrow(cd), 2)
  expect_true("age" %in% names(cd))
  expect_true("group" %in% names(cd))
})

test_that("gds_from_neurovols converts SE to variance", {
  skip_if_no_neuroim2()

  vdim <- c(3, 3, 3)
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))

  beta_vols <- list(
    `s1` = neuroim2::NeuroVol(array(1, dim = vdim), nspace)
  )

  # SE of 2 should become variance of 4
  se_vols <- list(
    `s1` = neuroim2::NeuroVol(array(2, dim = vdim), nspace)
  )

  gds_obj <- gds_from_neurovols(beta_vols, se = se_vols, mask = "none")

  var_vals <- assay(gds_obj, "var")
  expect_true(all(var_vals == 4))  # SE^2 = 2^2 = 4
})

test_that("gds_from_neurovol_nested handles nested list (subjects x contrasts)", {
  skip_if_no_neuroim2()

  vdim <- c(3, 3, 3)
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))

  # Create nested structure: subjects -> contrasts -> NeuroVol
  beta <- list(
    `s1` = list(
      faces = neuroim2::NeuroVol(array(1, dim = vdim), nspace),
      places = neuroim2::NeuroVol(array(2, dim = vdim), nspace)
    ),
    `s2` = list(
      faces = neuroim2::NeuroVol(array(3, dim = vdim), nspace),
      places = neuroim2::NeuroVol(array(4, dim = vdim), nspace)
    )
  )

  var <- list(
    `s1` = list(
      faces = neuroim2::NeuroVol(array(0.1, dim = vdim), nspace),
      places = neuroim2::NeuroVol(array(0.2, dim = vdim), nspace)
    ),
    `s2` = list(
      faces = neuroim2::NeuroVol(array(0.3, dim = vdim), nspace),
      places = neuroim2::NeuroVol(array(0.4, dim = vdim), nspace)
    )
  )

  gds_obj <- gds_from_neurovol_nested(beta, var = var, mask = "none")

  expect_s3_class(gds_obj, "gds")
  expect_equal(subjects(gds_obj), c("s1", "s2"))
  expect_equal(contrasts(gds_obj), c("faces", "places"))
  expect_equal(dim(assay(gds_obj, "beta")), c(prod(vdim), 2, 2))

  # Check values are correct
  # s1, faces should be all 1s
  expect_true(all(assay(gds_obj, "beta")[, 1, 1] == 1))
  # s2, places should be all 4s
  expect_true(all(assay(gds_obj, "beta")[, 2, 2] == 4))
})

test_that("gds_from_neurovol_nested handles list of NeuroVecs", {
  skip_if_no_neuroim2()

  vdim <- c(3, 3, 3)
  n_con <- 2

  # Create 4D volumes for each subject
  arr1 <- array(c(rep(1, prod(vdim)), rep(2, prod(vdim))), dim = c(vdim, n_con))
  arr2 <- array(c(rep(3, prod(vdim)), rep(4, prod(vdim))), dim = c(vdim, n_con))

  nspace <- neuroim2::NeuroSpace(dim = c(vdim, n_con), spacing = c(1, 1, 1), origin = c(0, 0, 0))

  beta <- list(
    `s1` = neuroim2::NeuroVec(arr1, nspace),
    `s2` = neuroim2::NeuroVec(arr2, nspace)
  )

  expect_warning(
    gds_obj <- gds_from_neurovol_nested(beta,
                                        contrasts = c("faces", "places"),
                                        mask = "none"),
    "No variance or SE provided; using unit variance",
    fixed = TRUE
  )

  expect_equal(subjects(gds_obj), c("s1", "s2"))
  expect_equal(contrasts(gds_obj), c("faces", "places"))
  expect_equal(dim(assay(gds_obj, "beta")), c(prod(vdim), 2, 2))

  # s1, con1 (faces) should be all 1s
  expect_true(all(assay(gds_obj, "beta")[, 1, 1] == 1))
  # s2, con2 (places) should be all 4s
  expect_true(all(assay(gds_obj, "beta")[, 2, 2] == 4))
})

test_that("gds_from_neurovol_nested attaches col_data and metadata", {
  skip_if_no_neuroim2()

  vdim <- c(3, 3, 3)
  nspace <- neuroim2::NeuroSpace(dim = vdim, spacing = c(1, 1, 1), origin = c(0, 0, 0))

  beta <- list(
    `s1` = list(con1 = neuroim2::NeuroVol(array(1, dim = vdim), nspace)),
    `s2` = list(con1 = neuroim2::NeuroVol(array(2, dim = vdim), nspace))
  )

  covars <- data.frame(
    age = c(25, 30),
    group = c("A", "B"),
    row.names = c("s1", "s2")
  )

  custom_meta <- list(
    study = "My Study",
    notes = "Test data"
  )

  expect_warning(
    gds_obj <- gds_from_neurovol_nested(beta,
                                        col_data = covars,
                                        metadata = custom_meta,
                                        mask = "none"),
    "No variance or SE provided; using unit variance",
    fixed = TRUE
  )

  cd <- col_data(gds_obj)
  expect_equal(nrow(cd), 2)
  expect_true("age" %in% names(cd))

  meta <- metadata(gds_obj)
  expect_equal(meta$study, "My Study")
  expect_equal(meta$notes, "Test data")
})

# =============================================================================
# SPLIT TESTS
# =============================================================================

test_that("split.gds splits by column name", {
  # Create a simple GDS with 4 subjects in 2 groups
  beta <- array(rnorm(10 * 4 * 1), dim = c(10, 4, 1))
  var <- array(abs(rnorm(10 * 4 * 1)) + 0.1, dim = c(10, 4, 1))

  sp <- space_sample_labels(paste0("s", 1:10))

  col_data <- data.frame(
    subject = paste0("subj", 1:4),
    group = c("A", "A", "B", "B"),
    age = c(25, 30, 28, 35),
    row.names = paste0("subj", 1:4)
  )

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = paste0("subj", 1:4),
    contrasts = "con1",
    col_data = col_data
  )

  # Split by group column
  split_result <- split(gds_obj, "group")

  expect_type(split_result, "list")
  expect_named(split_result, c("A", "B"))

  # Check group A
  expect_equal(length(subjects(split_result$A)), 2)
  expect_equal(subjects(split_result$A), c("subj1", "subj2"))

  # Check group B
  expect_equal(length(subjects(split_result$B)), 2)
  expect_equal(subjects(split_result$B), c("subj3", "subj4"))

  # Check assay dimensions
  expect_equal(dim(assay(split_result$A, "beta")), c(10, 2, 1))
  expect_equal(dim(assay(split_result$B, "beta")), c(10, 2, 1))

  # Check col_data is subset correctly
  expect_equal(nrow(col_data(split_result$A)), 2)
  expect_equal(col_data(split_result$A)$age, c(25, 30))
})

test_that("split.gds works with factor argument", {
  beta <- array(rnorm(5 * 3 * 1), dim = c(5, 3, 1))
  var <- array(abs(rnorm(5 * 3 * 1)) + 0.1, dim = c(5, 3, 1))

  sp <- space_sample_labels(paste0("s", 1:5))

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("s1", "s2", "s3"),
    contrasts = "con1"
  )

  # Split using a factor directly
  f <- factor(c("X", "Y", "X"))
  split_result <- split(gds_obj, f)

  expect_named(split_result, c("X", "Y"))
  expect_equal(length(subjects(split_result$X)), 2)
  expect_equal(length(subjects(split_result$Y)), 1)
})

test_that(".subset_gds_subjects preserves structure", {
  beta <- array(1:24, dim = c(4, 3, 2))
  var <- array(rep(0.1, 24), dim = c(4, 3, 2))

  sp <- space_sample_labels(paste0("s", 1:4))

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("A", "B", "C"),
    contrasts = c("c1", "c2")
  )

  # Subset to first two subjects
  sub_gds <- fmrigds:::.subset_gds_subjects(gds_obj, 1:2)

  expect_equal(subjects(sub_gds), c("A", "B"))
  expect_equal(contrasts(sub_gds), c("c1", "c2"))
  expect_equal(dim(assay(sub_gds, "beta")), c(4, 2, 2))

  # Check values are correct (should be columns 1,2 of original)
  expect_equal(assay(sub_gds, "beta")[, 1, 1], assay(gds_obj, "beta")[, 1, 1])
  expect_equal(assay(sub_gds, "beta")[, 2, 1], assay(gds_obj, "beta")[, 2, 1])
})

# =============================================================================
# EXPORT TESTS (neuroim2 dependent)
# =============================================================================

test_that("as_neurovol_list requires voxel space", {
  skip_if_no_neuroim2()

  # Create GDS with non-voxel space
  beta <- array(rnorm(10 * 2 * 1), dim = c(10, 2, 1))
  var <- array(abs(rnorm(10 * 2 * 1)) + 0.1, dim = c(10, 2, 1))
  sp <- space_sample_labels(paste0("s", 1:10))

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("s1", "s2"),
    contrasts = "con1"
  )

  expect_error(
    as_neurovol_list(gds_obj),
    "voxel space"
  )
})

test_that("as_neurovec requires voxel space", {
  skip_if_no_neuroim2()

  beta <- array(rnorm(10 * 2 * 1), dim = c(10, 2, 1))
  var <- array(abs(rnorm(10 * 2 * 1)) + 0.1, dim = c(10, 2, 1))
  sp <- space_sample_labels(paste0("s", 1:10))

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("s1", "s2"),
    contrasts = "con1"
  )

  expect_error(
    as_neurovec(gds_obj),
    "voxel space"
  )
})

test_that("extract_group validates inputs", {
  beta <- array(rnorm(10 * 4 * 1), dim = c(10, 4, 1))
  var <- array(abs(rnorm(10 * 4 * 1)) + 0.1, dim = c(10, 4, 1))
  sp <- space_sample_labels(paste0("s", 1:10))

  col_data <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = paste0("subj", 1:4)
  )

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = paste0("subj", 1:4),
    contrasts = "con1",
    col_data = col_data
  )

  # Non-existent column
  expect_error(
    extract_group(gds_obj, "nonexistent", "A"),
    "not found in col_data"
  )

  # Non-existent level
  expect_error(
    extract_group(gds_obj, "group", "C"),
    "not found in column"
  )
})

test_that("as_neurovol_list works with voxel space", {
  skip_if_no_neuroim2()

  # Create a small voxel-space GDS
  vdim <- c(4, 4, 4)
  n_vox <- prod(vdim)
  n_subj <- 3

  beta <- array(rnorm(n_vox * n_subj * 1), dim = c(n_vox, n_subj, 1))
  var <- array(abs(rnorm(n_vox * n_subj * 1)) + 0.1, dim = c(n_vox, n_subj, 1))

  affine <- diag(4)
  sp <- space_voxel(dim = vdim, affine = affine, storage = "dense")

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("s1", "s2", "s3"),
    contrasts = "con1"
  )

  vols <- as_neurovol_list(gds_obj, assay = "beta")

  expect_type(vols, "list")
  expect_length(vols, 3)
  expect_named(vols, c("s1", "s2", "s3"))

  # Check that each element is a NeuroVol
  expect_true(inherits(vols$s1, "NeuroVol"))
  expect_equal(dim(vols$s1), vdim)
})

test_that("as_neurovec works with voxel space", {
  skip_if_no_neuroim2()

  vdim <- c(4, 4, 4)
  n_vox <- prod(vdim)
  n_subj <- 3

  beta <- array(rnorm(n_vox * n_subj * 1), dim = c(n_vox, n_subj, 1))
  var <- array(abs(rnorm(n_vox * n_subj * 1)) + 0.1, dim = c(n_vox, n_subj, 1))

  affine <- diag(4)
  sp <- space_voxel(dim = vdim, affine = affine, storage = "dense")

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("s1", "s2", "s3"),
    contrasts = "con1"
  )

  vec4d <- as_neurovec(gds_obj, assay = "beta", along = "subject")

  expect_true(inherits(vec4d, "NeuroVec"))
  expect_equal(dim(vec4d), c(vdim, n_subj))
})

test_that("as_neurovol_list handles packed voxel space", {
  skip_if_no_neuroim2()

  # Create a packed voxel-space GDS (with mask)
  vdim <- c(4, 4, 4)
  n_vox <- prod(vdim)
  mask_idx <- c(1, 5, 10, 20, 30, 50, 60)  # 7 voxels in mask
  n_masked <- length(mask_idx)
  n_subj <- 2

  beta <- array(rnorm(n_masked * n_subj * 1), dim = c(n_masked, n_subj, 1))
  var <- array(abs(rnorm(n_masked * n_subj * 1)) + 0.1, dim = c(n_masked, n_subj, 1))

  affine <- diag(4)
  sp <- space_voxel(dim = vdim, affine = affine, mask_idx = mask_idx, storage = "packed")

  gds_obj <- new_gds(
    assays = list(beta = beta, var = var),
    space = sp,
    subjects = c("s1", "s2"),
    contrasts = "con1"
  )

  vols <- as_neurovol_list(gds_obj, assay = "beta")

  expect_length(vols, 2)
  expect_true(inherits(vols$s1, "NeuroVol"))
  expect_equal(dim(vols$s1), vdim)

  # Check that non-masked voxels are zero and masked voxels have values
  vol_arr <- as.array(vols$s1)
  expect_equal(sum(vol_arr != 0), n_masked)
})
