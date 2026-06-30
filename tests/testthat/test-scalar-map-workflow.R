test_that("gds_from_scalar_maps builds beta-only GDS and attaches col_data", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")

  td <- tempfile("scalar-maps-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  files <- file.path(td, paste0("sub-", 1:2, ".nii"))
  RNifti::writeNifti(array(1, c(2, 2, 2)), files[1])
  RNifti::writeNifti(array(2, c(2, 2, 2)), files[2])
  cd <- data.frame(subject = c("s1", "s2"), group = c("ctl", "pt"))

  fmrigds:::.reset_synthetic_variance_warning()
  expect_warning(
    g <- gds_from_scalar_maps(files, subject = cd$subject, contrast = "metric", col_data = cd),
    "No variance or SE provided; using unit variance",
    fixed = TRUE
  )

  expect_s3_class(g, "gds")
  expect_equal(subjects(g), c("s1", "s2"))
  expect_equal(contrasts(g), "metric")
  expect_equal(rownames(col_data(g)), c("s1", "s2"))
  expect_true(all(assay(g, "var") == 1))
})

test_that("group_ols and two_sample create voxelwise OLS plans", {
  beta <- array(c(
    1, 2, 3,
    2, 3, 4,
    4, 5, 6,
    5, 6, 7
  ), dim = c(3, 4, 1))
  var <- array(1, dim = dim(beta))
  cd <- data.frame(group = c("ctl", "ctl", "pt", "pt"),
                   row.names = paste0("s", 1:4))
  g <- new_gds(
    list(beta = beta, var = var),
    space_sample_labels(paste0("v", 1:3)),
    rownames(cd),
    "metric",
    col_data = cd
  )

  fit <- group_ols(g, ~ group) |> compute()
  expect_true("coef:grouppt" %in% names(assays(fit)))
  expect_true("t_coef:grouppt" %in% names(assays(fit)))

  p <- two_sample(g, group = "group", baseline = "ctl")
  expect_s3_class(p, "gds_plan")
  expect_equal(p$nodes[[1]]$method, "ols:voxelwise")
  expect_true("X" %in% names(p$nodes[[1]]$options))
})

test_that("write_nifti_assays writes selected image assays and returns manifest", {
  skip_if_not_installed("RNifti")

  td <- tempfile("nifti-assays-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  arr <- array(1:8, dim = c(8, 1, 1))
  g <- new_gds(
    list(
      "t_coef:grouppt" = arr,
      "p_coef:grouppt" = arr / 10,
      var = array(1, dim = dim(arr))
    ),
    space_voxel(c(2, 2, 2), diag(4)),
    "meta",
    "metric"
  )

  manifest <- write_nifti_assays(
    g,
    out_dir = td,
    prefix = "metric",
    assays = c("t_coef:grouppt", "p_coef:grouppt", "missing")
  )

  expect_equal(nrow(manifest), 3)
  expect_equal(sum(manifest$written), 2)
  expect_true(any(manifest$skipped_reason == "assay not found", na.rm = TRUE))
  expect_true(all(file.exists(manifest$path[manifest$written])))
  expect_true(all(grepl("metric_.*_coef_grouppt\\.nii\\.gz$", basename(manifest$path[manifest$written]))))

  second <- write_nifti_assays(g, out_dir = td, prefix = "metric", assays = "t_coef:grouppt")
  expect_false(second$written)
  expect_equal(second$skipped_reason, "file exists")
})

test_that("#13 write_nifti_assays embeds the contrast token in filenames", {
  skip_if_not_installed("RNifti")
  arr <- array(as.numeric(1:8), dim = c(8, 1, 1))
  g <- new_gds(
    list("t_coef:x" = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2, 2, 2), diag(4)), "meta", "condA"
  )
  td <- tempfile("ct-"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  m <- write_nifti_assays(g, out_dir = td, prefix = "study", assays = "t_coef:x")
  expect_true(all(m$written))
  expect_match(basename(m$path[m$written]), "contrast-condA")
})

test_that("#13 separate-contrast writes do not collide/overwrite", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
  arr <- array(0, dim = c(8, 1, 2))
  arr[, 1, 1] <- 1:8
  arr[, 1, 2] <- (1:8) * 100
  g <- new_gds(
    list("t_coef:x" = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2, 2, 2), diag(4)), "meta", c("condA", "condB")
  )
  td <- tempfile("loop-"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  mA <- write_nifti_assays(g, out_dir = td, prefix = "study", assays = "t_coef:x", contrasts = "condA")
  mB <- write_nifti_assays(g, out_dir = td, prefix = "study", assays = "t_coef:x", contrasts = "condB")
  expect_true(mA$written && mB$written)
  expect_false(identical(mA$path, mB$path))
  expect_equal(length(list.files(td)), 2L)
  vA <- as.numeric(as.array(neuroim2::read_vol(mA$path[mA$written])))
  vB <- as.numeric(as.array(neuroim2::read_vol(mB$path[mB$written])))
  expect_equal(sort(vA[vA != 0]), as.numeric(1:8))
  expect_equal(sort(vB[vB != 0]), as.numeric((1:8) * 100))
})

test_that("#13 write_nifti_assays errors on within-call filename collisions", {
  skip_if_not_installed("RNifti")
  arr <- array(0, dim = c(8, 1, 2))
  arr[, 1, 1] <- 1:8
  arr[, 1, 2] <- 9:16
  # 'a/b' and 'a_b' both sanitise to 'a_b' -> distinct contrasts, same filename
  g <- new_gds(
    list("t_coef:x" = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2, 2, 2), diag(4)), "meta", c("a/b", "a_b")
  )
  td <- tempfile("col-"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_error(
    write_nifti_assays(g, out_dir = td, prefix = "study", assays = "t_coef:x"),
    "collision"
  )
})

test_that("#13 multi-subject export: no double sub- prefix, no one-per-call collision", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
  arr <- array(0, dim = c(8, 2, 1))
  arr[, 1, 1] <- 1:8
  arr[, 2, 1] <- (1:8) * 10
  # BIDS-style subject labels must not be double-prefixed to sub-sub-01.
  g <- new_gds(
    list("beta" = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2, 2, 2), diag(4)), c("sub-01", "sub-02"), "cond"
  )
  td <- tempfile("multi-subj-"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # Export one subject per call into the same dir/prefix -> must not collide.
  mA <- write_nifti_assays(g, out_dir = td, prefix = "study", assays = "beta", subjects = "sub-01")
  mB <- write_nifti_assays(g, out_dir = td, prefix = "study", assays = "beta", subjects = "sub-02")
  expect_true(mA$written && mB$written)
  expect_false(identical(mA$path, mB$path))
  expect_match(basename(mA$path[mA$written]), "sub-01")
  expect_false(grepl("sub-sub-01", basename(mA$path[mA$written])))
  expect_equal(length(list.files(td)), 2L)
  vA <- as.numeric(as.array(neuroim2::read_vol(mA$path[mA$written])))
  vB <- as.numeric(as.array(neuroim2::read_vol(mB$path[mB$written])))
  expect_equal(sort(vA[vA != 0]), as.numeric(1:8))
  expect_equal(sort(vB[vB != 0]), as.numeric((1:8) * 10))
})
