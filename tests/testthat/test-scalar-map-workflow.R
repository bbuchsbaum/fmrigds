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
