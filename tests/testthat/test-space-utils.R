test_that("space_from_nifti creates voxel space when neuroim2 available", {
  testthat::skip_if_not_installed("neuroim2")
  # Create a tiny 3D volume via neuroim2 helpers
  sp <- neuroim2::NeuroSpace(c(2,2,2), spacing = c(1,1,1))
  vol <- neuroim2::NeuroVol(array(0, dim = c(2,2,2)), sp)
  tf <- tempfile(fileext = ".nii")
  on.exit(unlink(tf))
  neuroim2::write_vol(vol, tf)
  sp <- space_from_nifti(tf)
  expect_s3_class(sp, "gds_space")
  expect_equal(sp$type, "voxel")
  expect_equal(sp$dim, c(2L,2L,2L))
})

test_that("space_subset trims parcels and labels", {
  sp_par <- space_parcels(labels = letters[1:5])
  sub <- space_subset(sp_par, 2:4)
  expect_equal(sub$labels, letters[2:4])

  sp_lab <- space_sample_labels(labels = c("a","b","c"))
  sub2 <- space_subset(sp_lab, c(TRUE, FALSE, TRUE))
  expect_equal(sub2$labels, c("a","c"))
})

test_that("common_mask returns indices for both spaces", {
  sp1 <- space_voxel(dim = c(3,3,1), affine = diag(4), mask_bitmap = array(c(1,0,1, 1,1,0, 0,1,1) > 0, dim = c(3,3,1)), storage = "packed")
  sp2 <- space_voxel(dim = c(3,3,1), affine = diag(4), mask_bitmap = array(c(1,1,0, 1,0,1, 0,1,1) > 0, dim = c(3,3,1)), storage = "packed")
  cm <- common_mask(sp1, sp2, rule = "intersection")
  expect_true(is.list(cm) && all(c("idx1","idx2") %in% names(cm)))
  expect_true(length(cm$idx1) <= length(sp1$mask_idx))
  expect_true(length(cm$idx2) <= length(sp2$mask_idx))
})
