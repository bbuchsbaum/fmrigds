# Coverage tests for space-utils.R

# ===========================================================================
# space_subset
# ===========================================================================

test_that("space_subset parcels by integer index", {
  sp <- space_parcels(c("A", "B", "C", "D"))
  result <- space_subset(sp, 2:3)
  expect_s3_class(result, "space_parcels")
  expect_equal(result$labels, c("B", "C"))
})

test_that("space_subset parcels by logical index", {
  sp <- space_parcels(c("A", "B", "C"))
  result <- space_subset(sp, c(TRUE, FALSE, TRUE))
  expect_equal(result$labels, c("A", "C"))
})

test_that("space_subset sample_labels by index", {
  sp <- space_sample_labels(c("x", "y", "z"))
  result <- space_subset(sp, c(1, 3))
  expect_s3_class(result, "space_sample_labels")
  expect_equal(result$labels, c("x", "z"))
})

test_that("space_subset sample_labels by logical index", {
  sp <- space_sample_labels(c("x", "y", "z"))
  result <- space_subset(sp, c(FALSE, TRUE, TRUE))
  expect_equal(result$labels, c("y", "z"))
})

test_that("space_subset packed voxel space", {
  sp <- space_voxel(c(3, 3, 3), diag(4), mask_idx = c(1L, 5L, 10L, 15L), storage = "packed")
  result <- space_subset(sp, c(1, 3))
  expect_equal(result$mask_idx, c(1L, 10L))
  expect_equal(result$storage, "packed")
})

test_that("space_subset dense voxel with pack=TRUE", {
  sp <- space_voxel(c(2, 2, 2), diag(4))
  result <- space_subset(sp, c(1L, 4L, 7L), pack = TRUE)
  expect_equal(result$storage, "packed")
  expect_equal(result$mask_idx, c(1L, 4L, 7L))
})

test_that("space_subset dense voxel without pack warns", {
  sp <- space_voxel(c(2, 2, 2), diag(4))
  expect_warning(space_subset(sp, c(1L, 4L)), "cannot be reduced")
})

test_that("space_subset errors on unsupported space", {
  sp <- space_basis(5)
  expect_error(space_subset(sp, 1:3), "Unsupported")
})

# ===========================================================================
# common_mask
# ===========================================================================

test_that("common_mask intersection of two voxel spaces", {
  m1 <- array(c(TRUE, FALSE, TRUE, TRUE), c(2, 2, 1))
  m2 <- array(c(TRUE, TRUE, FALSE, TRUE), c(2, 2, 1))
  sp1 <- space_voxel(c(2, 2, 1), diag(4), mask_bitmap = m1, storage = "packed")
  sp2 <- space_voxel(c(2, 2, 1), diag(4), mask_bitmap = m2, storage = "packed")
  result <- common_mask(sp1, sp2, rule = "intersection")
  expect_true(is.list(result))
  expect_true(length(result$idx1) > 0)
  expect_true(length(result$idx2) > 0)
})

test_that("common_mask union of two voxel spaces", {
  m1 <- array(c(TRUE, FALSE, FALSE, FALSE), c(2, 2, 1))
  m2 <- array(c(FALSE, FALSE, TRUE, FALSE), c(2, 2, 1))
  sp1 <- space_voxel(c(2, 2, 1), diag(4), mask_bitmap = m1, storage = "packed")
  sp2 <- space_voxel(c(2, 2, 1), diag(4), mask_bitmap = m2, storage = "packed")
  result <- common_mask(sp1, sp2, rule = "union")
  expect_true(length(result$idx1) >= 1)
})

test_that("common_mask errors on non-voxel spaces", {
  sp <- space_parcels(c("A", "B"))
  expect_error(common_mask(sp, sp), "voxel spaces")
})

test_that("common_mask errors on dim mismatch", {
  sp1 <- space_voxel(c(2, 2, 2), diag(4))
  sp2 <- space_voxel(c(3, 3, 3), diag(4))
  expect_error(common_mask(sp1, sp2), "dims differ")
})

test_that("common_mask with dense voxel spaces (no mask)", {
  sp1 <- space_voxel(c(2, 2, 1), diag(4))
  sp2 <- space_voxel(c(2, 2, 1), diag(4))
  result <- common_mask(sp1, sp2, rule = "intersection")
  expect_equal(length(result$idx1), 4)
})

# ===========================================================================
# apply_common_mask
# ===========================================================================

test_that("apply_common_mask subsets two GDS objects", {
  m1 <- array(c(TRUE, FALSE, TRUE, TRUE), c(2, 2, 1))
  m2 <- array(c(TRUE, TRUE, FALSE, TRUE), c(2, 2, 1))
  sp1 <- space_voxel(c(2, 2, 1), diag(4), mask_bitmap = m1, storage = "packed")
  sp2 <- space_voxel(c(2, 2, 1), diag(4), mask_bitmap = m2, storage = "packed")

  g1 <- new_gds(
    assays = list(beta = array(1:3, dim = c(3, 1, 1)), var = array(0.1, dim = c(3, 1, 1))),
    space = sp1, subjects = "s1", contrasts = "c1"
  )
  g2 <- new_gds(
    assays = list(beta = array(4:6, dim = c(3, 1, 1)), var = array(0.1, dim = c(3, 1, 1))),
    space = sp2, subjects = "s1", contrasts = "c1"
  )

  result <- apply_common_mask(g1, g2)
  expect_true(is.list(result))
  expect_s3_class(result$g1, "gds")
  expect_s3_class(result$g2, "gds")
  # Both should have same number of samples after masking
  expect_equal(dim(result$g1$assays$beta)[1], dim(result$g2$assays$beta)[1])
})

# ===========================================================================
# assert_compatible_spaces
# ===========================================================================

test_that("assert_compatible_spaces passes for matching spaces", {
  sp <- space_sample_labels(c("a", "b"))
  expect_silent(assert_compatible_spaces(sp, sp))
})
