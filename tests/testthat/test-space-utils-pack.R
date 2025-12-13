test_that("space_subset packs dense voxel space when requested", {
  sp <- space_voxel(dim = c(2,2,1), affine = diag(4), storage = "dense")
  sub <- space_subset(sp, c(1L, 3L), pack = TRUE)
  expect_equal(sub$storage, "packed")
  expect_equal(as.integer(sub$mask_idx), c(1L,3L))
})

