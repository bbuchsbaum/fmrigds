# Coverage tests for space.R, space-utils.R, summary.R

# ===========================================================================
# space constructors
# ===========================================================================

test_that("space_parcels creates valid object", {
  sp <- space_parcels(c("ROI_A", "ROI_B", "ROI_C"))
  expect_s3_class(sp, "space_parcels")
  expect_s3_class(sp, "gds_space")
  expect_equal(sp$labels, c("ROI_A", "ROI_B", "ROI_C"))
  expect_null(sp$lookup)
})

test_that("space_parcels with lookup and membership", {
  sp <- space_parcels(c("R1", "R2"), lookup = data.frame(id = 1:2), membership = list(1:10, 11:20))
  expect_equal(sp$lookup$id, 1:2)
  expect_true(!is.null(sp$membership))
})

test_that("space_parcels errors on non-character labels", {
  expect_error(space_parcels(1:3), "non-empty character")
  expect_error(space_parcels(character(0)), "non-empty character")
})

test_that("space_surface creates valid object", {
  verts <- matrix(rnorm(9), ncol = 3)
  faces <- matrix(c(1L, 2L, 3L), ncol = 3)
  sp <- space_surface(verts, faces, "L")
  expect_s3_class(sp, "space_surface")
  expect_s3_class(sp, "gds_space")
  expect_equal(sp$hemi, "L")
})

test_that("space_surface with template_id", {
  verts <- matrix(rnorm(9), ncol = 3)
  faces <- matrix(c(1L, 2L, 3L), ncol = 3)
  sp <- space_surface(verts, faces, "R", template_id = "fsaverage")
  expect_equal(sp$template_id, "fsaverage")
})

test_that("space_surface validation errors", {
  expect_error(space_surface(1:3, matrix(1:3, ncol = 3), "L"), "matrix with 3 columns")
  expect_error(space_surface(matrix(1:6, ncol = 2), matrix(1:3, ncol = 3), "L"), "matrix with 3 columns")
  expect_error(space_surface(matrix(1:9, ncol = 3), 1:3, "L"), "matrix with 3 columns")
  expect_error(space_surface(matrix(1:9, ncol = 3), matrix(1:3, ncol = 3), "X"), "must be one of")
})

test_that("space_basis creates valid object", {
  sp <- space_basis(10, basis_name = "ICA")
  expect_s3_class(sp, "space_basis")
  expect_s3_class(sp, "gds_space")
  expect_equal(sp$k, 10L)
  expect_equal(sp$basis_name, "ICA")
})

test_that("space_basis with projector matrix", {
  proj <- matrix(rnorm(50), nrow = 5)
  sp <- space_basis(5, projector = proj)
  expect_true(is.matrix(sp$projector))
})

test_that("space_basis validation errors", {
  expect_error(space_basis(-1), "positive scalar")
  expect_error(space_basis(c(1, 2)), "positive scalar")
  expect_error(space_basis(5, projector = "not_matrix"), "matrix or function")
})

test_that("space_sample_labels creates valid object", {
  sp <- space_sample_labels(c("x", "y", "z"))
  expect_s3_class(sp, "space_sample_labels")
  expect_equal(sp$labels, c("x", "y", "z"))
})

test_that("space_sample_labels errors on invalid input", {
  expect_error(space_sample_labels(character(0)), "non-empty character")
  expect_error(space_sample_labels(1:3), "non-empty character")
})

test_that("space_voxel validation errors", {
  expect_error(space_voxel(c(2, 2), diag(4)), "length 3")
  expect_error(space_voxel(c(2, 2, 2), diag(3)), "4x4 matrix")
  expect_error(space_voxel(c(2, 2, 2), diag(4), mask_bitmap = "not_logical"), "must be logical")
  expect_error(space_voxel(c(2, 2, 2), diag(4), mask_bitmap = array(TRUE, dim = c(3, 3, 3))), "dimensions matching")
  expect_error(space_voxel(c(2, 2, 2), diag(4), mask_idx = c(0, 1)), "out of range")
  expect_error(space_voxel(c(2, 2, 2), diag(4), mask_idx = c(1, 100)), "out of range")
  expect_error(space_voxel(c(2, 2, 2), diag(4), mask_idx = "abc"), "must be numeric")
})

test_that("space_voxel derives mask_idx from mask_bitmap", {
  bm <- array(FALSE, dim = c(2, 2, 2))
  bm[1, 1, 1] <- TRUE
  bm[2, 1, 1] <- TRUE
  sp <- space_voxel(c(2, 2, 2), diag(4), mask_bitmap = bm)
  expect_equal(sp$mask_idx, c(1L, 2L))
})

# ===========================================================================
# summary.R: explain and validate
# ===========================================================================

test_that("explain.gds produces output", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  out <- capture.output(explain(g))
  expect_true(any(grepl("GDS", out)))
  expect_true(any(grepl("sample", out, ignore.case = TRUE)))
})

test_that("explain.gds with voxel space", {
  g <- new_gds(
    assays = list(beta = array(1:8, dim = c(8, 1, 1)), var = array(0.1, dim = c(8, 1, 1))),
    space = space_voxel(c(2, 2, 2), diag(4)),
    subjects = "s1", contrasts = "c1"
  )
  out <- capture.output(explain(g))
  expect_true(any(grepl("[Vv]oxel", out)))
})

test_that("explain.gds with parcels space", {
  g <- new_gds(
    assays = list(beta = array(1:2, dim = c(2, 1, 1)), var = array(0.1, dim = c(2, 1, 1))),
    space = space_parcels(c("P1", "P2")),
    subjects = "s1", contrasts = "c1"
  )
  out <- capture.output(explain(g))
  expect_true(any(grepl("[Pp]arcel", out)))
})

test_that("explain.gds_plan produces output", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  p <- as_plan(g)
  out <- capture.output(explain(p))
  expect_true(length(out) > 0)
})

test_that("validate.gds_plan works", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  p <- as_plan(g) |> derive("se")
  result <- validate(p)
  expect_true(is.list(result) || is.logical(result) || is.character(result))
})

test_that("explain_plan returns tidy table", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  p <- as_plan(g) |> derive("se") |> reduce("fixed")
  tbl <- explain_plan(p)
  expect_s3_class(tbl, "data.frame")
  expect_true(nrow(tbl) >= 2)
})

# ===========================================================================
# gds-class validators and helpers
# ===========================================================================

test_that(".validate_gds_assays errors on unnamed list", {
  expect_error(.validate_gds_assays(list(array(1, dim = c(1,1,1)))), "named")
})

test_that(".validate_gds_assays upcasts 2D to 3D", {
  result <- .validate_gds_assays(list(signal = matrix(1:4, 2, 2)))
  expect_equal(length(dim(result[[1]])), 3L)
})

test_that(".validate_gds_dims errors on mismatched dimensions", {
  a1 <- array(1:8, dim = c(2, 2, 2))
  expect_error(
    .validate_gds_dims(list(beta = a1), c("s1", "s2", "s3"), c("c1", "c2")),
    "subject|length"
  )
})

test_that(".normalise_col_data creates empty frame for NULL", {
  result <- .normalise_col_data(NULL, c("s1", "s2"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that(".normalise_row_data creates empty frame for NULL", {
  result <- .normalise_row_data(NULL, 5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
})

test_that("add_provenance_node appends to metadata list", {
  meta <- gds_metadata()
  meta2 <- add_provenance_node(meta, "test_op", list(note = "hello"))
  expect_true(length(meta2$provenance$graph) > 0)
})

test_that("subjects.gds_plan returns subject vector", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  p <- as_plan(g)
  expect_equal(subjects(p), c("s1", "s2"))
})

test_that("contrasts.gds_plan returns contrast vector", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  p <- as_plan(g)
  expect_equal(contrasts(p), "c1")
})
