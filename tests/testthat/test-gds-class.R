test_that("new_gds validates assays and dimensions", {
  beta <- array(rnorm(24), dim = c(4, 3, 2))
  var <- array(runif(24, 0.1, 1), dim = c(4, 3, 2))
  space <- structure(list(type = "dummy", storage = "dense", dim = c(2, 2, 1)), class = "gds_space")

  expect_error(new_gds(list(beta = beta), space, letters[1:3], letters[1:2]),
               "either 'var' or 'se'", fixed = TRUE)

  gds <- new_gds(
    assays = list(beta = beta, var = var),
    space = space,
    subjects = letters[1:3],
    contrasts = paste0("c", 1:2)
  )

  expect_s3_class(gds, "gds")
  expect_equal(dim(assay(gds, "beta")), c(4, 3, 2))
  expect_equal(subjects(gds), letters[1:3])
  expect_equal(contrasts(gds), paste0("c", 1:2))
  expect_equal(nrow(col_data(gds)), 3)
  expect_equal(nrow(row_data(gds)), 4)
})

test_that("provenance builder appends entries", {
  meta <- gds_metadata()
  meta <- add_provenance_node(meta, "subset", list(subject = "s01"))

  expect_length(meta$provenance$graph, 1)
  expect_match(meta$provenance$log[[1]], "subset(")
})

test_that("space constructors validate inputs", {
  affine <- diag(4)
  vox <- space_voxel(c(2, 2, 1), affine)
  expect_s3_class(vox, "space_voxel")
  expect_equal(vox$storage, "dense")

  mask <- array(TRUE, dim = c(2, 2, 1))
  vox_packed <- space_voxel(c(2, 2, 1), affine, mask_bitmap = mask, storage = "packed")
  expect_equal(length(vox_packed$mask_idx), sum(mask))

  parcels <- space_parcels(c("A", "B"))
  expect_s3_class(parcels, "space_parcels")

  surf <- space_surface(matrix(0, nrow = 3, ncol = 3), matrix(c(1, 2, 3), nrow = 1), "L")
  expect_s3_class(surf, "space_surface")

  basis <- space_basis(5)
  expect_s3_class(basis, "space_basis")
})

test_that("map constructors enforce inputs", {
  from <- space_basis(3)
  to <- space_basis(3)

  m <- matrix(diag(3), nrow = 3)
  lin <- map_linear(from, to, operator = m)
  expect_s3_class(lin, "map_linear")

  fam <- MapFamily("test", from, to, type = "orthogonal", by_subject = list(sub01 = m))
  expect_equal(fam$type, "orthogonal")

  ortho <- OrthogonalFamily("o", from, to, list(sub01 = m))
  expect_true(ortho$traits$orthogonal)

  expect_error(UncertaintyRule("cov_provider"), "cov_provider", fixed = TRUE)
})

test_that("plan infrastructure records operations", {
  src <- gds_source("tabular", list(path = "dummy"))
  plan <- gds_plan(src)
  plan <- subset(plan, subject = "s01")
  expect_s3_class(plan, "gds_plan")
  expect_length(plan$nodes, 1)
  expect_equal(plan$nodes[[1]]$subject, "s01")
})

test_that("assay registry defaults are installed", {
  expect_true(can_map_linear("beta"))
  info <- assay_info("var")
  expect_equal(info$variance_of, "beta")
})

test_that("derive verb appends op", {
  src <- gds_source("tabular", list())
  plan <- gds_plan(src)
  plan <- derive(plan, what = c("var", "t"))
  expect_equal(plan$nodes[[1]]$what, c("var", "t"))
})

test_that("tabular adapter roundtrip", {
  testthat::skip_if_not_installed("data.table")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  df <- data.frame(
    sample = rep(c("roi1", "roi2"), each = 4),
    subject = rep(c("s1", "s2"), times = 4),
    contrast = rep(c("c1", "c2"), each = 2, times = 2),
    beta = rnorm(8),
    var = runif(8, 0.1, 1)
  )
  data.table::fwrite(df, tmp)

  plan <- gds(tmp)
  g <- compute(plan)

  expect_s3_class(g, "gds")
  expect_equal(dim(assay(g, "beta")), c(2, 2, 2))
  expect_equal(subjects(g), c("s1", "s2"))
  expect_equal(contrasts(g), c("c1", "c2"))
})
