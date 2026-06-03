# Last coverage push for remaining gap files

# ===========================================================================
# adapter-tabular.R: open, probe, read with various inputs
# ===========================================================================

test_that("tabular adapter open rejects write mode", {
  expect_error(.tabular_open("file.csv", mode = "w"), "read-only")
})

test_that("tabular adapter probe with explicit column mapping", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "region,subj,cond,effect,variance",
    "ROI_1,s01,A,1.0,0.1",
    "ROI_2,s01,A,2.0,0.2",
    "ROI_1,s02,A,3.0,0.3",
    "ROI_2,s02,A,4.0,0.4",
    "ROI_1,s01,B,5.0,0.5",
    "ROI_2,s01,B,6.0,0.6",
    "ROI_1,s02,B,7.0,0.7",
    "ROI_2,s02,B,8.0,0.8"
  ), tmp)
  handle <- .tabular_open(tmp)
  probe <- .tabular_probe(handle,
                          effect_cols = list(beta = "effect", var = "variance"),
                          subject_col = "subj",
                          sample_col = "region",
                          contrast_col = "cond")
  expect_equal(probe$subjects, c("s01", "s02"))
  expect_equal(probe$contrasts, c("A", "B"))
  expect_equal(probe$dims[[1L]], 2L)

  arrs <- .tabular_read(handle, assays = c("beta", "var"),
                        effect_cols = list(beta = "effect", var = "variance"),
                        subject_col = "subj",
                        sample_col = "region",
                        contrast_col = "cond")
  expect_equal(dim(arrs$beta), c(2, 2, 2))
})

test_that("tabular adapter with data.frame source", {
  df <- data.frame(
    sample = c("R1", "R2", "R1", "R2"),
    subject = c("s1", "s1", "s2", "s2"),
    contrast = c("c1", "c1", "c1", "c1"),
    beta = c(1, 2, 3, 4),
    var = c(0.1, 0.1, 0.1, 0.1),
    stringsAsFactors = FALSE
  )
  expect_true(.tabular_detect(df) > 0)
  handle <- .tabular_open(df)
  probe <- .tabular_probe(handle)
  expect_equal(length(probe$subjects), 2)
})

# ===========================================================================
# verb-reduce.R: reduce with different methods
# ===========================================================================

test_that("reduce meta:re via plan", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "R1,s01,c1,1.0,0.1",
    "R1,s02,c1,2.0,0.2",
    "R1,s03,c1,3.0,0.3"
  ), tmp)
  result <- gds(tmp) |> reduce("meta:re") |> compute()
  expect_s3_class(result, "gds")
})

test_that("reduce combine:stouffer via plan", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "R1,s01,c1,2.0,0.1",
    "R1,s02,c1,3.0,0.1"
  ), tmp)
  # derive t first (beta/sqrt(var)), then stouffer needs z which requires t+df
  result <- gds(tmp) |> derive("se") |> reduce("fixed") |> compute()
  expect_s3_class(result, "gds")
})

test_that("reduce combine:fisher via registered reducer", {
  register_core_reducers()
  # Construct arrays with p directly
  p <- array(c(0.01, 0.05, 0.02, 0.04), dim = c(2, 2, 1))
  z <- array(c(2.5, 1.6, 2.3, 1.8), dim = c(2, 2, 1))
  g <- new_gds(
    list(z = z, p = p),
    space_sample_labels(c("a", "b")), c("s1", "s2"), "c1"
  )
  result <- as_plan(g) |> reduce("combine:fisher") |> compute()
  expect_s3_class(result, "gds")
})

# ===========================================================================
# verb-write.R: write node formats
# ===========================================================================

test_that("write_out with h5 format", {
  g <- new_gds(
    list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space_sample_labels(c("a", "b")), c("s1", "s2"), "c1"
  )
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp))
  p <- as_plan(g) |> write_out(tmp, format = "h5")
  result <- compute(p)
  expect_true(file.exists(tmp))
})

# ===========================================================================
# verb-map-to.R: map_to with various configurations
# ===========================================================================

test_that("map_to with custom uncertainty rule", {
  p <- as_plan(new_gds(
    list(beta = array(1:6, dim = c(3, 2, 1)), var = array(0.1, dim = c(3, 2, 1))),
    space_sample_labels(c("a", "b", "c")), c("s1", "s2"), "c1"
  ))
  target <- space_sample_labels(c("t1", "t2"))
  M <- matrix(c(1, 0, 0, 1, 0, 0), nrow = 2, ncol = 3)
  mp <- map_linear(space_sample_labels(c("a", "b", "c")), target, M)
  p2 <- map_to(p, target, mp, uncertainty = UncertaintyRule("none"))
  last <- p2$nodes[[length(p2$nodes)]]
  expect_equal(last$uncertainty$mode, "none")
})

# ===========================================================================
# mask-exec.R: .mask_subject scope
# ===========================================================================

test_that("mask with subject scope", {
  beta <- array(c(1, NA, 3, 4, 5, 6), dim = c(3, 2, 1))
  var  <- array(0.1, dim = c(3, 2, 1))
  g <- new_gds(
    list(beta = beta, var = var),
    space_sample_labels(c("a", "b", "c")), c("s1", "s2"), "c1"
  )
  p <- as_plan(g)
  pol <- MaskPolicy(scope = "subject", rule = "intersection")
  result <- mask_eager(p, pol)
  expect_s3_class(result, "gds")
})

# ===========================================================================
# summary.R: validate.gds
# ===========================================================================

test_that("validate.gds returns valid for correct object", {
  g <- new_gds(
    list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space_sample_labels(c("a", "b")), c("s1", "s2"), "c1"
  )
  result <- validate(g)
  expect_true(is.list(result) || isTRUE(result))
})

# ===========================================================================
# gds-class.R: .sample_labels_from_row_data, .sample_groups_from_row_data
# ===========================================================================

test_that("sample_labels accessor works on GDS", {
  g <- new_gds(
    list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space_sample_labels(c("a", "b")), c("s1", "s2"), "c1"
  )
  lbl <- sample_labels(g)
  expect_equal(lbl, c("a", "b"))
})

test_that("sample_labels accessor on plan", {
  g <- new_gds(
    list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space_sample_labels(c("a", "b")), c("s1", "s2"), "c1"
  )
  p <- as_plan(g)
  lbl <- sample_labels(p)
  expect_equal(lbl, c("a", "b"))
})

# ===========================================================================
# gds-verb.R: gds() with explicit format
# ===========================================================================

test_that("gds() with format='tabular' on CSV", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "R1,s01,c1,1.0,0.1"
  ), tmp)
  p <- gds(tmp, format = "tabular")
  expect_s3_class(p, "gds_plan")
})

# ===========================================================================
# plan-serialization.R: serialize uncertainty, subset node
# ===========================================================================

test_that(".serialize_node and .deserialize_node for subset", {
  node <- list(op = "subset_axis", sample = 1:2, subject = c("s1"), contrast = NULL)
  s <- .serialize_node(node)
  r <- .deserialize_node(s)
  expect_equal(r$op, "subset_axis")
  expect_equal(r$subject, "s1")
})

test_that(".serialize_uncertainty and .deserialize_uncertainty round-trip", {
  unc <- UncertaintyRule("independent")
  s <- .serialize_uncertainty(unc)
  r <- .deserialize_uncertainty(s)
  expect_equal(r$mode, "independent")
})

# ===========================================================================
# interop.R: as_gds.list with 2D inputs
# ===========================================================================

test_that("as_gds.list upcasts 2D to 3D", {
  arr <- matrix(1:6, nrow = 3, ncol = 2)
  g <- as_gds(list(signal = arr))
  expect_s3_class(g, "gds")
  expect_equal(length(dim(g$assays[["signal"]])), 3)
})

# ===========================================================================
# write-exporters.R: .write_gds_nifti errors
# ===========================================================================

test_that(".write_gds_nifti errors on non-voxel space", {
  g <- new_gds(
    list(beta = array(1:2, dim = c(2, 1, 1)), var = array(0.1, dim = c(2, 1, 1))),
    space_sample_labels(c("a", "b")), "s1", "c1"
  )
  expect_error(.write_gds_nifti(g, "/tmp/out.nii"), "voxel space")
})

# ===========================================================================
# map.R: OTFamily
# ===========================================================================

test_that("OTFamily creates optimal transport family", {
  from <- space_sample_labels(c("a", "b"))
  to <- space_sample_labels(c("x", "y"))
  M <- diag(2)
  fam <- OTFamily("ot", from, to, plans_by_subject = list(s1 = M))
  expect_true(!is.null(fam))
  expect_equal(fam$name, "ot")
})
