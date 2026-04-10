# Final coverage push: target files in the 58-85% range

# ===========================================================================
# map-exec.R: .map_with_effect_scale, .map_without_effect_scale, combine paths
# ===========================================================================

test_that("map_to with effect_scale map applies M to beta and propagates variance", {
  beta <- array(c(1, 2, 3, 4, 5, 6), dim = c(3, 2, 1))
  var  <- array(0.1, dim = c(3, 2, 1))
  se   <- array(sqrt(0.1), dim = c(3, 2, 1))
  arrays <- list(beta = beta, var = var, se = se)

  M <- matrix(c(0.5, 0.5, 0, 0, 0.5, 0.5), nrow = 2, ncol = 3, byrow = TRUE)
  target <- space_sample_labels(c("t1", "t2"))
  mp <- map_linear(space_sample_labels(c("a", "b", "c")), target, M)
  unc <- UncertaintyRule("independent")
  node <- list(op = "map", target_space = target, map = mp, uncertainty = unc)
  result <- apply_map_to(node, arrays)
  expect_equal(dim(result$arrays$beta), c(2, 2, 1))
  # Variance should be propagated: M diag(v) M'
  expect_true("var" %in% names(result$arrays))
})

test_that("map_to with z-scores uses combine stouffer path", {
  z <- array(c(2, 1, 0.5, 3, 2, 1), dim = c(3, 2, 1))
  p <- array(2 * pnorm(-abs(z)), dim = c(3, 2, 1))
  arrays <- list(z = z, p = p)

  M <- matrix(c(0.5, 0.5, 0, 0, 0.5, 0.5), nrow = 2, ncol = 3, byrow = TRUE)
  target <- space_sample_labels(c("t1", "t2"))
  from_sp <- space_sample_labels(c("a", "b", "c"))
  mp <- map_linear(from_sp, target, M)
  unc <- UncertaintyRule("none")
  node <- list(op = "map", target_space = target, map = mp, uncertainty = unc,
               combine = "stouffer")
  result <- apply_map_to(node, arrays)
  expect_equal(dim(result$arrays$z), c(2, 2, 1))
})

# ===========================================================================
# mask-exec.R: .finite_mask, threshold rule
# ===========================================================================

test_that("mask with threshold rule", {
  beta <- array(c(1, NA, 3, 4, NA, 6, 1, 2, 3, 4, 5, 6), dim = c(3, 2, 2))
  var  <- array(0.1, dim = c(3, 2, 2))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = c("s1", "s2"), contrasts = c("c1", "c2")
  )
  p <- as_plan(g)
  pol <- MaskPolicy(scope = "group", rule = "threshold", threshold = 0.5)
  result <- mask_eager(p, pol)
  expect_s3_class(result, "gds")
})

# ===========================================================================
# adapter-tabular.R: tabular adapter with different column patterns
# ===========================================================================

test_that("tabular adapter detects CSV files", {
  expect_true(.tabular_detect("data.csv") > 0)
  expect_true(.tabular_detect("data.tsv") > 0)
  expect_true(.tabular_detect(data.frame(x = 1)) > 0)
  expect_false(.tabular_detect(123))
})

test_that("tabular adapter probe and read with minimal CSV", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "R1,s01,c1,1.0,0.1",
    "R2,s01,c1,2.0,0.2",
    "R1,s02,c1,3.0,0.3",
    "R2,s02,c1,4.0,0.4"
  ), tmp)
  handle <- .tabular_open(tmp)
  probe <- .tabular_probe(handle)
  expect_equal(probe$subjects, c("s01", "s02"))
  expect_equal(probe$contrasts, "c1")

  arrs <- .tabular_read(handle, assays = "beta",
                        effect_cols = probe$columns$effect_cols,
                        subject_col = probe$columns$subject_col,
                        sample_col = probe$columns$sample_col,
                        contrast_col = probe$columns$contrast_col)
  expect_equal(dim(arrs$beta), c(2, 2, 1))
})

test_that("tabular adapter close is no-op", {
  expect_null(.tabular_close(list()))
})

# ===========================================================================
# verb-reduce.R: reduce with options
# ===========================================================================

test_that("reduce with options passes through", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "R1,s01,c1,1.0,0.1",
    "R1,s02,c1,2.0,0.2"
  ), tmp)
  p <- gds(tmp) |> reduce("fixed", weights = "equal")
  result <- compute(p)
  expect_s3_class(result, "gds")
})

# ===========================================================================
# verb-write.R: write_out with CSV sink
# ===========================================================================

test_that("write_out node produces file when computed", {
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv))
  beta <- array(1:4, dim = c(2, 2, 1))
  var  <- array(0.1, dim = c(2, 2, 1))
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("a","b")),
               c("s1","s2"), "c1")
  p <- as_plan(g) |> write_out(tmp_csv, format = "csv")
  result <- compute(p)
  expect_true(file.exists(tmp_csv))
})

# ===========================================================================
# summary.R: .space_brief for all space types
# ===========================================================================

test_that(".space_brief describes voxel space", {
  sp <- space_voxel(c(2, 3, 4), diag(4))
  brief <- .space_brief(sp)
  expect_true(grepl("[Vv]oxel", brief))
})

test_that(".space_brief describes parcels space", {
  sp <- space_parcels(c("A", "B", "C"))
  brief <- .space_brief(sp)
  expect_true(grepl("[Pp]arcel", brief))
})

test_that(".space_brief describes sample_labels space", {
  sp <- space_sample_labels(c("x", "y"))
  brief <- .space_brief(sp)
  expect_true(nzchar(brief))
})

# ===========================================================================
# adapter-nifti.R: detection
# ===========================================================================

test_that("nifti_source constructs a source spec", {
  src <- nifti_source(beta = "/tmp/test.nii")
  expect_true(is.list(src))
  expect_true("beta" %in% names(src))
})

# ===========================================================================
# gds-class.R: .validate_gds_space
# ===========================================================================

test_that(".validate_gds_space validates sample count", {
  sp <- space_sample_labels(c("a", "b", "c"))
  expect_silent(.validate_gds_space(sp, 3L))
})

# ===========================================================================
# plan-serialization.R: serialize/deserialize spaces
# ===========================================================================

test_that(".serialize_space and .deserialize_space round-trip for sample_labels", {
  sp <- space_sample_labels(c("a", "b", "c"))
  serialized <- .serialize_space(sp)
  restored <- .deserialize_space(serialized)
  expect_s3_class(restored, "space_sample_labels")
  expect_equal(restored$labels, c("a", "b", "c"))
})

test_that(".serialize_space and .deserialize_space round-trip for parcels", {
  sp <- space_parcels(c("P1", "P2"))
  serialized <- .serialize_space(sp)
  restored <- .deserialize_space(serialized)
  expect_s3_class(restored, "space_parcels")
  expect_equal(restored$labels, c("P1", "P2"))
})

test_that(".serialize_space and .deserialize_space round-trip for voxel", {
  sp <- space_voxel(c(2, 2, 2), diag(4))
  serialized <- .serialize_space(sp)
  restored <- .deserialize_space(serialized)
  expect_s3_class(restored, "space_voxel")
  expect_equal(restored$dim, c(2L, 2L, 2L))
})

# ===========================================================================
# plan-serialization.R: serialize/deserialize nodes
# ===========================================================================

test_that(".serialize_node and .deserialize_node round-trip for derive", {
  node <- list(op = "derive", what = "se", options = list())
  s <- .serialize_node(node)
  r <- .deserialize_node(s)
  expect_equal(r$op, "derive")
  expect_equal(r$what, "se")
})

test_that(".serialize_node and .deserialize_node for reduce", {
  node <- list(op = "reduce", method = "fixed", weights = "1/var", options = list())
  s <- .serialize_node(node)
  r <- .deserialize_node(s)
  expect_equal(r$op, "reduce")
  expect_equal(r$method, "fixed")
})

# ===========================================================================
# propagate-variance.R
# ===========================================================================

test_that("propagate_variance_independent handles identity map", {
  M <- diag(3)
  beta_arr <- array(c(1, 2, 3), dim = c(3, 1, 1))
  var_arr  <- array(c(0.1, 0.2, 0.3), dim = c(3, 1, 1))
  result <- propagate_variance_independent(M, beta_arr, var_arr)
  expect_equal(dim(result$beta), c(3, 1, 1))
  expect_equal(dim(result$var), c(3, 1, 1))
  # Identity map should preserve values
  expect_equal(as.vector(result$beta), c(1, 2, 3))
})

test_that("propagate_variance_independent with dimension reduction", {
  M <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)
  beta_arr <- array(c(1, 2, 3, 4, 5, 6), dim = c(3, 2, 1))
  var_arr  <- array(0.1, dim = c(3, 2, 1))
  result <- propagate_variance_independent(M, beta_arr, var_arr)
  expect_equal(dim(result$beta), c(2, 2, 1))
})

# ===========================================================================
# compat.R: relabel_subjects
# ===========================================================================

test_that("relabel_subjects renames subjects in GDS", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("old1", "old2"), contrasts = "c1"
  )
  g2 <- relabel_subjects(g, c(old1 = "new1", old2 = "new2"))
  expect_equal(g2$subjects, c("new1", "new2"))
})

# ===========================================================================
# export-public.R: gds_to_tibble edge cases
# ===========================================================================

test_that("gds_to_tibble with drop_na filters incomplete rows", {
  beta <- array(c(1, NA), dim = c(2, 1, 1))
  var  <- array(c(0.1, NA), dim = c(2, 1, 1))
  g <- new_gds(list(beta = beta, var = var), space_sample_labels(c("a","b")),
               "s1", "c1")
  df <- gds_to_tibble(g, drop_na = FALSE)
  expect_equal(nrow(df), 2)
  df2 <- gds_to_tibble(g, drop_na = TRUE)
  expect_true(nrow(df2) < 2)
})
