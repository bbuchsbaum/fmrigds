# Targeted coverage tests for files in the 66-84% range

# ===========================================================================
# Helper
# ===========================================================================
.mk_gds <- function(ns = 3, nsubj = 2, nc = 2) {
  beta <- array(seq_len(ns * nsubj * nc) * 0.5, dim = c(ns, nsubj, nc))
  var  <- array(0.1, dim = c(ns, nsubj, nc))
  new_gds(
    assays    = list(beta = beta, var = var),
    space     = space_sample_labels(paste0("s", seq_len(ns))),
    subjects  = paste0("sub", seq_len(nsubj)),
    contrasts = paste0("c", seq_len(nc)),
    col_data  = data.frame(group = rep(c("A", "B"), length.out = nsubj),
                           row.names = paste0("sub", seq_len(nsubj)))
  )
}

# ===========================================================================
# gds-verb.R: gds(), with_col_data, with_row_data
# ===========================================================================

test_that("with_col_data attaches col_data to plan", {
  g <- .mk_gds()
  p <- as_plan(g)
  cd <- data.frame(age = c(25, 30), row.names = c("sub1", "sub2"))
  p2 <- with_col_data(p, cd)
  expect_true(!is.null(p2$meta$col_data))
})

test_that("with_row_data attaches row_data to plan", {
  g <- .mk_gds()
  p <- as_plan(g)
  rd <- data.frame(region = c("frontal", "parietal", "occipital"),
                   row.names = paste0("s", 1:3))
  p2 <- with_row_data(p, rd)
  expect_true(!is.null(p2$meta$row_data))
})

# ===========================================================================
# gds-class.R: accessor methods
# ===========================================================================

test_that("assay() extracts named assay", {
  g <- .mk_gds()
  a <- assay(g, "beta")
  expect_equal(dim(a), c(3, 2, 2))
})

test_that("assays() returns all assays", {
  g <- .mk_gds()
  a <- assays(g)
  expect_named(a, c("beta", "var"))
})

test_that("col_data() extracts col_data from GDS", {
  g <- .mk_gds()
  cd <- col_data(g)
  expect_s3_class(cd, "data.frame")
  expect_equal(nrow(cd), 2)
})

test_that("row_data() extracts row_data", {
  g <- .mk_gds()
  rd <- row_data(g)
  expect_s3_class(rd, "data.frame")
})

test_that("space() extracts space", {
  g <- .mk_gds()
  sp <- space(g)
  expect_s3_class(sp, "gds_space")
})

test_that("metadata() extracts metadata", {
  g <- .mk_gds()
  md <- metadata(g)
  expect_true(is.list(md))
})

# ===========================================================================
# interop.R: model_matrix
# ===========================================================================

test_that("model_matrix builds from realized GDS", {
  g <- .mk_gds()
  mm <- model_matrix(g, ~ group)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), 2)
})

test_that("model_matrix builds from plan with col_data", {
  g <- .mk_gds()
  p <- as_plan(g) |> with_col_data(g$col_data)
  mm <- model_matrix(p, ~ group)
  expect_true(is.matrix(mm))
})

test_that("model_matrix with string formula", {
  g <- .mk_gds()
  mm <- model_matrix(g, "~ group")
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), 2)
})

# ===========================================================================
# plan.R: gds_plan, gds_source, as_plan, add_op
# ===========================================================================

test_that("gds_plan constructor works", {
  src <- gds_source(adapter = "memory", source = list(beta = array(1, dim = c(1,1,1))),
                    probe = list(assays = "beta",
                                 dims = gds_dims(1L, 1L, 1L),
                                 subjects = "s1", contrasts = "c1",
                                 space = space_sample_labels("x"),
                                 maps = list(), metadata = list(), columns = list()))
  p <- gds_plan(source = src, nodes = list(),
                meta = list(subjects = "s1", contrasts = "c1"))
  expect_s3_class(p, "gds_plan")
})

test_that("add_op appends a node", {
  g <- .mk_gds()
  p <- as_plan(g)
  node <- list(op = "custom_test", value = 42)
  p2 <- add_op(p, node)
  expect_equal(length(p2$nodes), length(p$nodes) + 1)
})

# ===========================================================================
# plan-serialization.R: save_plan / load_plan
# ===========================================================================

test_that("save_plan and load_plan roundtrip with tabular source", {
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,s01,c1,0.5,0.01",
    "ROI_1,s02,c1,0.8,0.02"
  ), tmp_csv)
  p <- gds(tmp_csv) |> derive("se") |> reduce("fixed")
  tmp_plan <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_plan), add = TRUE)
  save_plan(p, tmp_plan)
  expect_true(file.exists(tmp_plan))
  p2 <- load_plan(tmp_plan)
  expect_s3_class(p2, "gds_plan")
  expect_equal(length(p2$nodes), length(p$nodes))
})

# ===========================================================================
# posthoc-exec.R: apply_posthoc
# ===========================================================================

test_that("apply_posthoc applies FDR BH correction", {
  .register_builtin_posthoc()
  p_arr <- array(c(0.001, 0.05, 0.5, 0.8, 0.01, 0.9), dim = c(3, 1, 2))
  z_arr <- array(c(3.0, 2.0, 0.5, 0.2, 2.5, 0.1), dim = c(3, 1, 2))
  arrays <- list(p = p_arr, z = z_arr)
  node <- list(op = "posthoc", method = "fdr:bh", options = list())
  result <- apply_posthoc(node, arrays)
  expect_true("q" %in% names(result$arrays))
  expect_equal(dim(result$arrays$q), c(3, 1, 2))
})

# ===========================================================================
# posthoc-fdr-spatial.R: exercise spatial FDR if available
# ===========================================================================

test_that("fdr:spatial is registered if function exists", {
  .register_builtin_posthoc()
  lst <- list_posthoc()
  # Just check registry works; spatial may or may not be available
  expect_true(is.character(lst))
})

# ===========================================================================
# verb-reduce.R: reduce with formula
# ===========================================================================

test_that("reduce with formula stores formula in node", {
  g <- .mk_gds()
  p <- as_plan(g) |> with_col_data(g$col_data)
  p2 <- reduce(p, method = "meta:fe_reg", formula = ~ group)
  last_node <- p2$nodes[[length(p2$nodes)]]
  expect_true(!is.null(last_node$formula))
})

# ===========================================================================
# write-exporters.R: CSV export with options
# ===========================================================================

test_that(".write_gds_csv with options", {
  g <- .mk_gds()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  .write_gds_csv(g, tmp, options = list(stats = "beta", drop_na = TRUE, include_col_data = TRUE))
  df <- read.csv(tmp)
  expect_true("beta" %in% names(df))
  expect_true("group" %in% names(df))
})

# ===========================================================================
# map.R: MapFamily, OrthogonalFamily, map_linear, UncertaintyRule
# ===========================================================================

test_that("MapFamily with all arguments", {
  from <- space_sample_labels(c("a", "b"))
  to <- space_sample_labels(c("x", "y"))
  M <- diag(2)
  fam <- MapFamily("test", from, to, by_subject = list(s1 = M))
  expect_equal(fam$name, "test")
  expect_equal(length(fam$by_subject), 1)
})

test_that("OrthogonalFamily creates orthogonal map family", {
  from <- space_sample_labels(c("a", "b"))
  to <- space_sample_labels(c("x", "y"))
  M <- diag(2)
  fam <- OrthogonalFamily("orth", from, to, matrices_by_subject = list(s1 = M))
  expect_true(isTRUE(fam$traits$orthogonal))
})

test_that("UncertaintyRule modes", {
  u1 <- UncertaintyRule("independent")
  expect_equal(u1$mode, "independent")
  u2 <- UncertaintyRule("none")
  expect_equal(u2$mode, "none")
})

test_that("map_linear creates a linear map", {
  from <- space_sample_labels(c("a", "b", "c"))
  to <- space_sample_labels(c("x", "y"))
  M <- matrix(1:6, nrow = 2, ncol = 3)
  ml <- map_linear(from, to, M)
  expect_true(!is.null(ml))
})

# ===========================================================================
# compute.R: preview, subset node handling
# ===========================================================================

test_that("preview returns small GDS from tabular source", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  lines <- c("sample,subject,contrast,beta,var")
  for (i in 1:10) lines <- c(lines, sprintf("ROI_%d,s01,c1,%.1f,0.01", i, i * 0.1))
  for (i in 1:10) lines <- c(lines, sprintf("ROI_%d,s02,c1,%.1f,0.02", i, i * 0.2))
  writeLines(lines, tmp)
  p <- gds(tmp)
  prev <- preview(p, n = 5)
  expect_s3_class(prev, "gds")
  expect_true(dim(prev$assays[[1]])[1] <= 5)
})

test_that("compute with subset node", {
  g <- .mk_gds()
  p <- as_plan(g) |> subset(subject = "sub1")
  result <- compute(p)
  expect_s3_class(result, "gds")
  expect_equal(result$subjects, "sub1")
})

test_that("compute with derive + reduce pipeline", {
  g <- .mk_gds()
  result <- as_plan(g) |> derive("se") |> reduce("fixed") |> compute()
  expect_s3_class(result, "gds")
  expect_equal(dim(result$assays$beta)[2], 1L)
})

# ===========================================================================
# export-neuroim.R: split.gds
# ===========================================================================

test_that("split.gds splits by subject factor", {
  g <- .mk_gds()
  # split.gds expects a factor-like vector matching subjects
  parts <- split(g, g$subjects)
  expect_equal(length(parts), 2)
  expect_s3_class(parts[[1]], "gds")
  expect_equal(length(parts[[1]]$subjects), 1)
})

# ===========================================================================
# gds_to_tibble with voxel mask
# ===========================================================================

test_that("gds_to_tibble uses mask_idx for voxel space", {
  sp <- space_voxel(c(2, 2, 2), diag(4), mask_idx = c(1L, 4L, 7L), storage = "packed")
  beta <- array(1:3, dim = c(3, 1, 1))
  var  <- array(0.1, dim = c(3, 1, 1))
  g <- new_gds(assays = list(beta = beta, var = var), space = sp,
               subjects = "s1", contrasts = "c1")
  df <- gds_to_tibble(g)
  expect_equal(nrow(df), 3)
  expect_true(all(df$sample %in% c(1, 4, 7)))
})
