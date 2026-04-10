# Tests for export-public.R, write-exporters.R, plan-alias.R

# Helper: minimal GDS for export tests
.make_export_gds <- function() {
  beta <- array(1:12, dim = c(3, 2, 2))
  var  <- array(rep(0.1, 12), dim = c(3, 2, 2))
  new_gds(
    assays    = list(beta = beta, var = var),
    space     = space_sample_labels(c("a", "b", "c")),
    subjects  = c("s1", "s2"),
    contrasts = c("c1", "c2"),
    col_data  = data.frame(group = c("ctrl", "trt"), row.names = c("s1", "s2"))
  )
}

# ===========================================================================
# gds_to_tibble / .gds_tidy_frame
# ===========================================================================

test_that("gds_to_tibble returns long-form data.frame", {
  g <- .make_export_gds()
  df <- gds_to_tibble(g)
  expect_s3_class(df, "data.frame")
  # 3 samples x 2 subjects x 2 contrasts = 12 rows
  expect_equal(nrow(df), 12)
  expect_true("sample" %in% names(df))
  expect_true("subject" %in% names(df))
  expect_true("contrast" %in% names(df))
  expect_true("beta" %in% names(df))
  expect_true("var" %in% names(df))
})

test_that("gds_to_tibble respects assay selection", {
  g <- .make_export_gds()
  df <- gds_to_tibble(g, assays = "beta")
  expect_true("beta" %in% names(df))
  expect_false("var" %in% names(df))
})

test_that("gds_to_tibble include_col_data joins subject metadata", {
  g <- .make_export_gds()
  df <- gds_to_tibble(g, include_col_data = TRUE)
  expect_true("group" %in% names(df))
})

test_that("gds_to_tibble drop_na removes all-NA rows", {
  beta <- array(c(1, NA, 3, 4, NA, 6, 7, NA, 9, 10, NA, 12), dim = c(3, 2, 2))
  var  <- array(c(0.1, NA, 0.1, 0.1, NA, 0.1, 0.1, NA, 0.1, 0.1, NA, 0.1), dim = c(3, 2, 2))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = c("s1", "s2"), contrasts = c("c1", "c2")
  )
  df <- gds_to_tibble(g, drop_na = TRUE)
  expect_true(nrow(df) < 12)
})

test_that("gds_to_tibble errors with no matching assays", {
  g <- .make_export_gds()
  expect_error(gds_to_tibble(g, assays = "nonexistent"), "No matching assays")
})

test_that(".gds_tidy_frame uses parcel labels for sample column", {
  beta <- array(1:4, dim = c(2, 2, 1))
  var  <- array(0.1, dim = c(2, 2, 1))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_parcels(c("ROI_A", "ROI_B")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  df <- gds_to_tibble(g)
  expect_true(all(c("ROI_A", "ROI_B") %in% df$sample))
})

# ===========================================================================
# .write_gds_csv
# ===========================================================================

test_that(".write_gds_csv writes a valid CSV file", {
  g <- .make_export_gds()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  .write_gds_csv(g, tmp)
  expect_true(file.exists(tmp))
  df <- read.csv(tmp)
  expect_equal(nrow(df), 12)
})

# ===========================================================================
# .write_export dispatches correctly
# ===========================================================================

test_that(".write_export dispatches to csv", {
  g <- .make_export_gds()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  .write_export(g, "csv", tmp, list())
  expect_true(file.exists(tmp))
})

test_that(".write_export errors for unsupported format", {
  g <- .make_export_gds()
  expect_error(.write_export(g, "xyz", "out.xyz", list()), "Unsupported")
})

# ===========================================================================
# plan() alias
# ===========================================================================

test_that("plan() is an alias for as_plan()", {
  g <- .make_export_gds()
  p <- plan(g)
  expect_s3_class(p, "gds_plan")
})
