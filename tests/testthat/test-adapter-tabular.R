# Test tabular adapter (CSV, TSV, Parquet)

test_that("tabular adapter detects CSV files", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines("sample,subject,contrast,beta,var", tmp)

  score <- .tabular_detect(tmp)
  expect_true(is.numeric(score))
  expect_equal(score, 0.8)
})

test_that("tabular adapter detects TSV files", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp))

  writeLines("sample\tsubject\tcontrast\tbeta", tmp)

  score <- .tabular_detect(tmp)
  expect_true(is.numeric(score))
  expect_equal(score, 0.8)
})

test_that("tabular adapter detects Parquet files", {
  skip_if_not_installed("arrow")

  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp))

  df <- data.frame(
    sample = c("ROI_1", "ROI_2"),
    subject = c("sub-01", "sub-01"),
    contrast = c("cope1", "cope1"),
    beta = c(0.5, 0.8)
  )

  arrow::write_parquet(df, tmp)

  score <- .tabular_detect(tmp)
  expect_true(is.numeric(score))
  expect_equal(score, 0.8)
})

test_that("tabular adapter can open CSV", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_2,sub-01,cope1,0.8,0.02"
  ), tmp)

  handle <- .tabular_open(tmp, mode = "r")

  expect_type(handle, "list")
  expect_true("path" %in% names(handle))
  expect_equal(handle$path, tmp)
})

test_that("tabular adapter can probe CSV with default columns", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_2,sub-01,cope1,0.8,0.02",
    "ROI_1,sub-02,cope1,0.6,0.015",
    "ROI_2,sub-02,cope1,0.9,0.025"
  ), tmp)

  handle <- .tabular_open(tmp)
  probe_result <- .tabular_probe(handle)

  expect_type(probe_result, "list")
  expect_true("assays" %in% names(probe_result))
  expect_true("dims" %in% names(probe_result))
  expect_true("subjects" %in% names(probe_result))
  expect_true("contrasts" %in% names(probe_result))
  expect_true("space" %in% names(probe_result))

  # Check dimensions
  expect_equal(unname(probe_result$dims["sample"]), 2)    # ROI_1, ROI_2
  expect_equal(unname(probe_result$dims["subject"]), 2)   # sub-01, sub-02
  expect_equal(unname(probe_result$dims["contrast"]), 1)  # cope1

  # Check detected assays
  expect_true("beta" %in% probe_result$assays)
  expect_true("var" %in% probe_result$assays)

  # Check space type
  expect_s3_class(probe_result$space, "Space")
  expect_equal(probe_result$space$type, "sample_labels")
})

test_that("tabular adapter can probe with custom column mappings", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "roi,subj_id,task,effect,se",
    "ROI_1,s01,rest,0.5,0.1",
    "ROI_2,s01,rest,0.8,0.15"
  ), tmp)

  handle <- .tabular_open(tmp)
  probe_result <- .tabular_probe(
    handle,
    effect_cols = list(beta = "effect", se = "se"),
    subject_col = "subj_id",
    sample_col = "roi",
    contrast_col = "task"
  )

  expect_equal(unname(probe_result$dims["sample"]), 2)
  expect_equal(unname(probe_result$dims["subject"]), 1)
  expect_equal(probe_result$contrasts, "rest")
  expect_true("beta" %in% probe_result$assays)
  expect_true("se" %in% probe_result$assays)
})

test_that("tabular adapter can read full data", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_2,sub-01,cope1,0.8,0.02",
    "ROI_1,sub-02,cope1,0.6,0.015",
    "ROI_2,sub-02,cope1,0.9,0.025"
  ), tmp)

  handle <- .tabular_open(tmp)
  probe_result <- .tabular_probe(handle)

  data <- .tabular_read(
    handle = handle,
    assays = c("beta", "var"),
    block = NULL,
    effect_cols = NULL,
    subject_col = "subject",
    sample_col = "sample",
    contrast_col = "contrast"
  )

  expect_type(data, "list")
  expect_true("beta" %in% names(data))
  expect_true("var" %in% names(data))

  # Check dimensions [sample, subject, contrast]
  expect_equal(dim(data$beta), c(2, 2, 1))
  expect_equal(dim(data$var), c(2, 2, 1))

  # Check values
  expect_equal(data$beta[1, 1, 1], 0.5)   # ROI_1, sub-01, cope1
  expect_equal(data$beta[2, 1, 1], 0.8)   # ROI_2, sub-01, cope1
  expect_equal(data$beta[1, 2, 1], 0.6)   # ROI_1, sub-02, cope1
  expect_equal(data$beta[2, 2, 1], 0.9)   # ROI_2, sub-02, cope1
})

test_that("tabular adapter handles multiple contrasts", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta",
    "ROI_1,sub-01,cope1,0.5",
    "ROI_1,sub-01,cope2,0.6",
    "ROI_2,sub-01,cope1,0.8",
    "ROI_2,sub-01,cope2,0.9"
  ), tmp)

  handle <- .tabular_open(tmp)
  probe_result <- .tabular_probe(handle)

  expect_equal(unname(probe_result$dims["contrast"]), 2)
  expect_equal(probe_result$contrasts, c("cope1", "cope2"))

  data <- .tabular_read(
    handle = handle,
    assays = "beta",
    block = NULL,
    subject_col = "subject",
    sample_col = "sample",
    contrast_col = "contrast"
  )

  expect_equal(dim(data$beta), c(2, 1, 2))  # 2 samples, 1 subject, 2 contrasts
})

test_that("tabular adapter can read with block (contrast subset)", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta",
    "ROI_1,sub-01,cope1,0.5",
    "ROI_1,sub-01,cope2,0.6",
    "ROI_2,sub-01,cope1,0.8",
    "ROI_2,sub-01,cope2,0.9"
  ), tmp)

  handle <- .tabular_open(tmp)
  probe_result <- .tabular_probe(handle)

  # Read only first contrast
  data <- .tabular_read(
    handle = handle,
    assays = "beta",
    block = list(contrast = 1),
    subject_col = "subject",
    sample_col = "sample",
    contrast_col = "contrast"
  )

  expect_equal(dim(data$beta)[3], 1)  # Only 1 contrast
})

test_that("tabular adapter can read with block (sample subset)", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta",
    "ROI_1,sub-01,cope1,0.5",
    "ROI_2,sub-01,cope1,0.8",
    "ROI_3,sub-01,cope1,1.0"
  ), tmp)

  handle <- .tabular_open(tmp)

  # Read only first 2 samples
  data <- .tabular_read(
    handle = handle,
    assays = "beta",
    block = list(sample = 1:2),
    subject_col = "subject",
    sample_col = "sample",
    contrast_col = "contrast"
  )

  expect_equal(dim(data$beta)[1], 2)  # Only 2 samples
})

test_that("tabular adapter handles missing values", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_2,sub-01,cope1,,0.02",  # Missing beta
    "ROI_1,sub-02,cope1,0.6,",   # Missing var
    "ROI_2,sub-02,cope1,0.9,0.025"
  ), tmp)

  handle <- .tabular_open(tmp)

  data <- .tabular_read(
    handle = handle,
    assays = c("beta", "var"),
    block = NULL,
    subject_col = "subject",
    sample_col = "sample",
    contrast_col = "contrast"
  )

  # Should have NAs
  expect_true(is.na(data$beta[2, 1, 1]))  # ROI_2, sub-01, cope1
  expect_true(is.na(data$var[1, 2, 1]))   # ROI_1, sub-02, cope1
})

test_that("tabular adapter works with Parquet format", {
  skip_if_not_installed("arrow")

  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp))

  df <- data.frame(
    sample = c("ROI_1", "ROI_2", "ROI_1", "ROI_2"),
    subject = c("sub-01", "sub-01", "sub-02", "sub-02"),
    contrast = c("cope1", "cope1", "cope1", "cope1"),
    beta = c(0.5, 0.8, 0.6, 0.9),
    var = c(0.01, 0.02, 0.015, 0.025),
    stringsAsFactors = FALSE
  )

  arrow::write_parquet(df, tmp)

  handle <- .tabular_open(tmp)
  probe_result <- .tabular_probe(handle)

  expect_equal(unname(probe_result$dims["sample"]), 2)
  expect_equal(unname(probe_result$dims["subject"]), 2)

  data <- .tabular_read(
    handle = handle,
    assays = c("beta", "var"),
    block = NULL,
    subject_col = "subject",
    sample_col = "sample",
    contrast_col = "contrast"
  )

  expect_equal(dim(data$beta), c(2, 2, 1))
  expect_equal(data$beta[1, 1, 1], 0.5)
})

test_that("tabular adapter close is idempotent", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines("sample,subject,contrast,beta", tmp)

  handle <- .tabular_open(tmp)

  # Should succeed multiple times
  expect_silent(.tabular_close(handle))
  expect_silent(.tabular_close(handle))
})

test_that("tabular adapter rejects non-tabular files", {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))

  writeLines("This is not a tabular file", tmp)

  score <- .tabular_detect(tmp)
  expect_false(score)
})

test_that("tabular adapter works via gds() interface", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_2,sub-01,cope1,0.8,0.02"
  ), tmp)

  plan <- gds(tmp)
  expect_s3_class(plan, "gds_plan")

  # Check metadata
  meta <- plan$metadata
  expect_equal(unname(meta$dims["sample"]), 2)
  expect_equal(unname(meta$dims["subject"]), 1)
  expect_equal(unname(meta$dims["contrast"]), 1)

  # Should be able to compute
  result <- compute(plan, assays = "beta")
  expect_type(result, "list")
  expect_true("beta" %in% names(result))
})
