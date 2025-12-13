test_that("with_col_data aligns subjects on plan", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-02,cope1,0.8,0.09"
  ), tmp)

  plan <- gds(tmp)
  subjects <- plan$meta$subjects
  expect_equal(subjects, c("sub-01", "sub-02"))

  cd <- data.frame(age = c(25, 30), row.names = c("sub-01", "sub-02"))
  plan2 <- with_col_data(plan, cd)
  expect_equal(plan2$meta$col_data, cd[subjects, , drop = FALSE])

  # Extra rows trigger warning and are dropped
  cd_extra <- rbind(cd, data.frame(age = 40, row.names = "sub-03"))
  expect_warning(plan3 <- with_col_data(plan, cd_extra), "Dropping extra rows")
  expect_equal(plan3$meta$col_data, cd[subjects, , drop = FALSE])

  # Missing subject errors
  cd_missing <- cd[1, , drop = FALSE]
  expect_error(with_col_data(plan, cd_missing), "missing subjects")
})

test_that("col_data persists through HDF5 write/read", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-02,cope1,0.8,0.09"
  ), tmp)

  cd <- data.frame(group = c("A", "B"), row.names = c("sub-01", "sub-02"), stringsAsFactors = FALSE)
  plan <- gds(tmp, col_data = cd)
  g <- compute(plan)
  expect_equal(col_data(g), cd)

  out <- tempfile(fileext = ".h5")
  on.exit(unlink(out), add = TRUE)
  plan_write <- write_out(plan, out, format = "h5")
  compute(plan_write)

  plan2 <- gds(out)
  expect_equal(plan2$meta$col_data, cd)
  g2 <- compute(plan2)
  expect_equal(col_data(g2), cd)
})

test_that("design matrices recorded in metadata", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-01,cope1,0.8,0.09",
    "ROI_1,sub-02,cope1,1.0,0.04",
    "ROI_2,sub-02,cope1,1.2,0.09"
  ), tmp)

  cd <- data.frame(age = c(25, 30), row.names = c("sub-01", "sub-02"))
  plan <- gds(tmp, col_data = cd)
  plan_reg <- reduce(plan, method = "meta:fe_reg", formula = ~ age)
  g <- compute(plan_reg)
  dm <- metadata(g)$design_mats
  expect_true(length(dm) >= 1)
  expect_equal(dm[[1]]$columns, c("(Intercept)", "age"))
  expect_true(!is.null(dm[[1]]$hash))
})

test_that("exporters can join col_data", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-02,cope1,0.8,0.09"
  ), tmp)

  cd <- data.frame(group = c("A", "B"), row.names = c("sub-01", "sub-02"), stringsAsFactors = FALSE)
  plan <- gds(tmp, col_data = cd)
  g <- compute(plan)

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  plan_export <- write_out(plan, csv, format = "csv", options = list(include_col_data = TRUE))
  compute(plan_export)
  df <- utils::read.csv(csv, stringsAsFactors = FALSE)
  expect_true("group" %in% names(df))
  expect_setequal(unique(df$group), c("A", "B"))
})

test_that("with_col_data works on realised GDS", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-02,cope1,0.8,0.09"
  ), tmp)
  g <- compute(gds(tmp))
  cd <- data.frame(group = c("A","B"), row.names = c("sub-01","sub-02"))
  g2 <- with_col_data(g, cd)
  expect_equal(col_data(g2), cd)
})
