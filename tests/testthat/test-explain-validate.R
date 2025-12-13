test_that("explain() prints for gds and plan", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  g <- as_gds(df)
  out_g <- capture.output(explain(g))
  expect_true(any(grepl("GDS", out_g)))
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE)
  plan <- gds(tmp)
  out_p <- capture.output(explain(plan))
  expect_true(any(grepl("GDS plan", out_p)))
})

test_that("validate() checks basic structure", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  g <- as_gds(df)
  expect_true(validate(g))
  # Tamper with col_data to trigger error
  g$col_data <- data.frame(x = 1, row.names = "oops")
  expect_error(validate(g), "rownames matching subjects")
})

