test_that("explain_plan returns a tidy node table", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE)
  plan <- gds(tmp)
  plan2 <- derive(plan, c("z", "p"))
  plan3 <- reduce(plan2, method = "fixed")
  tab <- explain_plan(plan3)
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), length(plan3$nodes))
  expect_true(all(c("index", "op", "summary") %in% names(tab)))
})

test_that("preview executes a small block", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE)
  plan <- gds(tmp)
  arrs <- preview(plan, n = 1, assays = c("beta", "var"))
  expect_true(is.list(arrs))
  expect_equal(dim(arrs$beta)[1], 1)
  small <- preview(plan, n = 1)
  expect_s3_class(small, "gds")
  expect_equal(dim(assay(small, "beta"))[1], 1)
})

