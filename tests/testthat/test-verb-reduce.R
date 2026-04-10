test_that("reduce() returns a plan object (not error)", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_2,sub-01,cope1,0.8,0.02"
  ), tmp)

  plan <- gds(tmp)
  plan2 <- reduce(plan, method = "fixed")
  expect_s3_class(plan2, "gds_plan")
})

test_that("reduce() aligns formula data to current plan subject order", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.01",
    "ROI_1,sub-02,cope1,0.8,0.02"
  ), tmp)

  plan <- gds(tmp) |> subset(subject = c("sub-02", "sub-01"))
  cd <- data.frame(age = c(10, 20), row.names = c("sub-01", "sub-02"))

  plan2 <- reduce(plan, method = "meta:fe_reg", formula = ~ age, data = cd)
  expect_equal(unname(plan2$nodes[[2]]$options$X[, "age"]), c(20, 10))
})
