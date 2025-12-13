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

