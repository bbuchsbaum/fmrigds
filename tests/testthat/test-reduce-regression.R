test_that("meta-regression expands parameters into param-suffixed assays", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  # two samples, two subjects, single contrast with beta/var
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-01,cope1,0.8,0.09",
    "ROI_1,sub-02,cope1,1.0,0.04",
    "ROI_2,sub-02,cope1,1.2,0.09"
  ), tmp)

  # Subject-level covariate for regression (two subjects)
  cd <- data.frame(x = c(0, 1), row.names = c("sub-01", "sub-02"))

  plan <- gds(tmp, col_data = cd)
  # Regress beta on intercept + x
  plan2 <- reduce(plan, method = "meta:fe_reg", formula = ~ x)
  g <- compute(plan2)

  # Expect param-suffixed assays for coefficients and SEs
  a <- names(assays(g))
  expect_true("coef:(Intercept)" %in% a)
  expect_true("coef:x" %in% a)
  expect_true("se_coef:(Intercept)" %in% a)
  expect_true("se_coef:x" %in% a)

  # Shapes are [sample x 1 x contrast] = [2 x 1 x 1]
  expect_equal(dim(assay(g, "coef:(Intercept)")), c(2, 1, 1))
  expect_equal(dim(assay(g, "coef:x")), c(2, 1, 1))
  expect_equal(dim(assay(g, "se_coef:(Intercept)")), c(2, 1, 1))
  expect_equal(dim(assay(g, "se_coef:x")), c(2, 1, 1))

  # Group-level outputs have subjects = "meta"
  expect_equal(subjects(g), "meta")

  # col_data still exists but only has the subject column with "meta"
  cd <- col_data(g)
  expect_s3_class(cd, "data.frame")
  expect_equal(nrow(cd), 1)
  expect_equal(rownames(cd), "meta")
  expect_true("subject" %in% names(cd))
})
