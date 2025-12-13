test_that("model_matrix builds X from col_data on plan", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,s1,c1,0.5,0.04",
    "ROI_2,s2,c1,0.8,0.09"
  ), tmp)

  plan <- gds(tmp)
  cd <- data.frame(age = c(25, 30), row.names = c("s1", "s2"))
  plan2 <- with_col_data(plan, cd)
  X <- model_matrix(plan2, ~ 1 + age)
  expect_true(is.matrix(X))
  expect_equal(colnames(X), c("(Intercept)", "age"))
  expect_equal(nrow(X), 2)
})

test_that("model_matrix works with realised GDS", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  g <- as_gds(df)
  g$col_data <- data.frame(group = c(0, 1), row.names = c("s1", "s2"))
  X <- model_matrix(g, ~ 1 + group)
  expect_equal(colnames(X), c("(Intercept)", "group"))
  expect_equal(nrow(X), 2)
})

