test_that("as_gds.data.frame builds a valid GDS", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  g <- as_gds(df)
  expect_s3_class(g, "gds")
  expect_equal(dim(assay(g, "beta")), c(2, 2, 1))
  expect_equal(subjects(g), c("s1", "s2"))
  expect_equal(contrasts(g), "c1")
  expect_true(all(c("beta", "var") %in% names(assays(g))))
})

test_that("as_gds.list uses dimnames when available", {
  beta <- array(runif(2 * 3 * 1), dim = c(2, 3, 1),
                dimnames = list(c("A", "B"), paste0("sub", 1:3), "c1"))
  var <- array(runif(2 * 3 * 1), dim = c(2, 3, 1),
               dimnames = dimnames(beta))
  g <- as_gds(list(beta = beta, var = var))
  expect_equal(subjects(g), paste0("sub", 1:3))
  expect_equal(contrasts(g), "c1")
})

test_that("as_gds.array wraps a single array into an assay", {
  x <- array(runif(2 * 2 * 2), dim = c(2, 2, 2))
  g <- as_gds(x, assay_name = "z")
  expect_true("z" %in% names(assays(g)))
  expect_equal(dim(assay(g, "z")), c(2, 2, 2))
})
