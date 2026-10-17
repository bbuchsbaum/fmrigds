# reduce() must forward the weights scheme with exact name matching: an option
# whose name merely starts with "weights" must not stand in for it.

test_that("weights scheme reaches the reducer next to a weights_* option", {
  seen <- new.env()
  register_reducer(
    name = "test:weights_probe",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      seen$weights <- opts[["weights"]]
      seen$custom <- opts[["weights_custom"]]
      list(beta_g = colMeans(beta))
    },
    requires = "beta",
    provides = "beta_g",
    options_schema = list(weights = c("1/var", "equal", "custom"))
  )
  on.exit(rm(list = "test:weights_probe", envir = fmrigds:::.gds_reducers), add = TRUE)
  beta <- array(rnorm(24), c(4, 6, 1), dimnames = list(NULL, sprintf("s%d", 1:6), "c1"))
  g <- as_gds(list(beta = beta, var = array(1, dim(beta), dimnames(beta))))
  w <- c(1, 2, 1, 2, 1, 2)
  compute(reduce(g, method = "test:weights_probe", weights = "custom",
                 options = list(weights_custom = w)))
  expect_identical(seen$weights, "custom")
  expect_identical(seen$custom, w)
})
