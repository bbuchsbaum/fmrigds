test_that("apply_map_to maps beta/var with independent uncertainty", {
  beta <- array(c(1, 3), c(2, 1, 1))
  var <- array(c(0.25, 1), c(2, 1, 1))
  df <- array(c(10, 20), c(2, 1, 1))
  arrays <- list(
    beta = beta,
    var = var,
    se = sqrt(var),
    t = beta / sqrt(var),
    df = df
  )

  M <- matrix(c(0.6, 0.8), nrow = 1)
  node <- op_map(
    target_space = space_parcels("target"),
    map = M,
    uncertainty = UncertaintyRule("independent", df_rule = "satterthwaite"),
    combine = NULL
  )

  res <- fmrigds:::apply_map_to(node, arrays)
  arr <- res$arrays

  expect_equal(arr$beta[1, 1, 1], sum(M * beta[, 1, 1]))
  expect_equal(arr$var[1, 1, 1], sum((M^2) * var[, 1, 1]))
  expect_equal(arr$se[1, 1, 1], sqrt(arr$var[1, 1, 1]))
  expect_equal(arr$t[1, 1, 1], arr$beta[1, 1, 1] / arr$se[1, 1, 1])

  df_expected <- fmrigds:::aggregate_df_satterthwaite(M, var, df)
  expect_equal(arr$df[1, 1, 1], df_expected[1, 1, 1])
})

test_that("apply_map_to combines z via Stouffer", {
  z <- array(c(2, -1), c(2, 1, 1))
  arrays <- list(z = z)
  M <- matrix(c(0.5, 0.5), nrow = 1)
  node <- op_map(
    target_space = space_parcels("target"),
    map = M,
    uncertainty = UncertaintyRule("independent"),
    combine = "stouffer"
  )

  res <- fmrigds:::apply_map_to(node, arrays)
  expect_equal(res$arrays$z[1, 1, 1], sum(M * z[, 1, 1]) / sqrt(sum(M^2)))
  expect_equal(res$arrays$p[1, 1, 1], 2 * pnorm(-abs(res$arrays$z[1, 1, 1])))
})

test_that("apply_map_to requires combine when beta missing", {
  arrays <- list(z = array(1, c(2, 1, 1)))
  node <- op_map(
    target_space = space_parcels("target"),
    map = matrix(c(0.5, 0.5), nrow = 1),
    uncertainty = UncertaintyRule("independent"),
    combine = NULL
  )
  expect_error(fmrigds:::apply_map_to(node, arrays), "combiner")
})
