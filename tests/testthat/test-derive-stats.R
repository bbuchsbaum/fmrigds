test_that("execute_derive computes core statistics", {
  arrays <- list(
    se = array(0.5, c(3, 2, 1)),
    beta = array(c(-1, 0.5, 1.5), c(3, 2, 1)),
    var = array(0.25, c(3, 2, 1)),
    df = array(30, c(3, 2, 1))
  )

  arrays <- gdsfmri:::execute_derive(arrays, c("var", "se", "t", "z", "p"))

  expect_equal(arrays$var[1], 0.25)
  expect_equal(arrays$se[1], 0.5)
  expect_equal(arrays$t[1], -2)

  p_two <- 2 * pt(2, 30, lower.tail = FALSE)
  expect_equal(arrays$p[1], p_two, tolerance = 1e-8)
  expect_equal(arrays$z[1], qnorm(1 - p_two/2) * sign(arrays$t[1]), tolerance = 1e-8)
})

test_that("derive_z from p requires sign", {
  arrays <- list(p = array(0.05, c(1, 1, 1)))
  expect_error(gdsfmri:::derive_z(arrays), "sign information")

  arrays$beta <- array(-1, c(1, 1, 1))
  z <- gdsfmri:::derive_z(arrays)
  expect_lt(z[1], 0)
})

test_that("Satterthwaite aggregation matches formula", {
  M <- matrix(c(0.6, 0.4), nrow = 1)
  var_src <- array(c(0.5, 1.0), c(2, 1, 1))
  df_src <- array(c(20, 40), c(2, 1, 1))

  df_out <- gdsfmri:::aggregate_df_satterthwaite(M, var_src, df_src)

  expected <- (sum((M^2) * c(0.5, 1.0)))^2 /
    sum((M^4) * (c(0.5, 1.0)^2) / c(20, 40))
  expect_equal(df_out[1, 1, 1], expected)
})
