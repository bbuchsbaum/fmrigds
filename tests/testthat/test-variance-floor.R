# Pins the explicit regularization policy of the inverse-variance reducers:
# var <= 0 / non-finite excludes the subject; an opt-in eps floors positive var.

.floor_gds <- function(beta, var) {
  dn <- list(paste0("v", seq_len(nrow(beta))), paste0("s", seq_len(ncol(beta))), "c1")
  as_gds(list(beta = array(beta, c(dim(beta), 1), dn), var = array(var, c(dim(var), 1), dn)))
}

test_that("meta:fe honors an explicit floor and excludes non-positive variances", {
  beta <- rbind(c(1, 2, 3, 4), c(1, 2, 3, 4))
  var <- rbind(c(1e-20, 1e-12, 0.5, 1), c(0, -1, 0.5, 1))
  expect_warning(
    expect_warning(
      fit <- compute(reduce(.floor_gds(beta, var), method = "meta:fe",
                            options = list(eps = 1e-12))),
      "Explicit eps variance floor applied"
    ),
    "fewer than 4 finite subjects"
  )
  w1 <- 1 / pmax(var[1, ], 1e-12)
  expect_equal(assay(fit, "beta_g")[1, 1, 1], sum(w1 * beta[1, ]) / sum(w1), tolerance = 1e-12)
  expect_equal(assay(fit, "var_g")[1, 1, 1], 1 / sum(w1), tolerance = 1e-12)
  w2 <- 1 / var[2, 3:4]
  expect_equal(assay(fit, "beta_g")[2, 1, 1], sum(w2 * beta[2, 3:4]) / sum(w2), tolerance = 1e-12)
  expect_equal(assay(fit, "n_eff")[, 1, 1], c(4, 2))

  # options$eps overrides the floor
  fit2 <- suppressWarnings(
    compute(reduce(.floor_gds(beta, var), method = "meta:fe", options = list(eps = 0.25)))
  )
  w3 <- 1 / pmax(var[1, ], 0.25)
  expect_equal(assay(fit2, "beta_g")[1, 1, 1], sum(w3 * beta[1, ]) / sum(w3), tolerance = 1e-12)
})

test_that("meta:re and meta:fe_reg apply the same floor", {
  set.seed(3)
  beta <- matrix(rnorm(12), 2, 6)
  var <- matrix(runif(12, 0.2, 1), 2, 6)
  var[1, 2] <- 1e-30
  g <- .floor_gds(beta, var)
  ref <- .floor_gds(beta, pmax(var, 1e-12))
  expect_warning(
    re <- compute(reduce(g, method = "meta:re", options = list(eps = 1e-12))),
    "Explicit eps variance floor applied"
  )
  re_ref <- compute(reduce(ref, method = "meta:re"))
  expect_identical(assay(re, "beta_g"), assay(re_ref, "beta_g"))
  expect_identical(assay(re, "tau2"), assay(re_ref, "tau2"))
  x <- rnorm(6)
  cd <- data.frame(x = x, row.names = paste0("s", 1:6))
  expect_warning(
    fr <- compute(reduce(with_col_data(g, cd), method = "meta:fe_reg",
                         formula = ~ x, options = list(eps = 1e-12))),
    "Explicit eps variance floor applied"
  )
  fr_ref <- compute(reduce(with_col_data(ref, cd), method = "meta:fe_reg", formula = ~ x))
  expect_identical(assay(fr, "coef:x"), assay(fr_ref, "coef:x"))
})
