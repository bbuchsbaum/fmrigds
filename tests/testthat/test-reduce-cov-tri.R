# Regression: a reducer that lists "cov_tri" in `provides` must produce only
# the per-pair `cov:<a>:<b>` assays and the packed attachment, never a bare
# `cov_tri` base assay (previously the packed L x samples matrix was
# flattened and truncated to n_samples, or left all-NA when not requested).

.cov_tri_fixture <- function() {
  set.seed(11)
  N <- 9; B <- 6
  subj <- paste0("s", seq_len(N))
  x <- rnorm(N)
  beta <- array(rnorm(B * N), c(B, N, 1),
                dimnames = list(paste0("v", seq_len(B)), subj, "c1"))
  g <- as_gds(list(beta = beta, var = array(1, dim(beta), dimnames(beta))))
  list(g = with_col_data(g, data.frame(x = x, row.names = subj)), x = x, beta = beta)
}

test_that("ols:voxelwise with return_cov = 'tri' emits no bare cov_tri assay", {
  fx <- .cov_tri_fixture()
  out <- compute(reduce(fx$g, method = "ols:voxelwise", formula = ~ x,
                        options = list(return_cov = "tri")))
  nms <- names(assays(out))
  expect_false("cov_tri" %in% nms)
  expect_setequal(grep("^cov", nms, value = TRUE),
                  c("cov:(Intercept):(Intercept)", "cov:(Intercept):x", "cov:x:x"))

  # Packed attachment still available and consistent with the per-pair assays
  att <- coef_cov_tri(out, contrast = "c1")
  expect_equal(dim(att$cov_tri), c(3L, 6L))
  expect_equal(as.numeric(assay(out, "cov:(Intercept):x")[, 1, 1]),
               as.numeric(att$cov_tri[2, ]))

  # ... and equal to vcov(lm()) per sample
  for (b in 1:6) {
    V <- vcov(lm(fx$beta[b, , 1] ~ fx$x))
    expect_equal(as.numeric(att$cov_tri[, b]), c(V[1, 1], V[1, 2], V[2, 2]),
                 tolerance = 1e-10)
  }
})

test_that("ols:voxelwise without covariance emits no cov_tri / cov:* assays", {
  fx <- .cov_tri_fixture()
  out <- compute(reduce(fx$g, method = "ols:voxelwise", formula = ~ x))
  expect_false(any(grepl("^cov", names(assays(out)))))
})
