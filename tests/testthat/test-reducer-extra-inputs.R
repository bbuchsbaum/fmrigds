# Contrastwise reducers can declare per-feature assays beyond beta/var/z/p
# (`extra_inputs`); the engine slices them per contrast like beta and passes
# them in opts$inputs. Demonstrated with a voxelwise-covariate OLS reducer
# checked against lm() per voxel.

.register_voxel_covariate_reducer <- function() {
  register_reducer(
    name = "test:ols_voxel_covariate",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      cov <- opts$inputs$gm                       # [subjects x samples], like beta
      stopifnot(identical(dim(cov), dim(beta)))
      B <- ncol(beta)
      coef <- matrix(NA_real_, 2, B); se <- matrix(NA_real_, 2, B)
      for (b in seq_len(B)) {
        Xb <- cbind(1, cov[, b])
        A <- solve(crossprod(Xb))
        bh <- A %*% crossprod(Xb, beta[, b])
        r <- beta[, b] - Xb %*% bh
        s2 <- sum(r^2) / (nrow(Xb) - 2)
        coef[, b] <- bh
        se[, b] <- sqrt(diag(A) * s2)
      }
      rownames(coef) <- rownames(se) <- c("(Intercept)", "gm")
      list(coef = coef, se_coef = se)
    },
    requires = "beta",
    provides = c("coef", "se_coef"),
    extra_inputs = "gm"
  )
}

.unregister_reducer <- function(name) {
  rm(list = name, envir = fmrigds:::.gds_reducers)
}

test_that("extra_inputs assays are sliced per contrast/block and reach the kernel", {
  .register_voxel_covariate_reducer()
  on.exit(.unregister_reducer("test:ols_voxel_covariate"), add = TRUE)
  set.seed(8)
  N <- 12; B <- 10
  subj <- paste0("s", 1:N)
  dn <- list(NULL, subj, c("c1", "c2"))
  gm <- array(runif(B * N * 2), c(B, N, 2), dimnames = dn)
  beta <- array(rnorm(B * N * 2), c(B, N, 2), dimnames = dn) + 3 * gm
  g <- as_gds(list(beta = beta, var = array(1, dim(beta), dn), gm = gm))

  fit <- compute(reduce(g, method = "test:ols_voxel_covariate"))
  blk <- compute(reduce(g, method = "test:ols_voxel_covariate"), block = list(sample = 4:6))
  expect_equal(dim(assay(blk, "coef:gm")), c(3L, 1L, 2L))
  sub <- compute(reduce(subset(g, subject = subj[-c(2, 5)]), method = "test:ols_voxel_covariate"))
  expect_false("gm" %in% names(assays(fit)))
  for (k in 1:2) for (b in seq_len(B)) {
    m <- coef(summary(lm(beta[b, , k] ~ gm[b, , k])))
    expect_equal(assay(fit, "coef:gm")[b, 1, k], m[2, 1], tolerance = 1e-10)
    expect_equal(assay(fit, "se_coef:gm")[b, 1, k], m[2, 2], tolerance = 1e-10)
    if (b %in% 4:6) {                 # compute(block =) evaluates samples 4:6 only
      expect_equal(assay(blk, "coef:gm")[b - 3, 1, k], m[2, 1], tolerance = 1e-10)
    }
    ms <- coef(lm(beta[b, -c(2, 5), k] ~ gm[b, -c(2, 5), k]))
    expect_equal(assay(sub, "coef:gm")[b, 1, k], unname(ms[2]), tolerance = 1e-10)
  }

  g2 <- as_gds(list(beta = beta, var = array(1, dim(beta), dn)))
  expect_error(compute(reduce(g2, method = "test:ols_voxel_covariate")),
               "requires per-feature assay\\(s\\) 'gm'")
  expect_error(compute(reduce(g, method = "test:ols_voxel_covariate",
                              options = list(inputs = list()))), "reserved")
})

test_that("register_reducer validates extra_inputs", {
  f <- function(...) list()
  expect_error(register_reducer("test:bad", f, "beta", "x", extra_inputs = c("a", "a")),
               "unique")
  expect_error(register_reducer("test:bad", f, "beta", "x", input_shape = "joint_contrast",
                                extra_inputs = "a"), "contrastwise")
  expect_null(get_reducer("test:bad"))
})
