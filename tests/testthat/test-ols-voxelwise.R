test_that("ols:voxelwise recovers known coefficients and attaches covariance", {
  set.seed(42)
  # Subjects and design
  N <- 10
  x <- rnorm(N)
  X <- model.matrix(~ 1 + x)
  b_true <- c(0.5, 2.0)
  # Make 3 samples (voxels) with different intercepts/slopes
  B <- 3
  coefs <- rbind(
    c(0.5, 2.0),
    c(1.0, 1.5),
    c(-0.2, 2.5)
  )
  # Build beta array [samples × subjects × contrasts]
  beta <- array(NA_real_, dim = c(B, N, 1))
  for (i in seq_len(B)) {
    y <- as.numeric(X %*% coefs[i, ] + rnorm(N, sd = 0.1))
    beta[i, , 1] <- y
  }
  # Build CSV for tabular ingestion
  samples <- paste0("v", 1:B)
  subjects <- paste0("s", seq_len(N))
  df <- data.frame(
    sample = rep(samples, each = N),
    subject = rep(subjects, times = B),
    contrast = "c1",
    beta = as.numeric(c(beta[1, , 1], beta[2, , 1], beta[3, , 1])),
    var = 1
  )
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(df, tmp, row.names = FALSE)
  cd <- data.frame(x = x, row.names = subjects)
  plan <- gds(tmp, col_data = cd)
  # Reduce via OLS with covariance triangles
  plan2 <- reduce(plan, method = "ols:voxelwise", formula = ~ 1 + x, options = list(return_cov = "tri"))
  out <- compute(plan2)
  # Check param assays exist and roughly match truth
  expect_true(all(c("coef:(Intercept)", "coef:x", "se_coef:(Intercept)", "se_coef:x") %in% names(assays(out))))
  est_int <- assay(out, "coef:(Intercept)")
  est_x <- assay(out, "coef:x")
  expect_equal(as.numeric(est_int[, 1, 1]), coefs[, 1], tolerance = 0.2)
  expect_equal(as.numeric(est_x[, 1, 1]), coefs[, 2], tolerance = 0.2)
  # Covariance attachment
  att <- coef_cov_tri(out, contrast = "c1")
  expect_equal(att$type, "cov_tri")
  expect_equal(att$terms, c("(Intercept)", "x"))
  L <- length(att$terms) * (length(att$terms) + 1) / 2
  expect_equal(nrow(att$cov_tri), L)

  # Covariance assays promoted with clear naming
  expect_true(all(c(
    "cov:(Intercept):(Intercept)",
    "cov:(Intercept):x",
    "cov:x:x"
  ) %in% names(assays(out))))
  cov_ii <- assay(out, "cov:(Intercept):(Intercept)")[, 1, 1]
  cov_ix <- assay(out, "cov:(Intercept):x")[, 1, 1]
  cov_xx <- assay(out, "cov:x:x")[, 1, 1]
  # Compare to packed rows from attachments (row order: (1,1), (1,2), (2,2))
  expect_equal(as.numeric(cov_ii), as.numeric(att$cov_tri[1, ]))
  expect_equal(as.numeric(cov_ix), as.numeric(att$cov_tri[2, ]))
  expect_equal(as.numeric(cov_xx), as.numeric(att$cov_tri[3, ]))
})

test_that("coef_array stacks term assays into [samples × terms × contrasts]", {
  set.seed(1)
  N <- 5; B <- 2
  x <- rnorm(N); X <- model.matrix(~ 1 + x)
  # Write tabular beta for two samples
  samples <- paste0("v", 1:B)
  subjects <- paste0("s", 1:N)
  resp <- sapply(1:B, function(i) as.numeric(X %*% c(i, i)))
  df <- data.frame(
    sample = rep(samples, each = N),
    subject = rep(subjects, times = B),
    contrast = "c1",
    beta = as.numeric(resp),
    var = 1
  )
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(df, tmp, row.names = FALSE)
  cd <- data.frame(x = x, row.names = subjects)
  g2 <- compute(reduce(gds(tmp, col_data = cd), method = "ols:voxelwise", formula = ~ 1 + x))
  A <- coef_array(g2, prefix = "coef:")
  expect_equal(dim(A), c(B, 2, 1))
})
