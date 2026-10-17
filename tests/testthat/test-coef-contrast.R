# coef_contrast(): per-feature linear contrasts of reducer coefficients,
# checked against lm()/vcov()/anova() and the meta-regression z statistics.

.contrast_fixture <- function(na_voxel = TRUE) {
  set.seed(2026)
  N <- 14; B <- 7
  subj <- sprintf("s%02d", seq_len(N))
  cd <- data.frame(
    age = rnorm(N, 40, 10),
    grp = factor(rep(c("A", "B"), length.out = N)),
    dose = runif(N),
    row.names = subj
  )
  X <- model.matrix(~ age + grp + dose, cd)
  beta <- array(NA_real_, c(B, N, 2),
                dimnames = list(paste0("v", seq_len(B)), subj, c("c1", "c2")))
  for (k in 1:2) for (b in seq_len(B)) {
    beta[b, , k] <- drop(X %*% rnorm(4)) + rnorm(N)
  }
  if (na_voxel) beta[3, c(2, 9), 1] <- NA
  g <- as_gds(list(beta = beta, var = array(1, dim(beta), dimnames(beta))))
  list(g = with_col_data(g, cd), cd = cd, beta = beta, X = X)
}

test_that("t contrasts match lm() estimate, SE, t, df and p (incl. listwise NA)", {
  fx <- .contrast_fixture()
  fit <- suppressWarnings(compute(reduce(
    fx$g, method = "ols:voxelwise", formula = ~ age + grp + dose,
    options = list(return_cov = "tri")
  )))
  C <- rbind(age_minus_dose = c(age = 1, grpB = 0, dose = -1),
             grp_plus = c(age = 0.5, grpB = 2, dose = 0))
  res <- coef_contrast(fit, C)

  for (k in 1:2) for (b in seq_len(7)) {
    y <- fx$beta[b, , k]
    m <- lm(y ~ age + grp + dose, data = fx$cd)
    V <- vcov(m)
    for (nm in rownames(C)) {
      cv <- setNames(numeric(4), names(coef(m)))
      cv[colnames(C)] <- C[nm, ]
      est <- sum(cv * coef(m))
      se <- sqrt(drop(t(cv) %*% V %*% cv))
      tv <- est / se
      expect_equal(assay(res, paste0("contrast:", nm))[b, 1, k], est, tolerance = 1e-10)
      expect_equal(assay(res, paste0("se_contrast:", nm))[b, 1, k], se, tolerance = 1e-10)
      expect_equal(assay(res, paste0("t_contrast:", nm))[b, 1, k], tv, tolerance = 1e-10)
      expect_equal(assay(res, paste0("df_contrast:", nm))[b, 1, k], m$df.residual)
      expect_equal(assay(res, paste0("p_contrast:", nm))[b, 1, k],
                   2 * pt(-abs(tv), m$df.residual), tolerance = 1e-10)
    }
  }
  # original assays are retained
  expect_true(all(names(assays(fit)) %in% names(assays(res))))
})

test_that("joint F test matches anova() of nested lm fits", {
  fx <- .contrast_fixture()
  fit <- suppressWarnings(compute(reduce(
    fx$g, method = "ols:voxelwise", formula = ~ age + grp + dose,
    options = list(return_cov = "tri")
  )))
  C <- rbind(c(age = 1, grpB = 0, dose = 0), c(age = 0, grpB = 1, dose = 0),
             c(age = 0, grpB = 0, dose = 1))
  res <- coef_contrast(fit, C, joint = TRUE, name = "omnibus")
  for (k in 1:2) for (b in seq_len(7)) {
    y <- fx$beta[b, , k]
    ok <- is.finite(y)
    full <- lm(y ~ age + grp + dose, data = fx$cd, subset = ok)
    red <- lm(y ~ 1, data = fx$cd, subset = ok)
    a <- anova(red, full)
    expect_equal(assay(res, "F_contrast:omnibus")[b, 1, k], a$F[2], tolerance = 1e-9)
    expect_equal(assay(res, "p_contrast:omnibus")[b, 1, k], a$`Pr(>F)`[2], tolerance = 1e-9)
    expect_equal(assay(res, "df1_contrast:omnibus")[b, 1, k], 3)
    expect_equal(assay(res, "df2_contrast:omnibus")[b, 1, k], full$df.residual)
  }
})

test_that("packed attachment is used when cov:* assays are absent", {
  fx <- .contrast_fixture(na_voxel = FALSE)
  fit <- compute(reduce(fx$g, method = "ols:voxelwise", formula = ~ age + grp + dose,
                        options = list(return_cov = "tri")))
  ref <- coef_contrast(fit, c(age = 1, dose = -1), name = "d")
  stripped <- fit
  stripped$assays <- stripped$assays[!grepl("^cov:", names(stripped$assays))]
  alt <- coef_contrast(stripped, c(age = 1, dose = -1), name = "d")
  expect_equal(assay(alt, "se_contrast:d"), assay(ref, "se_contrast:d"), tolerance = 1e-12)
})

test_that("single-term contrasts reproduce reducer statistics; z for meta-regression", {
  fx <- .contrast_fixture(na_voxel = FALSE)
  fit <- compute(reduce(fx$g, method = "ols:voxelwise", formula = ~ age + grp + dose))
  res <- coef_contrast(fit, c(grpB = 1), name = "g")
  expect_equal(assay(res, "t_contrast:g"), assay(fit, "t_coef:grpB"), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(assay(res, "p_contrast:g"), assay(fit, "p_coef:grpB"), tolerance = 1e-12,
               ignore_attr = TRUE)
  one <- coef_contrast(fit, c(grpB = 1), name = "g", alternative = "greater")
  expect_equal(as.numeric(assay(one, "p_contrast:g")),
               pt(as.numeric(assay(fit, "t_coef:grpB")), as.numeric(assay(fit, "df_res")),
                  lower.tail = FALSE), tolerance = 1e-12)

  # meta:fe_reg reports normal-approximation z; stat = "auto" follows it
  var <- array(runif(length(fx$beta), 0.5, 2), dim(fx$beta), dimnames(fx$beta))
  gm <- with_col_data(as_gds(list(beta = fx$beta, var = var)), fx$cd)
  mfit <- compute(reduce(gm, method = "meta:fe_reg", formula = ~ age + grp + dose))
  mres <- coef_contrast(mfit, c(age = 2), name = "a2")
  expect_false("t_contrast:a2" %in% names(assays(mres)))
  expect_equal(assay(mres, "z_contrast:a2"), assay(mfit, "z_coef:age"), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(assay(mres, "contrast:a2"), 2 * assay(mfit, "coef:age"), ignore_attr = TRUE)
  expect_error(coef_contrast(mfit, c(age = 1, dose = 1)), "No covariance")
})

test_that("contrast p-values feed fdr:bh via options$source", {
  fx <- .contrast_fixture(na_voxel = FALSE)
  fit <- compute(reduce(fx$g, method = "ols:voxelwise", formula = ~ age + grp + dose,
                        options = list(return_cov = "tri")))
  res <- coef_contrast(fit, c(age = 1, dose = -1), name = "d")
  q <- compute(posthoc(res, "fdr:bh", options = list(source = "p_contrast:d")))
  p <- assay(res, "p_contrast:d")
  expect_equal(assay(q, "q_contrast:d")[, 1, 1], p.adjust(p[, 1, 1], "BH"), ignore_attr = TRUE)
})

test_that("coef_contrast validates its inputs", {
  fx <- .contrast_fixture(na_voxel = FALSE)
  fit <- compute(reduce(fx$g, method = "ols:voxelwise", formula = ~ age + grp + dose))
  expect_error(coef_contrast(fit, c(agee = 1)), "Unknown term\\(s\\) in `C`: agee")
  expect_error(coef_contrast(fit, c(1, -1)), "named by coefficient term")
  expect_error(coef_contrast(fit, c(age = 0)), "non-zero weight")
  expect_error(coef_contrast(fit, c(age = 1, dose = -1)), "return_cov")
  expect_error(coef_contrast(reduce(fx$g, method = "ols:voxelwise", formula = ~ age), c(age = 1)),
               "realised gds")
})

test_that("contrasts of interaction terms read the right covariance (keys with ':')", {
  set.seed(7)
  N <- 30; B <- 5
  subj <- sprintf("s%02d", seq_len(N))
  cd <- data.frame(x = rnorm(N), y = rnorm(N), z = rnorm(N), row.names = subj)
  X <- model.matrix(~ x + x:y + y:z + z, cd)
  beta <- array(NA_real_, c(B, N, 1), dimnames = list(paste0("v", seq_len(B)), subj, "c1"))
  for (b in seq_len(B)) beta[b, , 1] <- drop(X %*% rnorm(ncol(X))) + rnorm(N)
  g <- with_col_data(as_gds(list(beta = beta, var = array(1, dim(beta), dimnames(beta)))), cd)
  fit <- compute(reduce(g, method = "ols:voxelwise", formula = ~ x + x:y + y:z + z,
                        options = list(return_cov = "tri")))
  # Cov(x:y, z) and Cov(x, y:z) share the key "cov:x:y:z"
  out <- coef_contrast(fit, c(`x:y` = 1, z = 1), name = "xyz")
  for (b in seq_len(B)) {
    m <- lm(beta[b, , 1] ~ x + x:y + y:z + z, data = cd)
    L <- setNames(numeric(length(coef(m))), names(coef(m)))
    L[c("x:y", "z")] <- 1
    expect_equal(assay(out, "se_contrast:xyz")[b, 1, 1],
                 sqrt(drop(t(L) %*% vcov(m) %*% L)), tolerance = 1e-10)
  }
  # Without the attachment the ambiguous assay key is refused, not guessed
  stripped <- fit
  stripped$metadata$attachments <- list()
  expect_error(coef_contrast(stripped, c(`x:y` = 1, z = 1)), "ambiguous")
})
