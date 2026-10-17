# ttest:welch must reproduce stats::t.test(var.equal = FALSE) per feature.

.welch_fixture <- function() {
  set.seed(314)
  nA <- 7; nB <- 11; N <- nA + nB; B <- 9
  subj <- sprintf("s%02d", seq_len(N))
  grp <- factor(sample(rep(c("ctl", "pat"), c(nA, nB))), levels = c("ctl", "pat"))
  beta <- array(NA_real_, c(B, N, 2),
                dimnames = list(paste0("v", seq_len(B)), subj, c("c1", "c2")))
  for (k in 1:2) for (b in seq_len(B)) {
    beta[b, , k] <- rnorm(N, mean = ifelse(grp == "pat", 0.4 * b, 0),
                          sd = ifelse(grp == "pat", 2, 0.5))
  }
  beta[2, which(grp == "pat")[1:2], 1] <- NA      # listwise per voxel
  beta[4, which(grp == "ctl")[-1], 2] <- NA       # only one control left -> NA
  g <- as_gds(list(beta = beta, var = array(1, dim(beta), dimnames(beta))))
  list(g = with_col_data(g, data.frame(grp = grp, row.names = subj)),
       grp = grp, beta = beta)
}

test_that("ttest:welch matches t.test(var.equal = FALSE) for every feature", {
  fx <- .welch_fixture()
  for (alt in c("two.sided", "greater", "less")) {
    fit <- compute(reduce(fx$g, method = "ttest:welch", formula = ~ grp,
                          options = list(alternative = alt)))
    for (k in 1:2) for (b in seq_len(9)) {
      y <- fx$beta[b, , k]
      x1 <- y[fx$grp == "pat"]; x0 <- y[fx$grp == "ctl"]
      x1 <- x1[is.finite(x1)]; x0 <- x0[is.finite(x0)]
      if (length(x0) < 2 || length(x1) < 2) {
        expect_true(is.na(assay(fit, "t_g")[b, 1, k]))
        expect_true(is.na(assay(fit, "p_g")[b, 1, k]))
        next
      }
      tt <- t.test(x1, x0, var.equal = FALSE, alternative = alt)
      expect_equal(assay(fit, "beta_g")[b, 1, k], unname(diff(rev(tt$estimate))), tolerance = 1e-12)
      expect_equal(assay(fit, "se_g")[b, 1, k], tt$stderr, tolerance = 1e-12)
      expect_equal(assay(fit, "t_g")[b, 1, k], unname(tt$statistic), tolerance = 1e-12)
      expect_equal(assay(fit, "df")[b, 1, k], unname(tt$parameter), tolerance = 1e-12)
      expect_equal(assay(fit, "p_g")[b, 1, k], tt$p.value, tolerance = 1e-12)
      expect_equal(assay(fit, "p")[b, 1, k], tt$p.value, tolerance = 1e-12)
      # normal-equivalent z carries the same p-value
      expect_equal(fmrigds:::.p_from_z(assay(fit, "z")[b, 1, k], alt), tt$p.value,
                   tolerance = 1e-10)
      expect_equal(assay(fit, "n1")[b, 1, k], length(x1))
    }
  }
})

test_that("ttest:welch accepts options$group and agrees with the formula path", {
  fx <- .welch_fixture()
  a <- compute(reduce(fx$g, method = "ttest:welch", formula = ~ grp))
  b <- compute(reduce(fx$g, method = "ttest:welch", options = list(group = fx$grp)))
  expect_equal(assay(a, "t_g"), assay(b, "t_g"))
  expect_equal(assay(a, "df"), assay(b, "df"))
  expect_error(compute(reduce(fx$g, method = "ttest:welch")), "requires X or options\\$group")
  expect_error(compute(reduce(fx$g, method = "ttest:welch",
                              options = list(group = rep(1:3, 6)))), "exactly two groups")
})

test_that("ttest:welch matches the observed statistic of perm:twosample (welch)", {
  fx <- .welch_fixture()
  a <- compute(reduce(fx$g, method = "ttest:welch", formula = ~ grp))
  p <- compute(reduce(fx$g, method = "perm:twosample", formula = ~ grp,
                      options = list(n_perm = 20, seed = 1)))
  expect_equal(assay(a, "t_g"), assay(p, "t_g"), tolerance = 1e-10)
  expect_equal(assay(a, "df"), assay(p, "df"), tolerance = 1e-10)
})
