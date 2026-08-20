.weighted_onesample_reference <- function(x, weights) {
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[keep]
  weights <- weights[keep]
  if (length(x) < 2L) {
    return(c(beta_g = NA_real_, se_g = NA_real_, t_g = NA_real_, df = NA_real_))
  }
  sum_w <- sum(weights)
  sum_w2 <- sum(weights^2)
  beta_g <- sum(weights * x) / sum_w
  variance_denom <- sum_w - sum_w2 / sum_w
  variance <- sum(weights * (x - beta_g)^2) / variance_denom
  effective_n <- sum_w^2 / sum_w2
  se_g <- sqrt(variance / effective_n)
  c(beta_g = beta_g, se_g = se_g, t_g = beta_g / se_g, df = effective_n - 1)
}

.permutation_tail_score <- function(statistic, alternative) {
  if (identical(alternative, "less")) -statistic else if (identical(alternative, "greater")) statistic else abs(statistic)
}

.all_sign_flips <- function(n_subject) {
  out <- as.matrix(expand.grid(rep(list(c(-1L, 1L)), n_subject)))
  storage.mode(out) <- "integer"
  out
}

test_that("perm:onesample aligns with one-sample t statistics", {
  set.seed(101)
  n_subject <- 10L
  n_sample <- 4L
  y <- matrix(stats::rnorm(n_subject * n_sample, mean = 0.25), nrow = n_subject)
  beta <- array(t(y), dim = c(n_sample, n_subject, 1L))
  g <- as_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space = space_sample_labels(paste0("v", seq_len(n_sample))),
    subjects = paste0("s", seq_len(n_subject)),
    contrasts = "c1"
  )

  out <- compute(reduce(
    g,
    method = "perm:onesample",
    options = list(n_perm = 255L, seed = 99L)
  ))

  expected_t <- apply(y, 2L, function(x) unname(stats::t.test(x, mu = 0)$statistic))
  expected_p <- apply(y, 2L, function(x) stats::t.test(x, mu = 0)$p.value)
  expect_equal(as.numeric(assay(out, "t_g")[, 1, 1]), expected_t, tolerance = 1e-10)
  expect_equal(as.numeric(assay(out, "p_g")[, 1, 1]), expected_p, tolerance = 1e-10)
  expect_true(all(assay(out, "p_perm")[, 1, 1] >= 0 & assay(out, "p_perm")[, 1, 1] <= 1))
  expect_true(all(assay(out, "p_fwer")[, 1, 1] >= assay(out, "p_perm")[, 1, 1]))
})

test_that("perm:onesample tail-specific FWER matches an exhaustive R oracle", {
  beta <- matrix(c(
    -1.8, 0.2, 1.0,
    -1.2, 0.4, 0.7,
    -0.9, 0.1, 1.4,
    -1.5, 0.6, 0.8,
    -0.7, 0.3, 1.1
  ), nrow = 5L, byrow = TRUE)
  weights <- matrix(rep(c(1, 2, 4, 1, 3), ncol(beta)), nrow = nrow(beta))
  signs <- .all_sign_flips(nrow(beta))

  for (alternative in c("two.sided", "less", "greater")) {
    tail <- fmrigds:::.perm_tail_code(alternative)
    observed <- vapply(seq_len(ncol(beta)), function(j) {
      .weighted_onesample_reference(beta[, j], weights[, j])[["t_g"]]
    }, numeric(1))
    null_statistics <- vapply(seq_len(nrow(signs)), function(i) {
      vapply(seq_len(ncol(beta)), function(j) {
        .weighted_onesample_reference(signs[i, ] * beta[, j], weights[, j])[["t_g"]]
      }, numeric(1))
    }, numeric(ncol(beta)))
    null_statistics <- t(null_statistics)
    observed_score <- .permutation_tail_score(observed, alternative)
    null_score <- apply(null_statistics, 2L, .permutation_tail_score, alternative = alternative)
    max_null <- apply(null_score, 1L, max)
    expected_perm <- vapply(seq_along(observed), function(j) {
      (sum(null_score[, j] >= observed_score[j]) + 1) / (nrow(signs) + 1)
    }, numeric(1))
    expected_fwer <- vapply(seq_along(observed), function(j) {
      (sum(max_null >= observed_score[j]) + 1) / (nrow(signs) + 1)
    }, numeric(1))

    result <- perm_onesample_t_cpp(beta, signs, weights = weights, tail = tail)
    expect_equal(as.numeric(result$t_g), observed, tolerance = 1e-12)
    expect_equal(as.numeric(result$p_perm), expected_perm, tolerance = 1e-12)
    expect_equal(as.numeric(result$p_fwer), expected_fwer, tolerance = 1e-12)
    expect_true(all(result$p_fwer >= result$p_perm - 1e-12))
  }
})

test_that("perm:onesample uses fixed weights throughout inference", {
  beta <- array(c(0, 0, 0, 10), dim = c(1L, 4L, 1L))
  variance <- array(c(1, 1, 1, 100), dim = dim(beta))
  precision <- 1 / variance
  g <- new_gds(
    list(beta = beta, var = variance, n_eff = precision),
    space_sample_labels("v1"),
    paste0("s", seq_len(4L)),
    "c1"
  )
  permutation_options <- list(n_perm = 63L, seed = 7L)

  default <- compute(reduce(g, method = "perm:onesample", options = permutation_options))
  equal <- compute(reduce(g, method = "perm:onesample", weights = "equal", options = permutation_options))
  inverse_variance <- compute(reduce(g, method = "perm:onesample", weights = "1/var", options = permutation_options))
  effective_sample <- compute(reduce(g, method = "perm:onesample", weights = "n_eff", options = permutation_options))
  custom <- compute(reduce(
    g,
    method = "perm:onesample",
    weights = "custom",
    options = c(permutation_options, list(custom_weights = precision))
  ))
  large_scale_custom <- compute(reduce(
    g,
    method = "perm:onesample",
    weights = "custom",
    options = c(permutation_options, list(custom_weights = 1e200 * precision))
  ))
  small_scale_custom <- compute(reduce(
    g,
    method = "perm:onesample",
    weights = "custom",
    options = c(permutation_options, list(custom_weights = 1e-200 * precision))
  ))

  expect_equal(assays(default), assays(equal), tolerance = 1e-12)
  expect_equal(as.numeric(assay(equal, "beta_g")), 2.5, tolerance = 1e-12)
  expect_equal(
    as.numeric(assay(inverse_variance, "beta_g")),
    sum(precision * beta) / sum(precision),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(assays(equal), assays(inverse_variance))))
  expect_equal(assays(inverse_variance), assays(effective_sample), tolerance = 1e-12)
  expect_equal(assays(inverse_variance), assays(custom), tolerance = 1e-12)
  expect_equal(assays(custom), assays(large_scale_custom), tolerance = 1e-12)
  expect_equal(assays(custom), assays(small_scale_custom), tolerance = 1e-12)
  expect_true(all(assay(custom, "p_fwer") >= assay(custom, "p_perm") - 1e-12))
})

test_that("perm:onesample rejects unavailable or invalid requested weights", {
  beta <- array(seq_len(8L), dim = c(2L, 4L, 1L))
  variance <- array(1, dim = dim(beta))
  g <- new_gds(
    list(beta = beta, var = variance),
    space_sample_labels(c("v1", "v2")),
    paste0("s", seq_len(4L)),
    "c1"
  )

  expect_error(
    compute(reduce(g, method = "perm:onesample", weights = "n_eff", options = list(n_perm = 7L))),
    "n_eff"
  )
  expect_error(
    compute(reduce(
      g,
      method = "perm:onesample",
      weights = "custom",
      options = list(n_perm = 7L, custom_weights = array(0, dim = dim(beta)))
    )),
    "positive"
  )

  synthetic <- new_gds(
    list(beta = beta, var = variance),
    space_sample_labels(c("v1", "v2")),
    paste0("s", seq_len(4L)),
    "c1",
    metadata = list(synthetic_var = TRUE)
  )
  expect_error(
    reduce(synthetic, method = "perm:onesample", weights = "1/var"),
    "synthetic"
  )
  expect_s3_class(reduce(synthetic, method = "perm:onesample", weights = "equal"), "gds_plan")
})

test_that("perm:twosample aligns with Welch two-sample t statistics", {
  set.seed(202)
  n0 <- 6L
  n1 <- 7L
  group <- factor(rep(c("control", "case"), c(n0, n1)), levels = c("control", "case"))
  n_subject <- length(group)
  n_sample <- 3L
  y <- matrix(stats::rnorm(n_subject * n_sample), nrow = n_subject)
  y[group == "case", ] <- y[group == "case", ] + 0.6
  beta <- array(t(y), dim = c(n_sample, n_subject, 1L))
  subjects <- paste0("s", seq_len(n_subject))
  g <- as_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space = space_sample_labels(paste0("v", seq_len(n_sample))),
    subjects = subjects,
    contrasts = "c1",
    col_data = data.frame(group = group, row.names = subjects)
  )

  out <- compute(reduce(
    g,
    method = "perm:twosample",
    formula = ~ group,
    options = list(n_perm = 299L, seed = 44L)
  ))

  expected_t <- apply(y, 2L, function(x) {
    unname(stats::t.test(x[group == "case"], x[group == "control"], var.equal = FALSE)$statistic)
  })
  expected_p <- apply(y, 2L, function(x) {
    stats::t.test(x[group == "case"], x[group == "control"], var.equal = FALSE)$p.value
  })
  expect_equal(as.numeric(assay(out, "t_g")[, 1, 1]), expected_t, tolerance = 1e-10)
  expect_equal(as.numeric(assay(out, "p_g")[, 1, 1]), expected_p, tolerance = 1e-10)
  expect_true(all(assay(out, "p_perm")[, 1, 1] >= 0 & assay(out, "p_perm")[, 1, 1] <= 1))
  expect_true(all(assay(out, "p_fwer")[, 1, 1] >= assay(out, "p_perm")[, 1, 1]))
})

test_that("perm:twosample FWER uses the requested tail", {
  set.seed(222)
  n_sample <- 40L
  group <- factor(rep(c("control", "case"), each = 6L), levels = c("control", "case"))
  y <- matrix(stats::rnorm(length(group) * n_sample, sd = 0.1), nrow = length(group))
  y[group == "case", ] <- y[group == "case", ] - 2
  beta <- array(t(y), dim = c(n_sample, length(group), 1L))
  subjects <- paste0("s", seq_along(group))
  g <- as_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space = space_sample_labels(paste0("v", seq_len(n_sample))),
    subjects = subjects,
    contrasts = "c1",
    col_data = data.frame(group = group, row.names = subjects)
  )

  for (alternative in c("two.sided", "less", "greater")) {
    out <- compute(reduce(
      g,
      method = "perm:twosample",
      formula = ~ group,
      options = list(n_perm = 99L, seed = 22L, alternative = alternative)
    ))
    p_perm <- as.numeric(assay(out, "p_perm"))
    p_fwer <- as.numeric(assay(out, "p_fwer"))
    expect_true(all(is.finite(p_perm)))
    expect_true(all(is.finite(p_fwer)))
    expect_true(all(p_fwer >= p_perm - 1e-12), info = alternative)
  }
})

test_that("perm:twosample supports direct group options without formula", {
  set.seed(303)
  group <- rep(c(0L, 1L), each = 5L)
  beta <- array(stats::rnorm(20L), dim = c(2L, 10L, 1L))
  g <- as_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space = space_sample_labels(c("a", "b")),
    subjects = paste0("s", seq_len(10L)),
    contrasts = "c1"
  )

  out <- compute(reduce(
    g,
    method = "perm:twosample",
    options = list(group = group, n_perm = 49L, seed = 2L)
  ))

  expect_true(all(c("t_g", "p_g", "p_perm", "p_fwer") %in% names(assays(out))))
  expect_true(all(is.finite(assay(out, "p_perm")[, 1, 1])))

  expect_error(
    reduce(g, method = "perm:twosample", weights = "1/var", options = list(group = group)),
    "only weights"
  )
})
