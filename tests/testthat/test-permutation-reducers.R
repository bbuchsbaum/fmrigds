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
})
