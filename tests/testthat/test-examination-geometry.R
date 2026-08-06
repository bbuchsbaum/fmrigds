test_that("implicit projection depends on global feature IDs, not block partition", {
  ids <- paste0("contrast|feature-", seq_len(31))
  full <- fmrigds:::.geometry_projection(ids, dimension = 12L, seed = "digest")
  blocked <- rbind(
    fmrigds:::.geometry_projection(ids[1:7], 12L, "digest"),
    fmrigds:::.geometry_projection(ids[8:19], 12L, "digest"),
    fmrigds:::.geometry_projection(ids[20:31], 12L, "digest")
  )
  expect_identical(full, blocked)
  expect_true(all(rowSums(full != 0) == 3L))
  expect_gt(sum(colSums(abs(full)) > 0), 8L)
})

test_that("full subject-rank geometry agrees with explicit residual SVD", {
  g <- .group_examination_fixture()
  control <- examination_control(
    block_size = 6L,
    geometry = list(rank = 10L, oversample = 0L, cap = 8)
  )
  exam <- examine_group(reduce(as_plan(g), method = "meta:fe"), control = control)

  beta <- t(assay(g, "beta")[, , 1])
  var <- t(assay(g, "var")[, , 1])
  fit <- fmrigds:::core_meta_fe_kernel(beta, var)
  diagnostic <- fmrigds:::.diagnose_meta_fe_block(
    fit, beta, var, NULL,
    matrix(1, 1, 1, dimnames = list("pooled_effect", "pooled_effect")),
    list(), control$tolerance
  )
  E <- diagnostic$predictive_resid
  E[!diagnostic$surprise_eligible] <- 0
  E <- pmin(E, control$geometry$cap)
  E <- pmax(E, -control$geometry$cap)
  E <- E / sqrt(ncol(E))
  singular <- svd(E, nu = 0, nv = 0)$d
  expected_energy <- singular^2 / sum(E^2)
  expected_energy <- expected_energy[expected_energy > control$tolerance$degeneracy]

  expect_equal(exam$embedding$captured_energy, 1, tolerance = 1e-10)
  expect_equal(
    exam$embedding$explained_energy,
    head(expected_energy, length(exam$embedding$explained_energy)),
    tolerance = 1e-8
  )
  coordinates <- as.matrix(
    exam$embedding$coordinates[grep("^dimension", names(exam$embedding$coordinates))]
  )
  expect_equal(tcrossprod(coordinates), tcrossprod(E), tolerance = 1e-8)
})

test_that("captured residual energy remains a bounded fraction", {
  control <- examination_control(
    geometry = list(rank = 2L, oversample = 0L, stability_replicates = 0L)
  )
  state <- list(
    basis = list(Q = diag(2), sketch_rank = 2L, requested_rank = 2L),
    C = diag(c(0.6, 0.4 + 4 * .Machine$double.eps)),
    total_energy = 1,
    split_C = array(0, c(2L, 2L, 0L)),
    subjects = c("s1", "s2")
  )

  embedding <- fmrigds:::.finalize_residual_geometry(state, control)

  expect_equal(embedding$captured_energy, 1)
  expect_equal(sum(embedding$explained_energy), embedding$captured_energy)
})

test_that("geometry coordinates and fidelity are block invariant", {
  plan <- reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  a <- examine_group(
    plan,
    control = examination_control(
      block_size = 4L,
      geometry = list(rank = 6L, oversample = 3L)
    )
  )
  b <- examine_group(
    plan,
    control = examination_control(
      block_size = 19L,
      geometry = list(rank = 6L, oversample = 3L)
    )
  )
  expect_equal(a$embedding$coordinates, b$embedding$coordinates, tolerance = 1e-8)
  expect_equal(a$embedding$captured_energy, b$embedding$captured_energy, tolerance = 1e-10)
  expect_equal(a$embedding$explained_energy, b$embedding$explained_energy, tolerance = 1e-10)
})

test_that("contrast balancing prevents duplicated contrast features from dominating", {
  control <- examination_control(
    geometry = list(rank = 4L, oversample = 0L, cap = 2, balance_contrasts = TRUE)
  )
  diagnostic <- list(
    predictive_resid = matrix(c(10, -10, 1, -1), 2, 2),
    surprise_eligible = matrix(TRUE, 2, 2)
  )
  state_short <- list(n_sample = 2L)
  short <- fmrigds:::.geometry_residual_block(
    diagnostic, c("a", "b"), "short", state_short, control
  )
  diagnostic_long <- list(
    predictive_resid = diagnostic$predictive_resid[, rep(1:2, each = 4)],
    surprise_eligible = matrix(TRUE, 2, 8)
  )
  state_long <- list(n_sample = 8L)
  long <- fmrigds:::.geometry_residual_block(
    diagnostic_long, paste0("x", 1:8), "long", state_long, control
  )
  expect_equal(sum(short$E^2), sum(long$E^2), tolerance = 1e-12)
  expect_true(max(abs(short$E)) <= 2 / sqrt(2))
})

test_that("stable site structure disappears after model adjustment without clustering", {
  n <- 12L
  p <- 80L
  ids <- paste0("s", seq_len(n))
  site <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
  signal <- sin(seq(0, 2 * pi, length.out = p))
  beta <- array(NA_real_, c(p, n, 1L))
  for (i in seq_len(n)) {
    beta[, i, 1] <- signal + ifelse(site[i] == "B", 2, 0) +
      0.1 * sin(seq_len(p) * (0.11 + i / 50))
  }
  g <- new_gds(
    list(beta = beta, var = array(0.04, dim(beta))),
    space_sample_labels(paste0("v", seq_len(p))),
    ids,
    "task",
    col_data = data.frame(site = site, row.names = ids)
  )
  control <- examination_control(
    geometry = list(rank = 4L, oversample = 4L)
  )
  unadjusted <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = control
  )
  adjusted <- examine_group(
    reduce(as_plan(g), method = "meta:fe_reg", formula = ~ site),
    estimands = "siteB",
    control = control
  )
  site_code <- as.numeric(site == "B")
  expect_gt(abs(cor(unadjusted$embedding$coordinates$dimension1, site_code)), 0.99)
  adjusted_dimensions <- grep(
    "^dimension", names(adjusted$embedding$coordinates), value = TRUE
  )
  adjusted_correlations <- vapply(adjusted_dimensions, function(name) {
    abs(cor(adjusted$embedding$coordinates[[name]], site_code))
  }, numeric(1))
  expect_lt(max(adjusted_correlations), 0.1)
  expect_null(unadjusted$pairwise)
  expect_false("cluster" %in% names(unadjusted$embedding$coordinates))
})

test_that("single-feature apparent isolation is reported as unstable geometry", {
  n <- 10L
  p <- 80L
  beta <- array(0, c(p, n, 1L))
  for (i in seq_len(n)) beta[, i, 1] <- 0.02 * sin(seq_len(p) + i)
  beta[1, 10, 1] <- 10
  g <- new_gds(
    list(beta = beta, var = array(0.04, dim(beta))),
    space_sample_labels(paste0("v", seq_len(p))),
    paste0("s", seq_len(n)),
    "task"
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = examination_control(
      geometry = list(rank = 4L, oversample = 4L, stability_replicates = 2L)
    )
  )
  expect_lt(max(exam$embedding$coordinates$stability), 0.7)
  expect_identical(exam$embedding$stability_method, "deterministic_split_features")
})
