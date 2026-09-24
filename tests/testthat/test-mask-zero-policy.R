.zero_background_fixture <- function() {
  beta <- array(0, c(5, 20, 2))
  counts <- cbind(c(15, 5, 7, 6, 0), c(5, 15, 7, 0, 0))
  for (k in 1:2) {
    for (i in 1:5) {
      if (counts[i, k]) beta[i, seq_len(counts[i, k]), k] <- seq_len(counts[i, k])
    }
  }
  new_gds(
    list(beta = beta, var = array(1, dim(beta))),
    space_sample_labels(paste0("v", 1:5)), paste0("s", 1:20), c("a", "b"),
    row_data = data.frame(annotation = letters[1:5], row.names = paste0("v", 1:5))
  )
}

test_that("zero is a valid observation unless explicitly declared background", {
  beta <- array(c(0, 1, 0, 1), c(1, 4, 1))
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))),
               space_sample_labels("binary"), paste0("s", 1:4), "task")
  for (policy in list(MaskPolicy(), MaskPolicy(zero_is_missing = FALSE))) {
    masked <- compute(mask(g, policy))
    expect_equal(assay(masked, "beta"), assay(g, "beta"))
    fit <- compute(one_sample(masked))
    expect_equal(as.numeric(assay(fit, "coef:(Intercept)")), 0.5)
    expect_equal(as.numeric(assay(fit, "n_obs")), 4)
    expect_equal(as.numeric(assay(fit, "p_coef:(Intercept)")), t.test(c(0, 1, 0, 1))$p.value)
  }
})

test_that("background policy excludes observations and checks coverage per contrast", {
  g <- .zero_background_fixture()
  policy <- MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE)
  masked <- compute(mask(g, policy))
  expected <- assay(g, "beta")[1:3, , , drop = FALSE]
  expected[expected == 0] <- NA_real_
  expected[2, , 1] <- NA_real_
  expected[1, , 2] <- NA_real_
  expect_equal(unname(assay(masked, "beta")), unname(expected))
  expect_equal(sample_labels(masked), paste0("v", 1:3))
  expect_equal(row_data(masked)$annotation, letters[1:3])
  expect_equal(dim(assay(masked, "beta")), c(3, 20, 2))
  for (a in assays(masked)) expect_true(all(is.na(a[is.na(expected)])))
  expect_equal(sum(assay(g, "beta") == 0), 140) # original data are unchanged

  fit_warnings <- character()
  fit <- withCallingHandlers(compute(one_sample(mask(g, policy))), warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(fit_warnings, 2L)
  expect_true(all(grepl("per-voxel listwise deletion", fit_warnings)))
  expect_equal(as.numeric(assay(fit, "n_obs")), c(15, 0, 7, 0, 15, 7))
  expect_equal(as.numeric(assay(fit, "df_res")), c(14, NA, 6, NA, 14, 6))
  for (k in 1:2) {
    for (i in 1:3) {
      y <- expected[i, , k]
      if (all(is.na(y))) {
        expect_true(is.na(assay(fit, "p_coef:(Intercept)")[i, 1, k]))
      } else {
        reference <- stats::t.test(y)
        expect_equal(assay(fit, "coef:(Intercept)")[i, 1, k], unname(reference$estimate))
        expect_equal(assay(fit, "t_coef:(Intercept)")[i, 1, k], unname(reference$statistic))
        expect_equal(assay(fit, "p_coef:(Intercept)")[i, 1, k], reference$p.value)
      }
    }
  }
  corrected <- compute(posthoc(fit, "fdr:bh", options = list(source = "p_coef:(Intercept)")))
  for (k in 1:2) {
    p <- as.numeric(assay(fit, "p_coef:(Intercept)")[, 1, k])
    q <- as.numeric(assay(corrected, "q_coef:(Intercept)")[, 1, k])
    expect_equal(q[is.finite(p)], p.adjust(p[is.finite(p)], method = "BH"))
    expect_true(all(is.na(q[is.na(p)])))
  }
})

test_that("only exact zero effects trigger background exclusion in every assay", {
  beta <- array(c(0, -2, 1e-20, 4), c(1, 4, 1))
  arrays <- list(beta = beta, var = array(c(1, 1, 0, 1), dim(beta)),
                 z = array(7, dim(beta)), p = array(0, dim(beta)))
  result <- apply_mask_policy(
    op_mask_policy(MaskPolicy(rule = "union", zero_is_missing = TRUE)),
    arrays, space_sample_labels("v")
  )
  for (a in result$arrays) expect_true(is.na(a[1, 1, 1]))
  expect_equal(as.numeric(result$arrays$beta), c(NA, -2, 1e-20, 4))
  expect_equal(as.numeric(result$arrays$var), c(NA, 1, 0, 1))
  expect_equal(as.numeric(result$arrays$p), c(NA, 0, 0, 0))

  z_only <- apply_mask_policy(
    op_mask_policy(MaskPolicy(rule = "union", zero_is_missing = TRUE)),
    list(z = beta, p = array(0.5, dim(beta))), space_sample_labels("v")
  )
  expect_true(is.na(z_only$arrays$p[1, 1, 1]))
  expect_equal(as.numeric(z_only$arrays$z), c(NA, -2, 1e-20, 4))
})

test_that("coverage counts selected subjects and rejects zero and nonfinite background", {
  g <- .zero_background_fixture()
  selected <- as_plan(g) |> subset(subject = paste0("s", 1:10)) |>
    mask(MaskPolicy(rule = "threshold", threshold = 0.5, zero_is_missing = TRUE)) |>
    compute()
  expect_equal(sample_labels(selected), paste0("v", 1:4))
  expect_equal(sum(is.finite(assay(selected, "beta")[2, , 1])), 5)

  beta <- array(c(0, NA, NaN, Inf, -Inf), c(1, 5, 1))
  expect_error(apply_mask_policy(
    op_mask_policy(MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE)),
    list(beta = beta), space_sample_labels("v")
  ), "Mask removed all samples")
})

test_that("group intersection and union respect independent contrast support", {
  beta <- array(c(0, 1, 2, 3), c(1, 2, 2))
  arrays <- list(beta = beta)
  run <- function(rule) apply_mask_policy(
    op_mask_policy(MaskPolicy(rule = rule, zero_is_missing = TRUE)),
    arrays, space_sample_labels("v")
  )$arrays$beta
  expect_equal(as.numeric(run("intersection")), c(NA, NA, 2, 3))
  expect_equal(as.numeric(run("union")), c(NA, 1, 2, 3))
  # The same per-contrast coverage applies to missing values without the option.
  arrays$beta[1, 1, 1] <- NA_real_
  normal <- apply_mask_policy(op_mask_policy(MaskPolicy()), arrays, space_sample_labels("v"))
  expect_equal(as.numeric(normal$arrays$beta), c(NA, NA, 2, 3))
})

test_that("background masks update packed voxel coordinates", {
  beta <- array(c(0, 2, 0, 4), c(2, 2, 1))
  result <- apply_mask_policy(
    op_mask_policy(MaskPolicy(zero_is_missing = TRUE)), list(beta = beta),
    space_voxel(c(4, 1, 1), diag(4), mask_idx = c(2L, 4L), storage = "packed")
  )
  expect_equal(result$space$mask_idx, 4L)
  expect_equal(result$subset$samples, 2L)
})

test_that("variance-weighted fits use the non-background observations", {
  beta <- array(c(0, 2, 4, 6), c(1, 4, 1))
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))),
               space_sample_labels("v"), paste0("s", 1:4), "task")
  plan <- mask(g, MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE))
  expect_warning(fit <- compute(reduce(plan, method = "meta:fe")), "fewer than 4")
  expect_equal(as.numeric(assay(fit, "beta_g")), 4)
  expect_equal(as.numeric(assay(fit, "n_eff")), 3)
  expect_equal(as.numeric(assay(fit, "se_g")), sqrt(1 / 3))
})

test_that("group OLS retains design requirements after background removal", {
  beta <- array(c(0, 0, 1, 2, 3, 5), c(1, 6, 1))
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))),
               space_sample_labels("v"), paste0("s", 1:6), "task",
               col_data = data.frame(group = factor(c("a", "a", "b", "b", "b", "b")),
                                     row.names = paste0("s", 1:6)))
  expect_warning(fit <- g |>
    mask(MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE)) |>
    group_ols(~ group) |> compute(), "per-voxel listwise deletion")
  expect_equal(as.numeric(assay(fit, "n_obs")), 4)
  expect_true(all(is.na(assay(fit, "p_coef:groupb"))))
})

test_that("background exclusion reaches permutation statistics and null maxima", {
  beta <- array(c(0, 0, 2, 3, 4, 0, 6, 0), c(2, 4, 1))
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))),
               space_sample_labels(c("keep", "drop")), paste0("s", 1:4), "task")
  signs <- as.matrix(expand.grid(rep(list(c(-1L, 1L)), 4)))
  fit <- g |>
    mask(MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE)) |>
    reduce(method = "perm:onesample", options = list(signs = signs)) |> compute()
  y <- c(2, 4, 6)
  expected <- t.test(y)
  null <- apply(signs[, 2:4, drop = FALSE], 1, function(s) unname(t.test(s * y)$statistic))
  p <- (1 + sum(abs(null) >= abs(expected$statistic))) / (1 + nrow(signs))
  expect_equal(sample_labels(fit), "keep")
  expect_equal(as.numeric(assay(fit, "t_g")), unname(expected$statistic))
  expect_equal(as.numeric(assay(fit, "df")), 2)
  expect_equal(as.numeric(assay(fit, "p_perm")), p)
  expect_equal(as.numeric(assay(fit, "p_fwer")), p)
})

test_that("zero conventions survive saved plans and change portable provenance", {
  plan <- gds_plan(gds_source("tabular", list(path = "dummy"), list()))
  ordinary <- mask(plan, MaskPolicy(rule = "threshold", threshold = 1 / 3))
  background <- mask(plan, MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE))
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  save_plan(background, path)
  loaded <- load_plan(path)
  expect_identical(loaded$nodes[[1]]$policy$zero_is_missing, TRUE)
  expect_identical(loaded$nodes[[1]]$policy, background$nodes[[1]]$policy)
  ordinary_receipt <- .portable_plan_receipt(ordinary)
  background_receipt <- .portable_plan_receipt(background)
  expect_false(identical(ordinary_receipt, background_receipt))
  expect_identical(.portable_plan_node(loaded$nodes[[1]])$params$zero_is_missing, TRUE)
  old_node <- list(op = "mask_policy", scope = "group", rule = "threshold", threshold = 0.5)
  expect_identical(.deserialize_node(old_node)$policy$zero_is_missing, FALSE)
})

test_that("group examination sees the same missingness as the group fit", {
  g <- .zero_background_fixture()
  beta <- assay(g, "beta")[1:3, , , drop = FALSE]
  beta[beta == 0] <- NA_real_
  beta[2, , 1] <- NA_real_
  beta[1, , 2] <- NA_real_
  variance <- array(1, dim(beta))
  variance[is.na(beta)] <- NA_real_
  reference <- new_gds(list(beta = beta, var = variance),
                       space_sample_labels(paste0("v", 1:3)), subjects(g), contrasts(g))
  control <- examination_control(block_size = 2L)
  masked_exam <- g |>
    mask(MaskPolicy(rule = "threshold", threshold = 1 / 3, zero_is_missing = TRUE)) |>
    reduce(method = "meta:fe") |> examine_group(control = control)
  reference_exam <- reference |> reduce(method = "meta:fe") |> examine_group(control = control)
  expect_equal(assays(masked_exam$group_maps), assays(reference_exam$group_maps))
  expect_equal(masked_exam$subject_data, reference_exam$subject_data)
  expect_equal(masked_exam$contrast_data, reference_exam$contrast_data)
  expect_equal(masked_exam$estimand_data, reference_exam$estimand_data)
  # Each contrast has two supported locations, despite the shared three-row axis.
  coverage <- masked_exam$contrast_data
  expect_equal(coverage$coverage_fraction[coverage$subject == "s1"], c(1, 1))
  expect_equal(coverage$coverage_fraction[coverage$subject == "s8"], c(0.5, 0.5))
  expect_equal(coverage$coverage_fraction[coverage$subject == "s16"], c(0, 0))
  expect_true(masked_exam$provenance$staging$cleanup_succeeded)
})

test_that("an entirely unsupported contrast has unavailable examination coverage", {
  beta <- array(0, c(2, 4, 2))
  beta[, , 1] <- matrix(1:8, 2, 4)
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))),
               space_sample_labels(c("v1", "v2")), paste0("s", 1:4), c("a", "empty"))
  exam <- g |> mask(MaskPolicy(zero_is_missing = TRUE)) |>
    reduce(method = "meta:fe") |> examine_group()
  coverage <- exam$contrast_data
  expect_true(all(is.na(coverage$coverage_fraction[coverage$contrast == "empty"])))
  expect_equal(coverage$coverage_fraction[coverage$contrast == "a"], rep(1, 4))
  expect_equal(exam$subject_data$coverage_fraction, rep(1, 4))
})

test_that("invalid zero conventions and coverage thresholds fail clearly", {
  for (value in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE))) {
    expect_error(MaskPolicy(zero_is_missing = value), "zero_is_missing")
  }
  for (value in list(NA, Inf, -0.1, 1.1, "0.33", numeric(), c(0.2, 0.3))) {
    expect_error(MaskPolicy(threshold = value), "threshold")
  }
  arrays <- list(beta = array(c(0, 2), c(1, 2, 1)))
  policy <- MaskPolicy(rule = "custom", zero_is_missing = TRUE, custom = function(a) {
    expect_true(is.na(a$beta[1, 1, 1]))
    TRUE
  })
  expect_equal(as.numeric(apply_mask_policy(op_mask_policy(policy), arrays,
                                           space_sample_labels("v"))$arrays$beta), c(NA, 2))
  expect_error(apply_mask_policy(
    op_mask_policy(MaskPolicy(rule = "custom", custom = function(a) NA)),
    arrays, space_sample_labels("v")
  ), "non-missing logical")
})
