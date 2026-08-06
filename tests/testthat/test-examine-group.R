test_that("group examination separates surprise, precision leverage, and noise", {
  exam <- .group_examination_fixture() |>
    as_plan() |>
    reduce(method = "meta:fe") |>
    examine_group(control = examination_control(block_size = 7L))

  expect_s3_class(exam, "gds_examination")
  expect_s3_class(exam$group_maps, "gds")
  expect_identical(exam$cohort$model, "meta:fe")
  expect_identical(exam$config$weight_contract, "inverse_variance")

  subject <- exam$subject_data
  expect_identical(subject$review_status[subject$subject == "s10"], "review")
  expect_identical(subject$review_source[subject$subject == "s10"], "surprise")
  expect_match(subject$review_reason[subject$subject == "s10"], "negative map gain")
  expect_identical(subject$review_status[subject$subject == "s8"], "review")
  expect_identical(subject$review_source[subject$subject == "s8"], "influence")
  expect_identical(subject$review_status[subject$subject == "s9"], "none")

  contrast <- exam$contrast_data
  expect_lt(
    contrast$zero_intercept_gain[contrast$subject == "s10"],
    0
  )
  influence <- exam$estimand_data
  expect_lt(
    influence$influence_energy[influence$subject == "s9"],
    influence$influence_energy[influence$subject == "s8"]
  )
  expect_true(all(exam$availability$status == "available"))
})

test_that("homogeneous null-style cohort ranks without dramatic review labels", {
  beta <- array(1, c(60, 12, 1))
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_sample_labels(paste0("v", seq_len(60))),
    paste0("s", seq_len(12)),
    "task"
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = examination_control(block_size = 11L)
  )
  expect_true(all(exam$subject_data$review_status == "none"))
  expect_true(all(is.finite(exam$subject_data$review_priority)))
  expect_true(all(exam$contrast_data$correlation_status == "degenerate_expected"))
})

test_that("poor coverage is a validity concern rather than biological surprise", {
  beta <- array(1, c(30, 10, 1))
  beta[16:30, 1, 1] <- NA_real_
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_sample_labels(paste0("v", seq_len(30))),
    paste0("s", seq_len(10)),
    "task"
  )
  exam <- examine_group(reduce(as_plan(g), method = "meta:fe"))
  row <- exam$subject_data[exam$subject_data$subject == "s1", , drop = FALSE]
  expect_equal(row$coverage_fraction, 0.5)
  expect_identical(row$review_status, "review")
  expect_identical(row$review_source, "quality")
  expect_match(row$review_reason, "coverage_fraction")
})

test_that("external quality metrics require an explicit direction and threshold", {
  g <- .group_examination_fixture()
  ids <- subjects(g)
  g <- with_col_data(
    g,
    data.frame(mean_fd = c(rep(0.1, 9), 0.8), row.names = ids)
  )
  control <- examination_control(
    review = list(
      quality = list(
        mean_fd = list(direction = "high", threshold = 0.5)
      )
    )
  )
  expect_error(
    examine_group(reduce(as_plan(g), method = "meta:fe"), control = control),
    "columns in `quality`"
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    quality = "mean_fd",
    control = control
  )
  row <- exam$subject_data[exam$subject_data$subject == "s10", , drop = FALSE]
  expect_identical(row$review_status, "review")
  expect_identical(row$review_source, "quality")
  expect_match(row$review_reason, "mean_fd")
})

test_that("model adjustment explains a supplied site-defined cohort split", {
  n <- 12L
  p <- 40L
  ids <- paste0("s", seq_len(n))
  site <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
  signal <- sin(seq(0, 2 * pi, length.out = p))
  beta <- array(NA_real_, c(p, n, 1L))
  for (i in seq_len(n)) {
    beta[, i, 1] <- signal + ifelse(site[i] == "B", 2, 0) +
      0.05 * cos(seq_len(p) + i)
  }
  g <- new_gds(
    list(beta = beta, var = array(0.04, dim(beta))),
    space_sample_labels(paste0("v", seq_len(p))),
    ids,
    "task",
    col_data = data.frame(site = site, row.names = ids)
  )
  unadjusted <- examine_group(
    reduce(as_plan(g), method = "meta:fe_reg", formula = ~ 1),
    estimands = "(Intercept)"
  )
  adjusted <- examine_group(
    reduce(as_plan(g), method = "meta:fe_reg", formula = ~ site),
    estimands = "siteB"
  )
  expect_gt(
    min(unadjusted$contrast_data$surprise_energy),
    max(adjusted$contrast_data$surprise_energy) * 10
  )
  expect_true(all(adjusted$subject_data$review_status == "none"))
})

test_that("examination is invariant to valid sample block partitions", {
  plan <- reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  small <- examine_group(plan, control = examination_control(block_size = 3L))
  large <- examine_group(plan, control = examination_control(block_size = 17L))

  expect_equal(small$subject_data, large$subject_data, tolerance = 1e-12)
  expect_equal(small$contrast_data, large$contrast_data, tolerance = 1e-12)
  expect_equal(small$estimand_data, large$estimand_data, tolerance = 1e-12)
  expect_identical(
    small$config$retained_subjects,
    large$config$retained_subjects
  )
  expect_identical(
    small$provenance$examination_digest,
    large$provenance$examination_digest
  )
  for (name in names(assays(small$group_maps))) {
    expect_equal(
      assay(small$group_maps, name),
      assay(large$group_maps, name),
      tolerance = 1e-12,
      info = name
    )
  }
})

test_that("memory and native HDF5 examinations agree", {
  skip_if_not_installed("hdf5r")
  g <- .group_examination_fixture()
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  write_gds_h5(g, path)

  memory <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = examination_control(block_size = 9L)
  )
  h5 <- examine_group(
    reduce(gds(path, format = "h5"), method = "meta:fe"),
    control = examination_control(block_size = 9L)
  )
  expect_equal(memory$subject_data, h5$subject_data, tolerance = 1e-12)
  expect_equal(memory$contrast_data, h5$contrast_data, tolerance = 1e-12)
  expect_equal(memory$estimand_data, h5$estimand_data, tolerance = 1e-12)
  expect_equal(
    assay(memory$group_maps, "n_weight"),
    assay(h5$group_maps, "n_weight"),
    tolerance = 1e-12
  )
})

test_that("OLS effect-only examination is explicit and regression maps are method-sensitive", {
  n <- 10L
  p <- 30L
  subjects <- paste0("s", seq_len(n))
  group <- factor(rep(c("control", "patient"), each = n / 2),
                  levels = c("control", "patient"))
  age <- seq(20, 65, length.out = n)
  signal <- sin(seq(0, 2 * pi, length.out = p))
  beta <- array(NA_real_, c(p, n, 1L))
  for (i in seq_len(n)) {
    beta[, i, 1] <- signal + 0.5 * (group[i] == "patient") +
      0.01 * age[i] + 0.05 * cos(seq_len(p) + i)
  }
  g <- new_gds(
    list(beta = beta, var = array(1, dim(beta))),
    space_sample_labels(paste0("v", seq_len(p))),
    subjects,
    "task",
    col_data = data.frame(group = group, age = age, row.names = subjects),
    metadata = list(synthetic_var = TRUE)
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "ols:voxelwise", formula = ~ group + age),
    estimands = "grouppatient",
    control = examination_control(block_size = 8L)
  )
  expect_identical(exam$cohort$variance_mode, "effect_only_synthetic")
  expect_true("max_leverage" %in% names(assays(exam$group_maps)))
  expect_false("n_weight" %in% names(assays(exam$group_maps)))
  expect_false("max_weight_fraction" %in% names(assays(exam$group_maps)))
  expect_identical(unique(exam$estimand_data$estimand), "grouppatient")
})

test_that("missing covariate exclusion and write pruning are auditable", {
  g <- .group_examination_fixture()
  ids <- subjects(g)
  g <- with_col_data(
    g,
    data.frame(age = seq_along(ids), row.names = ids)
  )
  g$col_data$age[2] <- NA_real_
  path <- tempfile(fileext = ".h5")
  plan <- as_plan(g) |>
    reduce(method = "ols:voxelwise", formula = ~ age) |>
    write_out(path, format = "h5")
  exam <- examine_group(
    plan,
    estimands = "age",
    na_action = "exclude",
    control = examination_control(block_size = 13L)
  )
  expect_false(file.exists(path))
  expect_false("s2" %in% exam$subject_data$subject)
  expect_identical(exam$provenance$excluded_subjects$subject, "s2")
  expect_length(exam$provenance$discarded_write_node_ids, 1L)
})

test_that("responsible-subject maps carry categorical non-interpolation metadata", {
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  )
  name <- "argmax_delta_stat:pooled_effect"
  expect_true(name %in% names(assays(exam$group_maps)))
  expect_type(assay(exam$group_maps, name), "integer")
  info <- metadata(exam$group_maps)$examination
  expect_true(name %in% info$categorical_assays)
  expect_true(name %in% info$non_interpolable_assays)
  expect_identical(info$responsible_subject_lookup$subject, exam$subject_data$subject)
})

test_that("unsupported reducers and synthetic precision fail explicitly", {
  g <- .group_examination_fixture()
  expect_error(
    examine_group(as_plan(g), method = "combine:stouffer"),
    "no implemented"
  )
  synthetic <- g
  synthetic$metadata$synthetic_var <- TRUE
  expect_error(
    examine_group(as_plan(synthetic), method = "meta:fe"),
    "genuine first-level variance"
  )
})
