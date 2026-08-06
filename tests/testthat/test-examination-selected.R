.quiet_review_control <- function(block_size = 5L) {
  examination_control(
    block_size = block_size,
    retain_n = 0L,
    review = list(
      surprise = list(
        energy_threshold = 100,
        tail_threshold = 1,
        residual_threshold = 100
      ),
      influence = list(
        energy_threshold = 100,
        max_abs_threshold = 100
      )
    )
  )
}

test_that("second pass retains maps only for the frozen selected set", {
  g <- .group_examination_fixture()
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    retain = "s3",
    control = .quiet_review_control(block_size = 6L)
  )
  expect_identical(subjects(exam$subject_maps), "s3")
  expect_true(all(c(
    "observed", "expected", "predictive_residual",
    "fit_contribution:pooled_effect",
    "delta_effect:pooled_effect",
    "delta_stat:pooled_effect"
  ) %in% names(assays(exam$subject_maps))))
  expect_equal(
    drop(assay(exam$subject_maps, "observed")[, 1, 1]),
    drop(assay(g, "beta")[, 3, 1]),
    tolerance = 0
  )
  expect_identical(
    metadata(exam$subject_maps)$examination$selection_frozen_before_exact_refit,
    TRUE
  )
  expect_equal(exam$provenance$scan_receipt$retained_map_count, 6L)
})

test_that("selected random-effects exact refit matches brute-force deletion", {
  set.seed(42)
  beta <- array(rnorm(15 * 8, 0.8, 0.6), c(15, 8, 1))
  var <- array(runif(15 * 8, 0.08, 0.35), dim(beta))
  g <- new_gds(
    list(beta = beta, var = var),
    space_sample_labels(paste0("v", seq_len(15))),
    paste0("s", seq_len(8)),
    "task"
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:re"),
    retain = "s2",
    control = .quiet_review_control(block_size = 4L)
  )
  expect_identical(subjects(exam$subject_maps), "s2")
  expect_identical(
    sort(unique(exam$estimand_data$mode)),
    c("tau2_fixed_full", "tau2_refit_exact")
  )
  expect_equal(
    sum(exam$estimand_data$ranking_stage == "selected_refit"),
    1L
  )

  beta_matrix <- t(beta[, , 1])
  var_matrix <- t(var[, , 1])
  full <- fmrigds:::core_meta_re_dl_kernel(beta_matrix, var_matrix)
  deleted <- fmrigds:::core_meta_re_dl_kernel(
    beta_matrix[-2, , drop = FALSE],
    var_matrix[-2, , drop = FALSE]
  )
  expected_delta <- full$z_g - deleted$z_g
  expect_equal(
    drop(assay(exam$subject_maps, "delta_stat_exact:pooled_effect")[, 1, 1]),
    expected_delta,
    tolerance = 1e-10
  )
  expect_equal(
    drop(assay(exam$subject_maps, "tau2_deleted_exact")[, 1, 1]),
    deleted$tau2,
    tolerance = 1e-12
  )
  exact_subjects <- exam$estimand_data$subject[
    exam$estimand_data$ranking_stage == "selected_refit"
  ]
  expect_identical(exact_subjects, "s2")
  expect_true(all(exam$subject_data$review_status == "none"))
})
