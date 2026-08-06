test_that("BH conclusions are recomputed only for retained subjects", {
  g <- .group_examination_fixture()
  plan <- as_plan(g) |>
    reduce(method = "meta:fe") |>
    posthoc("fdr:bh")
  exam <- examine_group(
    plan,
    control = examination_control(retain_n = 3L)
  )
  expect_s3_class(exam$conclusion$full_maps, "gds")
  expect_s3_class(exam$conclusion$deleted_maps, "gds")
  expect_identical(
    subjects(exam$conclusion$deleted_maps),
    exam$config$retained_subjects
  )
  expect_identical(subjects(exam$conclusion$full_maps), "meta")
  expect_true(all(exam$conclusion$results$status == "available"))
  expect_true(all(exam$conclusion$results$mode == "exact"))
  expect_true(all(exam$conclusion$results$method == "fdr:bh"))
  expect_true(all(c(
    "full_significant_n", "deleted_significant_n", "transition_count",
    "gained_n", "lost_n"
  ) %in% names(exam$conclusion$results)))

  full <- compute(plan)
  expect_equal(
    unname(assay(exam$conclusion$full_maps, "adjusted_p:fdr_bh:pooled_effect")),
    unname(assay(full, "q")),
    tolerance = 1e-12
  )

  selected <- exam$config$retained_subjects[1L]
  brute <- as_plan(g) |>
    subset(subject = setdiff(subjects(g), selected)) |>
    reduce(method = "meta:fe") |>
    posthoc("fdr:bh") |>
    compute()
  expect_equal(
    unname(assay(exam$conclusion$deleted_maps, "adjusted_p:fdr_bh:pooled_effect")[
      , match(selected, subjects(exam$conclusion$deleted_maps)), , drop = FALSE
    ]),
    unname(assay(brute, "q")),
    tolerance = 1e-10
  )
})

test_that("random-effects conclusions use exact selected refits without reranking", {
  exam <- examine_group(
    as_plan(.group_examination_fixture()) |>
      reduce(method = "meta:re") |>
      posthoc("fdr:by"),
    control = examination_control(retain_n = 2L)
  )
  expect_true(all(exam$conclusion$results$mode == "tau2_refit_exact"))
  expect_identical(
    subjects(exam$conclusion$deleted_maps),
    exam$config$retained_subjects
  )
  expect_identical(
    exam$subject_data$subject[exam$subject_data$retained],
    exam$config$retained_subjects
  )
})

test_that("OLS conclusion recomputation uses deleted degrees of freedom", {
  g <- .group_examination_fixture()
  ids <- subjects(g)
  g <- with_col_data(
    g,
    data.frame(
      group = factor(rep(c("control", "patient"), each = 5)),
      row.names = ids
    )
  )
  plan <- as_plan(g) |>
    reduce(method = "ols:voxelwise", formula = ~ group) |>
    posthoc("fdr:bh", options = list(source = "p_coef:grouppatient"))
  exam <- examine_group(
    plan,
    estimands = "grouppatient",
    control = examination_control(retain_n = 2L)
  )
  full <- compute(plan)
  expect_equal(
    unname(assay(
      exam$conclusion$full_maps,
      "adjusted_p:fdr_bh:grouppatient"
    )),
    unname(assay(full, "q_coef:grouppatient")),
    tolerance = 1e-10
  )
  selected <- exam$config$retained_subjects[1L]
  brute <- as_plan(g) |>
    subset(subject = setdiff(ids, selected)) |>
    reduce(method = "ols:voxelwise", formula = ~ group) |>
    posthoc("fdr:bh", options = list(source = "p_coef:grouppatient")) |>
    compute()
  expect_equal(
    unname(assay(
      exam$conclusion$deleted_maps,
      "adjusted_p:fdr_bh:grouppatient"
    )[, match(selected, subjects(exam$conclusion$deleted_maps)), , drop = FALSE]),
    unname(assay(brute, "q_coef:grouppatient")),
    tolerance = 1e-8
  )
})

test_that("unsupported conclusion tails are declared and never executed", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  register_posthoc(
    "test:conclusion-unsupported",
    function(arrays, opts) {
      calls$n <- calls$n + 1L
      list(q = arrays$p)
    },
    requires = "p",
    provides = "q",
    case_deletion = NULL
  )
  on.exit(unregister_posthoc("test:conclusion-unsupported"), add = TRUE)
  exam <- examine_group(
    as_plan(.group_examination_fixture()) |>
      reduce(method = "meta:fe") |>
      posthoc("test:conclusion-unsupported"),
    control = examination_control(retain_n = 2L)
  )
  expect_identical(calls$n, 0L)
  expect_null(exam$conclusion$full_maps)
  expect_null(exam$conclusion$deleted_maps)
  expect_true(all(exam$conclusion$results$status == "unsupported_method"))
  expect_false(exam$conclusion$availability$supported)
})

test_that("stochastic selected recomputation requires an explicit seed", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  register_posthoc(
    "test:seeded-conclusion",
    function(arrays, opts) {
      calls$n <- calls$n + 1L
      list(q = arrays$p)
    },
    requires = "p",
    provides = "q",
    case_deletion = list(
      supported = TRUE,
      mode = "selected_recompute",
      requires_full_space = TRUE,
      deterministic = FALSE,
      seed_contract = "options$seed"
    )
  )
  on.exit(unregister_posthoc("test:seeded-conclusion"), add = TRUE)
  unseeded <- examine_group(
    as_plan(.group_examination_fixture()) |>
      reduce(method = "meta:fe") |>
      posthoc("test:seeded-conclusion"),
    control = examination_control(retain_n = 2L)
  )
  expect_identical(calls$n, 0L)
  expect_true(all(unseeded$conclusion$results$status == "seed_required"))

  seeded <- examine_group(
    as_plan(.group_examination_fixture()) |>
      reduce(method = "meta:fe") |>
      posthoc("test:seeded-conclusion", options = list(seed = 42L)),
    control = examination_control(retain_n = 2L)
  )
  expect_identical(calls$n, 1L + length(seeded$config$retained_subjects))
  expect_true(all(seeded$conclusion$results$status == "available"))
  expect_true(all(seeded$conclusion$results$seed == 42L))
})
