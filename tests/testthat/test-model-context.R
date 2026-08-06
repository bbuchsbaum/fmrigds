.model_context_fixture <- function() {
  subjects <- c("s1", "s2", "s3", "s4", "s5")
  beta <- array(seq_len(3 * 5), c(3, 5, 1))
  g <- new_gds(
    assays = list(beta = beta, var = array(0.25, dim(beta))),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = subjects,
    contrasts = "task",
    col_data = data.frame(
      group = factor(
        c("control", "patient", "control", "patient", "control"),
        levels = c("control", "patient")
      ),
      age = c(20, 30, 40, 50, 60),
      row.names = subjects
    )
  )
  g
}

test_that("model context is built after subject subset and reorder", {
  plan <- as_plan(.model_context_fixture()) |>
    subset(subject = c("s5", "s2", "s4", "s1")) |>
    reduce(method = "meta:fe_reg", formula = ~ group + age)
  compiled <- fmrigds:::compile_examination_plan(plan)
  context <- fmrigds:::.build_reducer_model_context(
    compiled,
    estimands = "grouppatient"
  )

  expect_s3_class(context, "reducer_model_context")
  expect_identical(context$subjects, c("s5", "s2", "s4", "s1"))
  expect_identical(rownames(context$X), context$subjects)
  expect_equal(unname(context$X[, "age"]), c(60, 30, 50, 20))
  expect_identical(colnames(context$X), c("(Intercept)", "grouppatient", "age"))
  expect_identical(context$factor_levels$group, c("control", "patient"))
  expect_identical(colnames(context$estimand_matrix), colnames(context$X))
  expect_identical(rownames(context$estimand_matrix), "grouppatient")
  expect_equal(unname(context$estimand_matrix[1, ]), c(0, 1, 0))
  expect_null(plan$nodes[[2L]]$options$X)
})

test_that("missing covariates fail by default and exclusion is recorded", {
  g <- .model_context_fixture()
  g$col_data$age[c(2, 5)] <- NA_real_
  compiled <- as_plan(g) |>
    reduce(method = "ols:voxelwise", formula = ~ group + age) |>
    fmrigds:::compile_examination_plan()

  expect_error(
    fmrigds:::.build_reducer_model_context(compiled),
    "s2, s5"
  )
  context <- fmrigds:::.build_reducer_model_context(
    compiled,
    na_action = "exclude"
  )
  expect_identical(context$subjects, c("s1", "s3", "s4"))
  expect_identical(context$excluded_subjects$subject, c("s2", "s5"))
  expect_true(all(context$excluded_subjects$reason == "missing_covariate"))
  expect_identical(context$included_subject_index, c(1L, 3L, 4L))
  expect_identical(rownames(context$X), context$subjects)
})

test_that("context digest includes covariates and realized subject order", {
  g <- .model_context_fixture()
  build <- function(object, order = c("s1", "s2", "s3", "s4", "s5")) {
    as_plan(object) |>
      subset(subject = order) |>
      reduce(method = "ols:voxelwise", formula = ~ group + age) |>
      fmrigds:::compile_examination_plan() |>
      fmrigds:::.build_reducer_model_context()
  }
  original <- build(g)
  changed <- g
  changed$col_data$age[1L] <- 21
  reordered <- build(g, rev(subjects(g)))

  expect_false(identical(original$digest, build(changed)$digest))
  expect_false(identical(original$digest, reordered$digest))
})

test_that("linear estimands retain exact design coordinates", {
  compiled <- as_plan(.model_context_fixture()) |>
    reduce(method = "ols:voxelwise", formula = ~ group + age) |>
    fmrigds:::compile_examination_plan()

  context <- fmrigds:::.build_reducer_model_context(
    compiled,
    estimands = list(
      age = c(age = 1),
      patient_at_40 = c("(Intercept)" = 1, grouppatient = 1, age = 40)
    )
  )
  expect_identical(rownames(context$estimand_matrix), c("age", "patient_at_40"))
  expect_equal(
    unname(context$estimand_matrix["age", ]),
    c(0, 0, 1)
  )
  expect_equal(
    unname(context$estimand_matrix["patient_at_40", ]),
    c(1, 1, 40)
  )

  expect_error(
    fmrigds:::.build_reducer_model_context(compiled, estimands = "unknown"),
    "Unknown estimand"
  )
})

test_that("nonestimable linear combinations are rejected", {
  g <- .model_context_fixture()
  g$col_data$twice_age <- 2 * g$col_data$age
  compiled <- as_plan(g) |>
    reduce(method = "ols:voxelwise", formula = ~ age + twice_age) |>
    fmrigds:::compile_examination_plan()

  expect_error(
    fmrigds:::.build_reducer_model_context(
      compiled,
      estimands = list(age_only = c(age = 1))
    ),
    "not estimable"
  )
})

test_that("reducer preflight reports semantics and implementation availability", {
  fe <- as_plan(.model_context_fixture()) |>
    reduce(method = "meta:fe") |>
    fmrigds:::compile_examination_plan() |>
    fmrigds:::.build_reducer_model_context()
  fe_preflight <- fmrigds:::.preflight_reducer_diagnostics(fe)
  expect_false(fe_preflight$model_contract$uses_X)
  expect_identical(fe_preflight$model_contract$weight_mode, "inverse_variance")
  expect_true("prediction" %in% fe_preflight$diagnostics$capabilities)
  expect_true(all(fe_preflight$availability$status == "available"))

  register_reducer(
    "test:no-diagnostics",
    function(beta, var, X, z, p, df, df1, df2, opts) {
      list(beta_g = colMeans(beta))
    },
    requires = "beta",
    provides = "beta_g"
  )
  compiled <- fmrigds:::compile_examination_plan(
    as_plan(.model_context_fixture()),
    method = "test:no-diagnostics"
  )
  custom <- fmrigds:::.build_reducer_model_context(compiled)
  preflight <- fmrigds:::.preflight_reducer_diagnostics(custom)
  expect_null(preflight$model_contract)
  expect_true(all(preflight$availability$status == "unsupported_reducer"))
})

test_that("reducers that ignore X reject covariate formulas", {
  expect_error(
    as_plan(.model_context_fixture()) |>
      reduce(method = "meta:re", formula = ~ age),
    "does not consume a design matrix"
  )
})

test_that("examination compiler rejects model overrides on an existing reducer", {
  plan <- as_plan(.model_context_fixture()) |>
    reduce(method = "meta:fe")
  expect_error(
    fmrigds:::compile_examination_plan(plan, method = "ols:voxelwise"),
    "inherits"
  )
  expect_error(
    fmrigds:::compile_examination_plan(plan, formula = ~ age),
    "inherits"
  )
})

test_that("examination_control validates named subsystem settings", {
  control <- examination_control(
    block_size = 64L,
    geometry = list(rank = 16L, oversample = 4L, cap = 6),
    review = list(min_stability = 0.8),
    exact_refit_n = 3L
  )
  expect_s3_class(control, "gds_examination_control")
  expect_identical(control$block_size, 64L)
  expect_identical(control$geometry$rank, 16L)
  expect_equal(control$review$min_stability, 0.8)
  expect_error(examination_control(block_size = 0), "positive")
  expect_error(
    examination_control(review = list(min_stability = 2)),
    "between 0 and 1"
  )
  expect_error(
    examination_control(geometry = list(mystery = 1)),
    "Unknown geometry"
  )
})
