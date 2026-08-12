.frame_fit_fixture <- function(instrument = FALSE) {
  set.seed(812L)
  subject_id <- factor(c("s1", "s1", "s2", "s2", "s2", "s3", "s3", "s4", "s4"))
  condition <- factor(c("A", "B", "A", "B", "A", "A", "B", "A", "B"))
  X <- stats::model.matrix(~condition)
  random_intercept <- c(s1 = -0.25, s2 = 0.18, s3 = 0.32, s4 = -0.12)
  feature_effects <- rbind(c(0.8, 0.4), c(1.2, -0.2), c(-0.3, 0.7))
  sampling_variance <- matrix(
    runif(length(subject_id) * 3L, 0.01, 0.05),
    nrow = length(subject_id),
    ncol = 3L
  )
  beta <- vapply(seq_len(3L), function(feature) {
    as.numeric(X %*% feature_effects[feature, ]) +
      random_intercept[as.character(subject_id)] +
      rnorm(length(subject_id), sd = sqrt(sampling_variance[, feature] + 0.03))
  }, numeric(length(subject_id)))
  sources <- list(
    beta = fmridataset::memory_source(beta),
    variance = fmridataset::memory_source(sampling_variance)
  )
  if (instrument) sources <- lapply(sources, fmridataset::counting_source)
  space <- fmridataset::volume_space(
    dim = c(2L, 2L, 1L),
    affine = diag(4),
    support = 1:3,
    template = "frame-fit"
  )
  frame <- fmridataset::fmri_frame(
    assays = sources,
    observations = data.frame(
      .obs_id = paste0("obs-", seq_along(subject_id)),
      subject_id = subject_id,
      condition = condition
    ),
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(space)),
      space = space
    )
  )
  list(
    frame = frame,
    beta = beta,
    variance = sampling_variance,
    subject_id = subject_id,
    condition = condition
  )
}

test_that("fit_group returns a spatially aligned result frame", {
  fixture <- .frame_fit_fixture(instrument = TRUE)
  spec <- multidesign::design_spec(
    fixed = ~condition,
    random = ~ 1 | subject_id
  )
  fit <- fit_group(
    fixture$frame,
    design = spec,
    block_size = 1L,
    memory_budget = 1024^2
  )

  expect_s3_class(fit, "fmri_group_fit")
  expect_s3_class(fit$result, "fmri_frame")
  expect_identical(fit$result$metadata$result_schema_version, 1L)
  expect_identical(fit$result$metadata$result_kind, "statistical")
  expect_s3_class(fit$result$provenance, "provenance_graph")
  expect_identical(
    fmridataset::feature_ids(fit$result),
    fmridataset::feature_ids(fixture$frame)
  )
  expect_identical(
    fmridataset::space_digest(fmridataset::space(fit$result)),
    fmridataset::space_digest(fmridataset::space(fixture$frame))
  )
  expect_identical(
    fmridataset::observation_ids(fit$result),
    c("(Intercept)", "conditionB", "vc_intercept", "vc_resid")
  )
  expect_true(all(fit$diagnostics$converged))
  expect_identical(
    fmridataset::source_counts(fmridataset::assay(fixture$frame, "beta")$source)$reads,
    3
  )
  expect_identical(
    fmridataset::source_counts(fmridataset::assay(fixture$frame, "variance")$source)$reads,
    3
  )
})

test_that("frame fits are invariant to feature block width", {
  fixture <- .frame_fit_fixture()
  spec <- multidesign::design_spec(~condition, ~ 1 | subject_id)
  one <- fit_group(fixture$frame, design = spec, block_size = 1L)
  all <- fit_group(fixture$frame, design = spec, block_size = 3L)

  for (assay_name in names(fmridataset::assays(one$result))) {
    expect_equal(
      fmridataset::collect_assay(one$result, assay_name),
      fmridataset::collect_assay(all$result, assay_name),
      tolerance = 1e-10,
      info = assay_name
    )
  }
})

test_that("variance-aware frame fit agrees with metafor", {
  skip_if_not_installed("metafor")
  fixture <- .frame_fit_fixture()
  fit <- fit_group(
    fixture$frame,
    design = multidesign::design_spec(~condition, ~ 1 | subject_id)
  )
  result <- fmridataset::collect_assay(fit$result, "estimate")
  standard_error <- fmridataset::collect_assay(fit$result, "std_error")

  for (feature in seq_len(ncol(fixture$beta))) {
    data <- data.frame(
      yi = fixture$beta[, feature],
      vi = fixture$variance[, feature],
      condition = fixture$condition,
      subject = fixture$subject_id,
      observation = factor(seq_len(nrow(fixture$beta)))
    )
    reference <- metafor::rma.mv(
      yi,
      V = vi,
      mods = ~condition,
      random = list(~ 1 | subject, ~ 1 | observation),
      data = data,
      method = "REML"
    )
    expect_equal(result[1:2, feature], unname(stats::coef(reference)), tolerance = 2e-3)
    expect_equal(standard_error[1:2, feature], reference$se, tolerance = 2e-3)
    expect_equal(result[3:4, feature], reference$sigma2, tolerance = 3e-3)
  }
})

test_that("fit_group enforces variance random-effect and budget contracts", {
  fixture <- .frame_fit_fixture()
  spec <- multidesign::design_spec(~condition, ~ 1 | subject_id)

  bad <- fixture$frame
  bad$assays$variance$source$data[1L, 1L] <- 0
  expect_error(fit_group(bad, design = spec), "strictly positive")
  expect_error(
    fit_group(
      fixture$frame,
      design = multidesign::design_spec(
        ~condition,
        ~ 1 + condition | subject_id
      )
    ),
    "random intercept only"
  )
  expect_error(
    fit_group(fixture$frame, design = spec, block_size = 3L, memory_budget = 16L),
    "memory_budget"
  )
})

test_that("fit_group honors synchronized observation and feature views", {
  fixture <- .frame_fit_fixture()
  view <- fixture$frame[c(1L, 2L, 3L, 4L, 6L, 7L, 8L, 9L), 1:2]
  fit <- fit_group(
    view,
    design = multidesign::design_spec(~condition, ~ 1 | subject_id),
    block_size = 1L
  )
  expect_identical(
    fmridataset::feature_ids(fit$result),
    fmridataset::feature_ids(view)
  )
  expect_identical(ncol(fit$result), 2L)
})

test_that("group_plan is metadata-only, serializable, and delegates through the registry", {
  fixture <- .frame_fit_fixture(instrument = TRUE)
  plan <- group_plan(
    fixture$frame,
    design = multidesign::design_spec(~condition, ~ 1 | subject_id),
    method = "lmm:ri_knownvar",
    block_size = 2L
  )

  expect_s3_class(plan, "fmri_group_plan")
  expect_identical(plan$method, "lmm:ri_knownvar")
  expect_true(is.function(get_reducer(plan$method)$frame_fun))
  expect_gt(length(serialize(plan, NULL)), 0L)
  expect_identical(digest_plan(unserialize(serialize(plan, NULL))), digest_plan(plan))
  explained <- explain(plan)
  expect_identical(explained$digest, digest_plan(plan))
  expect_identical(explained$method, "lmm:ri_knownvar")
  expect_identical(
    fmridataset::source_counts(fmridataset::assay(fixture$frame, "beta")$source)$reads,
    0
  )
  expect_identical(
    fmridataset::source_counts(fmridataset::assay(fixture$frame, "variance")$source)$reads,
    0
  )

  planned <- compute(plan)
  eager <- fit_group(
    fixture$frame,
    design = multidesign::design_spec(~condition, ~ 1 | subject_id),
    block_size = 2L
  )
  for (assay_name in names(fmridataset::assays(planned$result))) {
    expect_equal(
      fmridataset::collect_assay(planned$result, assay_name),
      fmridataset::collect_assay(eager$result, assay_name),
      tolerance = 1e-10
    )
  }
})

test_that("group plans map omitted design rows back to frame observations", {
  fixture <- .frame_fit_fixture()
  frame <- fixture$frame
  frame$observations$data$condition[3L] <- NA
  spec <- multidesign::design_spec(
    ~condition,
    ~ 1 | subject_id,
    na_action = "omit"
  )

  omitted <- fit_group(frame, design = spec, block_size = 2L)
  reference <- fit_group(
    frame[-3L, ],
    design = multidesign::design_spec(~condition, ~ 1 | subject_id),
    block_size = 2L
  )

  expect_identical(
    omitted$plan$observation_ids,
    fmridataset::observation_ids(frame)[-3L]
  )
  for (assay_name in names(fmridataset::assays(omitted$result))) {
    expect_equal(
      fmridataset::collect_assay(omitted$result, assay_name),
      fmridataset::collect_assay(reference$result, assay_name),
      tolerance = 1e-10
    )
  }
})

test_that("known-variance frame reduction handles missingness featurewise", {
  fixture <- .frame_fit_fixture()
  frame <- fixture$frame
  frame$assays$variance$source$data[1L, 1L] <- NA_real_
  fit <- fit_group(
    frame,
    design = multidesign::design_spec(~condition, ~ 1 | subject_id),
    block_size = 2L
  )

  expect_identical(fit$diagnostics$n_obs, c(8L, 9L, 9L))
  expect_true(all(is.finite(
    fmridataset::collect_assay(fit$result, "estimate")[1:2, ]
  )))

  bad <- fixture$frame
  bad$assays$variance$source$data[1L, 1L] <- -0.01
  expect_error(
    fit_group(bad, design = multidesign::design_spec(~condition, ~ 1 | subject_id)),
    "strictly positive"
  )
})

test_that("OLS frame reducer agrees with dense lm fits and block widths", {
  fixture <- .frame_fit_fixture()
  spec <- multidesign::design_spec(~condition)
  one <- compute(group_plan(
    fixture$frame,
    design = spec,
    method = "ols:voxelwise",
    variance = NULL,
    block_size = 1L
  ))
  all <- compute(group_plan(
    fixture$frame,
    design = spec,
    method = "ols:voxelwise",
    variance = NULL,
    block_size = 3L
  ))

  actual <- fmridataset::collect_assay(one$result, "estimate")
  reference <- vapply(seq_len(ncol(fixture$beta)), function(feature) {
    stats::coef(stats::lm(fixture$beta[, feature] ~ fixture$condition))
  }, numeric(2L))
  expect_equal(actual, reference, tolerance = 1e-10, ignore_attr = TRUE)
  for (assay_name in names(fmridataset::assays(one$result))) {
    expect_equal(
      fmridataset::collect_assay(one$result, assay_name),
      fmridataset::collect_assay(all$result, assay_name),
      tolerance = 1e-10
    )
  }
  expect_identical(one$diagnostics$n_obs, rep(9L, 3L))

  eager <- fit_group(
    fixture$frame,
    design = spec,
    method = "ols:voxelwise",
    variance = NULL,
    block_size = 2L
  )
  expect_equal(
    fmridataset::collect_assay(eager$result, "estimate"),
    reference,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("OLS frame plans ignore an unused variance default", {
  fixture <- .frame_fit_fixture(instrument = TRUE)
  frame <- fixture$frame
  variance_source <- fmridataset::assay(frame, "variance")$source
  frame$assays$variance <- NULL

  plan <- group_plan(
    frame,
    design = multidesign::design_spec(~condition),
    method = "ols:voxelwise",
    memory_budget = 8L * nrow(frame)
  )
  expect_null(plan$variance)
  expect_identical(plan$block_size, 1L)
  expect_s3_class(compute(plan), "fmri_group_fit")
  expect_identical(fmridataset::source_counts(variance_source)$reads, 0)
})

test_that("frame plans reject reducers without a frame execution contract", {
  fixture <- .frame_fit_fixture()
  expect_error(
    group_plan(
      fixture$frame,
      design = multidesign::design_spec(~ 1),
      method = "combine:fisher",
      variance = NULL
    ),
    "does not support fmri_frame"
  )
  expect_error(
    group_plan(
      fixture$frame,
      design = multidesign::design_spec(~condition, ~ 1 | subject_id),
      method = "lmm:ri_knownvar",
      variance = NULL
    ),
    "requires frame inputs: var"
  )
  expect_error(
    group_plan(
      fixture$frame[, integer()],
      design = multidesign::design_spec(~condition, ~ 1 | subject_id)
    ),
    "at least one feature"
  )
  expect_error(
    group_plan(
      fixture$frame,
      design = multidesign::design_spec(~condition, ~ 1 | subject_id),
      memory_budget = 1L
    ),
    "cannot hold one"
  )
})

test_that("frame plan preflight rejects source drift before reading", {
  fixture <- .frame_fit_fixture(instrument = TRUE)
  plan <- group_plan(
    fixture$frame,
    design = multidesign::design_spec(~condition, ~ 1 | subject_id)
  )
  plan$frame$assays$beta$source <- fmridataset::memory_source(
    fixture$beta + 1
  )

  expect_error(compute(plan), "source changed")
  expect_identical(
    fmridataset::source_counts(fmridataset::assay(fixture$frame, "variance")$source)$reads,
    0
  )
})
