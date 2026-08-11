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
