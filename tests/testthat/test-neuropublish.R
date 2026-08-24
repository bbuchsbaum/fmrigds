.neuropublish_fit <- function(n_subjects = 8L,
                              contrasts = c("faces", "places"),
                              method = "meta:re",
                              formula = NULL,
                              missing = FALSE,
                              weights = "1/var",
                              options = list()) {
  set.seed(1701)
  n_voxels <- 4L
  n_contrasts <- length(contrasts)
  beta <- array(
    stats::rnorm(n_voxels * n_subjects * n_contrasts),
    c(n_voxels, n_subjects, n_contrasts)
  )
  variance <- array(
    stats::runif(length(beta), min = 0.08, max = 0.35),
    dim(beta)
  )
  if (isTRUE(missing)) {
    beta[1L, 1L:2L, 1L] <- NA_real_
    variance[1L, 1L:2L, 1L] <- NA_real_
  }
  subject_ids <- paste0("private-subject-", seq_len(n_subjects))
  source <- new_gds(
    assays = list(beta = beta, var = variance),
    space = space_voxel(c(2L, 2L, 1L), diag(4), template_id = "test-grid"),
    subjects = subject_ids,
    contrasts = contrasts,
    col_data = data.frame(
      age = seq(20, 60, length.out = n_subjects),
      row.names = subject_ids
    )
  )
  plan <- as_plan(source)
  plan <- if (is.null(formula)) {
    reduce(plan, method = method, weights = weights, options = options)
  } else {
    reduce(plan, method = method, weights = weights, formula = formula, options = options)
  }
  suppressWarnings(compute(plan))
}

.field_values <- function(spec, member) {
  vapply(spec$fields, `[[`, character(1L), member)
}

test_that("pure mapper publishes cohort size, explicit measures, and stable ids", {
  fit <- .neuropublish_fit(n_subjects = 26L, method = "meta:re", missing = TRUE)
  spec <- .neuropublish_spec_gds(fit, diagnostics = "n_eff")
  repeated <- .neuropublish_spec_gds(fit, diagnostics = "n_eff")

  expect_identical(spec, repeated)
  expect_identical(spec$analysis$sampleSize, 26L)
  expect_identical(spec$analysis$method$reducerId, "meta:re")
  expect_identical(spec$analysis$method$nInputSubjects, 26L)
  expect_identical(length(spec$analysis$estimands), 2L)
  expect_true(all(c("faces", "places") %in% .field_values(spec, "contrast")))
  expect_true(all(c("beta_g", "se_g", "z_g", "tau2", "n_eff") %in% .field_values(spec, "assay")))

  tau <- spec$fields[vapply(spec$fields, function(field) identical(field$assay, "tau2"), logical(1L))]
  expect_length(tau, 2L)
  expect_true(all(vapply(tau, `[[`, character(1L), "label") ==
    "Between-study heterogeneity (\u03c4\u00b2)"))
  expect_true(all(vapply(tau, `[[`, character(1L), "measure") ==
    "org.bbuchsbaum.fmrigds.measure/between-study-heterogeneity-variance"))

  expect_lt(min(assay(fit, "n_eff"), na.rm = TRUE), 26)
  expect_identical(spec$analysis$sampleSize, 26L)
  all_ids <- c(
    spec$domain$id,
    spec$analysis$id,
    vapply(spec$analysis$estimands, `[[`, character(1L), "id"),
    .field_values(spec, "id"),
    .field_values(spec, "assetId")
  )
  expect_identical(anyDuplicated(all_ids), 0L)
  expect_true(all(grepl(.np_local_id_pattern, all_ids)))
})

test_that("fixed effects and permutation outputs use declared mappings only", {
  fixed <- .neuropublish_spec_gds(.neuropublish_fit(method = "meta:fe"))
  expect_identical(
    unique(.field_values(fixed, "assay")),
    c("beta_g", "se_g", "z_g")
  )

  permutation <- .neuropublish_fit(
    n_subjects = 5L,
    contrasts = "task",
    method = "perm:onesample",
    weights = "custom",
    options = list(
      n_perm = 15L,
      seed = 7L,
      custom_weights = c(1, 2, 3, 4, 5)
    )
  )
  spec <- .neuropublish_spec_gds(
    permutation,
    diagnostics = c("p_perm", "p_fwer")
  )
  measures <- stats::setNames(.field_values(spec, "measure"), .field_values(spec, "assay"))
  expect_identical(
    unname(measures[["p_perm"]]),
    "org.bbuchsbaum.fmrigds.measure/permutation-p-value"
  )
  expect_identical(
    unname(measures[["p_fwer"]]),
    "org.bbuchsbaum.fmrigds.measure/fwer-adjusted-p-value"
  )
  expect_identical(spec$analysis$method$weight$status, "verified")
  expect_identical(spec$analysis$method$weight$summary$length, 5L)
  encoded <- .canonical_portable_json(spec)
  expect_false(grepl("custom_weights", encoded, fixed = TRUE))
  expect_false(grepl("private-subject-", encoded, fixed = TRUE))
})

test_that("regression terms remain distinct estimands with a portable design receipt", {
  fit <- .neuropublish_fit(method = "meta:fe_reg", formula = ~ age)
  spec <- .neuropublish_spec_gds(fit)
  labels <- vapply(spec$analysis$estimands, `[[`, character(1L), "label")

  expect_true(all(c(
    "faces \u2014 (Intercept)", "places \u2014 (Intercept)",
    "faces \u2014 age", "places \u2014 age"
  ) %in% labels))
  expect_identical(spec$analysis$method$design$status, "verified")
  expect_identical(spec$analysis$method$design$columns, c("(Intercept)", "age"))
  expect_match(spec$analysis$method$design$digest$value, "^[0-9a-f]{64}$")
  expect_false(grepl("private-subject-", .canonical_portable_json(spec), fixed = TRUE))
})

test_that("hostile custom names never acquire trusted inferential semantics", {
  fit <- .neuropublish_fit(method = "meta:fe", contrasts = "task")
  hostile <- "vendor/t-statistic"
  fit$assays[[hostile]] <- assay(fit, "z_g")
  spec <- .neuropublish_spec_gds(fit, diagnostics = hostile)
  field <- spec$fields[[which(.field_values(spec, "assay") == hostile)]]

  expect_identical(field$measure, "org.bbuchsbaum.fmrigds.measure/custom-scalar-map")
  expect_false(field$trustedMeasure)
  expect_false(grepl("org.neuropublish.measure/t-statistic", field$measure, fixed = TRUE))
})

test_that("publication admission fails closed", {
  arbitrary <- new_gds(
    list(
      beta = array(1, c(4, 1, 1)),
      var = array(0.25, c(4, 1, 1))
    ),
    space_voxel(c(2, 2, 1), diag(4)),
    "meta",
    "task"
  )
  expect_error(.neuropublish_spec_gds(arbitrary), "reduction receipt")

  parcel <- .neuropublish_fit(method = "meta:fe", contrasts = "task")
  parcel$space <- space_parcels(paste0("parcel-", seq_len(4)))
  expect_error(.neuropublish_spec_gds(parcel), "voxel-space")

  categorical <- .neuropublish_fit(method = "meta:fe", contrasts = "task")
  categorical$metadata$categorical_assays <- "z_g"
  expect_error(.neuropublish_spec_gds(categorical, assays = "z_g"), "Categorical")
  expect_error(.neuropublish_spec_gds(categorical, assays = "absent"), "Unknown requested")

  subject_level <- arbitrary
  subject_level$subjects <- "s1"
  expect_error(.neuropublish_spec_gds(subject_level), "group output")

  invalid <- .neuropublish_fit(method = "meta:fe", contrasts = "task")
  invalid$metadata$provenance$graph[[1L]]$inputs <- "missing-source"
  expect_error(.neuropublish_spec_gds(invalid), "invalid provenance graph")
})

test_that("portable projection strips paths and emits honest warnings", {
  fit <- .neuropublish_fit(method = "meta:re", contrasts = "task")
  private_path <- "/Users/private/sub-001/secret-beta.nii.gz"
  fit$metadata$provenance$entities[[1L]]$private <- list(locator = private_path)
  fit$metadata$synthetic_var <- TRUE
  fit$metadata$sample_labels_synthetic <- TRUE
  spec <- .neuropublish_spec_gds(fit)
  encoded <- .canonical_portable_json(spec)
  warning_ids <- vapply(spec$warnings, `[[`, character(1L), "id")

  expect_false(grepl(private_path, encoded, fixed = TRUE))
  expect_false(grepl("secret-beta", encoded, fixed = TRUE))
  expect_true(all(c(
    "synthetic-variance", "synthetic-sample-labels",
    "source-identity-incomplete", "provenance-incomplete"
  ) %in% warning_ids))
})

test_that("Neuropublish dispatch builds an admitted bundle", {
  skip_if_not_installed("neuropublish")
  skip_if_not_installed("neuroim2")

  loadNamespace("neuropublish")
  expect_identical(
    getS3method("as_neuropublish", "gds", envir = asNamespace("neuropublish")),
    as_neuropublish.gds
  )
  fit <- .neuropublish_fit(method = "meta:re", contrasts = "task")
  assets <- tempfile("fmrigds-np-assets-")
  result <- neuropublish::as_neuropublish(fit, staging = assets)
  expect_s3_class(result, "neuropublish_result")
  expect_identical(result$manifest$analyses[[1L]]$sampleSize, 8L)

  staging <- tempfile("fmrigds-np-staging-")
  neuropublish::np_write_bundle(result$manifest, staging)
  expect_true(file.exists(file.path(staging, "manifest.json")))
  manifest_json <- paste(readLines(file.path(staging, "manifest.json"), warn = FALSE), collapse = "\n")
  expect_false(grepl("private-subject-", manifest_json, fixed = TRUE))
  expect_false(grepl(normalizePath(assets, winslash = "/", mustWork = FALSE), manifest_json, fixed = TRUE))

  if (!neuropublish::np_has_npub()) skip("npub CLI is unavailable")
  packed <- neuropublish::np_pack(staging, paste0(staging, ".npub"))
  expect_equal(nrow(neuropublish::np_validate(packed$dir)), 0L)
})
