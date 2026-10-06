# Regression tests for review findings in frame projection and examination.

.review_quiet_control <- function(block_size = 5L) {
  examination_control(
    block_size = block_size,
    retain_n = 0L,
    review = list(
      surprise = list(
        energy_threshold = 100,
        tail_threshold = 1,
        residual_threshold = 100
      ),
      influence = list(energy_threshold = 100, max_abs_threshold = 100)
    )
  )
}

.review_projection_frame <- function() {
  subjects <- c("sub-02", "sub-01", "sub-03")
  contrasts <- c("task", "baseline")
  grid <- expand.grid(
    subject_id = subjects,
    contrast_id = contrasts,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$.obs_id <- paste(grid$subject_id, grid$contrast_id, sep = "::")
  grid$condition <- factor(ifelse(grid$contrast_id == "task", "on", "off"))
  beta <- matrix(seq_len(nrow(grid) * 4L), nrow = nrow(grid)) / 10
  variance <- beta / 10 + 0.05
  spatial <- fmridataset::volume_space(
    dim = c(2L, 2L, 1L),
    affine = diag(4),
    support = 1:4,
    template = "review-fixture"
  )
  frame <- fmridataset::fmri_frame(
    assays = list(
      beta = fmridataset::memory_source(beta),
      variance = fmridataset::memory_source(variance)
    ),
    observations = grid,
    features = fmridataset::feature_axis(
      data.frame(.feature_id = fmridataset::feature_ids(spatial)),
      space = spatial
    )
  )
  list(frame = frame, beta = beta, grid = grid)
}

test_that("as_fmri_frame regenerates stale projection IDs after reduce()", {
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("multidesign")
  fixture <- .review_projection_frame()
  g <- as_gds(fixture$frame)
  reduced <- compute(one_sample(g))
  expect_false(identical(subjects(reduced), subjects(g)))

  restored <- fmridataset::as_fmri_frame(reduced)
  expect_s3_class(restored, "fmri_frame")
  observations <- fmridataset::observations(restored)
  expect_identical(nrow(observations), length(subjects(reduced)) * length(contrasts(reduced)))
  expect_identical(
    fmridataset::observation_ids(restored),
    paste(observations$subject_id, observations$contrast_id, sep = "::")
  )
  expect_s3_class(restored$provenance, "provenance_graph")
  expect_true(any(vapply(
    fmridataset::provenance_records(restored$provenance),
    function(record) identical(record$operation, "fmrigds::as_fmri_frame"),
    logical(1)
  )))
})

test_that("as_fmri_frame does not misattach stored IDs after a subject reorder", {
  skip_if_not_installed("fmridataset")
  skip_if_not_installed("multidesign")
  fixture <- .review_projection_frame()
  g <- as_gds(fixture$frame)
  reordered <- suppressWarnings(
    compute(subset(as_plan(g), subject = rev(subjects(g))))
  )
  expect_identical(subjects(reordered), rev(subjects(g)))

  restored <- fmridataset::as_fmri_frame(reordered)
  observations <- fmridataset::observations(restored)
  expect_identical(
    observations$.obs_id,
    paste(observations$subject_id, observations$contrast_id, sep = "::")
  )
  # Stored annotations travel with their subject-contrast cell.
  expected_condition <- ifelse(observations$contrast_id == "task", "on", "off")
  expect_identical(as.character(observations$condition), expected_condition)
  beta <- fmridataset::collect_assay(restored, "beta")
  source_row <- match(observations$.obs_id, fixture$grid$.obs_id)
  expect_equal(beta, fixture$beta[source_row, , drop = FALSE], ignore_attr = TRUE)
})

test_that("frame paths fail informatively without fmridataset/multidesign", {
  local_mocked_bindings(.frame_dep_available = function(package) FALSE)
  expect_error(.require_frame_deps("fit_group()"), "fmridataset")
  expect_error(.require_frame_deps("fit_group()"), "multidesign")
})

test_that("deletion of a non-contributing subject keeps the BH family intact", {
  g <- .group_examination_fixture()
  beta <- assay(g, "beta")
  beta[1:10, 3, 1] <- NA_real_
  g <- new_gds(
    list(beta = beta, var = assay(g, "var")),
    space(g), subjects(g), contrasts(g)
  )
  plan <- as_plan(g) |> reduce(method = "meta:fe") |> posthoc("fdr:bh")
  exam <- examine_group(plan, retain = "s3", control = .review_quiet_control())

  delta <- assay(exam$subject_maps, "delta_stat:pooled_effect")[, "s3", 1]
  expect_true(all(delta[1:10] == 0))
  brute <- as_plan(g) |>
    subset(subject = setdiff(subjects(g), "s3")) |>
    reduce(method = "meta:fe") |>
    posthoc("fdr:bh") |>
    compute()
  deleted_q <- assay(
    exam$conclusion$deleted_maps, "adjusted_p:fdr_bh:pooled_effect"
  )[, "s3", 1]
  expect_false(anyNA(deleted_q))
  expect_equal(unname(deleted_q), unname(assay(brute, "q")[, 1, 1]), tolerance = 1e-10)
})

test_that("stability_replicates = 0 disables stability gating without crashing", {
  g <- .group_examination_fixture()
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = examination_control(geometry = list(stability_replicates = 0L))
  )
  expect_s3_class(exam, "gds_examination")
  expect_true(all(is.na(exam$contrast_data$surprise_stability)))
  expect_true(any(exam$subject_data$review_status == "review"))
})

test_that("exact random-effects refits are capped by exact_refit_n", {
  g <- .group_examination_fixture()
  exam <- examine_group(
    as_plan(g) |> reduce(method = "meta:re") |> posthoc("fdr:bh"),
    control = examination_control(retain_n = 4L, exact_refit_n = 2L)
  )
  expect_gte(sum(exam$subject_data$retained), 3L)
  refit <- exam$estimand_data[exam$estimand_data$ranking_stage == "selected_refit", ]
  expect_lte(nrow(refit), 2L)
  priority <- exam$subject_data[exam$subject_data$retained, ]
  top <- head(priority$subject[order(-priority$review_priority, priority$subject)], 2L)
  expect_setequal(refit$subject, top)
  # Retained subjects outside the cap still receive screening conclusions.
  expect_true(all(exam$conclusion$results$status == "available"))
  others <- setdiff(exam$config$retained_subjects, top)
  expect_true(all(
    exam$conclusion$results$mode[exam$conclusion$results$subject %in% others] ==
      "tau2_fixed_full"
  ))
})

.reference_linear_deletion <- function(beta, X, estimands) {
  n_subject <- nrow(beta)
  n_sample <- ncol(beta)
  out <- array(NA_real_, c(n_subject, nrow(estimands), n_sample))
  stat <- function(y, Xm) {
    fit <- stats::lm.fit(Xm, y)
    df <- length(y) - ncol(Xm)
    s2 <- sum(fit$residuals^2) / df
    A <- solve(crossprod(Xm))
    drop(estimands %*% fit$coefficients) /
      sqrt(diag(estimands %*% A %*% t(estimands)) * s2)
  }
  for (b in seq_len(n_sample)) {
    ok <- is.finite(beta[, b])
    full <- stat(beta[ok, b], X[ok, , drop = FALSE])
    for (i in which(ok)) {
      keep <- ok & seq_len(n_subject) != i
      out[i, , b] <- full - stat(beta[keep, b], X[keep, , drop = FALSE])
    }
  }
  out
}

test_that("streamlined deletion diagnostics match brute-force refits", {
  set.seed(11)
  n <- 9L
  p <- 6L
  X <- cbind(1, rep(c(0, 1), length.out = n))
  colnames(X) <- c("(Intercept)", "group")
  estimands <- diag(2)
  dimnames(estimands) <- list(colnames(X), colnames(X))
  beta <- matrix(rnorm(n * p), n, p)
  beta[2, 3] <- NA_real_
  fit <- suppressWarnings(get_reducer("ols:voxelwise")$fun(
    beta, NULL, X, NULL, NULL, NULL, NULL, NULL, list()
  ))
  diagnostic <- .diagnose_ols_block(
    fit, beta, NULL, X, estimands, list(), examination_control()$tolerance
  )
  reference <- .reference_linear_deletion(beta, X, estimands)
  finite <- is.finite(reference)
  expect_equal(diagnostic$delta_stat[finite], reference[finite], tolerance = 1e-8)
  # The structurally absent subject has an unchanged statistic.
  expect_equal(diagnostic$delta_stat[2, , 3], c(0, 0))
})

test_that("influence histogram quantiles match a per-value reference", {
  g <- .group_examination_fixture()
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = examination_control(retain_n = 10L)
  )
  upper <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8, 16, Inf)
  screening <- exam$estimand_data[exam$estimand_data$ranking_stage == "screening", ]
  for (subject in subjects(exam$subject_maps)) {
    values <- abs(assay(exam$subject_maps, "delta_stat:pooled_effect")[, subject, 1])
    values <- values[is.finite(values)]
    counts <- numeric(length(upper))
    for (value in values) {
      bin <- which(value <= upper)[1L]
      counts[bin] <- counts[bin] + 1
    }
    row <- screening[screening$subject == subject, ]
    expect_identical(
      row$abs_delta_q95_approx,
      .histogram_quantile(counts, upper, 0.95, max(values))
    )
    expect_identical(row$eligible_n, as.numeric(length(values)))
  }
})

test_that("exact random-effects regression refits match brute-force deletion", {
  set.seed(5)
  n <- 8L
  p <- 5L
  X <- cbind(1, rep(c(0, 1), length.out = n))
  colnames(X) <- c("(Intercept)", "group")
  estimands <- matrix(c(0, 1), 1, dimnames = list("group", colnames(X)))
  beta <- matrix(rnorm(n * p, 0.5), n, p)
  var <- matrix(runif(n * p, 0.1, 0.3), n, p)
  reducer <- get_reducer("meta:re_reg")
  opts <- validate_reducer_options(reducer$options_schema %||% list(), list())
  full <- reducer$fun(beta, var, X, NULL, NULL, NULL, NULL, NULL, opts)
  exact <- .exact_random_deletion_block(
    beta, var, X, estimands, full, c(2L, 5L), reducer, list(),
    examination_control()$tolerance
  )
  z_of <- function(fit, rows) {
    vapply(seq_len(p), function(b) {
      w <- 1 / (var[rows, b] + fit$tau2[b])
      A <- solve(crossprod(X[rows, ] * sqrt(w)))
      drop(estimands %*% fit$coef[, b]) / sqrt(drop(estimands %*% A %*% t(estimands)))
    }, numeric(1))
  }
  full_z <- z_of(full, seq_len(n))
  for (position in 1:2) {
    i <- c(2L, 5L)[position]
    keep <- seq_len(n) != i
    deleted <- reducer$fun(beta[keep, ], var[keep, ], X[keep, ], NULL, NULL,
                           NULL, NULL, NULL, opts)
    deleted_z <- vapply(seq_len(p), function(b) {
      w <- 1 / (var[keep, b] + deleted$tau2[b])
      A <- solve(crossprod(X[keep, ] * sqrt(w)))
      drop(estimands %*% deleted$coef[, b]) /
        sqrt(drop(estimands %*% A %*% t(estimands)))
    }, numeric(1))
    expect_equal(exact$delta_stat[position, 1, ], full_z - deleted_z, tolerance = 1e-10)
  }
})

test_that("dense voxel storage with a mask bitmap uses the full grid", {
  bitmap <- array(c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), c(2, 2, 2))
  dense <- space_voxel(c(2L, 2L, 2L), diag(4), mask_bitmap = bitmap, storage = "dense")
  packed <- space_voxel(c(2L, 2L, 2L), diag(4), mask_bitmap = bitmap, storage = "packed")
  expect_true(.is_image_assay(array(0, c(8, 2, 1)), dense))
  expect_false(.is_image_assay(array(0, c(6, 2, 1)), dense))
  expect_true(.is_image_assay(array(0, c(6, 2, 1)), packed))

  set.seed(3)
  beta <- array(rnorm(8 * 6), c(8, 6, 1))
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    dense, paste0("s", 1:6), "task"
  )
  plan <- as_plan(g)
  expect_identical(.plan_sample_labels(plan), as.character(1:8))

  exam <- examine_group(
    reduce(plan, method = "meta:fe"),
    retain = "s1",
    control = .review_quiet_control()
  )
  map <- .examination_subject_map_data(
    exam, "s1", requested_assays = "observed", slice = 2L
  )
  expect_equal(map$value, as.vector(beta[, 1, 1])[5:8])
})

test_that("as_fmri_frame keeps dense voxel GDS support on the full grid", {
  skip_if_not_installed("fmridataset")
  bitmap <- array(c(TRUE, FALSE, TRUE, TRUE), c(2, 2, 1))
  g <- new_gds(
    list(beta = array(seq_len(4 * 2), c(4, 2, 1)), var = array(1, c(4, 2, 1))),
    space_voxel(c(2L, 2L, 1L), diag(4), mask_bitmap = bitmap, storage = "dense"),
    c("s1", "s2"), "task"
  )
  frame <- fmridataset::as_fmri_frame(g)
  expect_identical(fmridataset::space(frame)$support, 1:4)
  expect_identical(ncol(frame), 4L)
})

test_that("split stability rescales split energy so stable structure scores one", {
  control <- examination_control(
    geometry = list(rank = 2L, oversample = 0L, stability_replicates = 2L)
  )
  subjects <- paste0("s", 1:4)
  u <- c(2, -1, 0.5, -1.5)
  n_feature <- 40L
  E <- outer(u, rep(c(1, -1), length.out = n_feature))
  basis <- list(
    Q = qr.Q(qr(cbind(u, c(1, 1, -1, -1)))),
    sketch_rank = 2L,
    requested_rank = 2L,
    pass1_energy = sum(E^2),
    status = "available"
  )
  state <- .initialize_geometry_pass2(
    basis,
    list(n_sample = n_feature, subjects = subjects, contrasts = "task",
         estimands = "pooled_effect"),
    control
  )
  diagnostic <- list(
    predictive_resid = E,
    surprise_eligible = matrix(TRUE, 4, n_feature)
  )
  state <- .accumulate_geometry_pass2(
    state, diagnostic, paste0("v", seq_len(n_feature)), "task",
    rep(1:2, length.out = n_feature), control
  )
  embedding <- .finalize_residual_geometry(state, control)
  expect_equal(embedding$coordinates$stability, rep(1, 4), tolerance = 1e-8)
})

test_that("balance_contrasts equalizes contrasts with different eligible counts", {
  set.seed(9)
  n_sample <- 30L
  subjects <- paste0("s", 1:8)
  beta <- array(rnorm(n_sample * 8 * 2, sd = 0.3), c(n_sample, 8, 2))
  beta[, 1:2, 1] <- beta[, 1:2, 1] + 2
  beta[, 5:6, 2] <- beta[, 5:6, 2] - 2
  beta[6:n_sample, , 2] <- NA_real_
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_sample_labels(paste0("v", seq_len(n_sample))),
    subjects, c("a", "b")
  )
  plan <- reduce(as_plan(g), method = "meta:fe")
  run <- function(balance) {
    examine_group(plan, control = examination_control(
      retain_n = 0L,
      geometry = list(rank = 4L, oversample = 2L, balance_contrasts = balance)
    ))$embedding
  }
  balanced <- run(TRUE)
  unbalanced <- run(FALSE)
  expect_false(isTRUE(all.equal(
    balanced$explained_energy, unbalanced$explained_energy
  )))
  scale <- .geometry_contrast_scale(
    list(contrasts = c("a", "b"), subjects = subjects, n_sample = n_sample,
         geometry_eligible_contrast = c(30 * 8, 5 * 8)),
    examination_control()
  )
  expect_equal(scale, 1 / sqrt(c(30, 5)))
})

test_that("fixed-effect diagnostics honor min_subjects and the eps floor", {
  g <- .group_examination_fixture()
  beta <- assay(g, "beta")
  var <- assay(g, "var")
  beta[1, 5:10, 1] <- NA_real_
  var[2, 1, 1] <- 1e-14
  g <- new_gds(list(beta = beta, var = var), space(g), subjects(g), contrasts(g))
  plan <- reduce(as_plan(g), method = "meta:fe", options = list(min_subjects = 5L))
  exam <- suppressWarnings(examine_group(plan, control = .review_quiet_control()))
  reference <- suppressWarnings(compute(plan))
  stat <- assay(exam$group_maps, "stat:pooled_effect")[, 1, 1]
  expect_true(is.na(stat[1]))
  expect_equal(unname(stat), unname(assay(reference, "z")[, 1, 1]), tolerance = 1e-10)
})

test_that("exact refit influence counts only features the subject contributed to", {
  set.seed(42)
  beta <- array(rnorm(20 * 8, 0.8, 0.6), c(20, 8, 1))
  var <- array(runif(20 * 8, 0.08, 0.35), dim(beta))
  beta[1:10, 2, 1] <- NA_real_
  g <- new_gds(
    list(beta = beta, var = var),
    space_sample_labels(paste0("v", seq_len(20))),
    paste0("s", seq_len(8)),
    "task"
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:re"),
    retain = "s2",
    control = .review_quiet_control(block_size = 7L)
  )
  refit <- exam$estimand_data[
    exam$estimand_data$ranking_stage == "selected_refit" &
      exam$estimand_data$subject == "s2",
  ]
  expect_identical(refit$eligible_n, 10)
  delta <- assay(exam$subject_maps, "delta_stat_exact:pooled_effect")[, "s2", 1]
  expect_true(all(delta[1:10] == 0))
  expect_equal(
    refit$influence_energy,
    sqrt(mean(pmin(delta[11:20]^2, 8^2)))
  )
})

test_that("report covariates show factor labels and heatmap titles are attributes", {
  g <- .group_examination_fixture()
  ids <- subjects(g)
  g <- with_col_data(
    g,
    data.frame(
      group = factor(rep(c("control", "patient"), each = 5)),
      row.names = ids
    )
  )
  exam <- examine_group(
    as_plan(g) |> reduce(method = "ols:voxelwise", formula = ~ group),
    estimands = "grouppatient",
    retain = "s8",
    control = examination_control(retain_n = 0L)
  )
  file <- tempfile(fileext = ".html")
  on.exit(unlink(file), add = TRUE)
  write_report(exam, file, subjects = "s8")
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  expect_match(html, "<td>patient</td>", fixed = TRUE)
  expect_false(grepl("<td>2</td>", html, fixed = TRUE) &&
                 !grepl("<td>patient</td>", html, fixed = TRUE))
  expect_false(grepl("<title>[^<]*</title></td>", html))
  expect_match(html, "class=\"heat\" style=\"[^\"]*\" title=\"")
  expect_error(write_report(exam, file, overwrite = NA), "TRUE or FALSE")
})

test_that("one_sample subset drops subjects with missing covariates", {
  g <- .group_examination_fixture()
  ids <- subjects(g)
  g <- with_col_data(
    g,
    data.frame(age = c(20, NA, 35, 40, NA, 50, 60, 22, 30, 41), row.names = ids)
  )
  plan <- one_sample(g, subset = age > 25)
  expect_identical(
    .plan_subjects_for_model_matrix(plan),
    ids[c(3, 4, 6, 7, 9, 10)]
  )
  expect_s3_class(compute(plan), "gds")
})
