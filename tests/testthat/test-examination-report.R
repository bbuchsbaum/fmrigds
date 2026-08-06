test_that("review reasons are deterministic and component-specific", {
  plan <- reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  small <- examine_group(plan, control = examination_control(block_size = 3L))
  large <- examine_group(plan, control = examination_control(block_size = 17L))
  small_reason <- vapply(
    small$subject_data$subject,
    function(subject) fmrigds:::.examination_review_reason(small, subject),
    character(1)
  )
  large_reason <- vapply(
    large$subject_data$subject,
    function(subject) fmrigds:::.examination_review_reason(large, subject),
    character(1)
  )
  expect_identical(small_reason, large_reason)
  expect_match(small_reason[small$subject_data$subject == "s8"], "group influence")
  expect_match(small_reason[small$subject_data$subject == "s8"], "mode exact")
  expect_match(small_reason[small$subject_data$subject == "s10"], "model surprise")
  expect_match(small_reason[small$subject_data$subject == "s10"], "mode exact")
  expect_match(small_reason[small$subject_data$subject == "s1"], "No absolute")
})

test_that("self-contained report links cohort views to retained subjects", {
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  )
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  result <- write_report(exam, path)
  expect_identical(result, normalizePath(path))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "<!doctype html>", fixed = TRUE)
  expect_match(html, "Group examination", fixed = TRUE)
  expect_match(html, "Action plane", fixed = TRUE)
  expect_match(html, "Residual geometry", fixed = TRUE)
  expect_match(html, "Subject x contrast x estimand", fixed = TRUE)
  expect_match(html, "predictive_residual", fixed = TRUE)
  expect_match(html, "delta_stat", fixed = TRUE)
  for (subject in exam$config$retained_subjects) {
    id <- fmrigds:::.html_id(subject)
    expect_match(html, paste0("id=\"subject-", id, "\""), fixed = TRUE)
    expect_match(html, paste0("href=\"#subject-", id, "\""), fixed = TRUE)
  }
  expect_false(grepl("<script[^>]+src=", html))
  expect_false(grepl("<link[^>]+href=", html))
  expect_false(grepl("<img[^>]+src=", html))
  expect_false(grepl("https?://", html))
  expect_error(write_report(exam, path), "already exists")
})

test_that("no-review report uses restrained language", {
  beta <- array(1, c(24, 8, 1))
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_sample_labels(paste0("v", seq_len(24))),
    paste0("s", seq_len(8)),
    "task"
  )
  exam <- examine_group(reduce(as_plan(g), method = "meta:fe"))
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  write_report(exam, path)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "No subject met the absolute review criteria", fixed = TRUE)
  expect_match(html, "does not classify subjects", fixed = TRUE)
  expect_false(grepl("bad subject|outlier probability|automatically exclude", html, ignore.case = TRUE))
})

test_that("ranked localization uses supplied atlas labels", {
  g <- .group_examination_fixture()
  g$row_data <- data.frame(
    region = rep(c("visual", "temporal", "frontal", "motor"), each = 10),
    row.names = paste0("feature-", seq_len(40))
  )
  exam <- examine_group(reduce(as_plan(g), method = "meta:fe"))
  subject <- exam$config$retained_subjects[1L]
  ranked <- fmrigds:::.examination_ranked_regions(exam, subject)
  expect_true("region" %in% names(ranked))
  expect_true(all(ranked$region %in% unique(g$row_data$region)))
  expect_lte(nrow(ranked), 10L)
})

test_that("report validates retained subject scope and output directory", {
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  )
  expect_error(
    write_report(exam, tempfile(fileext = ".html"), subjects = "not-retained"),
    "retained maps"
  )
  missing <- file.path(tempfile("missing-report-dir-"), "report.html")
  expect_error(write_report(exam, missing), "does not exist")
  expect_error(write_report(exam, tempfile(fileext = ".html"), title = NA_character_), "title")
})

test_that("report keeps threshold-dependent conclusions structurally separate", {
  exam <- examine_group(
    as_plan(.group_examination_fixture()) |>
      reduce(method = "meta:fe") |>
      posthoc("fdr:bh"),
    control = examination_control(retain_n = 2L)
  )
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  write_report(exam, path)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Threshold-dependent conclusion sensitivity", fixed = TRUE)
  expect_match(html, "continuous influence", fixed = TRUE)
  expect_match(html, "full_significant_n", fixed = TRUE)
  expect_match(html, "transition_count", fixed = TRUE)
  expect_s3_class(exam$conclusion$full_maps, "gds")
  expect_s3_class(exam$conclusion$deleted_maps, "gds")
})

test_that("voxel-space selected views use non-interpolated axial rasters", {
  skip_if_not_installed("ggplot2")
  set.seed(19)
  voxel_dim <- c(4L, 4L, 3L)
  n_sample <- prod(voxel_dim)
  n_subject <- 8L
  beta <- array(rnorm(n_sample * n_subject, sd = 0.2), c(n_sample, n_subject, 1L))
  beta[, 8L, 1L] <- beta[, 8L, 1L] + rep(c(-1, 1), length.out = n_sample)
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_voxel(voxel_dim, diag(4), storage = "dense"),
    paste0("s", seq_len(n_subject)),
    "task"
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    retain = "s8"
  )
  data <- fmrigds:::.examination_subject_map_data(exam, "s8", slice = 2L)
  expect_identical(attr(data, "display"), "voxel")
  expect_identical(unique(data$slice), 2L)
  expect_equal(nrow(data), voxel_dim[1L] * voxel_dim[2L] * 4L)
  figure <- plot(exam, subject = "s8", slice = 2L)
  expect_s3_class(figure, "ggplot")
  expect_warning(ggplot2::ggplot_build(figure), NA)

  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  write_report(exam, path, subjects = "s8")
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "<rect", fixed = TRUE)
  expect_false(grepl("interpolat", html, ignore.case = TRUE))
})

test_that("insufficient presentation state is explicit and quiet", {
  skip_if_not_installed("ggplot2")
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  )
  exam$subject_data$review_status <- "insufficient"
  exam$subject_data$review_source <- NA_character_
  exam$subject_data$surprise_score <- NA_real_
  exam$subject_data$influence_score <- NA_real_
  exam$contrast_data$surprise_energy <- NA_real_
  exam$contrast_data$surprise_status <- "insufficient_samples"
  exam$estimand_data$influence_energy <- NA_real_
  exam$estimand_data$status <- "nonestimable"
  expect_warning(figure <- plot(exam), NA)
  expect_s3_class(figure, "ggplot")
  expect_warning(ggplot2::ggplot_build(figure), NA)
  expect_match(
    fmrigds:::.examination_review_reason(exam, exam$subject_data$subject[1L]),
    "Insufficient eligible data"
  )
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  write_report(exam, path)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Action-plane metrics are unavailable", fixed = TRUE)
})
