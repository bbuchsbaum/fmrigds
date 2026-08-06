test_that("action plane exposes the surprise-influence quadrants", {
  skip_if_not_installed("ggplot2")
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:fe"),
    control = examination_control(block_size = 7L)
  )
  data <- fmrigds:::.examination_action_data(exam)
  low_surprise_high_influence <- data[data$subject == "s8", , drop = FALSE]
  high_surprise_low_influence <- data[data$subject == "s10", , drop = FALSE]
  expect_lt(
    low_surprise_high_influence$surprise_energy,
    exam$config$control$review$surprise$energy_threshold
  )
  expect_gt(
    low_surprise_high_influence$influence_energy,
    exam$config$control$review$influence$energy_threshold
  )
  expect_gt(
    high_surprise_low_influence$surprise_energy,
    exam$config$control$review$surprise$energy_threshold
  )
  expect_lt(
    high_surprise_low_influence$influence_energy,
    exam$config$control$review$influence$energy_threshold
  )

  figure <- plot(exam)
  expect_s3_class(figure, "ggplot")
  expect_warning(ggplot2::ggplot_build(figure), NA)
  expect_identical(figure$labels$x, "Model surprise (capped residual energy)")
  expect_identical(
    figure$labels$y,
    "Group influence (capped deletion-statistic energy)"
  )
  expect_match(figure$labels$caption, "Influence mode: exact")
  expect_match(figure$labels$caption, "eligible features")
  expect_equal(length(figure$layers), 4L)
})

test_that("site and group columns survive into plotting data", {
  skip_if_not_installed("ggplot2")
  g <- .group_examination_fixture()
  ids <- subjects(g)
  g <- with_col_data(
    g,
    data.frame(
      group = factor(rep(c("control", "patient"), each = 5)),
      site = factor(rep(c("A", "B"), times = 5)),
      age = seq(20, 65, length.out = 10),
      row.names = ids
    )
  )
  exam <- examine_group(
    reduce(as_plan(g), method = "meta:fe_reg", formula = ~ group + age),
    estimands = "grouppatient"
  )
  expect_true(all(c("group", "site", "age") %in% names(exam$subject_data)))
  action <- fmrigds:::.examination_action_data(exam, group_by = "site")
  expect_identical(levels(action$display_group), c("A", "B"))
  expect_s3_class(plot(exam, group_by = "site"), "ggplot")
  expect_s3_class(plot(exam, type = "embedding", group_by = "group"), "ggplot")
})

test_that("heatmap, embedding, and selected-subject views are static ggplots", {
  skip_if_not_installed("ggplot2")
  exam <- examine_group(
    reduce(as_plan(.group_examination_fixture()), method = "meta:fe")
  )
  expect_s3_class(plot(exam, type = "heatmap"), "ggplot")
  embedding <- plot(exam, type = "embedding")
  expect_s3_class(embedding, "ggplot")
  expect_match(embedding$labels$caption, "Captured residual energy")

  retained <- exam$config$retained_subjects[1L]
  subject_view <- plot(exam, subject = retained)
  expect_s3_class(subject_view, "ggplot")
  expect_warning(ggplot2::ggplot_build(subject_view), NA)
  expect_match(subject_view$labels$caption, "map modes")
  expect_error(plot(exam, subject = "not-retained"), "retained subject")
})

test_that("degenerate summaries produce explicit quiet plot states", {
  skip_if_not_installed("ggplot2")
  beta <- array(1, c(30, 8, 1))
  g <- new_gds(
    list(beta = beta, var = array(0.1, dim(beta))),
    space_sample_labels(paste0("v", seq_len(30))),
    paste0("s", seq_len(8)),
    "task"
  )
  exam <- examine_group(reduce(as_plan(g), method = "meta:fe"))
  expect_true(all(exam$subject_data$review_status == "none"))
  expect_warning(action <- plot(exam), NA)
  expect_warning(heatmap <- plot(exam, type = "heatmap"), NA)
  expect_warning(embedding <- plot(exam, type = "embedding"), NA)
  expect_s3_class(action, "ggplot")
  expect_s3_class(heatmap, "ggplot")
  expect_s3_class(embedding, "ggplot")
  expect_false(any(!is.na(action$data$review_label)))
})

test_that("multi-contrast action views facet without changing availability", {
  skip_if_not_installed("ggplot2")
  g <- .group_examination_fixture()
  g$assays <- lapply(g$assays, function(value) {
    array(
      c(value, value * 0.5),
      c(dim(value)[1L], dim(value)[2L], 2L),
      dimnames = list(NULL, dimnames(value)[[2L]], c("task-a", "task-b"))
    )
  })
  g$contrasts <- c("task-a", "task-b")
  g$dims[["contrast"]] <- 2L
  exam <- examine_group(reduce(as_plan(g), method = "meta:fe"))
  data <- fmrigds:::.examination_action_data(exam)
  expect_identical(sort(unique(data$contrast)), c("task-a", "task-b"))
  expect_s3_class(plot(exam), "ggplot")
})
