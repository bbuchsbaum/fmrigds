.write_cli_fixture <- function() {
  path <- tempfile(fileext = ".csv")
  writeLines(c(
    "sample,subject,contrast,time,beta,var",
    "ROI_1,sub-01,baseline,0,0.5,0.04",
    "ROI_1,sub-01,task,1,0.8,0.04",
    "ROI_1,sub-02,baseline,0,0.4,0.05",
    "ROI_1,sub-02,task,1,0.9,0.05",
    "ROI_2,sub-01,baseline,0,0.6,0.03",
    "ROI_2,sub-01,task,1,0.7,0.03",
    "ROI_2,sub-02,baseline,0,0.55,0.02",
    "ROI_2,sub-02,task,1,0.95,0.02"
  ), path)
  path
}

.capture_cli <- function(args) {
  paste(capture.output(fmrigds:::.cli_main(args)), collapse = "\n")
}

test_that("cli plan command builds and saves plans", {
  csv <- .write_cli_fixture()
  json <- tempfile(fileext = ".json")
  on.exit(unlink(c(csv, json)), add = TRUE)

  txt <- .capture_cli(c(
    "plan",
    "--input", csv,
    "--contrast-data-cols", "time",
    "--subset", "contrast=task",
    "--derive", "t,p",
    "--reduce", "fixed",
    "--posthoc", "fdr:bh",
    "--save-plan", json,
    "--json"
  ))

  info <- jsonlite::fromJSON(txt)
  loaded <- load_plan(json)

  expect_equal(info$adapter, "tabular")
  expect_equal(vapply(loaded$nodes, `[[`, character(1), "op"), c("subset_axis", "derive", "reduce", "posthoc"))
  expect_equal(rownames(contrast_data(loaded)), c("baseline", "task"))
  expect_equal(contrast_data(loaded)$time, c(0, 1))
})

test_that("cli preview emits JSON preview tables", {
  csv <- .write_cli_fixture()
  on.exit(unlink(csv), add = TRUE)

  txt <- .capture_cli(c(
    "preview",
    "--input", csv,
    "--subset", "contrast=task",
    "--reduce", "fixed",
    "--n", "1",
    "--show-assay", "beta",
    "--json"
  ))

  out <- jsonlite::fromJSON(txt)

  expect_equal(out$summary$subjects, "meta")
  expect_equal(nrow(out$data), 1L)
  expect_equal(out$data$sample, "ROI_1")
  expect_equal(out$data$subject, "meta")
  expect_equal(out$data$contrast, "task")
  expect_true(is.numeric(out$data$beta))
})

test_that("cli preview raw mode returns arrays as JSON", {
  csv <- .write_cli_fixture()
  on.exit(unlink(csv), add = TRUE)

  txt <- .capture_cli(c(
    "preview",
    "--input", csv,
    "--n", "1",
    "--show-assay", "beta",
    "--raw",
    "--json"
  ))

  out <- jsonlite::fromJSON(txt)

  expect_true("beta" %in% names(out))
  expect_equal(unname(unlist(out$beta$dim)), c(1, 2, 2))
  expect_length(out$beta$values, 4L)
})

test_that("cli list reducers exposes reducer registry in JSON", {
  txt <- .capture_cli(c("list", "reducers", "--json"))
  rows <- jsonlite::fromJSON(txt)

  expect_true("meta:fe" %in% rows$name)
  expect_true("combine:stouffer" %in% rows$name)
})
