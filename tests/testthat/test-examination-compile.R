test_that("plan node identities are stable, unique, and serialization-safe", {
  source_path <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      sample = rep(c("a", "b"), each = 2),
      subject = rep(c("s1", "s2"), 2),
      contrast = "c1",
      beta = 1:4,
      var = 1
    ),
    source_path,
    row.names = FALSE
  )
  on.exit(unlink(source_path), add = TRUE)
  base <- gds(source_path)
  one <- base |> derive("se")
  first_id <- one$nodes[[1L]]$node_id
  two <- one |> reduce(method = "fixed")

  expect_match(first_id, "^node-[0-9]{4}-[0-9a-f]+$")
  expect_identical(two$nodes[[1L]]$node_id, first_id)
  expect_false(identical(two$nodes[[1L]]$node_id, two$nodes[[2L]]$node_id))
  expect_identical((base |> derive("se"))$nodes[[1L]]$node_id, first_id)

  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  save_plan(two, path)
  restored <- load_plan(path)
  expect_identical(
    vapply(restored$nodes, `[[`, character(1), "node_id"),
    vapply(two$nodes, `[[`, character(1), "node_id")
  )

  changed_identity <- two
  changed_identity$nodes[[1L]]$node_id <- "node-9999-deadbeef"
  expect_identical(digest_plan(changed_identity), digest_plan(two))
})

test_that("examination compiler splits grammar and discards writes", {
  g <- new_gds(
    assays = list(beta = array(1:12, c(3, 4, 1)), var = array(1, c(3, 4, 1))),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = paste0("s", 1:4),
    contrasts = "c1"
  )
  path_before <- tempfile(fileext = ".csv")
  path_after <- tempfile(fileext = ".csv")
  plan <- as_plan(g) |>
    subset(subject = c("s4", "s2")) |>
    write_out(path_before, format = "csv") |>
    derive("se") |>
    reduce(method = "fixed") |>
    posthoc("fdr:bh") |>
    write_out(path_after, format = "csv")

  compiled <- fmrigds:::compile_examination_plan(plan)
  expect_s3_class(compiled, "gds_examination_plan")
  expect_equal(vapply(compiled$prefix, `[[`, character(1), "op"),
               c("subset_axis", "derive"))
  expect_identical(compiled$reducer$op, "reduce")
  expect_equal(vapply(compiled$conclusion_tail, `[[`, character(1), "op"), "posthoc")
  expect_length(compiled$discarded_writes, 2L)
  expect_identical(compiled$scan_strategy, "direct")
  expect_identical(compiled$axis_selection$subject, c(4L, 2L))
  expect_false(file.exists(path_before))
  expect_false(file.exists(path_after))
  expect_true(all(nzchar(compiled$node_ids)))

  again <- fmrigds:::compile_examination_plan(plan)
  expect_identical(
    compiled[c("node_ids", "scan_strategy", "axis_selection", "model_spec")],
    again[c("node_ids", "scan_strategy", "axis_selection", "model_spec")]
  )
})

test_that("examination compiler composes leading subset positions", {
  g <- new_gds(
    assays = list(beta = array(1:30, c(5, 3, 2)), var = array(1, c(5, 3, 2))),
    space = space_sample_labels(letters[1:5]),
    subjects = c("s1", "s2", "s3"),
    contrasts = c("c1", "c2")
  )
  compiled <- as_plan(g) |>
    subset(sample = 2:5, subject = c("s3", "s1")) |>
    derive("se") |>
    subset(sample = 2:3, subject = "s1", contrast = "c2") |>
    reduce(method = "fixed") |>
    fmrigds:::compile_examination_plan()

  expect_identical(compiled$axis_selection$sample, c(3L, 4L))
  expect_identical(compiled$axis_selection$subject, 1L)
  expect_identical(compiled$axis_selection$contrast, 2L)
  expect_equal(vapply(compiled$execution_prefix, `[[`, character(1), "op"), "derive")
})

test_that("examination compiler stages materializing and unsafe prefixes", {
  g <- new_gds(
    assays = list(beta = array(1:12, c(3, 4, 1)), var = array(1, c(3, 4, 1))),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = paste0("s", 1:4),
    contrasts = "c1"
  )
  mapped <- as_plan(g) |>
    map_to(space_sample_labels("all"), matrix(1, 1, 3)) |>
    reduce(method = "fixed") |>
    fmrigds:::compile_examination_plan()
  expect_identical(mapped$scan_strategy, "stage")
  expect_match(paste(mapped$scan_reasons, collapse = " "), "map")

  masked_then_subset <- as_plan(g) |>
    mask(MaskPolicy(rule = "intersection")) |>
    subset(sample = 1L) |>
    reduce(method = "fixed") |>
    fmrigds:::compile_examination_plan()
  expect_identical(masked_then_subset$scan_strategy, "stage")
  expect_match(paste(masked_then_subset$scan_reasons, collapse = " "), "subset")
})

test_that("examination compiler rejects ambiguous plan grammar", {
  g <- new_gds(
    assays = list(beta = array(1:8, c(2, 4, 1)), var = array(1, c(2, 4, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = paste0("s", 1:4),
    contrasts = "c1"
  )

  expect_error(
    as_plan(g) |> reduce(method = "fixed") |> reduce(method = "fixed") |>
      fmrigds:::compile_examination_plan(),
    "exactly one reducer"
  )
  expect_error(
    as_plan(g) |> posthoc("fdr:bh") |> reduce(method = "fixed") |>
      fmrigds:::compile_examination_plan(),
    "post-hoc.*before"
  )
  expect_error(
    as_plan(g) |> reduce(method = "fixed") |> derive("se") |>
      fmrigds:::compile_examination_plan(),
    "after the reducer"
  )
})

test_that("adapter scan capabilities default conservatively", {
  expect_true(fmrigds:::.adapter_capabilities("memory")$sample_blocks)
  expect_true(fmrigds:::.adapter_capabilities("h5")$persistent_handle)
  expect_true(fmrigds:::.adapter_capabilities("nifti")$requires_staging)

  register_adapter(
    "legacy_capability_fixture",
    detect = function(source) 0,
    open = function(source, ...) list(),
    probe = function(handle, ...) NULL,
    read = function(handle, assays, block = NULL, ...) list(),
    close = function(handle) NULL
  )
  caps <- fmrigds:::.adapter_capabilities("legacy_capability_fixture")
  expect_false(caps$sample_blocks)
  expect_false(caps$persistent_handle)
  expect_false(caps$cheap_revisit)
})
