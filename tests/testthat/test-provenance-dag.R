.provenance_test_metadata <- function(n = 1L) {
  files <- lapply(seq_len(n), function(i) {
    path <- tempfile(paste0("prov-source-", i, "-"))
    writeLines(paste("source", i), path)
    path
  })
  entities <- .source_entities_from_files(lapply(seq_along(files), function(i) list(
    path = files[[i]], role = "input", ordinal = i
  )))
  metadata <- .metadata_with_source_entities(gds_metadata(), entities)
  attr(metadata, "test_files") <- unlist(files)
  metadata
}

test_that("plan activities form a closed source-rooted DAG", {
  beta <- array(seq_len(16), c(2, 4, 2))
  se <- array(0.5, c(2, 4, 2))
  source <- new_gds(
    list(beta = beta, se = se),
    space_sample_labels(c("x", "y")),
    paste0("s", 1:4),
    c("A", "B")
  )
  plan <- derive(source, "var") |>
    reduce(method = "meta:re")
  fit <- compute(plan)
  provenance <- metadata(fit)$provenance
  ids <- vapply(plan$nodes, `[[`, character(1L), "node_id")

  expect_length(provenance$graph, 2L)
  expect_identical(vapply(provenance$graph, `[[`, character(1L), "id"), ids)
  expect_identical(provenance$graph[[2L]]$inputs, ids[[1L]])
  expect_identical(
    provenance$graph[[1L]]$inputs,
    vapply(provenance$entities, `[[`, character(1L), "id")
  )
  expect_true(validate_provenance_graph(fit))
  expect_identical(metadata(fit)$reduction_receipts[[ids[[2L]]]]$nInputSubjects, 4L)
})

test_that("activity digest includes parents while occurrence id stays distinct", {
  metadata <- .provenance_test_metadata(2L)
  on.exit(unlink(attr(metadata, "test_files")), add = TRUE)
  entity_ids <- vapply(metadata$provenance$entities, `[[`, character(1L), "id")
  left <- add_provenance_node(metadata, "derive", list(what = "z"), entity_ids[[1L]], id = "left")
  right <- add_provenance_node(metadata, "derive", list(what = "z"), entity_ids[[2L]], id = "right")
  expect_false(identical(
    left$provenance$graph[[1L]]$digest$value,
    right$provenance$graph[[1L]]$digest$value
  ))

  repeated <- add_provenance_node(left, "derive", list(what = "z"), "left", id = "repeat")
  expect_false(identical(repeated$provenance$graph[[1L]]$id, repeated$provenance$graph[[2L]]$id))

  timestamp <- as.POSIXct("2026-08-24 12:00:00", tz = "UTC")
  occurrence_a <- provenance_node("derive", list(what = "z"), "left", timestamp = timestamp)
  occurrence_b <- provenance_node("derive", list(what = "z"), "left", timestamp = timestamp)
  expect_false(identical(occurrence_a$id, occurrence_b$id))
  expect_identical(occurrence_a$digest, occurrence_b$digest)
})

test_that("independent non-file lineages keep distinct source occurrences", {
  left <- .ensure_nonfile_source(gds_metadata(), "memory", "ingest-source")
  right <- .ensure_nonfile_source(gds_metadata(), "memory", "ingest-source")
  left_id <- left$provenance$entities[[1L]]$id
  right_id <- right$provenance$entities[[1L]]$id

  expect_false(identical(left_id, right_id))
  merged <- .merge_provenance_lineages(left, right)
  expect_length(merged$provenance$entities, 2L)
  expect_identical(length(unique(c(merged$parents_a, merged$parents_b))), 2L)
})

test_that("graph validator rejects missing, duplicate, self, and cyclic endpoints", {
  metadata <- .provenance_test_metadata(1L)
  on.exit(unlink(attr(metadata, "test_files")), add = TRUE)
  root <- metadata$provenance$entities[[1L]]$id
  valid <- add_provenance_node(metadata, "derive", list(what = "z"), root, id = "a")
  valid <- add_provenance_node(valid, "reduce", list(method = "meta:re"), "a", id = "b")

  missing <- valid$provenance
  missing$graph[[2L]]$inputs <- "absent"
  expect_match(validate_provenance_graph(missing, error = FALSE), "missing input", all = FALSE)

  duplicate <- valid$provenance
  duplicate$graph[[2L]]$id <- "a"
  expect_match(validate_provenance_graph(duplicate, error = FALSE), "globally unique", all = FALSE)

  self <- valid$provenance
  self$graph[[2L]]$inputs <- "b"
  expect_match(validate_provenance_graph(self, error = FALSE), "self-edge", all = FALSE)

  cycle <- valid$provenance
  cycle$graph[[1L]]$inputs <- "b"
  cycle$graph[[2L]]$inputs <- "a"
  expect_match(validate_provenance_graph(cycle, error = FALSE), "cycle", all = FALSE)

  wrong_heads <- valid$provenance
  wrong_heads$heads <- "a"
  expect_match(
    validate_provenance_graph(wrong_heads, error = FALSE),
    "declared provenance heads",
    all = FALSE
  )
})

test_that("legacy provenance remains readable and explicitly incomplete", {
  legacy <- gds_metadata(provenance = list(
    graph = list(list(op = "derive", params = list(what = "z"), inputs = list())),
    log = "legacy",
    digest = "old"
  ))
  expect_identical(legacy$provenance$status, "legacy-incomplete")
  expect_identical(legacy$provenance$graphReceipt$status, "unavailable")
  expect_match(
    validate_provenance_graph(legacy, error = FALSE),
    "legacy activities",
    all = FALSE
  )
  expect_error(
    add_provenance_node(legacy, "reduce", list(method = "meta:re"), inputs = "old"),
    "Cannot append to legacy provenance"
  )
})

test_that("serialized plans preserve occurrence ids and edge structure", {
  path <- tempfile(fileext = ".csv")
  json <- tempfile(fileext = ".json")
  on.exit(unlink(c(path, json)), add = TRUE)
  writeLines(c(
    "sample,subject,contrast,beta,se",
    "r1,s1,c1,1,0.5",
    "r1,s2,c1,2,0.5",
    "r1,s3,c1,3,0.5"
  ), path)
  plan <- derive(gds(path, effect_cols = list(beta = "beta", se = "se")), "var") |>
    reduce(method = "meta:re")
  save_plan(plan, json)
  restored <- load_plan(json)
  original_ids <- vapply(plan$nodes, `[[`, character(1L), "node_id")
  restored_ids <- vapply(restored$nodes, `[[`, character(1L), "node_id")
  expect_identical(restored_ids, original_ids)

  fit <- compute(restored)
  graph <- metadata(fit)$provenance$graph
  expect_identical(vapply(graph, `[[`, character(1L), "id"), original_ids)
  expect_identical(graph[[2L]]$inputs, original_ids[[1L]])
})

test_that("add_provenance_node requires explicit parents after the first activity", {
  metadata <- .provenance_test_metadata(1L)
  on.exit(unlink(attr(metadata, "test_files")), add = TRUE)
  root <- metadata$provenance$entities[[1L]]$id
  metadata <- add_provenance_node(metadata, "derive", list(what = "z"), root)
  expect_error(
    add_provenance_node(metadata, "reduce", list(method = "meta:re")),
    "`inputs` is required"
  )
})

test_that("register_map extends the current activity head", {
  beta <- array(seq_len(8), c(4, 2, 1))
  source <- new_gds(
    list(beta = beta, se = array(0.5, dim(beta))),
    space_parcels(paste0("roi", seq_len(4))),
    c("s1", "s2"),
    "task"
  )
  fit <- compute(derive(as_plan(source), "var"))
  previous_head <- metadata(fit)$provenance$heads
  family <- make_linear_family(
    "identity",
    space(fit),
    space(fit),
    list(s1 = diag(4), s2 = diag(4))
  )
  registered <- register_map(fit, family)
  graph <- metadata(registered)$provenance$graph

  expect_length(graph, 2L)
  expect_identical(graph[[2L]]$inputs, previous_head)
  expect_identical(metadata(registered)$provenance$heads, graph[[2L]]$id)
  expect_true(validate_provenance_graph(registered))
})

test_that("experimental cancellation retains both addressable input lineages", {
  n_subjects <- 8L
  make_split <- function(offset) {
    pattern <- rep(c(-1, 1), length.out = n_subjects) + offset
    beta <- array(rep(pattern, each = 2L), c(2, n_subjects, 1))
    new_gds(
      list(beta = beta, var = array(0.1, dim(beta))),
      space_voxel(c(2, 1, 1), diag(4)),
      paste0("s", seq_len(n_subjects)),
      "task"
    )
  }
  result <- experimental_cancellation(
    make_split(0),
    make_split(0.01),
    delta = 0.1,
    equivalence = 0.1,
    candidate_threshold = 1,
    min_subjects = 3L
  )
  provenance <- metadata(result)$provenance
  activity <- provenance$graph[[length(provenance$graph)]]

  expect_identical(activity$op, "experimental_cancellation")
  expect_length(activity$inputs, 2L)
  expect_true(all(activity$inputs %in% vapply(
    provenance$entities,
    `[[`,
    character(1L),
    "id"
  )))
  expect_true(validate_provenance_graph(result))
})
