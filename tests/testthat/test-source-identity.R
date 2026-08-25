test_that("file identity is byte-based and path-independent", {
  first <- tempfile("identity-a-")
  second <- tempfile("identity-b-")
  on.exit(unlink(c(first, second)), add = TRUE)
  writeBin(charToRaw("same bytes"), first)
  file.copy(first, second)

  a <- .capture_file_identity(first)
  b <- .capture_file_identity(second)
  expect_identical(a$identityStatus, "verified")
  expect_identical(a$digest, b$digest)
  expect_identical(a$byteSize, b$byteSize)

  entities_a <- .source_entities_from_files(list(list(path = first, role = "input")))
  entities_b <- .source_entities_from_files(list(list(path = second, role = "input")))
  expect_identical(entities_a[[1L]]$digest, entities_b[[1L]]$digest)
  expect_false(identical(entities_a[[1L]]$id, entities_b[[1L]]$id))

  writeBin(charToRaw("changed bytes"), first)
  changed <- .capture_file_identity(first)
  expect_false(identical(a$digest$value, changed$digest$value))
})

test_that("source identity states distinguish policy, absence, and mutation", {
  path <- tempfile("identity-state-")
  on.exit(unlink(path), add = TRUE)
  writeLines("one", path)

  skipped <- .capture_file_identity(path, policy = "none")
  expect_identical(skipped$identityStatus, "not-requested")
  expect_null(skipped$digest)

  missing <- .capture_file_identity(
    paste0(path, "-missing"),
    on_error = "record"
  )
  expect_identical(missing$identityStatus, "unavailable")
  expect_match(missing$reason, "missing")

  receipt <- .source_entities_from_files(list(list(
    path = path, role = "table", ordinal = 1L
  )))
  writeLines("two", path)
  mutated <- .verify_source_entities(receipt, "before-read", fail = FALSE)
  expect_length(mutated, 1L)
  expect_identical(mutated[[1L]]$identityStatus, "changed-during-read")
  expect_error(
    .verify_source_entities(receipt, "before-read"),
    "Source identity changed"
  )

  writeLines("probe version", path)
  pre_read <- .source_stat(path)
  writeLines("a different and longer version", path)
  during_probe <- .source_entities_from_files(
    list(list(
      path = path,
      role = "table",
      ordinal = 1L,
      expected_stat = pre_read
    )),
    on_error = "record"
  )
  expect_identical(during_probe[[1L]]$identityStatus, "changed-during-read")
  expect_null(during_probe[[1L]]$digest)
  expect_error(
    .source_entities_from_files(list(list(
      path = path,
      role = "table",
      ordinal = 1L,
      expected_stat = pre_read
    ))),
    "between adapter read and identity capture"
  )
})

test_that("unreadable hashing failures remain structured unavailable states", {
  path <- tempfile("identity-unreadable-")
  on.exit(unlink(path), add = TRUE)
  writeLines("secret", path)
  local_mocked_bindings(
    .hash_source_file = function(locator) stop("permission denied", call. = FALSE)
  )
  receipt <- .capture_file_identity(path, on_error = "record")
  expect_identical(receipt$identityStatus, "unavailable")
  expect_match(receipt$reason, "permission denied")
})

test_that("duplicate physical paths are hashed once without losing roles", {
  path <- tempfile("identity-role-")
  on.exit(unlink(path), add = TRUE)
  writeLines("payload", path)
  calls <- 0L
  original <- .hash_source_file
  local_mocked_bindings(
    .hash_source_file = function(locator) {
      calls <<- calls + 1L
      original(locator)
    }
  )
  entities <- .source_entities_from_files(list(
    list(path = path, role = "beta", ordinal = 1L, pair = 1L),
    list(path = path, role = "standard-error", ordinal = 2L, pair = 1L)
  ))
  expect_length(entities, 1L)
  expect_identical(calls, 1L)
  expect_equal(
    vapply(entities[[1L]]$roles, `[[`, character(1L), "role"),
    c("beta", "standard-error")
  )
})

test_that("compute verifies a successful probe without hashing again", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "r1,s1,c1,1,0.2",
    "r1,s2,c1,2,0.2"
  ), path)
  calls <- 0L
  original <- .hash_source_file
  local_mocked_bindings(
    .hash_source_file = function(locator) {
      calls <<- calls + 1L
      original(locator)
    }
  )
  plan <- gds(path)
  expect_identical(calls, 1L)
  compute(plan)
  expect_identical(calls, 1L)
})

test_that("portable source projection strips private locators and basenames", {
  dir <- tempfile("sub-01-secret-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "sub-01_beta.nii")
  writeBin(charToRaw("volume bytes"), path)
  entities <- .source_entities_from_files(list(list(
    path = path, role = "beta", ordinal = 1L, pair = 1L
  )))
  portable <- .portable_source_entities(entities)
  encoded <- .canonical_portable_json(portable)
  expect_false(grepl(normalizePath(dir, winslash = "/"), encoded, fixed = TRUE))
  expect_false(grepl("sub-01", encoded, fixed = TRUE))
  expect_null(portable[[1L]]$private)
})

test_that("NIfTI beta and SE inputs receive one receipt per physical file", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
  dir <- tempfile("source-nifti-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  beta <- file.path(dir, sprintf("sub-%02d_beta.nii", seq_len(26L)))
  se <- file.path(dir, sprintf("sub-%02d_se.nii", seq_len(26L)))
  for (i in seq_along(beta)) {
    RNifti::writeNifti(array(i, c(2, 2, 2)), beta[[i]])
    RNifti::writeNifti(array(i / 10, c(2, 2, 2)), se[[i]])
  }
  plan <- gds(nifti_source(beta, se), format = "nifti")
  entities <- plan$source$probe$metadata$provenance$entities
  expect_length(entities, 52L)
  expect_true(all(vapply(entities, function(x) x$identityStatus == "verified", logical(1L))))
  expect_equal(
    vapply(entities, function(x) x$roles[[1L]]$role, character(1L)),
    c(rep("beta", 26L), rep("standard-error", 26L))
  )
  expect_true(all(vapply(entities, function(x) x$byteSize > 0, logical(1L))))
  portable <- .canonical_portable_json(.portable_source_entities(entities))
  expect_false(grepl(normalizePath(dir, winslash = "/"), portable, fixed = TRUE))
  expect_false(grepl("sub-01", portable, fixed = TRUE))
})

test_that("HDF5 round-trip preserves lineage and connects its container receipt", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  beta <- array(seq_len(8), c(4, 2, 1))
  source <- new_gds(
    list(beta = beta, se = array(0.5, dim(beta))),
    space_voxel(c(2, 2, 1), diag(4)),
    c("s1", "s2"),
    "task"
  )
  first <- compute(derive(as_plan(source), "var"))
  write_gds_h5(first, path)

  reopened <- gds(path, format = "h5")
  restored <- reopened$source$probe$metadata$provenance
  expect_identical(
    vapply(restored$graph, `[[`, character(1L), "id"),
    vapply(metadata(first)$provenance$graph, `[[`, character(1L), "id")
  )
  expect_identical(
    vapply(restored$graph, function(node) node$digest$value, character(1L)),
    vapply(metadata(first)$provenance$graph, function(node) node$digest$value, character(1L))
  )
  container <- restored$entities[vapply(restored$entities, function(entity) {
    any(vapply(entity$roles, function(role) identical(role$role, "h5-container"), logical(1L)))
  }, logical(1L))]
  expect_length(container, 1L)
  expect_identical(container[[1L]]$identityStatus, "verified")

  second <- compute(derive(reopened, "se"))
  graph <- metadata(second)$provenance$graph
  expect_length(graph, 2L)
  expect_setequal(
    graph[[2L]]$inputs,
    c(graph[[1L]]$id, container[[1L]]$id)
  )
  expect_true(validate_provenance_graph(second))
})
