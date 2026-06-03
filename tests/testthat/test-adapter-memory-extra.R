# Extra tests for adapter-memory.R and stubs-core.R

# ===========================================================================
# adapter-memory: detect
# ===========================================================================

test_that("memory adapter detects GDS objects with score 1.0", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  expect_equal(.memory_detect(g), 1.0)
})

test_that("memory adapter detects list-of-3D-arrays with score 0.9", {
  src <- list(beta = array(1, dim = c(2, 2, 1)))
  expect_equal(.memory_detect(src), 0.9)
})

test_that("memory adapter returns 0 for non-matching sources", {
  expect_equal(.memory_detect("a string"), 0)
  expect_equal(.memory_detect(42), 0)
  expect_equal(.memory_detect(data.frame(x = 1)), 0)
  expect_equal(.memory_detect(list()), 0)
})

# ===========================================================================
# adapter-memory: open
# ===========================================================================

test_that("memory adapter opens a GDS object", {
  g <- new_gds(
    assays = list(beta = array(1:4, dim = c(2, 2, 1)), var = array(0.1, dim = c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = "c1"
  )
  handle <- .memory_open(g)
  expect_true(!is.null(handle$gds))
})

test_that("memory adapter opens a list of arrays", {
  src <- list(beta = array(1, dim = c(2, 2, 1)))
  handle <- .memory_open(src)
  expect_true(!is.null(handle$arrays))
})

test_that("memory adapter rejects write mode", {
  expect_error(.memory_open(list(), mode = "w"), "read-only")
})

test_that("memory adapter rejects unsupported source", {
  expect_error(.memory_open("not-a-gds-or-list"), "Unsupported")
})

# ===========================================================================
# adapter-memory: probe from GDS
# ===========================================================================

test_that("memory adapter probes a GDS object correctly", {
  g <- new_gds(
    assays = list(beta = array(1:8, dim = c(2, 2, 2)), var = array(0.1, dim = c(2, 2, 2))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"), contrasts = c("c1", "c2"),
    col_data = data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
  )
  handle <- .memory_open(g)
  probe <- .memory_probe(handle)

  expect_setequal(probe$assays, c("beta", "var"))
  expect_equal(probe$dims[[1L]], 2L)
  expect_equal(probe$dims[[2L]], 2L)
  expect_equal(probe$dims[[3L]], 2L)
  expect_equal(probe$subjects, c("s1", "s2"))
  expect_equal(probe$contrasts, c("c1", "c2"))
  expect_s3_class(probe$space, "gds_space")
  expect_true(!is.null(probe$col_data))
})

# ===========================================================================
# adapter-memory: probe from raw arrays
# ===========================================================================

test_that("memory adapter probes raw arrays with synthetic labels", {
  src <- list(beta = array(1:6, dim = c(3, 2, 1)))
  handle <- .memory_open(src)
  probe <- .memory_probe(handle)

  expect_equal(probe$assays, "beta")
  expect_equal(probe$dims[[1L]], 3L)
  expect_equal(probe$dims[[2L]], 2L)
  expect_equal(probe$dims[[3L]], 1L)
  expect_equal(probe$subjects, c("1", "2"))
  expect_equal(probe$contrasts, "contrast1")
  expect_s3_class(probe$space, "space_sample_labels")
})

test_that("memory adapter probe errors on invalid handle", {
  expect_error(.memory_probe(list()), "invalid handle")
})

# ===========================================================================
# adapter-memory: read
# ===========================================================================

test_that("memory adapter reads all assays", {
  src <- list(
    beta = array(1:12, dim = c(3, 2, 2)),
    var  = array(0.1, dim = c(3, 2, 2))
  )
  handle <- .memory_open(src)
  arrs <- .memory_read(handle, assays = c("beta", "var"))

  expect_named(arrs, c("beta", "var"))
  expect_equal(dim(arrs$beta), c(3, 2, 2))
})

test_that("memory adapter reads with block subsetting", {
  src <- list(beta = array(1:12, dim = c(3, 2, 2)))
  handle <- .memory_open(src)

  # Integer block
  arrs <- .memory_read(handle, assays = "beta",
                       block = list(sample = c(1L, 3L), subject = 1L))
  expect_equal(dim(arrs$beta), c(2, 1, 2))

  # Logical block
  arrs <- .memory_read(handle, assays = "beta",
                       block = list(contrast = c(TRUE, FALSE)))
  expect_equal(dim(arrs$beta), c(3, 2, 1))
})

test_that("memory adapter read errors on missing assay", {
  src <- list(beta = array(1, dim = c(2, 2, 1)))
  handle <- .memory_open(src)
  expect_error(.memory_read(handle, assays = "nonexistent"), "Assay not found")
})

# ===========================================================================
# adapter-memory: close
# ===========================================================================

test_that("memory adapter close is a no-op", {
  expect_null(.memory_close(list()))
})

# ===========================================================================
# stubs-core: subset_eager and derive_eager
# ===========================================================================

test_that("subset_eager subsets and computes in one step", {
  beta <- array(1:12, dim = c(3, 2, 2))
  var  <- array(0.1, dim = c(3, 2, 2))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("a", "b", "c")),
    subjects = c("s1", "s2"), contrasts = c("c1", "c2")
  )
  p <- as_plan(g)
  result <- subset_eager(p, subject = "s1")
  expect_s3_class(result, "gds")
  expect_equal(result$subjects, "s1")
})

test_that("derive_eager derives and computes in one step", {
  beta <- array(c(1, 2, 3, 4), dim = c(2, 1, 2))
  var  <- array(c(0.25, 0.16, 0.09, 0.04), dim = c(2, 1, 2))
  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("a", "b")),
    subjects = "s1", contrasts = c("c1", "c2")
  )
  p <- as_plan(g)
  result <- derive_eager(p, "se")
  expect_s3_class(result, "gds")
  expect_true("se" %in% names(result$assays))
})
