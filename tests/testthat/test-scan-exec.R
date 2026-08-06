test_that("direct scanner is partition-invariant and reports I/O", {
  beta <- array(seq_len(11 * 3 * 2), c(11, 3, 2))
  g <- new_gds(
    assays = list(beta = beta, var = array(4, dim(beta))),
    space = space_sample_labels(sprintf("f%02d", 1:11)),
    subjects = c("s1", "s2", "s3"),
    contrasts = c("c1", "c2")
  )
  compiled <- as_plan(g) |>
    subset(sample = c(2:6, 9:11), subject = c("s3", "s1"), contrast = "c2") |>
    derive("se") |>
    reduce(method = "fixed") |>
    fmrigds:::compile_examination_plan()

  run <- function(block_size) {
    fmrigds:::.scan_compiled_plan(
      compiled,
      assays = c("beta", "var"),
      block_size = block_size,
      initialize = function(context) list(total = 0, ids = integer(), finalized = 0L),
      update = function(state, arrays, block, context) {
        expect_true("se" %in% names(arrays))
        state$total <- state$total + sum(arrays$beta)
        state$ids <- c(state$ids, block$source_sample)
        state
      },
      finalize = function(state, receipt, context) {
        state$finalized <- state$finalized + 1L
        state
      }
    )
  }

  small <- run(2L)
  large <- run(5L)
  expected_idx <- c(2:6, 9:11)
  expect_equal(small$value$total, sum(beta[expected_idx, c(3, 1), 2]))
  expect_equal(small$value$total, large$value$total, tolerance = 1e-12)
  expect_identical(small$value$ids, expected_idx)
  expect_identical(large$value$ids, expected_idx)
  expect_identical(small$value$finalized, 1L)
  expect_equal(small$receipt$adapter_reads, ceiling(length(expected_idx) / 2))
  expect_gt(small$receipt$bytes_read, 0)
  expect_true(is.finite(small$receipt$elapsed_seconds))
  expect_true(is.finite(small$receipt$peak_rss_bytes))
  expect_gt(small$receipt$peak_rss_bytes, 0)
  expect_identical(small$receipt$scan_strategy, "direct")
  expect_identical(small$receipt$samples_emitted, length(expected_idx))
})

test_that("scanner closes and invokes error cleanup after update failure", {
  lifecycle <- new.env(parent = emptyenv())
  lifecycle$closed <- 0L
  lifecycle$errored <- 0L
  lifecycle$finalized <- 0L
  values <- array(1:8, c(4, 2, 1))

  register_adapter(
    "scan_lifecycle_fixture",
    detect = function(source) 0,
    open = function(source, ...) list(values = source),
    probe = function(handle, ...) NULL,
    read = function(handle, assays, block = NULL, ...) {
      idx <- block$sample
      list(beta = handle$values[idx, , , drop = FALSE])
    },
    close = function(handle) lifecycle$closed <- lifecycle$closed + 1L,
    capabilities = list(
      sample_blocks = TRUE,
      subject_blocks = TRUE,
      contrast_blocks = TRUE,
      persistent_handle = TRUE,
      preferred_axis = "sample",
      cheap_revisit = TRUE
    )
  )
  probe <- probe_contract(list(
    assays = "beta",
    dims = gds_dims(sample = 4, subject = 2, contrast = 1),
    subjects = c("s1", "s2"),
    contrasts = "c1",
    space = space_sample_labels(letters[1:4]),
    maps = list(),
    metadata = list(),
    columns = list()
  ))
  plan <- gds_plan(gds_source("scan_lifecycle_fixture", values, probe)) |>
    reduce(method = "ols:voxelwise", formula = ~ 1)
  compiled <- fmrigds:::compile_examination_plan(plan)

  expect_error(
    fmrigds:::.scan_compiled_plan(
      compiled,
      assays = "beta",
      block_size = 2L,
      initialize = function(context) list(),
      update = function(state, arrays, block, context) stop("injected update failure"),
      finalize = function(state, receipt, context) {
        lifecycle$finalized <- lifecycle$finalized + 1L
        state
      },
      on_error = function(state, error, receipt, context) {
        lifecycle$errored <- lifecycle$errored + 1L
      }
    ),
    "injected update failure"
  )
  expect_identical(lifecycle$closed, 1L)
  expect_identical(lifecycle$errored, 1L)
  expect_identical(lifecycle$finalized, 0L)
})

test_that("scanner never executes discarded write nodes", {
  beta <- array(1:8, c(2, 4, 1))
  g <- new_gds(
    assays = list(beta = beta, var = array(1, dim(beta))),
    space = space_sample_labels(c("a", "b")),
    subjects = paste0("s", 1:4),
    contrasts = "c1"
  )
  target <- tempfile(fileext = ".csv")
  compiled <- as_plan(g) |>
    write_out(target, "csv") |>
    reduce(method = "fixed") |>
    fmrigds:::compile_examination_plan()

  fmrigds:::.scan_compiled_plan(
    compiled,
    assays = c("beta", "var"),
    block_size = 1L,
    initialize = function(context) 0,
    update = function(state, arrays, block, context) state + length(arrays$beta),
    finalize = function(state, receipt, context) state
  )
  expect_false(file.exists(target))
})
