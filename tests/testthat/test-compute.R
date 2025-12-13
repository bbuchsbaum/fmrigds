test_that("compute derives statistics via plan", {
  skip_if_not_installed("data.table")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  df <- data.frame(
    sample = c("roi1", "roi2"),
    subject = "s1",
    contrast = "c1",
    beta = c(1, 2),
    se = c(0.5, 1)
  )
  data.table::fwrite(df, tmp)

  plan <- gds(tmp, effect_cols = list(beta = "beta", se = "se")) %>% derive(c("var", "t"))
  g <- compute(plan)

  expect_true("var" %in% names(assays(g)))
  expect_true("t" %in% names(assays(g)))
  expect_equal(assay(g, "var")[1, 1, 1], 0.25)
  expect_equal(assay(g, "t")[1, 1, 1], 2)
})

test_that("compute maps with variance propagation", {
  skip_if_not_installed("data.table")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  df <- data.frame(
    sample = c("roi1", "roi2"),
    subject = "s1",
    contrast = "c1",
    beta = c(1, 2),
    var = c(0.25, 1)
  )
  data.table::fwrite(df, tmp)

  target <- space_parcels("parcel1")
  M <- matrix(c(0.6, 0.8), nrow = 1)

  plan <- gds(tmp) %>%
    map_to(target_space = target, map = M, uncertainty = UncertaintyRule("independent"))

  g <- compute(plan)

  expect_equal(assay(g, "beta")[1, 1, 1], sum(M * c(1, 2)))
  expect_equal(space(g)$labels, "parcel1")
})

test_that("compute applies mask policies", {
  skip_if_not_installed("data.table")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  df <- data.frame(
    sample = c("roi1", "roi2"),
    subject = "s1",
    contrast = "c1",
    beta = c(1, NA),
    var = c(1, 1)
  )
  data.table::fwrite(df, tmp)

  plan <- gds(tmp) %>% mask(MaskPolicy(rule = "intersection"))
  g <- compute(plan)

  expect_equal(dim(assay(g, "beta"))[1], 1)
})

test_that("compute forwards block argument to adapter read", {
  env <- new.env(parent = emptyenv())

  register_adapter(
    name = "dummy_block",
    detect = function(source) 0,
    open = function(source, ...) list(),
    probe = function(handle, ...) {
      list(
        assays = c("beta", "var"),
        dims = c(sample = 4L, subject = 2L, contrast = 1L),
        subjects = c("s1", "s2"),
        contrasts = "c1",
        space = space_parcels(paste0("roi", 1:4)),
        maps = list(),
        metadata = gds_metadata(),
        columns = list()
      )
    },
    read = function(handle, assays, block, ...) {
      env$block <- block
      list(
        beta = array(0, c(4, 2, 1)),
        var = array(1, c(4, 2, 1))
      )
    },
    close = function(handle) NULL
  )

  src <- gds_source(
    "dummy_block",
    list(),
    probe_result = list(
      assays = c("beta", "var"),
      dims = c(sample = 4L, subject = 2L, contrast = 1L),
      subjects = c("s1", "s2"),
      contrasts = "c1",
      space = space_parcels(paste0("roi", 1:4)),
      maps = list(),
      metadata = gds_metadata(),
      columns = list()
    )
  )
  plan <- gds_plan(src)
  res <- compute(plan, block = list(sample = 1:2))
  expect_s3_class(res, "gds")
  expect_equal(env$block$sample, 1:2)
})
