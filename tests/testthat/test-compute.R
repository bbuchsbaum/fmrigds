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

  plan <- gds(tmp) %>% derive(c("var", "t"))
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
