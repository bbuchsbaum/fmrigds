test_that("gds_to_tibble returns expected long table", {
  df <- data.frame(
    sample = rep(c("ROI_1", "ROI_2"), each = 2),
    subject = rep(c("s1", "s2"), times = 2),
    contrast = "c1",
    beta = c(0.5, 0.6, 0.7, 0.8),
    var = c(0.04, 0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  g <- as_gds(df)
  tab <- gds_to_tibble(g, assays = c("beta", "var"))
  expect_true(all(c("sample", "subject", "contrast", "beta", "var") %in% names(tab)))
  expect_equal(nrow(tab), 2 * length(subjects(g)) * length(contrasts(g)))
})

