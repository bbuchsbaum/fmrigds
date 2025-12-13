test_that("posthoc registry has built-ins", {
  expect_true(all(c("fdr:bh", "fdr:by") %in% list_posthoc()))
  expect_true(!is.null(get_posthoc("fdr:bh")))
})

test_that("FDR BH matches p.adjust across samples", {
  set.seed(1)
  df <- data.frame(
    sample = rep(paste0("ROI_", 1:5), each = 3),
    subject = rep(c("s1", "s2", "s3"), times = 5),
    contrast = "c1",
    p = runif(15),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE)
  plan <- gds(tmp, effect_cols = list(p = "p"))
  plan <- posthoc(plan, method = "fdr:bh")
  g <- compute(plan)
  q <- assay(g, "q")
  p <- compute(gds(tmp, effect_cols = list(p = "p")), assays = "p")$p
  for (j in seq_len(dim(q)[2])) {
    for (k in seq_len(dim(q)[3])) {
      expect_equal(as.numeric(q[, j, k]), as.numeric(stats::p.adjust(p[, j, k], method = "BH")))
    }
  }
})

test_that("FDR BY matches p.adjust across samples", {
  set.seed(2)
  df <- data.frame(
    sample = rep(paste0("ROI_", 1:4), each = 2),
    subject = rep(c("s1", "s2"), times = 4),
    contrast = "c1",
    p = runif(8),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE)
  p_plan <- gds(tmp, effect_cols = list(p = "p"))
  by_plan <- posthoc(p_plan, method = "fdr:by")
  out <- compute(by_plan)
  q <- assay(out, "q")
  p <- compute(gds(tmp, effect_cols = list(p = "p")), assays = "p")$p
  for (j in seq_len(dim(q)[2])) {
    for (k in seq_len(dim(q)[3])) {
      expect_equal(as.numeric(q[, j, k]), as.numeric(stats::p.adjust(p[, j, k], method = "BY")))
    }
  }
})
