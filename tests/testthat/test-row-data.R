test_that("with_row_data aligns samples on plans and survives compute/subset", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,cope1,0.5,0.04",
    "ROI_2,sub-01,cope1,0.8,0.09",
    "ROI_1,sub-02,cope1,0.4,0.05",
    "ROI_2,sub-02,cope1,0.7,0.08"
  ), tmp)

  plan <- gds(tmp)
  rd <- data.frame(
    feature_group = c("g2", "g1"),
    sample = c("ROI_2", "ROI_1"),
    row.names = c("ROI_2", "ROI_1"),
    stringsAsFactors = FALSE
  )

  plan2 <- with_row_data(plan, rd)
  expect_equal(rownames(row_data(plan2)), c("ROI_1", "ROI_2"))
  expect_equal(sample_labels(plan2), c("ROI_1", "ROI_2"))
  expect_equal(as.character(sample_groups(plan2)), c("g1", "g2"))

  g <- subset(plan2, sample = 2) |> compute()
  expect_equal(sample_labels(g), "ROI_2")
  expect_equal(as.character(sample_groups(g)), "g2")
  expect_equal(row_data(g)$sample, "ROI_2")
})

test_that("spatial FDR resolves group metadata from row_data", {
  skip_if(!("fdr:spatial" %in% list_posthoc()))
  set.seed(24)
  n_subj <- 2
  n_samp <- 6
  p <- array(runif(n_samp * n_subj), dim = c(n_samp, n_subj, 1))
  labels <- paste0("ROI_", seq_len(n_samp))
  grp <- c("A", "A", "B", "B", "C", "C")

  g <- as_gds(
    list(p = p),
    space = space_sample_labels(labels),
    subjects = paste0("s", seq_len(n_subj)),
    contrasts = "cope1",
    row_data = data.frame(sample = labels, feature_group = grp, row.names = labels, stringsAsFactors = FALSE)
  )

  gout <- posthoc(g, method = "fdr:spatial") |> compute()
  q <- assay(gout, "q")[, 1, 1]
  expect_equal(q[1], q[2])
  expect_equal(q[3], q[4])
  expect_equal(q[5], q[6])
})

test_that("row_data persists through HDF5 write and read", {
  skip_if_not_installed("hdf5r")
  beta <- array(rnorm(4), dim = c(2, 2, 1))
  var <- array(runif(4, 0.1, 0.2), dim = c(2, 2, 1))
  labels <- c("ROI_1", "ROI_2")
  rd <- data.frame(sample = labels, feature_group = c("g1", "g2"), row.names = labels, stringsAsFactors = FALSE)

  g <- as_gds(
    list(beta = beta, var = var),
    space = space_sample_labels(labels),
    subjects = c("s1", "s2"),
    contrasts = "cope1",
    row_data = rd
  )

  out <- tempfile(fileext = ".h5")
  on.exit(unlink(out), add = TRUE)
  write_out(g, out, format = "h5") |> compute()

  plan2 <- gds(out)
  expect_equal(row_data(plan2), rd)
  g2 <- compute(plan2)
  expect_equal(row_data(g2), rd)
  expect_equal(as.character(sample_groups(g2)), c("g1", "g2"))
})
