test_that("with_contrast_data aligns contrasts on plans and survives compute/subset", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "ROI_1,sub-01,baseline,0.5,0.04",
    "ROI_1,sub-01,task,0.8,0.04",
    "ROI_1,sub-02,baseline,0.4,0.05",
    "ROI_1,sub-02,task,0.9,0.05"
  ), tmp)

  plan <- gds(tmp)
  kd <- data.frame(
    condition = c(1, 0),
    label = c("task", "baseline"),
    row.names = c("task", "baseline"),
    stringsAsFactors = FALSE
  )

  plan2 <- with_contrast_data(plan, kd)
  expect_equal(rownames(contrast_data(plan2)), c("baseline", "task"))
  expect_equal(as.numeric(contrast_data(plan2)$condition), c(0, 1))

  g <- subset(plan2, contrast = "task") |> compute()
  expect_equal(contrasts(g), "task")
  expect_equal(rownames(contrast_data(g)), "task")
  expect_equal(as.numeric(contrast_data(g)$condition), 1)
})

test_that("contrast_data persists through HDF5 write and read", {
  skip_if_not_installed("hdf5r")

  beta <- array(rnorm(8), dim = c(2, 2, 2))
  var <- array(runif(8, 0.1, 0.2), dim = c(2, 2, 2))
  kd <- data.frame(
    condition = c(0, 1),
    row.names = c("baseline", "task"),
    stringsAsFactors = FALSE
  )

  g <- new_gds(
    assays = list(beta = beta, var = var),
    space = space_sample_labels(c("ROI_1", "ROI_2")),
    subjects = c("s1", "s2"),
    contrasts = c("baseline", "task")
  )
  g <- with_contrast_data(g, kd)

  out <- tempfile(fileext = ".h5")
  on.exit(unlink(out), add = TRUE)
  write_out(g, out, format = "h5") |> compute()

  plan2 <- gds(out)
  expect_equal(contrast_data(plan2), kd)

  g2 <- compute(plan2)
  expect_equal(contrast_data(g2), kd)
})
