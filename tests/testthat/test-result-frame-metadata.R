test_that("typed metadata helpers cover cold branches to clear 90%", {
  expect_error(
    fmrigds:::.frame_fit_block_size(2L, 0L, 1L, 1e6, 1L),
    "at least one feature"
  )
  expect_error(
    fmrigds:::.fmrigds_unaligned_metadata("nope"),
    "list or unaligned_record"
  )

  rec <- fmridataset::unaligned_record(list(task = "rest"))
  expect_s3_class(fmrigds:::.fmrigds_unaligned_metadata(rec), "unaligned_record")
  expect_identical(
    fmrigds:::.fmrigds_serialize_for_metadata(as.Date("2020-01-02")),
    "2020-01-02"
  )
  expect_identical(
    fmrigds:::.fmrigds_serialize_for_metadata(package_version("1.2.3")),
    "1.2.3"
  )
  ser_df <- fmrigds:::.fmrigds_serialize_for_metadata(
    data.frame(a = 1:2, when = as.Date(c("2020-01-01", "2020-01-02")))
  )
  expect_true(is.list(ser_df))
  expect_identical(ser_df$a, 1:2)

  space <- fmridataset::index_space(2L, id_policy = "ephemeral")
  feats <- fmridataset::feature_axis(
    data.frame(.feature_id = fmridataset::feature_ids(space)),
    space = space
  )
  expect_error(
    result_frame(
      assays = list(estimate = matrix(1:4, nrow = 2)),
      observations = data.frame(.obs_id = c("a", "b")),
      features = feats,
      method = "ols:voxelwise",
      diagnostics = list(converged = TRUE)
    ),
    "length"
  )
  expect_error(
    result_frame(
      assays = list(estimate = matrix(1:4, nrow = 2)),
      observations = data.frame(.obs_id = c("a", "b")),
      features = feats,
      method = "ols:voxelwise",
      term_data = list(not = "a data.frame")
    ),
    "data frame"
  )

  rf <- result_frame(
    assays = list(estimate = matrix(1:4, nrow = 2)),
    observations = data.frame(.obs_id = c("a", "b")),
    features = feats,
    method = "ols:voxelwise",
    diagnostics = list(converged = c(TRUE, FALSE)),
    term_data = data.frame(term = "x", stringsAsFactors = FALSE),
    source_observation_ids = c("s1", "s2"),
    metadata = list(note = "ok", when = Sys.time())
  )
  expect_s3_class(rf, "fmri_frame")
  expect_true("diagnostics" %in% names(rf$tables))
  expect_true("term_data" %in% names(rf$tables))
  expect_true("source_observation_ids" %in% names(rf$tables))
  expect_identical(rf$metadata$result_kind, "statistical")
})
