.rectangular_frame_fixture <- function(instrument = FALSE) {
  subjects <- c("sub-02", "sub-01")
  contrasts <- c("task", "baseline", "followup")
  grid <- expand.grid(
    subject_id = subjects,
    contrast_id = contrasts,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[c(4L, 1L, 6L, 3L, 2L, 5L), , drop = FALSE]
  grid$.obs_id <- paste(grid$subject_id, grid$contrast_id, sep = "::")
  grid$condition <- factor(
    ifelse(grid$contrast_id == "baseline", "rest", "active"),
    levels = c("rest", "active")
  )

  beta <- matrix(seq_len(nrow(grid) * 3L), nrow = nrow(grid)) / 10
  variance <- beta / 10 + 0.05
  sources <- list(
    beta = fmridataset::memory_source(beta),
    variance = fmridataset::memory_source(variance)
  )
  if (instrument) sources <- lapply(sources, fmridataset::counting_source)
  spatial <- fmridataset::volume_space(
    dim = c(2L, 2L, 1L),
    affine = diag(4),
    support = c(1L, 2L, 4L),
    template = "projection-fixture"
  )
  subject <- fmridataset::entity_frame(
    data = data.frame(subject_id = subjects, age = c(70, 62)),
    key = "subject_id",
    entity_type = "subject"
  )
  contrast <- fmridataset::entity_frame(
    data = data.frame(
      contrast_id = contrasts,
      time = c(1, 0, 2),
      condition = factor(c("active", "rest", "active"),
                         levels = c("rest", "active"))
    ),
    key = "contrast_id",
    entity_type = "contrast"
  )
  frame <- fmridataset::fmri_frame(
    assays = sources,
    observations = grid,
    features = fmridataset::feature_axis(
      data.frame(
        .feature_id = fmridataset::feature_ids(spatial),
        region = c("left", "left", "right")
      ),
      space = spatial
    ),
    entities = list(subject = subject, contrast = contrast),
    relations = list(
      observation_subject = fmridataset::key_relation(
        "subject_id", target = "subject"
      ),
      observation_contrast = fmridataset::key_relation(
        "contrast_id", target = "contrast"
      )
    )
  )
  list(
    frame = frame,
    beta = beta,
    variance = variance,
    subjects = subjects,
    contrasts = contrasts
  )
}

test_that("rectangular frame projection is explicit and axis exact", {
  fixture <- .rectangular_frame_fixture()
  g <- as_gds(
    fixture$frame,
    subject = "subject_id",
    contrast = "contrast_id"
  )

  expect_s3_class(g, "gds")
  expect_identical(subjects(g), fixture$subjects)
  expect_identical(contrasts(g), fixture$contrasts)
  expect_identical(dim(assay(g, "beta")), c(3L, 2L, 3L))
  expect_identical(space(g)$mask_idx, c(1L, 2L, 4L))
  expect_identical(space(g)$template_id, "projection-fixture")
  expect_identical(rownames(row_data(g)), fmridataset::feature_ids(fixture$frame))
  expect_identical(rownames(col_data(g)), fixture$subjects)
  expect_identical(rownames(contrast_data(g)), fixture$contrasts)

  expected <- array(NA_real_, c(3L, 2L, 3L))
  observations <- fmridataset::observations(fixture$frame)
  for (row in seq_len(nrow(observations))) {
    s <- match(observations$subject_id[[row]], fixture$subjects)
    c <- match(observations$contrast_id[[row]], fixture$contrasts)
    expected[, s, c] <- fixture$beta[row, ]
  }
  expect_equal(unname(assay(g, "beta")), expected)
})

test_that("GDS round trip preserves assays, key order, factors, and spatial identity", {
  fixture <- .rectangular_frame_fixture()
  g <- as_gds(fixture$frame, subject = "subject_id", contrast = "contrast_id")
  restored <- fmridataset::as_fmri_frame(g)
  projected <- as_gds(restored, subject = "subject_id", contrast = "contrast_id")

  expect_identical(fmridataset::observation_ids(restored),
                   as.vector(outer(subjects(g), contrasts(g), paste, sep = "::")))
  expect_identical(levels(fmridataset::observations(restored)$condition),
                   c("rest", "active"))
  expect_identical(fmridataset::feature_ids(restored),
                   fmridataset::feature_ids(fixture$frame))
  expect_identical(
    fmridataset::space_digest(fmridataset::space(restored)),
    fmridataset::space_digest(fmridataset::space(fixture$frame))
  )
  expect_identical(subjects(projected), subjects(g))
  expect_identical(contrasts(projected), contrasts(g))
  for (name in names(assays(g))) {
    expect_equal(assay(projected, name), assay(g, name))
  }
  expect_identical(
    fmridataset::space_digest(fmridataset::space(
      unserialize(serialize(restored, NULL))
    )),
    fmridataset::space_digest(fmridataset::space(restored))
  )
})

test_that("rectangular projection is invariant to observation order", {
  fixture <- .rectangular_frame_fixture()
  original <- as_gds(
    fixture$frame, subject = "subject_id", contrast = "contrast_id"
  )
  reordered <- as_gds(
    fixture$frame[c(6L, 2L, 4L, 1L, 5L, 3L), ],
    subject = "subject_id",
    contrast = "contrast_id"
  )

  expect_identical(subjects(reordered), subjects(original))
  expect_identical(contrasts(reordered), contrasts(original))
  for (name in names(assays(original))) {
    expect_equal(assay(reordered, name), assay(original, name))
  }
})

test_that("projection rejects ragged, duplicate, and missing key domains before reads", {
  fixture <- .rectangular_frame_fixture(instrument = TRUE)
  frame <- fixture$frame
  beta_source <- fmridataset::assay(frame, "beta")$source
  variance_source <- fmridataset::assay(frame, "variance")$source

  expect_error(
    as_gds(frame[-1L, ], subject = "subject_id", contrast = "contrast_id"),
    "complete Cartesian grid"
  )
  duplicate <- fmridataset::observations(frame)
  duplicate$contrast_id[[2L]] <- duplicate$contrast_id[[1L]]
  frame$observations$data <- duplicate
  expect_error(
    as_gds(frame, subject = "subject_id", contrast = "contrast_id"),
    "unique subject-contrast pairs"
  )
  expect_error(
    as_gds(fixture$frame, subject = "missing", contrast = "contrast_id"),
    "subject column"
  )
  expect_identical(fmridataset::source_counts(beta_source)$bytes, 0)
  expect_identical(fmridataset::source_counts(variance_source)$bytes, 0)
})

test_that("projection honors frame views and explicit realization budgets", {
  fixture <- .rectangular_frame_fixture(instrument = TRUE)
  view <- fixture$frame[c(2L, 5L, 1L, 4L), c(3L, 1L)]
  beta_source <- fmridataset::assay(fixture$frame, "beta")$source

  expect_error(
    as_gds(view, subject = "subject_id", contrast = "contrast_id",
           memory_budget = 1L),
    "memory_budget"
  )
  expect_identical(fmridataset::source_counts(beta_source)$bytes, 0)

  g <- as_gds(
    view,
    subject = "subject_id",
    contrast = "contrast_id",
    memory_budget = 1024L
  )
  expect_identical(subjects(g), c("sub-02", "sub-01"))
  expect_identical(contrasts(g), c("task", "baseline"))
  expect_identical(rownames(row_data(g)), fmridataset::feature_ids(view))
  expect_gt(fmridataset::source_counts(beta_source)$bytes, 0)
})

test_that("group result frames project through an explicit singleton GDS axis", {
  fixture <- .rectangular_frame_fixture()
  fit <- fit_group(
    fixture$frame,
    design = multidesign::design_spec(~condition),
    method = "ols:voxelwise",
    variance = NULL,
    block_size = 2L
  )
  g <- as_gds(
    fit$result,
    subject = NULL,
    contrast = ".obs_id",
    subjects = "meta"
  )

  expect_identical(subjects(g), "meta")
  expect_identical(contrasts(g), fmridataset::observation_ids(fit$result))
  expect_equal(
    assay(g, "estimate")[, 1L, ],
    t(fmridataset::collect_assay(fit$result, "estimate")),
    ignore_attr = TRUE
  )
})

test_that("plain legacy GDS objects convert to canonical frames", {
  beta <- array(
    seq_len(3L * 2L * 2L) / 10,
    dim = c(3L, 2L, 2L),
    dimnames = list(c("voxel-1", "voxel-2", "voxel-4"),
                    c("s2", "s1"), c("b", "a"))
  )
  g <- new_gds(
    assays = list(beta = beta, var = beta / 10 + 0.1),
    space = space_voxel(
      c(2L, 2L, 1L), diag(4), mask_idx = c(1L, 2L, 4L),
      storage = "packed", template_id = "legacy"
    ),
    subjects = c("s2", "s1"),
    contrasts = c("b", "a"),
    col_data = data.frame(age = c(70, 62), row.names = c("s2", "s1")),
    row_data = data.frame(region = c("L", "L", "R"),
                          row.names = c("voxel-1", "voxel-2", "voxel-4"))
  )
  frame <- fmridataset::as_fmri_frame(g)

  expect_identical(dim(frame), c(4L, 3L))
  expect_identical(names(fmridataset::assays(frame)), c("beta", "variance"))
  expect_identical(fmridataset::observations(frame)$subject_id,
                   c("s2", "s1", "s2", "s1"))
  expect_identical(fmridataset::observations(frame)$contrast_id,
                   c("b", "b", "a", "a"))
  expect_identical(fmridataset::feature_ids(frame),
                   c("voxel-1", "voxel-2", "voxel-4"))
  expect_identical(fmridataset::space(frame)$support, c(1L, 2L, 4L))
  expect_equal(
    fmridataset::collect_assay(frame, "beta"),
    matrix(aperm(beta, c(2L, 3L, 1L)), nrow = 4L, ncol = 3L)
  )
})

test_that("projection validates key and assay mapping contracts", {
  fixture <- .rectangular_frame_fixture()
  expect_error(
    as_gds(fixture$frame, subject = "subject_id", contrast = "subject_id"),
    "distinct"
  )
  expect_error(
    as_gds(
      fixture$frame,
      subject = "subject_id",
      contrast = "contrast_id",
      assay_map = c(beta = "same", variance = "same")
    ),
    "duplicate"
  )
})
