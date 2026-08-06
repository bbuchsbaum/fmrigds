.register_counted_stage_adapter <- function(counter) {
  register_adapter(
    "test:stage",
    detect = function(source) FALSE,
    open = function(source, ...) source,
    probe = function(handle, ...) {
      arrays <- assays(handle)
      dimensions <- dim(arrays[[1L]])
      probe_contract(list(
        assays = names(arrays),
        dims = gds_dims(
          sample = dimensions[1L],
          subject = dimensions[2L],
          contrast = dimensions[3L]
        ),
        subjects = subjects(handle),
        contrasts = contrasts(handle),
        space = space(handle),
        maps = list(),
        columns = list(),
        metadata = metadata(handle),
        col_data = col_data(handle),
        row_data = row_data(handle)
      ))
    },
    read = function(handle, assays, block = NULL, ...) {
      counter$reads <- counter$reads + 1L
      selected <- handle$assays[assays]
      if (is.null(block)) return(selected)
      lapply(selected, function(value) {
        value[
          block$sample %||% seq_len(dim(value)[1L]),
          block$subject %||% seq_len(dim(value)[2L]),
          block$contrast %||% seq_len(dim(value)[3L]),
          drop = FALSE
        ]
      })
    },
    close = function(handle) invisible(NULL),
    capabilities = list(requires_staging = TRUE)
  )
}

test_that("staged and direct execution agree and source is read once", {
  skip_if_not_installed("hdf5r")
  g <- .group_examination_fixture()
  counter <- new.env(parent = emptyenv())
  counter$reads <- 0L
  .register_counted_stage_adapter(counter)
  source <- gds_source(
    "test:stage",
    g,
    probe_result = get_adapter("test:stage")$probe(g)
  )
  stage_dir <- tempfile("exam-stage-")
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)
  control <- examination_control(
    block_size = 7L,
    staging = list(tempdir = stage_dir)
  )
  staged <- examine_group(
    reduce(as_plan(source), method = "meta:fe"),
    control = control
  )
  direct <- examine_group(
    reduce(as_plan(g), method = "meta:fe"),
    control = control
  )

  expect_identical(counter$reads, 1L)
  expect_equal(staged$subject_data, direct$subject_data, tolerance = 1e-12)
  expect_equal(staged$contrast_data, direct$contrast_data, tolerance = 1e-12)
  expect_equal(staged$estimand_data, direct$estimand_data, tolerance = 1e-12)
  dimensions <- grep("^dimension", names(staged$embedding$coordinates), value = TRUE)
  expect_equal(
    staged$embedding$coordinates[, c("subject", dimensions), drop = FALSE],
    direct$embedding$coordinates[, c("subject", dimensions), drop = FALSE],
    tolerance = 1e-8
  )
  expect_equal(
    staged$embedding$captured_energy,
    direct$embedding$captured_energy,
    tolerance = 1e-10
  )
  expect_true(staged$provenance$staging$cleanup_succeeded)
  expect_identical(list.files(stage_dir), character())
  expect_identical(staged$provenance$scan_receipt$staging$adapter_reads, 1L)
  expect_gt(staged$provenance$scan_receipt$stage_size_bytes, 0)
  expect_true(is.finite(staged$provenance$scan_receipt$elapsed_seconds))
  expect_true(is.finite(staged$provenance$scan_receipt$peak_rss_bytes))
  expect_gt(staged$provenance$scan_receipt$peak_rss_bytes, 0)
  expect_gt(staged$provenance$scan_receipt$bytes_read, 0)
  expect_identical(
    staged$provenance$scan_receipt$retained_map_count,
    length(assays(staged$subject_maps))
  )
  expect_identical(
    staged$provenance$scan_receipt$staging$cleanup_policy,
    "always"
  )
})

test_that("support-unaware map prefixes stage to the same canonical data", {
  skip_if_not_installed("hdf5r")
  g <- .group_examination_fixture()
  map <- matrix(0, 4, 40)
  for (i in seq_len(4)) map[i, ((i - 1) * 10 + 1):(i * 10)] <- 0.1
  mapped_plan <- as_plan(g) |>
    map_to(
      space_sample_labels(paste0("parcel-", seq_len(4))),
      map,
      uncertainty = UncertaintyRule("independent")
    )
  staged <- examine_group(
    reduce(mapped_plan, method = "meta:fe"),
    control = examination_control(block_size = 2L)
  )
  canonical <- compute(mapped_plan)
  direct <- examine_group(
    reduce(as_plan(canonical), method = "meta:fe"),
    control = examination_control(block_size = 2L)
  )
  expect_equal(staged$subject_data, direct$subject_data, tolerance = 1e-12)
  expect_equal(staged$contrast_data, direct$contrast_data, tolerance = 1e-12)
  expect_equal(
    assay(staged$group_maps, "stat:pooled_effect"),
    assay(direct$group_maps, "stat:pooled_effect"),
    tolerance = 1e-12
  )
  expect_identical(space(staged$group_maps)$labels, paste0("parcel-", seq_len(4)))
})

test_that("staging cleanup runs after a post-stage preflight failure", {
  skip_if_not_installed("hdf5r")
  g <- .group_examination_fixture()
  counter <- new.env(parent = emptyenv())
  counter$reads <- 0L
  .register_counted_stage_adapter(counter)
  source <- gds_source(
    "test:stage",
    g,
    probe_result = get_adapter("test:stage")$probe(g)
  )
  register_reducer(
    "test:unsupported-examination",
    function(beta, var, X, z, p, df, df1, df2, opts) {
      list(beta_g = colMeans(beta))
    },
    requires = "beta",
    provides = "beta_g"
  )
  stage_dir <- tempfile("exam-stage-failure-")
  dir.create(stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)
  expect_error(
    examine_group(
      as_plan(source),
      method = "test:unsupported-examination",
      control = examination_control(staging = list(tempdir = stage_dir))
    ),
    "no implemented"
  )
  expect_identical(counter$reads, 1L)
  expect_identical(list.files(stage_dir), character())
})

test_that("NIfTI examination stages one source read before both scans", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("RNifti")
  skip_if_not_installed("neuroim2")
  directory <- tempfile("nifti-examination-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  subjects <- paste0("s", seq_len(6))
  files <- file.path(directory, paste0(subjects, "_beta.nii.gz"))
  for (i in seq_along(files)) {
    RNifti::writeNifti(array(i + seq_len(8) / 10, c(2, 2, 2)), files[i])
  }
  source <- nifti_source(
    beta = files,
    subjects = subjects,
    contrast = "metric"
  )
  stage_dir <- file.path(directory, "stage")
  dir.create(stage_dir)
  exam <- suppressWarnings(examine_group(
    reduce(
      gds(source, format = "nifti"),
      method = "ols:voxelwise",
      formula = ~ 1
    ),
    estimands = "(Intercept)",
    control = examination_control(
      block_size = 3L,
      staging = list(tempdir = stage_dir)
    )
  ))
  expect_identical(exam$provenance$staging$source_adapter, "nifti")
  expect_identical(exam$provenance$staging$adapter_reads, 1L)
  expect_identical(exam$cohort$variance_mode, "effect_only_synthetic")
  expect_identical(list.files(stage_dir), character())
})

test_that("source fingerprints change with file identity metadata", {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  writeLines("first", path)
  source <- gds_source("test", path)
  first <- fmrigds:::.source_fingerprint(source)
  writeLines(c("first", "second line"), path)
  second <- fmrigds:::.source_fingerprint(source)
  expect_false(identical(first$digest, second$digest))
})

test_that("space portability distinguishes function-valued bases", {
  expect_identical(
    space_portability(space_sample_labels(c("a", "b"))),
    "portable"
  )
  expect_identical(
    space_portability(space_basis(2, projector = diag(2))),
    "portable"
  )
  expect_identical(
    space_portability(space_basis(2, projector = function(x) x)),
    "nonportable_function"
  )
})
