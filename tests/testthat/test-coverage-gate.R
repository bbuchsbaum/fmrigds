# Focused coverage for the portfolio coverage-badge gate (>= 90%).
# Exercises real API/internal branches that existing suites leave cold.

# ---------------------------------------------------------------------------
# plan-optimizer.R: rewrite helpers remain callable even though optimize_plan
# currently preserves authored order for semantic safety.
# ---------------------------------------------------------------------------

test_that("plan optimizer rewrite helpers coalesce subsets/derives/masks", {
  a <- list(op = "subset_axis", sample = 1:5, subject = c("s1", "s2"), contrast = NULL)
  b <- list(op = "subset_axis", sample = 3:6, subject = "s1", contrast = "c1")
  merged <- fmrigds:::.merge_subset(a, b)
  expect_equal(merged$sample, 3:5)
  expect_equal(merged$subject, "s1")
  expect_equal(merged$contrast, "c1")

  nodes <- list(
    list(op = "subset_axis", sample = 1:3, subject = NULL, contrast = NULL),
    list(op = "derive", what = "t", options = list(a = 1)),
    list(op = "derive", what = c("z"), options = list(b = 2)),
    list(op = "subset_axis", sample = 2:4, subject = "s1", contrast = NULL),
    list(op = "mask_policy", policy = MaskPolicy(rule = "intersection")),
    list(op = "mask_policy", policy = MaskPolicy(rule = "threshold", threshold = 0.5)),
    list(op = "reduce", method = "fixed")
  )
  combined <- fmrigds:::.combine_subsets(nodes)
  expect_true(any(vapply(combined, function(n) identical(n$op, "subset_axis"), logical(1))))

  coalesced <- fmrigds:::.coalesce_derives(nodes)
  derive_nodes <- Filter(function(n) identical(n$op, "derive"), coalesced)
  expect_equal(length(derive_nodes), 1L)
  expect_true(all(c("t", "z") %in% derive_nodes[[1]]$what))

  reordered <- fmrigds:::.reorder_nodes(nodes)
  expect_equal(reordered[[1]]$op, "subset_axis")
  expect_equal(reordered[[length(reordered)]]$op, "reduce")

  fused <- fmrigds:::.fuse_masks(nodes)
  mask_nodes <- Filter(function(n) identical(n$op, "mask_policy"), fused)
  expect_equal(length(mask_nodes), 1L)
  expect_true(is.list(mask_nodes[[1]]$policy))
})

# ---------------------------------------------------------------------------
# map-exec.R: fisher combiner, error paths, and zero-weight rows
# ---------------------------------------------------------------------------

test_that("apply_map_to covers fisher/stouffer edges and guards", {
  target <- space_sample_labels(c("t1", "t2"))
  from_sp <- space_sample_labels(c("a", "b", "c"))
  M <- matrix(c(0.5, 0.5, 0, 0, 0, 0), nrow = 2, ncol = 3, byrow = TRUE)
  mp <- map_linear(from_sp, target, M)

  z <- array(c(1.2, -0.8, 0.5, 2.1, 0.4, -1.1), dim = c(3, 2, 1))
  p <- 2 * pnorm(-abs(z))

  fisher_z <- apply_map_to(
    list(op = "map", target_space = target, map = mp, uncertainty = UncertaintyRule("none"),
         combine = "fisher"),
    list(z = z)
  )
  expect_equal(dim(fisher_z$arrays$p), c(2, 2, 1))
  expect_true(all(is.na(fisher_z$arrays$p[2, , ])))

  fisher_p <- apply_map_to(
    list(op = "map", target_space = target, map = as.matrix(M), uncertainty = UncertaintyRule("none"),
         combine = "fisher"),
    list(p = p)
  )
  expect_true("chi2" %in% names(fisher_p$arrays))

  expect_error(
    apply_map_to(
      list(op = "map", target_space = target, map = "nope", uncertainty = UncertaintyRule("none")),
      list(z = z)
    ),
    "map must be matrix"
  )

  bad_mp <- mp
  bad_mp$by_subject <- list(s1 = M)
  expect_error(
    apply_map_to(
      list(op = "map", target_space = target, map = bad_mp, uncertainty = UncertaintyRule("none"),
           combine = "stouffer"),
      list(z = z)
    ),
    "align\\(\\)"
  )

  expect_error(
    apply_map_to(
      list(op = "map", target_space = target, map = M, uncertainty = UncertaintyRule("none")),
      list(z = z)
    ),
    "explicit combiner"
  )
  expect_error(
    apply_map_to(
      list(op = "map", target_space = target, map = M, uncertainty = UncertaintyRule("none"),
           combine = "stouffer"),
      list(p = p)
    ),
    "z-scores"
  )
  expect_error(
    apply_map_to(
      list(op = "map", target_space = target, map = M, uncertainty = UncertaintyRule("none"),
           combine = "fisher"),
      list(t = z)
    ),
    "Fisher combine"
  )

  beta <- array(1:6, dim = c(3, 2, 1))
  var <- array(0.25, dim = c(3, 2, 1))
  tvals <- beta / sqrt(var)
  mapped <- apply_map_to(
    list(
      op = "map",
      target_space = target,
      map = M,
      uncertainty = UncertaintyRule("independent")
    ),
    list(beta = beta, var = var, se = sqrt(var), t = tvals, z = tvals, df = array(20, dim = c(3, 2, 1)))
  )
  expect_equal(dim(mapped$arrays$beta), c(2, 2, 1))
  expect_true("z" %in% names(mapped$arrays))
})

# ---------------------------------------------------------------------------
# scalar-map-workflow.R: guards and metadata helpers
# ---------------------------------------------------------------------------

test_that("scalar-map helpers and two_sample/one_sample guards", {
  expect_error(
    fmrigds:::.scalar_maps_frame(data.frame(subject = "s1"), subject = "s1"),
    "`file` column"
  )
  expect_error(
    fmrigds:::.scalar_maps_frame(
      data.frame(file = "a.nii", stringsAsFactors = FALSE),
      subject = NULL,
      subject_supplied = TRUE
    ),
    "`subject` is required"
  )
  expect_error(
    fmrigds:::.scalar_maps_frame(
      c("a.nii", "b.nii"),
      subject = "only-one",
      contrast = "m"
    ),
    "Length of `subject`"
  )
  expect_error(
    fmrigds:::.scalar_maps_frame(
      c(s1 = "a.nii", s2 = "b.nii"),
      contrast = c("m1", "m2", "m3")
    ),
    "Length of `contrast`"
  )

  paths <- c(s01 = "s01_metric.nii.gz", s02 = "s02_metric.nii.gz")
  maps <- fmrigds:::.scalar_maps_frame(paths, contrast = "metric")
  expect_equal(maps$subject, c("s01", "s02"))

  expect_null(fmrigds:::.scalar_maps_col_data(paths))
  expect_null(fmrigds:::.normalise_keyed_col_data(NULL, c("s1")))
  expect_error(fmrigds:::.normalise_keyed_col_data(list(a = 1), c("s1")), "data.frame")

  cd <- data.frame(subject = c("s1", "s2"), group = c("A", "B"), stringsAsFactors = FALSE)
  normed <- fmrigds:::.normalise_keyed_col_data(cd, c("s1", "s2"))
  expect_equal(rownames(normed), c("s1", "s2"))

  g <- new_gds(
    assays = list(beta = array(1:4, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_sample_labels(c("r1", "r2")),
    subjects = c("s1", "s2"),
    contrasts = "c1",
    col_data = data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
  )
  expect_error(two_sample(g, group = "missing"), "Grouping column")
  expect_error(two_sample(g, group = "group", level = "B"), "`baseline` is required")
  expect_error(
    two_sample(g, group = "group", baseline = "X", level = "Z"),
    "No subjects match"
  )
  expect_error(
    two_sample(g, group = "group", baseline = "Z"),
    "Baseline level not present"
  )
  bare <- new_gds(
    assays = list(beta = array(1:2, c(2, 1, 1)), var = array(1, c(2, 1, 1))),
    space = space_sample_labels(c("r1", "r2")),
    subjects = "s1",
    contrasts = "c1"
  )
  bare$col_data <- NULL
  expect_error(two_sample(bare, group = "group"), "requires col_data")

  plan_ok <- one_sample(g, subset = group == "A")
  expect_s3_class(plan_ok, "gds_plan")
  bad_g <- g
  bad_g$col_data <- NULL
  expect_error(one_sample(bad_g, subset = TRUE), "`subset` requires col_data")
})

test_that("write_nifti_assays reports empty/non-image assay rows", {
  skip_if_not_installed("RNifti")
  sp <- space_voxel(dim = c(2, 2, 1), affine = diag(4), mask_idx = 1:4, storage = "packed")
  g <- new_gds(
    assays = list(
      beta = array(as.numeric(1:8), c(4, 2, 1)),
      var = array(1, c(4, 2, 1))
    ),
    space = sp,
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  td <- tempfile("nifti-cov-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  man <- write_nifti_assays(g, out_dir = td, assays = c("missing", "beta"), subjects = "s1")
  expect_true(any(man$skipped_reason == "assay not found", na.rm = TRUE))
  expect_true(any(man$written, na.rm = TRUE))
  # Non-image assay path: surgically shrink dim[1] after construction.
  g2 <- g
  g2$assays$beta <- array(1:2, c(2, 2, 1))
  man2 <- write_nifti_assays(g2, out_dir = td, assays = "beta", subjects = "s1", overwrite = TRUE)
  expect_true(any(man2$skipped_reason == "not an image assay", na.rm = TRUE))

  expect_error(
    write_nifti_assays(
      new_gds(
        assays = list(beta = array(1:2, c(2, 1, 1)), var = array(1, c(2, 1, 1))),
        space = space_sample_labels(c("a", "b")),
        subjects = "s1",
        contrasts = "c1"
      ),
      out_dir = td
    ),
    "voxel space"
  )
})

# ---------------------------------------------------------------------------
# summary / validate / explain branches
# ---------------------------------------------------------------------------

test_that("explain and validate cover parcel/sample/plan node summaries", {
  g_parcel <- new_gds(
    assays = list(beta = array(1:6, c(3, 2, 1)), var = array(1, c(3, 2, 1))),
    space = space_parcels(c("p1", "p2", "p3")),
    subjects = c("s1", "s2", "s3", "s4")[1:2],
    contrasts = "c1",
    metadata = list(
      map_families = list(id = "x"),
      design_mats = list(list(columns = c("(Intercept)", "group")))
    )
  )
  out <- capture.output(explain(g_parcel))
  expect_true(any(grepl("parcels", out)))
  expect_true(any(grepl("design\\[1\\]", out)))

  g_lab <- new_gds(
    assays = list(beta = array(1:4, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  expect_true(any(grepl("sample_labels", capture.output(explain(g_lab)))))

  plan <- as_plan(g_lab) |>
    derive(c("t", "z")) |>
    reduce("fixed") |>
    mask(MaskPolicy()) |>
    write_out(tempfile(fileext = ".h5"), format = "h5")
  pout <- capture.output(explain(plan))
  expect_true(any(grepl("reduce:fixed", pout)))
  expect_true(any(grepl("write:h5", pout)))
  expect_true(any(grepl("mask:", pout)))
  capture.output(explain_plan(plan))

  expect_true(validate(plan))
  bad_plan <- plan
  bad_plan$nodes[[length(bad_plan$nodes)]]$format <- "nope"
  expect_error(validate(bad_plan), "Unsupported write format")
  bad_map <- plan
  bad_map$nodes <- list(list(op = "map", target_space = NULL, map = matrix(1)))
  expect_error(validate(bad_map), "target_space")
  bad_ph <- plan
  bad_ph$nodes <- list(list(op = "posthoc", method = "not-a-method"))
  expect_error(validate(bad_ph), "Unknown post-hoc")
})
# ---------------------------------------------------------------------------
# gds-class validation helpers and space.default
# ---------------------------------------------------------------------------

test_that("new_gds validation and space.default error paths", {
  expect_error(fmrigds:::.validate_gds_assays(list()), "non-empty named list")
  expect_error(fmrigds:::.validate_gds_assays(list(1)), "named list")
  expect_error(
    fmrigds:::.validate_gds_assays(list(beta = array(1, c(2, 2, 1)), var = array(1, c(2, 1, 1)))),
    "identical dimensions"
  )
  expect_error(
    fmrigds:::.validate_gds_dims(
      list(beta = array(1, c(2, 2, 1))),
      subjects = "s1",
      contrasts = "c1"
    ),
    "subjects"
  )
  expect_error(
    fmrigds:::.validate_gds_dims(
      list(beta = array(1, c(2, 2, 1))),
      subjects = c("s1", "s2"),
      contrasts = c("c1", "c2")
    ),
    "contrasts"
  )
  expect_error(fmrigds:::.validate_gds_space(list(), 2), "gds_space")
  packed <- space_voxel(dim = c(2, 2, 1), affine = diag(4), mask_idx = 1:2, storage = "packed")
  expect_error(fmrigds:::.validate_gds_space(packed, 3), "mask_idx length")

  expect_error(fmrigds:::.normalise_col_data("x", "s1"), "data.frame")
  expect_error(
    fmrigds:::.normalise_col_data(data.frame(a = 1, row.names = "z"), "s1"),
    "rownames matching"
  )
  expect_error(fmrigds:::.normalise_row_data("x", 2), "data.frame")
  expect_error(fmrigds:::.normalise_row_data(data.frame(a = 1), 2), "one row per sample")

  expect_null(fmrigds:::.normalise_contrast_data(NULL, "c1"))
  expect_error(fmrigds:::.normalise_contrast_data("x", "c1"), "data.frame")
  expect_error(
    fmrigds:::.normalise_contrast_data(data.frame(a = 1), character()),
    "no contrasts"
  )
  expect_error(
    fmrigds:::.normalise_contrast_data(data.frame(a = 1, row.names = "z"), "c1"),
    "rownames matching"
  )
  # Duplicated rownames force the positional assignment path.
  positional <- data.frame(a = 1:2, row.names = c("a", "b"))
  attr(positional, "row.names") <- c("dup", "dup")
  cd <- fmrigds:::.normalise_contrast_data(positional, c("c1", "c2"))
  expect_equal(rownames(cd), c("c1", "c2"))

  expect_error(space(1L), "No applicable 'space\\(\\)' method|neuroim2")
})

# ---------------------------------------------------------------------------
# gds-verb alignment helpers
# ---------------------------------------------------------------------------

test_that("col/row/contrast data alignment helpers cover guards", {
  expect_error(fmrigds:::.align_col_data_for_subjects(1, "s1"), "data.frame")
  expect_error(
    fmrigds:::.align_col_data_for_subjects(data.frame(a = 1), "s1"),
    "missing subjects"
  )
  na_rows <- data.frame(a = 1)
  attr(na_rows, "row.names") <- NA_integer_
  expect_error(
    fmrigds:::.align_col_data_for_subjects(na_rows, "s1"),
    "rownames matching subjects"
  )
  dup_rows <- data.frame(a = 1:2, row.names = c("s1", "s2"))
  attr(dup_rows, "row.names") <- c("s1", "s1")
  expect_error(
    fmrigds:::.align_col_data_for_subjects(dup_rows, "s1"),
    "unique"
  )
  expect_error(
    fmrigds:::.align_col_data_for_subjects(
      data.frame(a = 1, row.names = "s2"),
      "s1"
    ),
    "missing subjects"
  )
  expect_warning(
    fmrigds:::.align_col_data_for_subjects(
      data.frame(a = 1:2, row.names = c("s1", "extra")),
      "s1",
      warn_extra = TRUE
    ),
    "Dropping extra rows"
  )

  expect_error(fmrigds:::.align_row_data_for_samples(1, "a"), "data.frame")
  expect_error(
    fmrigds:::.align_row_data_for_samples(data.frame(a = 1:2), NULL, n_samples = 3),
    "one row per sample"
  )
  expect_error(
    fmrigds:::.align_row_data_for_samples(
      data.frame(a = 1, row.names = "b"),
      "a"
    ),
    "missing samples"
  )
  expect_warning(
    fmrigds:::.align_row_data_for_samples(
      data.frame(a = 1:2, row.names = c("a", "extra")),
      "a",
      warn_extra = TRUE
    ),
    "Dropping extra rows"
  )
  expect_error(
    fmrigds:::.align_row_data_for_samples(data.frame(a = 1:2), c("a", "b", "c")),
    "one row per sample"
  )

  expect_error(fmrigds:::.align_contrast_data_for_contrasts(1, "c1"), "data.frame")
  expect_error(
    fmrigds:::.align_contrast_data_for_contrasts(data.frame(a = 1), "c1"),
    "missing contrasts"
  )
  na_con <- data.frame(a = 1)
  attr(na_con, "row.names") <- NA_integer_
  expect_error(
    fmrigds:::.align_contrast_data_for_contrasts(na_con, "c1"),
    "rownames matching contrasts"
  )
  dup_con <- data.frame(a = 1:2, row.names = c("c1", "c2"))
  attr(dup_con, "row.names") <- c("c1", "c1")
  expect_error(
    fmrigds:::.align_contrast_data_for_contrasts(dup_con, "c1"),
    "unique"
  )
  expect_error(
    fmrigds:::.align_contrast_data_for_contrasts(
      data.frame(a = 1, row.names = "c2"),
      "c1"
    ),
    "missing contrasts"
  )
  expect_warning(
    fmrigds:::.align_contrast_data_for_contrasts(
      data.frame(a = 1:2, row.names = c("c1", "extra")),
      "c1",
      warn_extra = TRUE
    ),
    "Dropping extra"
  )
})

# ---------------------------------------------------------------------------
# catalog print/summary/validation issue reporting
# ---------------------------------------------------------------------------

test_that("image_catalog print/summary cover large subject and assay-map branches", {
  subjects <- sprintf("sub-%02d", 1:8)
  files <- paste0("/tmp/", subjects, "_cope1.nii.gz")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = subjects,
    contrast = rep(c("A", "B"), length.out = 8),
    stringsAsFactors = FALSE
  )
  cat_obj <- new_image_catalog(
    files,
    metadata = meta,
    root_dir = "/tmp",
    pattern = "*.nii.gz",
    assay_map = list(
      beta = "cope",
      var = list(type = "varcope")
    )
  )
  pout <- capture.output(print(cat_obj))
  expect_true(any(grepl("Subjects: 8", pout)))
  expect_true(any(grepl("Assay mappings", pout)))
  sout <- capture.output(summary(cat_obj))
  expect_true(any(grepl("Assay mappings", sout)))
  expect_true(any(grepl("pattern 'cope'", sout)))

  # Discovery fallback via list.files when Sys.glob misses
  root <- tempfile("catalog-cov-")
  dir.create(file.path(root, "nested"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  file.create(file.path(root, "nested", "stat.nii.gz"))
  discovered <- fmrigds:::.catalog_discover_files(root, pattern = "nested/*.nii.gz", recursive = TRUE)
  expect_true(length(discovered) >= 1)
})

# ---------------------------------------------------------------------------
# adapter-tabular guards and LMM registration re-entry
# ---------------------------------------------------------------------------

test_that("tabular adapter rejects incomplete/missing inputs", {
  expect_error(fmrigds:::.tabular_open("/no/such/file.csv"), "does not exist")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("a,b\n1,2", tmp)
  handle <- fmrigds:::.tabular_open(tmp)
  expect_error(fmrigds:::.tabular_probe(handle), "Missing required columns|No effect columns")

  register_tabular_adapter()
  register_lmm_reducers()
  expect_true(!is.null(get_reducer("lmm:ri")))
  expect_true(!is.null(get_reducer("lmm:ri_slope1_knownvar")))
})

# ---------------------------------------------------------------------------
# neuroim mask helper without requiring neuroim2 objects
# ---------------------------------------------------------------------------

test_that(".resolve_neuroim_mask covers dense/packed/logical/invalid paths", {
  arr <- array(c(0, 1, 0, 2), dim = c(2, 2, 1))
  packed <- fmrigds:::.resolve_neuroim_mask(arr, mask = NULL, vdim = dim(arr))
  expect_equal(packed$storage, "packed")
  zeros <- fmrigds:::.resolve_neuroim_mask(array(0, dim(arr)), mask = NULL, vdim = dim(arr))
  expect_equal(zeros$storage, "dense")
  none <- fmrigds:::.resolve_neuroim_mask(arr, mask = "none", vdim = dim(arr))
  expect_null(none$mask_idx)
  logical_mask <- array(c(TRUE, FALSE, TRUE, FALSE), dim = dim(arr))
  lm <- fmrigds:::.resolve_neuroim_mask(arr, mask = logical_mask, vdim = dim(arr))
  expect_equal(lm$storage, "packed")
  expect_error(
    fmrigds:::.resolve_neuroim_mask(arr, mask = array(TRUE, c(1, 1, 1)), vdim = dim(arr)),
    "Mask dimensions"
  )
  expect_error(fmrigds:::.resolve_neuroim_mask(arr, mask = 1L, vdim = dim(arr)), "Invalid mask")
})
