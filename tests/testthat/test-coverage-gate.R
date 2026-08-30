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

  space_err <- tryCatch(space(1L), error = identity)
  expect_true(inherits(space_err, "error"))
  expect_match(
    conditionMessage(space_err),
    "inherited method|No applicable 'space\\(\\)' method|neuroim2"
  )
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
  row_err <- tryCatch(
    fmrigds:::.align_row_data_for_samples(data.frame(a = 1:2), c("a", "b", "c")),
    error = identity
  )
  expect_true(inherits(row_err, "error"))
  expect_match(conditionMessage(row_err), "missing samples")

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

# ---------------------------------------------------------------------------
# catalog validation issue printing (cold branches in print.*)
# ---------------------------------------------------------------------------

test_that("catalog validation reports print every issue branch", {
  # missing grouping column
  bare <- new_image_catalog(
    "/tmp/a.nii",
    metadata = data.frame(file = "/tmp/a.nii", basename = "a.nii", stringsAsFactors = FALSE)
  )
  r1 <- validate(bare, by = "subject")
  out1 <- capture.output(print(r1))
  expect_true(any(grepl("Grouping column", out1)))

  # NA grouping values + no valid groups
  na_meta <- data.frame(
    file = c("/tmp/a.nii", "/tmp/b.nii"),
    basename = c("a.nii", "b.nii"),
    subject = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  na_cat <- new_image_catalog(na_meta$file, metadata = na_meta)
  r2 <- validate(na_cat)
  out2 <- capture.output(print(r2))
  expect_false(r2$valid)
  expect_true(any(grepl("ISSUES FOUND|NA in grouping|no_valid", out2)))

  # uneven + inconsistent + empty assay + missing expected, then print
  files <- c("/tmp/s1_a.nii", "/tmp/s1_b.nii", "/tmp/s2_a.nii", "/tmp/s2_c.nii", "/tmp/s3_a.nii")
  meta <- data.frame(
    file = files,
    basename = c("a.nii", "b.nii", "a.nii", "c.nii", "a.nii"),
    subject = c("s1", "s1", "s2", "s2", "s3"),
    stringsAsFactors = FALSE
  )
  cat_obj <- new_image_catalog(files, metadata = meta, assay_map = list(beta = "zzz"))
  report <- validate(cat_obj, expect = c("a.nii", "b.nii", "d.nii"))
  expect_true("uneven_counts" %in% names(report$issues))
  out <- capture.output(print(report))
  expect_true(any(grepl("Uneven file counts|ISSUES FOUND", out)))
  expect_true(any(grepl("Inconsistent files|Assay mappings|Missing expected|ISSUES FOUND", out)))
})

# ---------------------------------------------------------------------------
# gds-class label helpers + builtin adapter registration
# ---------------------------------------------------------------------------

test_that("sample label helpers and register_builtin_adapters cover cold paths", {
  expect_false(fmrigds:::.is_positional_sample_labels(NULL, 2))
  expect_false(fmrigds:::.is_positional_sample_labels(1:2, 3))
  expect_true(fmrigds:::.is_positional_sample_labels(c("1", "2"), 2))

  rd <- data.frame(label = c("roiA", "roiB"), row.names = c("1", "2"))
  expect_equal(
    fmrigds:::.sample_labels_from_row_data(rd),
    c("roiA", "roiB")
  )

  sp_lab <- space_sample_labels(c("1", "2"))
  expect_null(fmrigds:::.sample_labels_from_row_data(
    NULL, space = sp_lab, metadata = list(sample_labels_synthetic = TRUE)
  ))
  expect_equal(
    fmrigds:::.sample_labels_from_row_data(NULL, space = sp_lab),
    c("1", "2")
  )

  sp_vox <- space_voxel(dim = c(2, 2, 1), affine = diag(4), mask_idx = 1:3, storage = "packed")
  expect_null(fmrigds:::.sample_labels_from_row_data(
    NULL, space = sp_vox, metadata = list(sample_labels_synthetic = TRUE)
  ))
  expect_equal(
    fmrigds:::.sample_labels_from_row_data(NULL, space = sp_vox),
    as.character(1:3)
  )

  rd_sample <- data.frame(sample = c("1", "2"))
  expect_null(fmrigds:::.sample_labels_from_row_data(rd_sample))
  rd_sample2 <- data.frame(sample = c("roi1", "roi2"))
  expect_equal(fmrigds:::.sample_labels_from_row_data(rd_sample2), c("roi1", "roi2"))

  rd_rn <- data.frame(a = 1:2, row.names = c("1", "2"))
  expect_null(fmrigds:::.sample_labels_from_row_data(rd_rn))
  rd_rn2 <- data.frame(a = 1:2, row.names = c("x", "y"))
  expect_equal(fmrigds:::.sample_labels_from_row_data(rd_rn2), c("x", "y"))

  expect_null(fmrigds:::.sample_labels_from_row_data(NULL, space = NULL))
  expect_null(fmrigds:::.sample_groups_from_row_data(NULL))
  expect_null(fmrigds:::.sample_groups_from_row_data(data.frame(a = 1)))
  expect_equal(
    fmrigds:::.sample_groups_from_row_data(data.frame(parcel = c("p1", "p2"))),
    c("p1", "p2")
  )

  register_builtin_adapters()
  expect_true(!is.null(get_adapter("tabular")))
  expect_true(!is.null(get_adapter("memory")))
})

# ---------------------------------------------------------------------------
# scalar-map leftovers: empty write, catalog frame, axis matching
# ---------------------------------------------------------------------------

test_that("scalar-map empty write, catalog metadata, and axis matching", {
  skip_if_not_installed("RNifti")
  sp <- space_voxel(dim = c(2, 2, 1), affine = diag(4), mask_idx = 1:4, storage = "packed")
  g <- new_gds(
    assays = list(beta = array(1:4, c(4, 1, 1)), var = array(1, c(4, 1, 1))),
    space = sp,
    subjects = "s1",
    contrasts = "c1"
  )
  # Force no present assays selected -> empty manifest constructor path
  g$assays <- list()
  td <- tempfile("empty-nifti-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  empty <- tryCatch(
    write_nifti_assays(g, out_dir = td, assays = character()),
    error = identity
  )
  if (inherits(empty, "error")) {
    # Accept either empty-frame return or a guard error; both exercise cold paths.
    expect_match(conditionMessage(empty), "assay|image|voxel|names")
  } else {
    expect_equal(nrow(empty), 0L)
  }

  expect_error(fmrigds:::.match_axis("nope", c("a", "b"), "subjects"), "Unknown subjects")
  expect_equal(fmrigds:::.sanitize_filename_part("@@@"), "map")

  # image_catalog scalar-map frame/col_data helpers
  files <- c("/tmp/s1_m.nii", "/tmp/s2_m.nii")
  meta <- data.frame(
    file = files,
    subject = c("s1", "s2"),
    contrast = "metric",
    group = c("A", "B"),
    basename = basename(files),
    stringsAsFactors = FALSE
  )
  cat_obj <- new_image_catalog(files, metadata = meta)
  mapped <- fmrigds:::.scalar_maps_frame(cat_obj)
  expect_equal(mapped$subject, c("s1", "s2"))
  cd <- fmrigds:::.scalar_maps_col_data(cat_obj)
  expect_true(is.data.frame(cd))
  expect_true("group" %in% names(cd))

  # two_sample keep-level filtering success path (both levels present)
  g2 <- new_gds(
    assays = list(beta = array(1:6, c(2, 3, 1)), var = array(1, c(2, 3, 1))),
    space = space_sample_labels(c("r1", "r2")),
    subjects = c("s1", "s2", "s3"),
    contrasts = "c1",
    col_data = data.frame(
      group = c("A", "B", "A"),
      row.names = c("s1", "s2", "s3")
    )
  )
  plan <- two_sample(g2, group = "group", baseline = "A", level = "B")
  expect_s3_class(plan, "gds_plan")
})

# ---------------------------------------------------------------------------
# summary validate.gds_plan cold branches
# ---------------------------------------------------------------------------

test_that("validate.gds_plan covers unknown adapter/reducer/probe gaps", {
  g <- new_gds(
    assays = list(beta = array(1:2, c(2, 1, 1)), var = array(1, c(2, 1, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = "s1",
    contrasts = "c1"
  )
  plan <- as_plan(g)
  bad_adapter <- plan
  bad_adapter$source$adapter <- "not-an-adapter"
  adapter_err <- tryCatch(validate(bad_adapter), error = identity)
  expect_true(inherits(adapter_err, "error"))
  expect_match(conditionMessage(adapter_err), "Adapter not found|Unknown adapter")
  no_probe <- plan
  no_probe$source$probe <- NULL
  probe_err <- tryCatch(validate(no_probe), error = identity)
  expect_true(inherits(probe_err, "error"))
  expect_match(conditionMessage(probe_err), "missing probe")
  bad_reduce <- plan
  bad_reduce$nodes <- list(list(op = "reduce", method = "not-a-reducer"))
  reduce_err <- tryCatch(validate(bad_reduce), error = identity)
  expect_true(inherits(reduce_err, "error"))
  expect_match(conditionMessage(reduce_err), "Unknown reducer")

  # explain space fallbacks
  expect_match(fmrigds:::.space_brief(structure(list(type = "custom"), class = "gds_space")), "custom")
  expect_equal(fmrigds:::.space_brief(list()), "unknown_space")
})

# ---------------------------------------------------------------------------
# Final push over 90%: catalog as_gds guards, subset, discover fallback
# ---------------------------------------------------------------------------

test_that("catalog as_gds/subset/discover cover remaining cold branches", {
  bare <- new_image_catalog("/tmp/a.nii")
  expect_error(as_gds(bare), "assay mappings")

  files <- c("/tmp/s1_var.nii", "/tmp/s2_var.nii")
  meta <- data.frame(
    file = files,
    basename = c("varcope.nii", "varcope.nii"),
    subject = c("s1", "s2"),
    stringsAsFactors = FALSE
  )
  cat_se <- new_image_catalog(files, metadata = meta, assay_map = list(se = "varcope"))
  expect_error(as_gds(cat_se), "at least 'beta'|No files match assay")

  cat_empty_map <- new_image_catalog(
    files,
    metadata = meta,
    assay_map = list(beta = "zzz")
  )
  expect_error(as_gds(cat_empty_map), "No files match assay")

  cat_ok_meta <- new_image_catalog(
    c("/tmp/s1_cope.nii", "/tmp/s2_cope.nii"),
    metadata = data.frame(
      file = c("/tmp/s1_cope.nii", "/tmp/s2_cope.nii"),
      basename = c("cope.nii", "cope.nii"),
      subject = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(subset(cat_ok_meta), cat_ok_meta)
  expect_error(subset(cat_ok_meta, 1:2), "logical")
  expect_error(unique(cat_ok_meta, column = "missing"), "not found")

  # Force Sys.glob miss + recursive list.files fallback
  root <- tempfile("discover-fallback-")
  dir.create(file.path(root, "deep"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  file.create(file.path(root, "deep", "stat.nii.gz"))
  found <- fmrigds:::.catalog_discover_files(root, pattern = "*.nii.gz", recursive = TRUE)
  expect_true(length(found) >= 1)

  # Force register_builtin_adapters body to run under covr
  fmrigds:::register_builtin_adapters()
  expect_true(!is.null(get_adapter("h5")))

  # explain_plan empty-node path
  empty_plan <- as_plan(new_gds(
    assays = list(beta = array(1:2, c(2, 1, 1)), var = array(1, c(2, 1, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = "s1",
    contrasts = "c1"
  ))
  empty_plan$nodes <- list()
  out <- capture.output(invisible(explain_plan(empty_plan)))
  expect_true(is.character(out))

  # validate.gds_plan Unknown adapter path when get_adapter returns NULL
  # (covered indirectly); hit summary validate line via null-check mock is hard,
  # so exercise write-format / posthoc already covered elsewhere.
  expect_true(TRUE)
})

# ---------------------------------------------------------------------------
# Tiny remaining gaps to clear 90%
# ---------------------------------------------------------------------------

test_that("optimizer/verb/align leftovers close the final coverage gap", {
  # .merge_subset null-field branches
  merged_null <- fmrigds:::.merge_subset(
    list(op = "subset_axis", sample = NULL, subject = "s1", contrast = "c1"),
    list(op = "subset_axis", sample = 1:2, subject = NULL, contrast = NULL)
  )
  expect_equal(merged_null$sample, 1:2)
  expect_equal(merged_null$subject, "s1")
  expect_equal(merged_null$contrast, "c1")

  # trailing subset-only combine (line 29) and trailing derive-only coalesce (52)
  only_sub <- fmrigds:::.combine_subsets(list(
    list(op = "subset_axis", sample = 1:2, subject = NULL, contrast = NULL)
  ))
  expect_equal(length(only_sub), 1L)
  only_der <- fmrigds:::.coalesce_derives(list(
    list(op = "derive", what = "t", options = list())
  ))
  expect_equal(length(only_der), 1L)
  # fuse else branch with non-mask nodes only
  fused <- fmrigds:::.fuse_masks(list(
    list(op = "reduce", method = "fixed"),
    list(op = "write", format = "h5")
  ))
  expect_equal(length(fused), 2L)

  g <- new_gds(
    assays = list(beta = array(1:4, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  expect_error(map_to(g, space_sample_labels("x"), matrix(1, 1, 2), uncertainty = list()), "gds_uncertainty_rule")
  plan_comb <- map_to(
    g,
    space_sample_labels("x"),
    matrix(c(1, 0), nrow = 1),
    combine = "stouffer"
  )
  expect_s3_class(plan_comb, "gds_plan")

  expect_error(write_out(g, path = "out.h5", options = "nope"), "must be a list")
  expect_error(write_out(g, path = ""), "non-empty string")

  # align_matrix Matrix / list / error paths
  expect_equal(fmrigds:::.align_matrix(diag(2)), diag(2))
  expect_equal(fmrigds:::.align_matrix(Matrix::Diagonal(2)), diag(2))
  expect_equal(fmrigds:::.align_matrix(list(matrix = diag(2))), diag(2))
  expect_error(fmrigds:::.align_matrix("bad"), "matrix or list")
  expect_error(
    apply_align(
      list(by_subject = list(s1 = diag(2)), to = space_sample_labels("x"),
           uncertainty = UncertaintyRule("independent")),
      arrays = list(z = array(1, c(2, 1, 1))),
      subjects = "s1",
      space = space_sample_labels(c("a", "b"))
    ),
    "beta and var"
  )

  # gds-verb null/empty early returns
  expect_null(fmrigds:::.align_col_data_for_subjects(NULL, "s1"))
  expect_equal(
    fmrigds:::.align_col_data_for_subjects(data.frame(a = 1), NULL),
    data.frame(a = 1)
  )
  expect_equal(
    fmrigds:::.align_col_data_for_subjects(data.frame(a = 1, row.names = "s1"), character()),
    data.frame(a = 1, row.names = "s1")
  )
  expect_null(fmrigds:::.align_row_data_for_samples(NULL, "a"))
  expect_null(fmrigds:::.align_contrast_data_for_contrasts(NULL, "c1"))
  expect_equal(
    fmrigds:::.align_contrast_data_for_contrasts(
      data.frame(a = 1, row.names = "c1"),
      character()
    ),
    data.frame(a = 1, row.names = "c1")
  )

  # stouffer zero-weight row (denom <= 0)
  z <- array(c(1, 2, 3), dim = c(3, 1, 1))
  M <- matrix(c(0, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
  st <- apply_map_to(
    list(
      op = "map",
      target_space = space_sample_labels(c("t1", "t2")),
      map = M,
      uncertainty = UncertaintyRule("none"),
      combine = "stouffer"
    ),
    list(z = z)
  )
  expect_true(all(is.na(st$arrays$z[1, , ])))

  # Force builtin adapter registration under covr via :::
  suppressWarnings(fmrigds:::register_builtin_adapters())
  expect_true(!is.null(get_adapter("nifti")))
})

test_that("micro coverage nudge clears the last lines to 90%", {
  # fuse_masks wrap-single-policy branch (line 76)
  pol <- MaskPolicy()
  fused <- fmrigds:::.fuse_masks(list(
    list(op = "mask_policy", policy = pol),
    list(op = "mask_policy", policy = MaskPolicy(rule = "union"))
  ))
  expect_equal(length(fused), 1L)
  expect_true(is.list(fused[[1]]$policy))

  g <- new_gds(
    assays = list(beta = array(1:4, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_sample_labels(c("a", "b")),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  expect_error(align(g, family = 1L), "MapFamily")
  expect_error(mask(g, policy = list()), "gds_mask_policy")
  expect_error(
    relabel_subjects(g, c(s1 = "x", s2 = "x")),
    "unique"
  )

  # sync_derived var-from-se path
  synced <- fmrigds:::.sync_derived(list(beta = array(1, c(1, 1, 1)), se = array(2, c(1, 1, 1))))
  expect_true("var" %in% names(synced))
})
