# Tests for image_catalog functionality

# Helper to create test directory structure
create_test_files <- function(root) {
  # Create directory structure: sub-01/analysis/, sub-02/analysis/, sub-03/analysis/
  for (subj in c("sub-01", "sub-02", "sub-03")) {
    dir.create(file.path(root, subj, "analysis"), recursive = TRUE)
    # Create mock NIfTI files (just empty files for testing)
    file.create(file.path(root, subj, "analysis", "cope1.nii.gz"))
    file.create(file.path(root, subj, "analysis", "cope2.nii.gz"))
    file.create(file.path(root, subj, "analysis", "varcope1.nii.gz"))
    file.create(file.path(root, subj, "analysis", "varcope2.nii.gz"))
  }
}

# =============================================================================
# Constructor Tests
# =============================================================================

test_that("new_image_catalog creates valid object", {
  files <- c("/tmp/a.nii", "/tmp/b.nii")
  cat <- new_image_catalog(files)

  expect_s3_class(cat, "image_catalog")
  expect_equal(cat$files, files)
  expect_true(is.data.frame(cat$metadata))
  expect_equal(nrow(cat$metadata), 2L)
  expect_true("file" %in% names(cat$metadata))
})

test_that("new_image_catalog accepts custom metadata", {
  files <- c("/tmp/a.nii", "/tmp/b.nii")
  meta <- data.frame(
    file = files,
    subject = c("01", "02"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  expect_equal(cat$metadata$subject, c("01", "02"))
})

test_that("print.image_catalog works", {
  files <- c("/tmp/a.nii", "/tmp/b.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta, root_dir = "/tmp", pattern = "*.nii")

  output <- capture.output(print(cat))
  expect_true(any(grepl("Image Catalog", output)))
  expect_true(any(grepl("Files: 2", output)))
  expect_true(any(grepl("Subjects: 2", output)))
})

test_that("summary.image_catalog works", {
  files <- c("/tmp/a.nii", "/tmp/b.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02"),
    contrast = c("A", "B"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  output <- capture.output(summary(cat))
  expect_true(any(grepl("Image Catalog Summary", output)))
  expect_true(any(grepl("Total files: 2", output)))
  expect_true(any(grepl("subject:", output)))
  expect_true(any(grepl("contrast:", output)))
})

# =============================================================================
# Discovery Tests
# =============================================================================

test_that("image_catalog discovers files with glob pattern", {
  skip_on_cran()

  root <- tempfile("catalog_test_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  create_test_files(root)

  cat <- image_catalog(root, pattern = "sub-*/analysis/*.nii.gz")

  expect_s3_class(cat, "image_catalog")
  expect_equal(length(cat$files), 12L)  # 3 subjects * 4 files
  expect_true(all(grepl("\\.nii\\.gz$", cat$files)))
})

test_that("image_catalog extracts metadata via path_regex", {
  skip_on_cran()

  root <- tempfile("catalog_test_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  create_test_files(root)

  cat <- image_catalog(
    root,
    pattern = "sub-*/analysis/*.nii.gz",
    path_regex = "sub-(?<subject>[0-9]+)"
  )

  expect_true("subject" %in% names(cat$metadata))
  expect_setequal(unique(cat$metadata$subject), c("01", "02", "03"))
})

test_that("image_catalog warns when no files found", {
  skip_on_cran()

  root <- tempfile("catalog_test_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_warning(
    cat <- image_catalog(root, pattern = "nonexistent/*.nii"),
    "No files found"
  )
  expect_equal(length(cat$files), 0L)
})

test_that("image_catalog errors on missing root", {
  expect_error(
    image_catalog("/nonexistent/path"),
    "Root directory does not exist"
  )
})

# =============================================================================
# Metadata Assignment Tests
# =============================================================================

test_that("assign_meta extracts values via regex", {
  files <- c("/tmp/cope1.nii", "/tmp/cope2.nii", "/tmp/varcope1.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  cat2 <- assign_meta(cat, "contrast_num", "cope([0-9]+)")

  expect_true("contrast_num" %in% names(cat2$metadata))
  expect_equal(cat2$metadata$contrast_num, c("1", "2", "1"))
})

test_that("assign_meta handles non-matches as NA", {
  files <- c("/tmp/cope1.nii", "/tmp/other.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  cat2 <- assign_meta(cat, "cope_num", "cope([0-9]+)")

  expect_equal(cat2$metadata$cope_num[1], "1")
  expect_true(is.na(cat2$metadata$cope_num[2]))
})

test_that("assign_meta errors on missing source column", {
  cat <- new_image_catalog("/tmp/a.nii")

  expect_error(
    assign_meta(cat, "x", "pattern", source = "nonexistent"),
    "not found in catalog metadata"
  )
})

test_that("join_meta merges external data", {
  files <- c("/tmp/a.nii", "/tmp/b.nii", "/tmp/c.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02", "01"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  external <- data.frame(
    subject = c("01", "02"),
    age = c(25, 30),
    group = c("control", "patient"),
    stringsAsFactors = FALSE
  )

  cat2 <- join_meta(cat, external, by = "subject")

  expect_true("age" %in% names(cat2$metadata))
  expect_true("group" %in% names(cat2$metadata))
  expect_equal(cat2$metadata$age, c(25, 30, 25))
  expect_equal(cat2$metadata$group, c("control", "patient", "control"))
})

test_that("join_meta errors on missing key", {
  cat <- new_image_catalog("/tmp/a.nii")
  external <- data.frame(subject = "01", age = 25)

  expect_error(
    join_meta(cat, external, by = "subject"),
    "not found in catalog"
  )
})

# =============================================================================
# Assay Mapping Tests
# =============================================================================

test_that("map_assays stores pattern-based mappings", {
  cat <- new_image_catalog("/tmp/a.nii")

  cat2 <- map_assays(cat, beta = "cope", se = "varcope")

  expect_equal(names(cat2$assay_map), c("beta", "se"))
  expect_equal(cat2$assay_map$beta, "cope")
  expect_equal(cat2$assay_map$se, "varcope")
})

test_that("map_assays stores column-based mappings", {
  cat <- new_image_catalog("/tmp/a.nii")

  cat2 <- map_assays(cat,
    beta = list(stat_type = "effect"),
    se = list(stat_type = "stderr")
  )

  expect_true(is.list(cat2$assay_map$beta))
  expect_equal(cat2$assay_map$beta$stat_type, "effect")
})

test_that("map_assays warns on unknown assay types", {
  cat <- new_image_catalog("/tmp/a.nii")

  expect_warning(
    map_assays(cat, custom_assay = "pattern"),
    "Non-standard assay types"
  )
})

test_that(".resolve_assay_map works with pattern", {
  files <- c("/tmp/cope1.nii", "/tmp/cope2.nii", "/tmp/varcope1.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)
  cat <- map_assays(cat, beta = "^cope")

  idx <- fmrigds:::.resolve_assay_map(cat, "beta")
  expect_equal(idx, c(1L, 2L))
})

test_that(".resolve_assay_map works with column match", {
  files <- c("/tmp/a.nii", "/tmp/b.nii", "/tmp/c.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    stat_type = c("cope", "cope", "varcope"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)
  cat <- map_assays(cat, beta = list(stat_type = "cope"))

  idx <- fmrigds:::.resolve_assay_map(cat, "beta")
  expect_equal(idx, c(1L, 2L))
})

# =============================================================================
# Validation Tests
# =============================================================================

test_that("validate.image_catalog returns valid for consistent catalog", {
  files <- c("/tmp/s1_cope.nii", "/tmp/s1_var.nii", "/tmp/s2_cope.nii", "/tmp/s2_var.nii")
  meta <- data.frame(
    file = files,
    basename = c("cope.nii", "var.nii", "cope.nii", "var.nii"),
    subject = c("s1", "s1", "s2", "s2"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  report <- validate(cat)

  expect_true(report$valid)
  expect_equal(length(report$issues), 0L)
})

test_that("validate.image_catalog detects uneven counts", {
  files <- c("/tmp/s1_a.nii", "/tmp/s1_b.nii", "/tmp/s2_a.nii")
  meta <- data.frame(
    file = files,
    basename = c("a.nii", "b.nii", "a.nii"),
    subject = c("s1", "s1", "s2"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  report <- validate(cat)

  expect_false(report$valid)
  expect_true("uneven_counts" %in% names(report$issues))
})

test_that("validate.image_catalog detects inconsistent files", {
  files <- c("/tmp/s1_a.nii", "/tmp/s1_b.nii", "/tmp/s2_a.nii", "/tmp/s2_c.nii")
  meta <- data.frame(
    file = files,
    basename = c("a.nii", "b.nii", "a.nii", "c.nii"),
    subject = c("s1", "s1", "s2", "s2"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  report <- validate(cat)

  expect_false(report$valid)
  expect_true("inconsistent_files" %in% names(report$issues))
  expect_true("s2" %in% names(report$issues$inconsistent_files))
})

test_that("validate.image_catalog detects empty assay mappings", {
  files <- c("/tmp/cope.nii")
  meta <- data.frame(
    file = files,
    basename = "cope.nii",
    subject = "s1",
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)
  cat <- map_assays(cat, beta = "cope", se = "varcope")  # se won't match

  report <- validate(cat)

  expect_false(report$valid)
  expect_true("empty_assay_map" %in% names(report$issues))
  expect_true("se" %in% report$issues$empty_assay_map)
})

test_that("validate.image_catalog checks expected files", {
  files <- c("/tmp/s1_a.nii", "/tmp/s2_a.nii")
  meta <- data.frame(
    file = files,
    basename = c("a.nii", "a.nii"),
    subject = c("s1", "s2"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  report <- validate(cat, expect = c("a.nii", "b.nii"))

  expect_false(report$valid)
  expect_true("missing_expected" %in% names(report$issues))
})

test_that("print.catalog_validation_report works", {
  files <- c("/tmp/s1_a.nii", "/tmp/s2_a.nii")
  meta <- data.frame(
    file = files,
    basename = c("a.nii", "a.nii"),
    subject = c("s1", "s2"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)
  report <- validate(cat)

  output <- capture.output(print(report))
  expect_true(any(grepl("Catalog Validation Report", output)))
  expect_true(any(grepl("Status: VALID", output)))
})

# =============================================================================
# Serialization Tests
# =============================================================================

test_that("write_catalog and read_catalog roundtrip", {
  skip_on_cran()

  files <- c("/tmp/a.nii", "/tmp/b.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02"),
    age = c(25, 30),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta, root_dir = "/tmp", pattern = "*.nii")
  cat <- map_assays(cat, beta = "a", se = "b")

  tmp_file <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_file), add = TRUE)

  write_catalog(cat, tmp_file)
  expect_true(file.exists(tmp_file))

  cat2 <- read_catalog(tmp_file)

  expect_s3_class(cat2, "image_catalog")
  expect_equal(cat2$files, cat$files)
  expect_equal(cat2$root_dir, cat$root_dir)
  expect_equal(cat2$pattern, cat$pattern)
  expect_equal(cat2$metadata$subject, cat$metadata$subject)
  expect_equal(cat2$metadata$age, cat$metadata$age)
  expect_equal(names(cat2$assay_map), names(cat$assay_map))
})

test_that("read_catalog errors on missing file", {
  expect_error(
    read_catalog("/nonexistent/file.json"),
    "does not exist"
  )
})

# =============================================================================
# Subsetting Tests
# =============================================================================

test_that("subset.image_catalog filters by condition", {
  files <- c("/tmp/a.nii", "/tmp/b.nii", "/tmp/c.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02", "01"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  sub_cat <- subset(cat, subject == "01")

  expect_equal(length(sub_cat$files), 2L)
  expect_true(all(sub_cat$metadata$subject == "01"))
})

test_that("subset.image_catalog handles NA correctly", {
  files <- c("/tmp/a.nii", "/tmp/b.nii", "/tmp/c.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", NA, "02"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  sub_cat <- subset(cat, subject == "01")

  expect_equal(length(sub_cat$files), 1L)
})

test_that("subset.image_catalog preserves assay_map", {
  files <- c("/tmp/a.nii", "/tmp/b.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)
  cat <- map_assays(cat, beta = "a")

  sub_cat <- subset(cat, subject == "01")

  expect_equal(sub_cat$assay_map, cat$assay_map)
})

test_that("unique.image_catalog returns unique values", {
  files <- c("/tmp/a.nii", "/tmp/b.nii", "/tmp/c.nii")
  meta <- data.frame(
    file = files,
    basename = basename(files),
    subject = c("01", "02", "01"),
    stringsAsFactors = FALSE
  )
  cat <- new_image_catalog(files, metadata = meta)

  subjects <- unique(cat, column = "subject")
  expect_setequal(subjects, c("01", "02"))
})

test_that("unique.image_catalog errors on missing column", {
  cat <- new_image_catalog("/tmp/a.nii")

  expect_error(
    unique(cat, column = "nonexistent"),
    "not found in catalog metadata"
  )
})

# =============================================================================
# Internal Helper Tests
# =============================================================================

test_that(".extract_capture_group_names extracts names", {
  pattern <- "sub-(?<subject>[^/]+)/ses-(?<session>[^/]+)"
  names <- fmrigds:::.extract_capture_group_names(pattern)

  expect_equal(names, c("subject", "session"))
})

test_that(".extract_capture_group_names returns empty for no groups", {
  pattern <- "sub-[0-9]+"
  names <- fmrigds:::.extract_capture_group_names(pattern)

  expect_equal(names, character())
})

test_that(".glob_to_regex converts patterns correctly", {
  expect_equal(fmrigds:::.glob_to_regex("*.nii"), "^[^/]*\\.nii$")
  expect_equal(fmrigds:::.glob_to_regex("**.nii"), "^.*\\.nii$")
  expect_equal(fmrigds:::.glob_to_regex("file?.nii"), "^file.\\.nii$")
})
