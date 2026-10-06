# Regression tests for I/O review findings (HDF5, adapters, catalog, exporters,
# NIfTI/neuroim2 import, CLI). Each test fails on the pre-fix code.

.review_affine <- function() {
  matrix(c(
    2, 0, 0, 0,
    0, 2, 0, 0,
    0, 0, 2, 0,
    -90, -126, -72, 1
  ), nrow = 4)
}

.review_voxel_gds <- function(space, n) {
  new_gds(
    assays = list(beta = array(seq_len(n * 2), c(n, 2, 1)), var = array(1, c(n, 2, 1))),
    space = space,
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
}

.review_write_nifti <- function(path, dims, spacing = c(2, 2, 2), origin = c(0, 0, 0), fill = 1) {
  sp <- neuroim2::NeuroSpace(dims, spacing = spacing, origin = origin)
  neuroim2::write_vol(neuroim2::NeuroVol(array(fill, dims), sp), path)
  path
}

# 1. HDF5 affine ---------------------------------------------------------------

test_that("write_gds_h5 round-trips a non-symmetric voxel affine unchanged", {
  skip_if_not_installed("hdf5r")
  aff <- .review_affine()
  sp <- space_voxel(dim = c(2, 2, 2), affine = aff, mask_idx = c(1L, 3L, 8L), storage = "packed")
  g <- .review_voxel_gds(sp, 3L)
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  write_gds_h5(g, f)

  plan <- gds(f)
  expect_equal(plan$source$probe$space$affine, aff)
  expect_equal(plan$source$probe$space$affine[1:3, 4], c(-90, -126, -72))
  expect_equal(plan$source$probe$space$mask_idx, c(1L, 3L, 8L))
  expect_identical(plan$source$probe$space$storage, "packed")
})

# 2. Dense voxel storage ---------------------------------------------------------

test_that("dense voxel spaces read back from HDF5 as dense without mask_idx", {
  skip_if_not_installed("hdf5r")
  sp <- space_voxel(dim = c(2, 2, 2), affine = .review_affine(), storage = "dense")
  g <- .review_voxel_gds(sp, 8L)
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  write_gds_h5(g, f)

  plan <- gds(f)
  space <- plan$source$probe$space
  expect_identical(space$storage, "dense")
  expect_null(space$mask_idx)
  res <- compute(plan)
  expect_equal(unname(assay(res, "beta")), unname(assay(g, "beta")))
})

# 3. Unsupported space types ------------------------------------------------------

test_that("write_gds_h5 refuses space types it cannot serialise", {
  skip_if_not_installed("hdf5r")
  g <- new_gds(
    assays = list(beta = array(1, c(3, 2, 1)), var = array(1, c(3, 2, 1))),
    space = space_basis(k = 3, basis_name = "pca"),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  expect_error(write_gds_h5(g, f), "cannot serialise a 'basis' space")
  expect_false(file.exists(f))
})

test_that("an empty provenance log is written in a form that reads back", {
  skip_if_not_installed("hdf5r")
  g <- new_gds(
    assays = list(beta = array(1, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_parcels(c("a", "b")),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  write_gds_h5(g, f)
  h5 <- hdf5r::H5File$new(f, mode = "r")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)
  log <- h5[["gds/provenance/log"]]$read()
  expect_identical(log, character())
})

# 4. Adapter auto-detection -------------------------------------------------------

test_that("native write_gds_h5 files auto-detect as the h5 adapter", {
  skip_if_not_installed("hdf5r")
  g <- new_gds(
    assays = list(beta = array(1, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_parcels(c("a", "b")),
    subjects = c("s1", "s2"),
    contrasts = "c1"
  )
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  write_gds_h5(g, f)
  expect_identical(detect_adapter(f), "h5")
  expect_identical(gds(f)$source$adapter, "h5")
})

test_that("detect_adapter breaks score ties by registration order, not by name", {
  reg <- fmrigds:::.adapter_registry
  on.exit(rm(list = intersect(c("zz_review_first", "aa_review_second"), ls(reg)), envir = reg), add = TRUE)
  probe <- function(handle, ...) stop("unused")
  is_marker <- function(source) if (inherits(source, "review_tie_marker")) 1 else FALSE
  register_adapter("zz_review_first", detect = is_marker, open = identity,
                   probe = probe, read = probe, close = function(h) NULL)
  register_adapter("aa_review_second", detect = is_marker, open = identity,
                   probe = probe, read = probe, close = function(h) NULL)
  src <- structure(list(), class = "review_tie_marker")
  expect_identical(detect_adapter(src), "zz_review_first")
})

test_that("fmristore refuses multiple native /gds files instead of using the first", {
  skip_if_not_installed("hdf5r")
  g <- new_gds(
    assays = list(beta = array(1, c(2, 1, 1)), var = array(1, c(2, 1, 1))),
    space = space_parcels(c("a", "b")),
    subjects = "s1",
    contrasts = "c1"
  )
  f1 <- tempfile(fileext = ".h5")
  f2 <- tempfile(fileext = ".h5")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  write_gds_h5(g, f1)
  write_gds_h5(g, f2)
  handle <- fmrigds:::.fmri_open(c(f1, f2))
  on.exit(fmrigds:::.fmri_close(handle), add = TRUE, after = FALSE)
  expect_error(fmrigds:::.fmri_probe(handle), "complete GDS")
})

# 5. Handles closed on probe failure ----------------------------------------------

test_that("gds() closes the adapter handle when the probe fails", {
  reg <- fmrigds:::.adapter_registry
  on.exit(rm(list = intersect("review_failprobe", ls(reg)), envir = reg), add = TRUE)
  counter <- new.env()
  counter$closed <- 0L
  register_adapter(
    "review_failprobe",
    detect = function(source) FALSE,
    open = function(source, ...) list(src = source),
    probe = function(handle, ...) stop("boom in probe"),
    read = function(handle, ...) NULL,
    close = function(handle) {
      counter$closed <- counter$closed + 1L
      invisible(NULL)
    }
  )
  expect_error(gds("anything", format = "review_failprobe"), "boom in probe")
  expect_identical(counter$closed, 1L)
})

# 6. Image catalog -----------------------------------------------------------------

.review_fsl_tree <- function() {
  root <- tempfile("fsl_")
  for (s in c("sub-01", "sub-02")) {
    d <- file.path(root, s, "stats")
    dir.create(d, recursive = TRUE)
    .review_write_nifti(file.path(d, "cope1.nii.gz"), c(3, 3, 3), fill = 1)
    .review_write_nifti(file.path(d, "varcope1.nii.gz"), c(3, 3, 3), fill = 4)
  }
  root
}

test_that("image_catalog default '**' pattern finds nested files", {
  skip_if_not_installed("neuroim2")
  root <- .review_fsl_tree()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cat <- image_catalog(root, path_regex = "(?<subject>sub-[0-9]+)")
  expect_length(cat$files, 4L)
  expect_setequal(basename(cat$files), c("cope1.nii.gz", "varcope1.nii.gz"))
  expect_false(any(grepl("\\\\", cat$files)))
  expect_setequal(unique(cat$metadata$subject), c("sub-01", "sub-02"))

  # **/ also matches zero directory levels
  writeLines("x", file.path(root, "top.nii"))
  cat2 <- image_catalog(root, pattern = "**/*.nii")
  expect_true("top.nii" %in% basename(cat2$files))
})

test_that("as_gds(<image_catalog>) uses catalog subjects and pairs var by subject", {
  skip_if_not_installed("neuroim2")
  root <- .review_fsl_tree()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cat <- image_catalog(root, path_regex = "(?<subject>sub-[0-9]+)") |>
    map_assays(beta = "^cope", var = "^varcope")
  plan <- as_gds(cat)
  expect_identical(plan$source$probe$subjects, c("sub-01", "sub-02"))
  expect_true("var" %in% plan$source$probe$assays)
  expect_false(isTRUE(plan$source$probe$metadata$synthetic_var))
  res <- compute(plan)
  expect_true(all(assay(res, "var") == 4))
  expect_true(all(assay(res, "beta") == 1))
})

test_that("as_gds(<image_catalog>) errors on unsupported or overlapping mappings", {
  skip_if_not_installed("neuroim2")
  root <- .review_fsl_tree()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cat <- image_catalog(root, path_regex = "(?<subject>sub-[0-9]+)")
  expect_error(as_gds(map_assays(cat, beta = "^cope", t = "^varcope")), "unsupported")
  # Unanchored "cope" also matches varcope files
  expect_error(as_gds(map_assays(cat, beta = "cope", var = "varcope")), "overlap")
})

test_that("assign_meta scopes top-level alternations", {
  files <- c("/x/cope1.nii.gz", "/x/tstat3.nii.gz")
  cat <- new_image_catalog(files, metadata = data.frame(
    file = files, basename = basename(files), stringsAsFactors = FALSE
  ))
  cat <- assign_meta(cat, "num", "cope([0-9]+)|tstat([0-9]+)", replacement = "\\1\\2")
  expect_identical(cat$metadata$num, c("1", "3"))
})

test_that("join_meta errors on duplicate keys in the joined data", {
  files <- c("/x/a.nii", "/x/b.nii")
  cat <- new_image_catalog(files, metadata = data.frame(
    file = files, basename = basename(files), subject = c("01", "02"),
    stringsAsFactors = FALSE
  ))
  ext <- data.frame(subject = c("01", "01", "02"), age = c(20, 21, 30))
  expect_error(join_meta(cat, ext, by = "subject"), "duplicate")
})

# 7. Tidy export keeps column types ------------------------------------------------

test_that("factor col_data survives gds_to_tibble and CSV export", {
  g <- new_gds(
    assays = list(beta = array(1:4, c(2, 2, 1)), var = array(1, c(2, 2, 1))),
    space = space_parcels(c("r1", "r2")),
    subjects = c("s1", "s2"),
    contrasts = "c1",
    col_data = data.frame(
      group = factor(c("ctl", "pat")), age = c(30.5, 40),
      row.names = c("s1", "s2")
    )
  )
  tb <- gds_to_tibble(g, include_col_data = TRUE)
  expect_true(is.factor(tb$group))
  expect_identical(as.character(tb$group[tb$subject == "s2"]), c("pat", "pat"))
  expect_equal(tb$age[tb$subject == "s1"], c(30.5, 30.5))

  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  fmrigds:::.write_gds_csv(g, f, list(include_col_data = TRUE))
  back <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_setequal(unique(back$group), c("ctl", "pat"))
})

# 8. Per-file grid checks -------------------------------------------------------------

test_that("NIfTI sources with mismatched grids are rejected", {
  skip_if_not_installed("neuroim2")
  d <- tempfile("grid_")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f1 <- .review_write_nifti(file.path(d, "sub-01_beta.nii.gz"), c(3, 3, 3))
  f2 <- .review_write_nifti(file.path(d, "sub-02_beta.nii.gz"), c(4, 3, 3))
  expect_error(gds(c(f1, f2), format = "nifti"), "grid mismatch")

  f3 <- .review_write_nifti(file.path(d, "sub-03_beta.nii.gz"), c(3, 3, 3), spacing = c(3, 3, 3))
  expect_error(gds(c(f1, f3), format = "nifti"), "affine")
})

test_that("beta-only NIfTI sources require unique subject keys", {
  skip_if_not_installed("neuroim2")
  d1 <- tempfile("a_"); d2 <- tempfile("b_")
  dir.create(d1); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  f1 <- .review_write_nifti(file.path(d1, "stat.nii.gz"), c(3, 3, 3))
  f2 <- .review_write_nifti(file.path(d2, "stat.nii.gz"), c(3, 3, 3))
  expect_error(gds(c(f1, f2), format = "nifti"), "duplicate subject")
})

test_that("gds_from_neurovols / nested reject volumes on different grids", {
  skip_if_not_installed("neuroim2")
  v1 <- neuroim2::NeuroVol(array(1, c(3, 3, 3)), neuroim2::NeuroSpace(c(3, 3, 3)))
  v2 <- neuroim2::NeuroVol(array(1, c(4, 3, 3)), neuroim2::NeuroSpace(c(4, 3, 3)))
  v3 <- neuroim2::NeuroVol(array(1, c(3, 3, 3)), neuroim2::NeuroSpace(c(3, 3, 3), spacing = c(2, 2, 2)))
  expect_error(
    suppressWarnings(gds_from_neurovols(list(a = v1, b = v2))),
    "Grid mismatch"
  )
  expect_error(
    suppressWarnings(gds_from_neurovols(list(a = v1, b = v3))),
    "affine"
  )
  expect_error(
    suppressWarnings(gds_from_neurovol_nested(list(
      a = list(c1 = v1, c2 = v1),
      b = list(c1 = v1, c2 = v2)
    ))),
    "Grid mismatch"
  )
})

# 9. Synthetic variance tagging --------------------------------------------------------

test_that("as_gds(NeuroVol/NeuroVec) tags the unit variance as synthetic", {
  skip_if_not_installed("neuroim2")
  sp <- neuroim2::NeuroSpace(c(3, 3, 3, 3))
  vec <- neuroim2::NeuroVec(array(rnorm(81) + 5, c(3, 3, 3, 3)), sp)
  g <- as_gds(vec, along = "subject")
  expect_true(isTRUE(g$metadata$synthetic_var))
  expect_true(isTRUE(attr(assay(g, "var"), "synthetic_unit_variance")))
  expect_error(compute(reduce(g, method = "fixed")), "synthetic")

  vol <- neuroim2::NeuroVol(array(1, c(3, 3, 3)), neuroim2::NeuroSpace(c(3, 3, 3)))
  g1 <- as_gds(vol)
  expect_true(isTRUE(g1$metadata$synthetic_var))
})

# 10. fmristore latent read uses the basis k ------------------------------------------

test_that("fmristore latent read orients embeddings with the basis k", {
  skip_if_not_installed("hdf5r")
  expect_equal(dim(fmrigds:::.fmri_orient_embedding(matrix(1:10, 2, 5), k = 2L)), c(2L, 5L))
  expect_equal(dim(fmrigds:::.fmri_orient_embedding(matrix(1:10, 5, 2), k = 2L)), c(2L, 5L))
  expect_identical(fmrigds:::.fmri_latent_k(c(2L, 10L)), 2L)
  expect_identical(fmrigds:::.fmri_latent_k(c(10L, 2L)), 2L)

  f <- tempfile(fileext = ".lv.h5")
  on.exit(unlink(f), add = TRUE)
  h5 <- hdf5r::H5File$new(f, mode = "w")
  h5$create_group("basis")
  h5$create_group("scans")
  h5[["scans"]]$create_group("s1")
  h5[["basis/basis_matrix"]] <- matrix(rnorm(20), 2, 10) # [k = 2, V = 10]
  h5[["scans/s1/embedding"]] <- matrix(as.numeric(1:10), 2, 5) # [k = 2, C = 5 > k]
  h5$close_all()

  handle <- fmrigds:::.fmri_open(f)
  on.exit(fmrigds:::.fmri_close(handle), add = TRUE, after = FALSE)
  pr <- fmrigds:::.fmri_probe(handle)
  out <- fmrigds:::.fmri_read(handle, assays = "beta")
  expect_equal(dim(out$beta), c(2L, 1L, 5L))
  expect_equal(as.integer(unclass(pr$dims)), dim(out$beta))
  expect_equal(out$beta[, 1, ], matrix(as.numeric(1:10), 2, 5))
})

# 11. NeuroSpace spacing from rotated affines -------------------------------------------

test_that(".gds_space_to_neuroim derives spacing from affine column norms", {
  skip_if_not_installed("neuroim2")
  th <- pi / 6
  R <- matrix(c(cos(th), sin(th), 0, -sin(th), cos(th), 0, 0, 0, 1), 3, 3)
  aff <- diag(4)
  aff[1:3, 1:3] <- R %*% diag(c(2, 2, 3))
  aff[1:3, 4] <- c(-90, -126, -72)
  sp <- space_voxel(dim = c(4, 4, 4), affine = aff)
  ns <- fmrigds:::.gds_space_to_neuroim(sp)
  expect_equal(as.numeric(neuroim2::spacing(ns)), c(2, 2, 3), tolerance = 1e-8)

  # 90-degree rotation: the diagonal is zero, which the diag()-based spacing
  # turned into an invalid zero spacing.
  aff90 <- diag(4)
  aff90[1:3, 1:3] <- matrix(c(0, 2, 0, -2, 0, 0, 0, 0, 2), 3, 3)
  ns90 <- fmrigds:::.gds_space_to_neuroim(space_voxel(dim = c(4, 4, 4), affine = aff90))
  expect_equal(as.numeric(neuroim2::spacing(ns90)), c(2, 2, 2), tolerance = 1e-8)
  expect_equal(unname(neuroim2::trans(ns90)), aff90, tolerance = 1e-6)
})

# 12. CLI parsing ----------------------------------------------------------------------

test_that("CLI subset values stay labels unless marked as indices", {
  p <- fmrigds:::.cli_parse_subset_spec("subject=001,002")
  expect_identical(p$values, c("001", "002"))
  p <- fmrigds:::.cli_parse_subset_spec("contrast=12")
  expect_identical(p$values, "12")
  p <- fmrigds:::.cli_parse_subset_spec("subject=idx:1,3")
  expect_identical(p$values, c(1L, 3L))
})

test_that("CLI treats negative numbers as option values, not flags", {
  expect_false(fmrigds:::.is_flag("-1"))
  expect_false(fmrigds:::.is_flag("-0.5"))
  expect_true(fmrigds:::.is_flag("--json"))
  parsed <- fmrigds:::.cli_parse_args(c("run", "--mask-threshold", "-1", "--json"))
  expect_identical(parsed$opts[["mask-threshold"]], "-1")
})

test_that("CLI metadata ids keep leading zeros and explicit id columns are required", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("id,age", "001,30", "010,40"), f)
  df <- fmrigds:::.cli_read_keyed_data(f, "id", kind = "col-data", required = TRUE)
  expect_identical(rownames(df), c("001", "010"))
  expect_error(
    fmrigds:::.cli_read_keyed_data(f, "subject_id", kind = "col-data", required = TRUE),
    "not found"
  )
})

test_that("tabular sources keep leading-zero subject ids", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(
    "sample,subject,contrast,beta,var",
    "r1,001,c1,0.5,0.1",
    "r1,002,c1,0.7,0.1"
  ), f)
  plan <- gds(f)
  expect_identical(plan$source$probe$subjects, c("001", "002"))
  res <- compute(plan)
  expect_equal(unname(assay(res, "beta")[1, , 1]), c(0.5, 0.7))
})

test_that("fmrigds_cli_exec returns an exit status instead of quitting", {
  out <- utils::capture.output(st <- fmrigds:::fmrigds_cli_exec("--version"))
  expect_identical(st, 0L)
  err <- utils::capture.output(
    st2 <- fmrigds:::fmrigds_cli_exec(c("no-such-command")),
    type = "message"
  )
  expect_identical(st2, 1L)
})
