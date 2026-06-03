skip_if_not_installed("neurotabs")

# -- Helper: build a minimal ROI nftab inline --------------------------------
.make_test_nftab <- function() {
  obs_cols <- list(
    row_id    = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "row_id"),
    subject   = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "subject"),
    condition = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "condition"),
    group     = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "group"),
    roi_1     = neurotabs::nf_col_schema("float32"),
    roi_2     = neurotabs::nf_col_schema("float32"),
    roi_3     = neurotabs::nf_col_schema("float32")
  )

  feat <- neurotabs::nf_feature(
    logical   = neurotabs::nf_logical_schema("vector", "roi", "float32", shape = 3L),
    encodings = list(neurotabs::nf_columns_encoding(c("roi_1", "roi_2", "roi_3")))
  )

  m <- neurotabs::nf_manifest(
    dataset_id          = "test-roi",
    row_id              = "row_id",
    observation_axes    = c("subject", "condition"),
    observation_columns = obs_cols,
    features            = list(roi_beta = feat)
  )

  obs <- data.frame(
    row_id    = c("r1", "r2", "r3", "r4"),
    subject   = c("s01", "s01", "s02", "s02"),
    condition = c("faces", "houses", "faces", "houses"),
    group     = c("ctrl", "ctrl", "pt", "pt"),
    roi_1     = c(0.3, 0.1, 0.4, 0.2),
    roi_2     = c(0.4, 0.2, 0.5, 0.3),
    roi_3     = c(0.3, 0.1, 0.4, 0.1),
    stringsAsFactors = FALSE
  )

  neurotabs::nftab(manifest = m, observations = obs)
}

# -- Helper: single-contrast nftab -------------------------------------------
.make_single_contrast_nftab <- function() {
  obs_cols <- list(
    row_id  = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "row_id"),
    subject = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "subject"),
    roi_1   = neurotabs::nf_col_schema("float32"),
    roi_2   = neurotabs::nf_col_schema("float32")
  )

  feat <- neurotabs::nf_feature(
    logical   = neurotabs::nf_logical_schema("vector", "roi", "float32", shape = 2L),
    encodings = list(neurotabs::nf_columns_encoding(c("roi_1", "roi_2")))
  )

  m <- neurotabs::nf_manifest(
    dataset_id          = "test-single",
    row_id              = "row_id",
    observation_axes    = "subject",
    observation_columns = obs_cols,
    features            = list(roi_beta = feat)
  )

  obs <- data.frame(
    row_id  = c("r1", "r2", "r3"),
    subject = c("s01", "s02", "s03"),
    roi_1   = c(1.0, 2.0, 3.0),
    roi_2   = c(4.0, 5.0, 6.0),
    stringsAsFactors = FALSE
  )

  neurotabs::nftab(manifest = m, observations = obs)
}

.make_primary_feature_nftab <- function(primary_feature = "roi_var") {
  obs_cols <- list(
    row_id    = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "row_id"),
    subject   = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "subject"),
    condition = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "condition"),
    beta_1    = neurotabs::nf_col_schema("float32"),
    beta_2    = neurotabs::nf_col_schema("float32"),
    var_1     = neurotabs::nf_col_schema("float32"),
    var_2     = neurotabs::nf_col_schema("float32")
  )

  beta_feat <- neurotabs::nf_feature(
    logical = neurotabs::nf_logical_schema("vector", "roi", "float32", shape = 2L),
    encodings = list(neurotabs::nf_columns_encoding(c("beta_1", "beta_2")))
  )
  var_feat <- neurotabs::nf_feature(
    logical = neurotabs::nf_logical_schema("vector", "roi", "float32", shape = 2L),
    encodings = list(neurotabs::nf_columns_encoding(c("var_1", "var_2")))
  )

  m <- neurotabs::nf_manifest(
    dataset_id = "test-primary-feature",
    row_id = "row_id",
    observation_axes = c("subject", "condition"),
    observation_columns = obs_cols,
    features = list(roi_beta = beta_feat, roi_var = var_feat),
    primary_feature = primary_feature
  )

  obs <- data.frame(
    row_id = c("r1", "r2"),
    subject = c("s01", "s01"),
    condition = c("faces", "houses"),
    beta_1 = c(1, 2),
    beta_2 = c(10, 20),
    var_1 = c(3, 4),
    var_2 = c(30, 40),
    stringsAsFactors = FALSE
  )

  neurotabs::nftab(manifest = m, observations = obs)
}

.make_supported_parcel_nftab <- function(support_id = "atlas-v1") {
  obs_cols <- list(
    row_id    = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "row_id"),
    subject   = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "subject"),
    condition = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "condition"),
    roi_1     = neurotabs::nf_col_schema("float32"),
    roi_2     = neurotabs::nf_col_schema("float32")
  )

  feat <- neurotabs::nf_feature(
    logical = neurotabs::nf_logical_schema(
      "vector",
      "roi",
      "float32",
      shape = 2L,
      support_ref = "atlas"
    ),
    encodings = list(neurotabs::nf_columns_encoding(c("roi_1", "roi_2")))
  )

  m <- neurotabs::nf_manifest(
    dataset_id = paste0("test-parcel-", support_id),
    row_id = "row_id",
    observation_axes = c("subject", "condition"),
    observation_columns = obs_cols,
    features = list(roi_beta = feat),
    supports = list(
      atlas = neurotabs::nf_support_parcel(
        support_id = support_id,
        space = "MNI152NLin2009cAsym",
        n_parcels = 2L
      )
    )
  )

  obs <- data.frame(
    row_id = c("r1", "r2"),
    subject = c("s01", "s01"),
    condition = c("faces", "houses"),
    roi_1 = c(1, 2),
    roi_2 = c(3, 4),
    stringsAsFactors = FALSE
  )

  neurotabs::nftab(manifest = m, observations = obs)
}

.make_surface_nftab <- function() {
  obs_cols <- list(
    row_id    = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "row_id"),
    subject   = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "subject"),
    condition = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "condition"),
    v1        = neurotabs::nf_col_schema("float32"),
    v2        = neurotabs::nf_col_schema("float32"),
    v3        = neurotabs::nf_col_schema("float32")
  )

  feat <- neurotabs::nf_feature(
    logical = neurotabs::nf_logical_schema(
      "surface",
      "vertex",
      "float32",
      shape = 3L,
      support_ref = "surf"
    ),
    encodings = list(neurotabs::nf_columns_encoding(c("v1", "v2", "v3")))
  )

  m <- neurotabs::nf_manifest(
    dataset_id = "test-surface",
    row_id = "row_id",
    observation_axes = c("subject", "condition"),
    observation_columns = obs_cols,
    features = list(surf_beta = feat),
    supports = list(
      surf = neurotabs::nf_support_surface(
        support_id = "fsaverage-left-demo",
        template = "fsaverage",
        mesh_id = "fsaverage-3k-left",
        topology_id = "fsaverage-3k-topo-left",
        hemisphere = "left"
      )
    )
  )

  obs <- data.frame(
    row_id = c("r1", "r2"),
    subject = c("s01", "s01"),
    condition = c("faces", "houses"),
    v1 = c(1, 2),
    v2 = c(3, 4),
    v3 = c(5, 6),
    stringsAsFactors = FALSE
  )

  neurotabs::nftab(manifest = m, observations = obs)
}

# ===========================================================================
# Tests
# ===========================================================================

test_that("nftab adapter: detect scores nftab objects highly", {
  tab <- .make_test_nftab()
  expect_equal(.nftab_detect(tab), 0.95)
  expect_equal(.nftab_detect(data.frame(x = 1)), 0)
  expect_equal(.nftab_detect("not-an-nftab"), 0)
})

test_that("nftab adapter: probe returns valid contract", {

  tab <- .make_test_nftab()
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle, feature = "roi_beta")

  expect_equal(probe$assays, "beta")
  expect_equal(probe$dims[[1L]], 3L)  # 3 ROIs

  expect_equal(probe$dims[[2L]], 2L)  # 2 subjects
  expect_equal(probe$dims[[3L]], 2L)  # 2 conditions
  expect_equal(probe$subjects, c("s01", "s02"))
  expect_equal(probe$contrasts, c("faces", "houses"))
  expect_s3_class(probe$space, "gds_space")
  expect_s3_class(probe$space, "space_sample_labels")
})

test_that("nftab adapter: probe with named features", {
  tab <- .make_test_nftab()
  handle <- .nftab_open(tab, features = c(beta = "roi_beta"))
  probe <- .nftab_probe(handle)

  expect_equal(probe$assays, "beta")
})

test_that("nftab adapter: probe defaults to first feature", {
  tab <- .make_test_nftab()
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle)

  expect_equal(probe$assays, "beta")
})

test_that("nftab adapter: probe respects primary_feature when no feature is supplied", {
  tab <- .make_primary_feature_nftab()
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle)

  expect_equal(probe$assays, "beta")
  expect_equal(probe$columns$effect_cols$beta, "roi_var")
})

test_that("nftab adapter: read produces correct 3D array", {
  tab <- .make_test_nftab()
  handle <- .nftab_open(tab, feature = "roi_beta")
  probe <- .nftab_probe(handle)

  arrs <- .nftab_read(
    handle,
    assays = "beta",
    effect_cols  = probe$columns$effect_cols,
    subject_col  = probe$columns$subject_col,
    contrast_col = probe$columns$contrast_col
  )

  expect_named(arrs, "beta")
  arr <- arrs$beta
  expect_equal(dim(arr), c(3L, 2L, 2L))  # [sample, subject, contrast]

  # s01, faces → roi values 0.3, 0.4, 0.3
  expect_equal(arr[, 1, 1], c(0.3, 0.4, 0.3))
  # s02, houses → roi values 0.2, 0.3, 0.1
  expect_equal(arr[, 2, 2], c(0.2, 0.3, 0.1))
})

test_that("nftab adapter: read with block subsetting", {
  tab <- .make_test_nftab()
  handle <- .nftab_open(tab, feature = "roi_beta")
  probe <- .nftab_probe(handle)

  arrs <- .nftab_read(
    handle,
    assays = "beta",
    block = list(subject = 1L, sample = c(1L, 3L)),
    effect_cols  = probe$columns$effect_cols,
    subject_col  = probe$columns$subject_col,
    contrast_col = probe$columns$contrast_col
  )

  arr <- arrs$beta
  expect_equal(dim(arr), c(2L, 1L, 2L))  # 2 samples, 1 subject, 2 contrasts
})

test_that("nftab adapter: col_data extracted from design", {
  tab <- .make_test_nftab()
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle, feature = "roi_beta")

  # "group" is constant per subject: s01→ctrl, s02→pt
  expect_false(is.null(probe$col_data))
  expect_true("group" %in% names(probe$col_data))
  expect_equal(rownames(probe$col_data), c("s01", "s02"))
  expect_equal(probe$col_data$group, c("ctrl", "pt"))
})

test_that("nftab adapter: parcel support metadata is preserved in probe space", {
  tab <- .make_supported_parcel_nftab("atlas-v1")
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle, feature = "roi_beta")

  expect_s3_class(probe$space, "space_parcels")
  expect_equal(probe$space$support_id, "atlas-v1")
  expect_equal(probe$space$reference_space, "MNI152NLin2009cAsym")
})

test_that("nftab adapter: support-aware compatibility rejects differing parcel supports", {
  probe1 <- .nftab_probe(.nftab_open(.make_supported_parcel_nftab("atlas-a")), feature = "roi_beta")
  probe2 <- .nftab_probe(.nftab_open(.make_supported_parcel_nftab("atlas-b")), feature = "roi_beta")

  expect_error(assert_compatible_spaces(probe1$space, probe2$space), "support_id")
})

test_that("nftab adapter: surface support probes as a surface space", {
  tab <- .make_surface_nftab()
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle, feature = "surf_beta")

  expect_s3_class(probe$space, "space_surface")
  expect_equal(probe$space$template_id, "fsaverage")
  expect_equal(probe$space$mesh_id, "fsaverage-3k-left")
  expect_equal(probe$space$topology_id, "fsaverage-3k-topo-left")
  expect_equal(probe$space$support_id, "fsaverage-left-demo")
  expect_equal(probe$space$hemi, "L")
})

test_that("nftab adapter: single contrast (no contrast axis)", {
  tab <- .make_single_contrast_nftab()
  handle <- .nftab_open(tab)
  probe <- .nftab_probe(handle, feature = "roi_beta")

  expect_equal(probe$contrasts, "all")
  expect_equal(probe$dims[[3L]], 1L)
  expect_equal(probe$subjects, c("s01", "s02", "s03"))
  expect_equal(probe$dims[[2L]], 3L)
})

test_that("nftab adapter: single contrast read", {
  tab <- .make_single_contrast_nftab()
  handle <- .nftab_open(tab, feature = "roi_beta")
  probe <- .nftab_probe(handle)

  arrs <- .nftab_read(
    handle,
    assays = "beta",
    effect_cols  = probe$columns$effect_cols,
    subject_col  = probe$columns$subject_col,
    contrast_col = probe$columns$contrast_col
  )

  arr <- arrs$beta
  expect_equal(dim(arr), c(2L, 3L, 1L))
  expect_equal(arr[, 1, 1], c(1.0, 4.0))
  expect_equal(arr[, 3, 1], c(3.0, 6.0))
})

test_that("nftab adapter: end-to-end gds() |> compute()", {
  tab <- .make_test_nftab()
  plan <- gds(tab, features = c(signal = "roi_beta"))

  expect_s3_class(plan, "gds_plan")

  result <- plan |> compute()
  expect_s3_class(result, "gds")
  expect_equal(result$subjects, c("s01", "s02"))
  expect_equal(result$contrasts, c("faces", "houses"))
  expect_equal(dim(result$assays[["signal"]]), c(3L, 2L, 2L))
})

test_that("nftab adapter: end-to-end with reduce (beta + var)", {
  # Build nftab with two features so we have beta + var for FE
  obs_cols <- list(
    row_id    = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "row_id"),
    subject   = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "subject"),
    condition = neurotabs::nf_col_schema("string", nullable = FALSE, semantic_role = "condition"),
    b1 = neurotabs::nf_col_schema("float64"),
    b2 = neurotabs::nf_col_schema("float64"),
    v1 = neurotabs::nf_col_schema("float64"),
    v2 = neurotabs::nf_col_schema("float64")
  )
  beta_feat <- neurotabs::nf_feature(
    logical   = neurotabs::nf_logical_schema("vector", "roi", "float64", shape = 2L),
    encodings = list(neurotabs::nf_columns_encoding(c("b1", "b2")))
  )
  var_feat <- neurotabs::nf_feature(
    logical   = neurotabs::nf_logical_schema("vector", "roi", "float64", shape = 2L),
    encodings = list(neurotabs::nf_columns_encoding(c("v1", "v2")))
  )
  m <- neurotabs::nf_manifest(
    dataset_id          = "test-bv",
    row_id              = "row_id",
    observation_axes    = c("subject", "condition"),
    observation_columns = obs_cols,
    features            = list(betas = beta_feat, vars = var_feat)
  )
  obs <- data.frame(
    row_id    = paste0("r", 1:4),
    subject   = c("s01", "s01", "s02", "s02"),
    condition = c("A", "B", "A", "B"),
    b1 = c(1.0, 0.5, 1.2, 0.8),
    b2 = c(0.6, 0.3, 0.9, 0.4),
    v1 = c(0.1, 0.1, 0.2, 0.1),
    v2 = c(0.1, 0.2, 0.1, 0.1),
    stringsAsFactors = FALSE
  )
  tab <- neurotabs::nftab(manifest = m, observations = obs)

  result <- gds(tab, features = c(beta = "betas", var = "vars")) |>
    reduce("fixed") |>
    compute()

  expect_s3_class(result, "gds")
  expect_equal(length(result$subjects), 1L)
  expect_equal(result$contrasts, c("A", "B"))
})

test_that("nftab adapter: as_gds.nftab convenience method", {
  tab <- .make_test_nftab()
  result <- as_gds(tab, features = c(signal = "roi_beta")) |> compute()

  expect_s3_class(result, "gds")
  expect_equal(dim(result$assays[["signal"]]), c(3L, 2L, 2L))
})

test_that("nftab adapter: user-specified axis overrides", {
  tab <- .make_test_nftab()
  handle <- .nftab_open(tab, subject_axis = "subject", contrast_axis = "condition")
  probe <- .nftab_probe(handle, feature = "roi_beta")

  expect_equal(probe$subjects, c("s01", "s02"))
  expect_equal(probe$contrasts, c("faces", "houses"))
})

test_that("register_nftab_adapter registers in the adapter registry", {
  register_nftab_adapter()
  adapter <- get_adapter("nftab")
  expect_true(!is.null(adapter))
  expect_true(is.function(adapter$detect))
})
