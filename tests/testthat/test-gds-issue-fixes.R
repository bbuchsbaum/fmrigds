# Regression tests for issues reported at github.com/bbuchsbaum/fmrigds
# #1 intercept-only OLS without col_data, #6 affine preservation on NIfTI
# export, #7 per-voxel NaN handling in the ols:voxelwise reducer.

test_that("#1 one_sample()/group_ols(~1) run without col_data", {
  set.seed(1)
  n_samples <- 3L; n_subj <- 8L
  beta <- array(rnorm(n_samples * n_subj, mean = 1.5), dim = c(n_samples, n_subj, 1L))
  g <- new_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space_sample_labels(paste0("v", seq_len(n_samples))),
    paste0("s", seq_len(n_subj)),
    "c1"
  )

  fit <- one_sample(g) |> compute()
  expect_true("coef:(Intercept)" %in% names(assays(fit)))
  # Intercept estimate equals the across-subject mean per sample.
  expect_equal(
    as.numeric(assay(fit, "coef:(Intercept)")[, 1, 1]),
    as.numeric(apply(beta[, , 1], 1, mean)),
    tolerance = 1e-8
  )

  # group_ols(~1) should likewise build the trivial design itself.
  fit2 <- group_ols(g, ~ 1) |> compute()
  expect_true("t_coef:(Intercept)" %in% names(assays(fit2)))
})

test_that("#7 ols:voxelwise applies per-voxel listwise deletion for NaN subjects", {
  set.seed(2)
  n_samples <- 4L; n_subj <- 6L
  beta <- array(rnorm(n_samples * n_subj, mean = 2), dim = c(n_samples, n_subj, 1L))
  beta[1L, 1:5, 1L] <- NaN          # sample 1: only 1 finite subject -> unestimable (NA)
  beta[2L, 1L, 1L]  <- NaN          # sample 2: 5 finite subjects -> estimated from subset
  # sample 3 + 4: all finite
  poison <- beta
  poison[, 6L, ] <- NaN             # subject 6 entirely NaN -> must NOT poison whole map

  g <- new_gds(
    list(beta = poison, var = array(1, dim = dim(poison))),
    space_sample_labels(paste0("v", seq_len(n_samples))),
    paste0("s", seq_len(n_subj)),
    "c1"
  )

  expect_warning(fit <- one_sample(g) |> compute(), "listwise deletion")
  co <- as.numeric(assay(fit, "coef:(Intercept)")[, 1, 1])
  n_obs <- as.numeric(assay(fit, "n_obs")[, 1, 1])

  # An all-NaN subject no longer poisons every voxel.
  expect_true(all(is.finite(co[3:4])))
  # Sample 1 has <2 finite obs (subjects 1-5 NaN, subject 6 NaN) -> NA.
  expect_true(is.na(co[1L]))
  # Sample 2: subjects 1 and 6 NaN -> 4 finite -> estimated.
  expect_true(is.finite(co[2L]))
  expect_equal(n_obs, c(0, 4, 5, 5))
  # Estimated intercept equals the mean of the finite subjects.
  expect_equal(co[2L], mean(poison[2L, , 1L], na.rm = TRUE), tolerance = 1e-8)
})

test_that("#11 space() does not mask neuroim2::space() for NeuroVol", {
  skip_if_not_installed("neuroim2")
  v <- neuroim2::read_vol(
    system.file("extdata", "global_mask_v4.nii", package = "neuroim2")
  )
  # fmrigds's space generic is on the search path, yet space() on a NeuroVol
  # must delegate to neuroim2's S4 generic and return its NeuroSpace.
  expect_identical(space(v), neuroim2::space(v))
  expect_s4_class(space(v), "NeuroSpace")
})

test_that("#11 space() still dispatches for gds and gds_plan", {
  g <- structure(list(space = "GDS_SPACE"), class = "gds")
  expect_identical(space(g), "GDS_SPACE")

  p <- structure(list(source = list(probe = list(space = "PLAN_SPACE"))),
                 class = "gds_plan")
  expect_identical(space(p), "PLAN_SPACE")
})

test_that("#6 write_nifti_assays preserves the spatial affine", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("RNifti")

  aff <- matrix(c(-2, 0, 0, 90,
                   0, 2, 0, -126,
                   0, 0, 2, -72,
                   0, 0, 0, 1), nrow = 4L, byrow = TRUE)
  arr <- array(as.numeric(1:8), dim = c(8L, 1L, 1L))
  g <- new_gds(
    list("t_coef:x" = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2L, 2L, 2L), aff),
    "meta",
    "c1"
  )

  td <- tempfile("nifti-affine-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  manifest <- write_nifti_assays(g, out_dir = td, assays = "t_coef:x")
  path <- manifest$path[manifest$written][1L]
  expect_true(file.exists(path))

  v <- neuroim2::read_vol(path)
  sp <- neuroim2::space(v)
  expect_equal(as.numeric(neuroim2::spacing(sp)), c(2, 2, 2), tolerance = 1e-4)
  expect_equal(as.numeric(neuroim2::origin(sp)), c(90, -126, -72), tolerance = 1e-4)
  expect_equal(unname(neuroim2::trans(sp)), unname(aff), tolerance = 1e-4)
  # scl_slope must be a usable scaling (1), not 0.
  expect_equal(RNifti::niftiHeader(path)$scl_slope, 1)
  # Data values survive the round-trip.
  expect_equal(as.numeric(as.array(v)), as.numeric(1:8), tolerance = 1e-4)
})

test_that("#5 variance-weighted reducers refuse a synthetic unit-variance placeholder", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("RNifti")

  td <- tempfile("syn-var-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  files <- file.path(td, paste0("sub-", 1:3, ".nii"))
  for (i in 1:3) RNifti::writeNifti(array(as.numeric(i), c(2, 2, 2)), files[i])

  suppressWarnings(
    g <- gds_from_scalar_maps(files, subject = paste0("s", 1:3), contrast = "metric")
  )
  expect_true(isTRUE(g$metadata$synthetic_var))

  # Variance-weighted reducers are refused with an actionable error.
  expect_error(reduce(g, method = "random"), "synthetic")
  expect_error(reduce(g, method = "fixed"), "synthetic")
  expect_error(reduce(g, method = "meta:fe"), "synthetic")

  # The unweighted OLS path is allowed.
  expect_s3_class(group_ols(g, ~ 1), "gds_plan")
})

test_that("#5 beta-only NIfTI OLS does not materialize synthetic variance", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("RNifti")

  td <- tempfile("beta-only-ols-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  files <- file.path(td, paste0("sub-", 1:3, "_beta.nii"))
  for (i in 1:3) RNifti::writeNifti(array(as.numeric(i), c(2, 2, 2)), files[i])

  plan <- gds(
    nifti_source(beta = files, subject = paste0("s", 1:3), contrast = "metric"),
    format = "nifti"
  )

  fmrigds:::.reset_synthetic_variance_warning()
  expect_warning(compute(plan), "No variance or SE provided; using unit variance", fixed = TRUE)
  expect_silent(fit <- group_ols(plan, ~ 1) |> compute())
  expect_true("coef:(Intercept)" %in% names(assays(fit)))
})

test_that("#5 synthetic unit-variance warning fires at most once per session (lazy)", {
  skip_if_not_installed("RNifti")

  td <- tempfile("once-warn-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  files <- file.path(td, paste0("sub-", 1:3, "_beta.nii"))
  for (i in 1:3) RNifti::writeNifti(array(as.numeric(i), c(2, 2, 2)), files[i])

  plan <- gds(
    nifti_source(beta = files, subject = paste0("s", 1:3), contrast = "metric"),
    format = "nifti"
  )

  fmrigds:::.reset_synthetic_variance_warning()
  expect_warning(compute(plan), "No variance or SE provided; using unit variance", fixed = TRUE)
  # Second materialisation in the same session is de-noised.
  expect_silent(compute(plan))
})

test_that("#5 eager beta-only constructors tag synthetic variance and reducers refuse it", {
  skip_if_not_installed("neuroim2")
  mk <- function(v) neuroim2::NeuroVol(array(v, c(3, 3, 3)), neuroim2::NeuroSpace(c(3, 3, 3)))
  beta <- list(`sub-01` = mk(1), `sub-02` = mk(2), `sub-03` = mk(3))

  g <- suppressWarnings(gds_from_neurovols(beta))
  expect_true(isTRUE(g$metadata$synthetic_var))

  expect_error(reduce(g, method = "fixed"),  "synthetic")
  expect_error(reduce(g, method = "random"), "synthetic")
  expect_error(reduce(g, method = "meta:fe"), "synthetic")
  # Unweighted path stays open, including a direct reduce() with default weights
  # (ols:voxelwise ignores variance and must not be blocked).
  expect_s3_class(group_ols(g, ~ 1), "gds_plan")
  expect_s3_class(reduce(g, method = "ols:voxelwise", formula = ~ 1), "gds_plan")
})

test_that("#5 reduce() guard catches a lazy beta-only NIfTI plan", {
  skip_if_not_installed("RNifti")

  td <- tempfile("lazy-guard-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  files <- file.path(td, paste0("sub-", 1:3, "_beta.nii"))
  for (i in 1:3) RNifti::writeNifti(array(as.numeric(i), c(2, 2, 2)), files[i])

  plan <- gds(
    nifti_source(beta = files, subject = paste0("s", 1:3), contrast = "metric"),
    format = "nifti"
  )

  expect_error(reduce(plan, method = "fixed"),  "synthetic")
  expect_error(reduce(plan, method = "random"), "synthetic")
})

test_that("#5 apply_reduce backstop refuses attr-tagged synthetic var, not a plain array(1)", {
  beta <- array(c(1, 2, 3, 4, 5, 6), dim = c(2L, 3L, 1L))
  node <- op_reduce(method = "fixed", weights = "1/var", by = "contrast", options = list())

  syn <- array(1, dim = dim(beta))
  attr(syn, "synthetic_unit_variance") <- TRUE
  expect_error(
    apply_reduce(node, list(beta = beta, var = syn), weights = "1/var",
                 subjects = paste0("s", 1:3)),
    "synthetic"
  )

  # Untagged unit variance (a legitimate real var that happens to be 1) is allowed.
  plain <- array(1, dim = dim(beta))
  expect_silent(suppressWarnings(
    apply_reduce(node, list(beta = beta, var = plain), weights = "1/var",
                 subjects = paste0("s", 1:3))
  ))
})

test_that("#4 list_assays() reports computed assays and previews reducer outputs", {
  set.seed(3)
  beta <- array(rnorm(4 * 6), dim = c(4L, 6L, 1L))
  g <- new_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space_sample_labels(paste0("v", 1:4)),
    paste0("s", 1:6),
    "c1"
  )
  fit <- one_sample(g) |> compute()

  tbl <- list_assays(fit)
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("assay", "role", "units") %in% names(tbl)))
  expect_true("coef:(Intercept)" %in% tbl$assay)

  # info = FALSE returns a plain character vector.
  nm <- list_assays(fit, info = FALSE)
  expect_type(nm, "character")
  expect_true("t_coef:(Intercept)" %in% nm)

  # reducer preview works without computing; alias resolves.
  prov <- list_assays(reducer = "random", info = FALSE)
  expect_true(all(c("beta_g", "se_g", "z_g", "p_g", "tau2") %in% prov))

  # input assays on an unrealised plan.
  plan_assays <- list_assays(group_ols(g, ~ 1), info = FALSE)
  expect_true("beta" %in% plan_assays)

  expect_error(list_assays(reducer = "not_a_reducer"), "Unknown reducer")
})

test_that("#14 multi-term FDR produces one q_coef:<term> per term", {
  set.seed(7)
  n_s <- 5L; n_subj <- 12L
  beta <- array(rnorm(n_s * n_subj, 1), dim = c(n_s, n_subj, 1L))
  g <- new_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space_sample_labels(paste0("v", seq_len(n_s))),
    paste0("s", seq_len(n_subj)),
    "c1"
  )
  cd <- data.frame(age = rnorm(n_subj), row.names = paste0("s", seq_len(n_subj)))

  fit <- group_ols(g, ~ age, col_data = cd) |> posthoc("fdr:bh") |> compute()
  nm <- names(assays(fit))
  expect_true(all(c("q_coef:(Intercept)", "q_coef:age") %in% nm))
  expect_false("q" %in% nm) # no ambiguous bare q for a regression model
  expect_equal(
    as.numeric(assay(fit, "q_coef:age")[, 1, 1]),
    stats::p.adjust(as.numeric(assay(fit, "p_coef:age")[, 1, 1]), "BH"),
    tolerance = 1e-8
  )

  # A genuine bare p (meta path) still yields a single q (back-compat).
  gv <- new_gds(
    list(beta = beta, var = array(0.5, dim = dim(beta))),
    space_sample_labels(paste0("v", seq_len(n_s))),
    paste0("s", seq_len(n_subj)),
    "c1"
  )
  fit2 <- reduce(gv, method = "fixed") |> posthoc("fdr:bh") |> compute()
  expect_true("q" %in% names(assays(fit2)))
})

test_that("#12 fdr:spatial gives an actionable error on a bare voxelwise GDS", {
  n_s <- 12L; n_subj <- 8L
  beta <- array(rnorm(n_s * n_subj), dim = c(n_s, n_subj, 1L))
  g <- new_gds(
    list(beta = beta, var = array(1, dim = dim(beta))),
    space_voxel(c(2L, 2L, 3L), diag(4)),
    paste0("s", seq_len(n_subj)),
    "shape_old"
  )
  g2 <- reduce(g, method = "ols:voxelwise", formula = ~ 1)
  err <- tryCatch(compute(posthoc(g2, "fdr:spatial")), error = function(e) conditionMessage(e))
  expect_match(err, "parcel/cluster id per sample")
  expect_match(err, "fdr:bh", fixed = TRUE)

  # With an explicit grouping it runs and produces a q_coef assay (single term).
  grp <- rep(1:3, length.out = n_s)
  fit <- compute(posthoc(g2, "fdr:spatial", options = list(group = grp)))
  expect_true(any(grepl("^q", names(assays(fit)))))
})

test_that("#7 meta reducers expose effective subject counts for finite inputs", {
  beta <- array(c(1, 2, 3, 4, NaN, 6), dim = c(2L, 3L, 1L))
  var <- array(1, dim = dim(beta))

  fixed_node <- op_reduce(method = "fixed", weights = "1/var", by = "contrast", options = list())
  expect_warning(
    fixed <- apply_reduce(fixed_node, list(beta = beta, var = var), weights = "1/var", subjects = paste0("s", 1:3)),
    "n_eff"
  )
  expect_true("n_eff" %in% names(fixed$arrays))
  expect_equal(as.numeric(fixed$arrays$n_eff[, 1, 1]), c(2, 3))
  expect_equal(as.numeric(fixed$arrays$beta_g[, 1, 1]), c(2, 4), tolerance = 1e-8)

  random_node <- op_reduce(method = "random", weights = "1/var", by = "contrast", options = list())
  expect_warning(
    random <- apply_reduce(random_node, list(beta = beta, var = var), weights = "1/var", subjects = paste0("s", 1:3)),
    "n_eff"
  )
  expect_true("n_eff" %in% names(random$arrays))
  expect_equal(as.numeric(random$arrays$n_eff[, 1, 1]), c(2, 3))
})

test_that("#10 synthetic positional sample labels are not exposed as real labels", {
  beta <- array(seq_len(15), dim = c(5L, 3L, 1L))
  se <- array(1, dim = dim(beta))

  g <- as_gds(
    list(beta = beta, se = se),
    space = NULL,
    subjects = paste0("sub-", 1:3),
    contrasts = "c1"
  )

  expect_null(sample_labels(g))
  expect_false("sample" %in% names(row_data(g)))
  expect_equal(nrow(row_data(g)), 5L)
})

test_that("#10 explicit sample labels are preserved", {
  beta <- array(seq_len(15), dim = c(5L, 3L, 1L))
  se <- array(1, dim = dim(beta))
  dimnames(beta) <- list(paste0("roi-", 1:5), paste0("sub-", 1:3), "c1")
  dimnames(se) <- dimnames(beta)

  g <- as_gds(list(beta = beta, se = se))
  expect_equal(sample_labels(g), paste0("roi-", 1:5))

  numeric_labels <- new_gds(
    list(beta = beta, se = se),
    space_sample_labels(as.character(seq_len(5L))),
    paste0("sub-", 1:3),
    "c1"
  )
  expect_equal(sample_labels(numeric_labels), as.character(seq_len(5L)))
})

test_that("#6 write_out(format = 'nifti') preserves the spatial affine", {
  skip_if_not_installed("neuroim2")

  aff <- matrix(c(3, 0, 0, 10,
                  0, 3, 0, 20,
                  0, 0, 3, 30,
                  0, 0, 0, 1), nrow = 4L, byrow = TRUE)
  arr <- array(as.numeric(1:8), dim = c(8L, 1L, 1L))
  g <- new_gds(
    list(beta = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(2L, 2L, 2L), aff),
    "meta",
    "c1"
  )

  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path), add = TRUE)
  .write_gds_nifti(g, path, options = list(stat = "beta"))

  v <- neuroim2::read_vol(path)
  expect_equal(as.numeric(neuroim2::spacing(neuroim2::space(v))), c(3, 3, 3), tolerance = 1e-4)
  expect_equal(as.numeric(neuroim2::origin(neuroim2::space(v))), c(10, 20, 30), tolerance = 1e-4)
})

test_that("#6 write_out(format = 'nifti') preserves the affine for >3D output (RNifti fallback)", {
  skip_if_not_installed("RNifti")

  # Multi-subject and/or multi-contrast output stays 4D/5D and is written via
  # the RNifti fallback (not the 3D neuroim2 path). The naive fallback wrote a
  # correct sform but a qform that lost voxel scale, so qform-preferring readers
  # (including RNifti's default xform()) saw unit spacing. Exercise that branch.
  aff <- matrix(c(-2, 0, 0, 90,
                   0, 2, 0, -126,
                   0, 0, 2, -72,
                   0, 0, 0, 1), nrow = 4L, byrow = TRUE)
  vdim <- c(4L, 5L, 6L)
  nvox <- prod(vdim)

  cases <- list(
    c(nsub = 2L, ncon = 1L), # 4D
    c(nsub = 1L, ncon = 2L), # 4D
    c(nsub = 2L, ncon = 2L)  # 5D
  )
  for (cs in cases) {
    arr <- array(seq_len(nvox * cs[["nsub"]] * cs[["ncon"]]) * 1.0,
                 dim = c(nvox, cs[["nsub"]], cs[["ncon"]]))
    g <- new_gds(
      list(beta = arr, var = array(1, dim = dim(arr))),
      space_voxel(vdim, aff, storage = "dense"),
      paste0("s", seq_len(cs[["nsub"]])),
      paste0("c", seq_len(cs[["ncon"]]))
    )
    path <- tempfile(fileext = ".nii.gz")
    .write_gds_nifti(g, path, options = list(stat = "beta"))

    img <- RNifti::readNifti(path)
    # The affine must read back correctly regardless of qform/sform preference.
    xq <- RNifti::xform(img, useQuaternionFirst = TRUE)
    xs <- RNifti::xform(img, useQuaternionFirst = FALSE)
    expect_equal(unname(xq)[1:3, ], unname(aff)[1:3, ], tolerance = 1e-3)
    expect_equal(unname(xs)[1:3, ], unname(aff)[1:3, ], tolerance = 1e-3)
    hdr <- RNifti::niftiHeader(path)
    expect_equal(as.numeric(hdr$pixdim[2:4]), c(2, 2, 2), tolerance = 1e-3)

    # Data volumes must land in the right place (first subject x first contrast).
    got <- array(as.array(img), dim = c(vdim, cs[["nsub"]], cs[["ncon"]]))
    expect_equal(as.numeric(got[, , , 1, 1]),
                 as.numeric(array(arr[, 1, 1], dim = vdim)), tolerance = 1e-4)
    unlink(path)
  }
})

test_that("#6 RNifti fallback leaves qform unset for a non-orthogonal (sheared) affine", {
  skip_if_not_installed("RNifti")

  # A sheared affine cannot be represented by the quaternion qform; stamping it
  # anyway would silently approximate the geometry. The writer must keep the
  # exact sform and disable the qform (code 0) so readers use the sform.
  shear <- matrix(c(2, 0.7, 0, 10,
                    0, 2,   0, 20,
                    0, 0,   2, 30,
                    0, 0,   0, 1), nrow = 4L, byrow = TRUE)
  arr <- array(seq_len(prod(c(60, 2, 1))) * 1.0, dim = c(60L, 2L, 1L)) # 4D fallback
  g <- new_gds(
    list(beta = arr, var = array(1, dim = dim(arr))),
    space_voxel(c(3L, 4L, 5L), shear, storage = "dense"),
    c("s1", "s2"),
    "c1"
  )
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path), add = TRUE)
  .write_gds_nifti(g, path, options = list(stat = "beta"))

  hdr <- RNifti::niftiHeader(path)
  expect_identical(as.integer(hdr$qform_code), 0L)
  expect_gt(hdr$sform_code, 0)
  # The sform (read with sform preference) keeps the exact sheared affine.
  xs <- RNifti::xform(RNifti::readNifti(path), useQuaternionFirst = FALSE)
  expect_equal(unname(xs)[1:3, ], unname(shear)[1:3, ], tolerance = 1e-3)
})
