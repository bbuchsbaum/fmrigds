test_that("nt:tfce_fwer registration follows neurothresh availability", {
  has_nt <- requireNamespace("neurothresh", quietly = TRUE)
  expect_equal("nt:tfce_fwer" %in% list_posthoc(), has_nt)
  expect_equal("nt:cluster_fdr_perm" %in% list_posthoc(), has_nt)
})

.make_nt_tfce_fixture <- function(seed = 1L, dims3 = c(5L, 5L, 5L), n_subj = 8L) {
  set.seed(seed)
  n_vox <- prod(dims3)
  subj <- replicate(n_subj, array(stats::rnorm(n_vox), dim = dims3), simplify = FALSE)
  mask <- array(stats::runif(n_vox) > 0.25, dim = dims3)
  if (sum(mask) < 10L) mask[] <- TRUE

  sf <- neurothresh::make_null_fun_subject_signflip(subj, mask = mask, seed = seed + 1L)
  mask_idx <- which(as.vector(mask))
  z_vec <- as.numeric(sf$z_vol)[mask_idx]
  z <- array(z_vec, dim = c(length(mask_idx), 1L, 1L))

  g <- as_gds(
    list(z = z),
    space = space_voxel(dim = dims3, affine = diag(4), mask_idx = mask_idx, storage = "packed"),
    subjects = "meta",
    contrasts = "c1"
  )

  list(g = g, null_fun = sf$null_fun)
}

test_that("nt:tfce_fwer computes q, sig_mask, and tfce on packed voxel data", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:tfce_fwer" %in% list_posthoc()))

  fx <- .make_nt_tfce_fixture(seed = 10L)
  out <- posthoc(
    fx$g,
    method = "nt:tfce_fwer",
    options = list(
      n_perm = 32L,
      alpha = 0.10,
      tail = "pos",
      seed = 123L,
      null_fun = fx$null_fun
    )
  ) |> compute()

  expect_true(all(c("q", "sig_mask", "tfce") %in% names(assays(out))))
  q <- assay(out, "q")[, 1, 1]
  expect_true(all(is.finite(q)))
  expect_true(all(q >= 0 & q <= 1))
  sig <- assay(out, "sig_mask")[, 1, 1]
  expect_true(all(sig %in% c(0, 1)))

  att <- metadata(out)$attachments
  key <- "posthoc/nt:tfce_fwer/subject=1/contrast=1"
  expect_true(key %in% names(att))
  expect_equal(att[[key]]$null_model, "explicit_null_fun")
})

test_that("nt:tfce_fwer requires explicit null_fun unless fallback is enabled", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:tfce_fwer" %in% list_posthoc()))

  fx <- .make_nt_tfce_fixture(seed = 20L)
  expect_error(
    posthoc(
      fx$g,
      method = "nt:tfce_fwer",
      options = list(n_perm = 8L, alpha = 0.10, tail = "pos")
    ) |> compute(),
    "requires options\\$null_fun unless allow_voxel_signflip = TRUE"
  )

  out <- posthoc(
    fx$g,
    method = "nt:tfce_fwer",
    options = list(
      n_perm = 12L,
      alpha = 0.10,
      tail = "pos",
      seed = 99L,
      allow_voxel_signflip = TRUE
    )
  ) |> compute()

  key <- "posthoc/nt:tfce_fwer/subject=1/contrast=1"
  expect_equal(metadata(out)$attachments[[key]]$null_model, "voxel_signflip_fallback")
})

test_that("nt:tfce_fwer rejects non-voxel spaces", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:tfce_fwer" %in% list_posthoc()))

  z <- array(stats::rnorm(6), dim = c(6, 1, 1))
  g <- as_gds(
    list(z = z),
    space = space_sample_labels(paste0("R", seq_len(6))),
    subjects = "meta",
    contrasts = "c1"
  )

  expect_error(
    posthoc(g, method = "nt:tfce_fwer", options = list(allow_voxel_signflip = TRUE, n_perm = 6L)) |> compute(),
    "requires voxel space"
  )
})

test_that("nt:tfce_fwer null rejection rate stays near target under global null", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:tfce_fwer" %in% list_posthoc()))
  skip_on_cran()

  alpha <- 0.10
  n_sims <- 8L
  reject <- logical(n_sims)

  for (s in seq_len(n_sims)) {
    fx <- .make_nt_tfce_fixture(seed = 1000L + s, dims3 = c(4L, 4L, 4L), n_subj = 8L)
    out <- posthoc(
      fx$g,
      method = "nt:tfce_fwer",
      options = list(
        n_perm = 24L,
        alpha = alpha,
        tail = "pos",
        seed = 2000L + s,
        null_fun = fx$null_fun
      )
    ) |> compute()
    reject[s] <- any(assay(out, "sig_mask")[, 1, 1] > 0.5, na.rm = TRUE)
  }

  reject_rate <- mean(reject)
  expect_lte(reject_rate, alpha + 0.30)
})

test_that("nt:cluster_fdr_perm computes q and sig_mask on packed voxel data", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:cluster_fdr_perm" %in% list_posthoc()))

  fx <- .make_nt_tfce_fixture(seed = 300L)
  out <- posthoc(
    fx$g,
    method = "nt:cluster_fdr_perm",
    options = list(
      n_perm = 28L,
      q = 0.10,
      cluster_thresh = 2.5,
      tail = "pos",
      seed = 33L,
      null_fun = fx$null_fun
    )
  ) |> compute()

  expect_true(all(c("q", "sig_mask") %in% names(assays(out))))
  q <- assay(out, "q")[, 1, 1]
  expect_true(all(is.finite(q)))
  expect_true(all(q >= 0 & q <= 1))
  sig <- assay(out, "sig_mask")[, 1, 1]
  expect_true(all(sig %in% c(0, 1)))

  att <- metadata(out)$attachments
  key <- "posthoc/nt:cluster_fdr_perm/subject=1/contrast=1"
  expect_true(key %in% names(att))
  expect_equal(att[[key]]$null_model, "explicit_null_fun")
})

test_that("nt:cluster_fdr_perm requires explicit null_fun unless fallback is enabled", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:cluster_fdr_perm" %in% list_posthoc()))

  fx <- .make_nt_tfce_fixture(seed = 350L)
  expect_error(
    posthoc(
      fx$g,
      method = "nt:cluster_fdr_perm",
      options = list(n_perm = 10L, q = 0.10, cluster_thresh = 2.5, tail = "pos")
    ) |> compute(),
    "requires options\\$null_fun unless allow_voxel_signflip = TRUE"
  )

  out <- posthoc(
    fx$g,
    method = "nt:cluster_fdr_perm",
    options = list(
      n_perm = 12L,
      q = 0.10,
      cluster_thresh = 2.5,
      tail = "pos",
      seed = 77L,
      allow_voxel_signflip = TRUE
    )
  ) |> compute()

  key <- "posthoc/nt:cluster_fdr_perm/subject=1/contrast=1"
  expect_equal(metadata(out)$attachments[[key]]$null_model, "voxel_signflip_fallback")
})

test_that("nt:cluster_fdr_perm supports two-sided mode with bounded q", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:cluster_fdr_perm" %in% list_posthoc()))

  fx <- .make_nt_tfce_fixture(seed = 375L)
  out <- posthoc(
    fx$g,
    method = "nt:cluster_fdr_perm",
    options = list(
      n_perm = 24L,
      q = 0.10,
      cluster_thresh = 2.5,
      tail = "two",
      two_sided_policy = "BH_all",
      seed = 88L,
      null_fun = fx$null_fun
    )
  ) |> compute()

  q <- assay(out, "q")[, 1, 1]
  expect_true(all(is.finite(q)))
  expect_true(all(q >= 0 & q <= 1))
})

test_that("nt:cluster_fdr_perm rejects non-voxel spaces", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:cluster_fdr_perm" %in% list_posthoc()))

  z <- array(stats::rnorm(6), dim = c(6, 1, 1))
  g <- as_gds(
    list(z = z),
    space = space_sample_labels(paste0("R", seq_len(6))),
    subjects = "meta",
    contrasts = "c1"
  )

  expect_error(
    posthoc(
      g,
      method = "nt:cluster_fdr_perm",
      options = list(allow_voxel_signflip = TRUE, n_perm = 6L, cluster_thresh = 2.5)
    ) |> compute(),
    "requires voxel space"
  )
})

test_that("nt:cluster_fdr_perm null rejection rate stays near target under global null", {
  skip_if_not_installed("neurothresh")
  skip_if(!("nt:cluster_fdr_perm" %in% list_posthoc()))
  skip_on_cran()

  q_level <- 0.10
  n_sims <- 8L
  reject <- logical(n_sims)

  for (s in seq_len(n_sims)) {
    fx <- .make_nt_tfce_fixture(seed = 5000L + s, dims3 = c(4L, 4L, 4L), n_subj = 8L)
    out <- posthoc(
      fx$g,
      method = "nt:cluster_fdr_perm",
      options = list(
        n_perm = 24L,
        q = q_level,
        cluster_thresh = 2.5,
        tail = "pos",
        seed = 6000L + s,
        null_fun = fx$null_fun
      )
    ) |> compute()
    reject[s] <- any(assay(out, "sig_mask")[, 1, 1] > 0.5, na.rm = TRUE)
  }

  reject_rate <- mean(reject)
  expect_lte(reject_rate, q_level + 0.35)
})
