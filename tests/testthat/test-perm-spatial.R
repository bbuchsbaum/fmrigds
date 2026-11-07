# Subject-level sign-flip TFCE / cluster inference via perm:onesample.

.sp_smooth <- function(a, iter = 2L) {
  # Repeated 6-neighbour averaging (zero padding): a cheap smooth field.
  d <- dim(a)
  for (it in seq_len(iter)) {
    b <- a
    s <- array(1, d)
    sh <- function(x, dim_i, k) {
      out <- array(0, d)
      idx <- lapply(d, seq_len)
      src <- idx
      if (k > 0) {
        idx[[dim_i]] <- (1 + k):d[dim_i]
        src[[dim_i]] <- 1:(d[dim_i] - k)
      } else {
        idx[[dim_i]] <- 1:(d[dim_i] + k)
        src[[dim_i]] <- (1 - k):d[dim_i]
      }
      out <- do.call(`[<-`, c(list(out), idx, list(value = do.call(`[`, c(list(x), src)))))
      out
    }
    for (dd in 1:3) for (k in c(-1L, 1L)) {
      b <- b + sh(a, dd, k)
      s <- s + sh(array(1, d), dd, k)
    }
    a <- b / s
  }
  a
}

.sp_fixture <- function(seed = 1L, dims = c(9L, 8L, 6L), n_subj = 10L,
                        effect = 0, blob = NULL, mask_frac = 0.85, smooth = 2L) {
  set.seed(seed)
  mask <- array(stats::runif(prod(dims)) < mask_frac, dims)
  mask_idx <- which(mask)
  maps <- lapply(seq_len(n_subj), function(s) {
    a <- .sp_smooth(array(stats::rnorm(prod(dims)), dims), smooth)
    a <- a / stats::sd(a)
    if (!is.null(blob)) a[blob] <- a[blob] + effect
    a
  })
  beta <- vapply(maps, function(a) a[mask_idx], numeric(length(mask_idx)))
  g <- new_gds(
    list(beta = array(beta, c(length(mask_idx), n_subj, 1L)),
         var = array(1, c(length(mask_idx), n_subj, 1L))),
    space_voxel(dims, diag(c(2, 2, 2, 1)), mask_idx = mask_idx, storage = "packed"),
    paste0("s", seq_len(n_subj)), "c1"
  )
  list(g = g, maps = maps, mask = mask, mask_idx = mask_idx, dims = dims)
}

.sp_run <- function(g, ...) {
  compute(reduce(g, method = "perm:onesample", options = list(...)))
}

.sp_att <- function(r) metadata(r)$attachments[["reduce/perm:onesample/contrast=c1"]]

test_that("spatial run keeps t, p_perm and p_fwer of the plain kernel", {
  fx <- .sp_fixture(seed = 2L)
  for (alt in c("two.sided", "greater", "less")) {
    plain <- .sp_run(fx$g, n_perm = 99L, seed = 3L, alternative = alt)
    sp <- .sp_run(fx$g, n_perm = 99L, seed = 3L, alternative = alt,
                  spatial = c("tfce", "cluster"))
    for (nm in c("beta_g", "se_g", "t_g", "df", "p_g", "p_perm", "p_fwer")) {
      expect_equal(assay(sp, nm), assay(plain, nm), info = paste(alt, nm))
    }
    # FWER p can never be below the uncorrected p (#22).
    expect_true(all(assay(sp, "p_fwer") >= assay(sp, "p_perm") - 1e-12), info = alt)
    expect_true(all(assay(plain, "p_fwer") >= assay(plain, "p_perm") - 1e-12), info = alt)
    att <- .sp_att(sp)
    expect_equal(att$n_null, 99L)
    expect_length(att$null$max_t, 99L)
    expect_length(att$null$tfce, 99L)
    expect_length(att$null$cluster_mass, 99L)
  }
  expect_false(any(c("tfce", "p_tfce_fwer", "cluster_id") %in% names(assays(plain))))
})

test_that("directional FWER follows alternative on a negative effect (#22)", {
  set.seed(42)
  nv <- 400; ns <- 12
  beta <- array(stats::rnorm(nv * ns, mean = -0.5, sd = 0.3), dim = c(nv, ns, 1))
  g <- new_gds(list(beta = beta, var = array(1, dim(beta))),
               space_sample_labels(paste0("v", 1:nv)), paste0("s", 1:ns), "c1")
  r <- .sp_run(g, n_perm = 199L, seed = 1L, alternative = "greater")
  expect_true(all(assay(r, "p_fwer") >= assay(r, "p_perm") - 1e-12))
  expect_true(min(assay(r, "p_fwer")) > 0.5)
})

test_that("one sign matrix drives every null; observed row excluded", {
  fx <- .sp_fixture(seed = 4L)
  set.seed(9)
  signs <- matrix(sample(c(-1L, 1L), 40 * 10, TRUE), 40, 10)
  signs <- rbind(1L, signs)
  r1 <- .sp_run(fx$g, signs = signs, spatial = c("tfce", "cluster"))
  att <- .sp_att(r1)
  expect_equal(att$n_null, 40L)
  # Running the null rows alone (no identity row) gives the same nulls.
  r2 <- .sp_run(fx$g, signs = signs[-1, ], spatial = c("tfce", "cluster"))
  expect_equal(.sp_att(r2)$null, att$null)
  expect_equal(assay(r2, "p_tfce_fwer"), assay(r1, "p_tfce_fwer"))
  # p = (1 + #{null >= obs}) / (1 + n_null)
  score <- abs(as.vector(assay(r1, "tfce")))
  p_ref <- vapply(score, function(s) (1 + sum(att$null$tfce >= s)) / 41, numeric(1))
  expect_equal(as.vector(assay(r1, "p_tfce_fwer")), p_ref)
})

test_that("results do not depend on the thread count", {
  fx <- .sp_fixture(seed = 5L)
  a <- .sp_run(fx$g, n_perm = 60L, seed = 2L, spatial = c("tfce", "cluster"), threads = 1L)
  b <- .sp_run(fx$g, n_perm = 60L, seed = 2L, spatial = c("tfce", "cluster"), threads = 3L)
  for (nm in names(assays(a))) expect_equal(assay(a, nm), assay(b, nm), info = nm)
  expect_equal(.sp_att(a)$null, .sp_att(b)$null)
})

test_that("TFCE and clusters equal neurothresh given identical sign flips", {
  skip_if_not_installed("neurothresh")
  nt <- asNamespace("neurothresh")
  skip_if(!exists("label_clusters", nt), "neurothresh without compiled TFCE")
  fx <- .sp_fixture(seed = 6L, dims = c(10L, 9L, 7L), n_subj = 9L)
  B <- 40L
  spec <- neurothresh::null_one_sample(fx$maps, mask = fx$mask, n_perm = B, seed = 7L)
  signs <- rbind(1L, matrix(as.integer(spec$meta$signs), B))
  for (alt in c("two.sided", "greater", "less")) {
    tail <- c(two.sided = "two", greater = "pos", less = "neg")[[alt]]
    ref <- neurothresh::tfce_fwer(spec$stat_vol, mask = fx$mask, null_spec = spec, tail = tail)
    r <- .sp_run(fx$g, signs = signs, alternative = alt, spatial = c("tfce", "cluster"),
                 cluster_threshold = 1.5)
    att <- .sp_att(r)
    expect_equal(att$null$tfce, ref$null_max, tolerance = 1e-10, info = alt)
    tf_ref <- as.numeric(if (tail == "two") ref$extras$tfce_two else ref$extras$tfce)
    expect_equal(abs(as.vector(assay(r, "tfce"))), tf_ref[fx$mask_idx],
                 tolerance = 1e-10, info = alt)
    expect_equal(as.vector(assay(r, "p_tfce_fwer")),
                 as.numeric(ref$pvals)[fx$mask_idx], info = alt)
    expect_equal(as.vector(assay(r, "z_g")), as.numeric(spec$stat_vol)[fx$mask_idx],
                 tolerance = 1e-10, info = alt)
    # Cluster null: largest excess mass over the tested tail(s) per flip.
    tails <- switch(tail, two = c("pos", "neg"), tail)
    mass_ref <- vapply(seq_len(B), function(b) {
      d <- array(NA_real_, fx$dims)
      d[fx$mask_idx] <- spec$draw(b)
      max(0, unlist(lapply(tails, function(tl) {
        neurothresh::label_clusters(d, 1.5, mask = fx$mask, tail = tl)$table$mass
      })))
    }, numeric(1))
    expect_equal(att$null$cluster_mass, mass_ref, tolerance = 1e-10, info = alt)
  }
})

test_that("vendored spatial core matches the neurothresh header", {
  hdr <- system.file("include", "neurothresh", "spatial_core.h", package = "neurothresh")
  skip_if(!nzchar(hdr), "neurothresh header not installed")
  ours <- test_path("..", "..", "src", "nt_spatial_core.h")
  skip_if(!file.exists(ours), "package sources not available")
  body <- function(f) {
    x <- readLines(f)
    x[seq.int(grep("^#ifndef NEUROTHRESH_SPATIAL_CORE_H", x), length(x))]
  }
  expect_identical(body(ours), body(hdr))
})

test_that("cluster table, cluster_id and p_cluster_fwer agree", {
  dims <- c(12L, 10L, 8L)
  blob <- array(FALSE, dims)
  blob[3:6, 3:6, 2:4] <- TRUE
  fx <- .sp_fixture(seed = 8L, dims = dims, n_subj = 12L, effect = 1.5, blob = blob,
                    mask_frac = 1)
  r <- .sp_run(fx$g, n_perm = 199L, seed = 1L, alternative = "greater",
               spatial = c("tfce", "cluster"), cluster_p = 0.01)
  att <- .sp_att(r)
  tab <- att$clusters
  cid <- as.vector(assay(r, "cluster_id"))
  pc <- as.vector(assay(r, "p_cluster_fwer"))
  expect_gt(nrow(tab), 0L)
  expect_true(all(tab$sign == 1L))
  expect_equal(tab$cluster_id, seq_len(nrow(tab)))
  expect_true(all(diff(tab$mass) <= 0))
  expect_equal(as.integer(table(factor(cid[cid > 0], levels = tab$cluster_id))), tab$size)
  for (k in tab$cluster_id) expect_true(all(pc[cid == k] == tab$p_fwer_mass[k]))
  expect_true(all(pc[cid == 0] == 1))
  # The planted blob is found and survives FWER by TFCE and cluster mass.
  inb <- blob[fx$mask_idx]
  expect_lt(min(as.vector(assay(r, "p_tfce_fwer"))[inb]), 0.05)
  expect_lt(tab$p_fwer_mass[1], 0.05)
  expect_gt(mean(inb[cid == 1]), 0.5)
  # Peak coordinates: grid indices and world coordinates via the affine.
  pk <- tab[1, ]
  expect_equal(fx$mask_idx[pk$peak_sample],
               pk$peak_i + dims[1] * ((pk$peak_j - 1) + dims[2] * (pk$peak_k - 1)))
  expect_equal(c(pk$peak_x, pk$peak_y, pk$peak_z), 2 * (c(pk$peak_i, pk$peak_j, pk$peak_k) - 1))
  # Peak stat is the cluster maximum of z_g.
  z <- as.vector(assay(r, "z_g"))
  expect_equal(pk$peak_stat, max(z[cid == 1]))
  # FWER p is never below uncorrected p; corrected TFCE p is in (0, 1].
  expect_true(all(assay(r, "p_fwer") >= assay(r, "p_perm")))
  expect_true(all(assay(r, "p_tfce_fwer") > 0 & assay(r, "p_tfce_fwer") <= 1))
})

test_that("two-sided TFCE is signed like t and less-tail clusters are negative", {
  dims <- c(12L, 10L, 8L)
  blob <- array(FALSE, dims)
  blob[3:6, 3:6, 2:4] <- TRUE
  fx <- .sp_fixture(seed = 9L, dims = dims, n_subj = 12L, effect = -1.5, blob = blob,
                    mask_frac = 1)
  r <- .sp_run(fx$g, n_perm = 99L, seed = 1L, alternative = "two.sided",
               spatial = c("tfce", "cluster"), cluster_p = 0.01)
  tf <- as.vector(assay(r, "tfce"))
  z <- as.vector(assay(r, "z_g"))
  expect_true(all(sign(tf[tf != 0]) == sign(z[tf != 0])))
  r2 <- .sp_run(fx$g, n_perm = 99L, seed = 1L, alternative = "less",
                spatial = "cluster", cluster_p = 0.01)
  tab <- .sp_att(r2)$clusters
  expect_true(all(tab$sign == -1L))
  expect_lt(tab$p_fwer_mass[1], 0.05)
  r3 <- .sp_run(fx$g, n_perm = 99L, seed = 1L, alternative = "greater",
                spatial = c("tfce", "cluster"), cluster_p = 0.01)
  inb <- blob[fx$mask_idx]
  expect_true(all(as.vector(assay(r3, "p_tfce_fwer"))[inb] > 0.5))
})

test_that("voxels below min_subjects are excluded and do not bridge clusters", {
  dims <- c(9L, 3L, 3L)
  n_subj <- 8L
  set.seed(10)
  mask_idx <- seq_len(prod(dims))
  # A bar along x with a strong effect; the middle voxel column x = 5 has
  # only one finite subject.
  bar <- array(FALSE, dims)
  bar[, 2, 2] <- TRUE
  beta <- matrix(stats::rnorm(prod(dims) * n_subj, sd = 0.2), prod(dims), n_subj)
  beta[which(bar), ] <- beta[which(bar), ] + 3
  gap <- which(arrayInd(mask_idx, dims)[, 1] == 5)
  beta[gap, -1] <- NA
  g <- new_gds(
    list(beta = array(beta, c(prod(dims), n_subj, 1L)),
         var = array(1, c(prod(dims), n_subj, 1L))),
    space_voxel(dims, diag(4), mask_idx = mask_idx, storage = "packed"),
    paste0("s", seq_len(n_subj)), "c1"
  )
  r <- .sp_run(g, n_perm = 49L, seed = 1L, alternative = "greater",
               spatial = c("tfce", "cluster"), cluster_threshold = 5)
  for (nm in c("t_g", "tfce", "p_tfce_fwer", "cluster_id", "p_cluster_fwer")) {
    expect_true(all(is.na(as.vector(assay(r, nm))[gap])), info = nm)
  }
  tab <- .sp_att(r)$clusters
  expect_equal(nrow(tab), 2L)
  expect_equal(sort(tab$size), c(4L, 4L))
  expect_equal(.sp_att(r)$n_excluded, length(gap))
  # With min_subjects = 1 the gap still has no t (n = 1), so it is excluded too.
  r1 <- .sp_run(g, n_perm = 19L, seed = 1L, alternative = "greater",
                spatial = "cluster", cluster_threshold = 5, min_subjects = 1L)
  expect_equal(nrow(.sp_att(r1)$clusters), 2L)
})

test_that("spatial statistic t uses t thresholds; option validation", {
  fx <- .sp_fixture(seed = 11L)
  r <- .sp_run(fx$g, n_perm = 19L, seed = 1L, spatial = "cluster", spatial_stat = "t",
               cluster_p = 0.01)
  expect_equal(.sp_att(r)$cluster$threshold, stats::qt(0.01, 9, lower.tail = FALSE))
  expect_equal(.sp_att(r)$statistic, "t")
  expect_error(.sp_run(fx$g, n_perm = 9L, spatial = "voxel"), "unknown options\\$spatial")
  expect_error(.sp_run(fx$g, n_perm = 9L, spatial = "tfce", connectivity = 8), "6, 18 or 26")
  expect_error(.sp_run(fx$g, n_perm = 9L, spatial = "tfce", tfce_dh = 0), "tfce_dh")
  expect_error(.sp_run(fx$g, n_perm = 9L, spatial = "cluster", cluster_p = 0.01,
                       cluster_threshold = 3), "not both")
  expect_error(.sp_run(fx$g, n_perm = 9L, spatial = "tfce", spatial_stat = "p"), "spatial_stat")
  g_lab <- new_gds(list(beta = array(stats::rnorm(40), c(4, 10, 1)), var = array(1, c(4, 10, 1))),
                   space_sample_labels(paste0("v", 1:4)), paste0("s", 1:10), "c1")
  expect_error(.sp_run(g_lab, n_perm = 9L, spatial = "tfce"), "voxel space")
})

test_that("connectivity 26 joins diagonal neighbours", {
  dims <- c(4L, 4L, 4L)
  n_subj <- 8L
  set.seed(12)
  beta <- matrix(stats::rnorm(64 * n_subj, sd = 0.1), 64, n_subj)
  beta[1, ] <- beta[1, ] + 3          # (1,1,1)
  beta[1 + 1 + 4 + 16, ] <- beta[22, ] + 3  # (2,2,2), corner neighbour
  g <- new_gds(list(beta = array(beta, c(64, n_subj, 1)), var = array(1, c(64, n_subj, 1))),
               space_voxel(dims, diag(4), mask_idx = 1:64, storage = "packed"),
               paste0("s", 1:n_subj), "c1")
  c6 <- .sp_att(.sp_run(g, n_perm = 9L, seed = 1L, alternative = "greater",
                        spatial = "cluster", cluster_threshold = 5))$clusters
  c26 <- .sp_att(.sp_run(g, n_perm = 9L, seed = 1L, alternative = "greater",
                         spatial = "cluster", cluster_threshold = 5, connectivity = 26))$clusters
  expect_equal(nrow(c6), 2L)
  expect_equal(nrow(c26), 1L)
  expect_equal(c26$size, 2L)
})

test_that("explicit inverse-variance weights are refused by perm:onesample (#21)", {
  fx <- .sp_fixture(seed = 13L)
  expect_error(reduce(fx$g, method = "perm:onesample", weights = "1/var"), "unweighted")
  expect_error(reduce(fx$g, method = "perm:onesample", weights = "custom"), "not supported")
  expect_s3_class(reduce(fx$g, method = "perm:onesample", weights = "equal"), "gds_plan")
  expect_s3_class(reduce(fx$g, method = "perm:onesample"), "gds_plan")
})

test_that("spatial assays export to NIfTI", {
  skip_if_not_installed("RNifti")
  fx <- .sp_fixture(seed = 14L)
  r <- .sp_run(fx$g, n_perm = 19L, seed = 1L, spatial = c("tfce", "cluster"), cluster_p = 0.05)
  out <- tempfile("perm_spatial_nifti")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  man <- write_nifti_assays(r, out, assays = c("tfce", "p_tfce_fwer", "cluster_id", "p_cluster_fwer"))
  expect_true(all(man$written))
  img <- RNifti::readNifti(man$path[man$assay == "cluster_id"])
  expect_equal(dim(img), fx$dims)
  vals <- as.numeric(img)[fx$mask_idx]
  expect_equal(vals, as.vector(assay(r, "cluster_id")))
})
