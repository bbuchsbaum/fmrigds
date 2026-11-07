# Calibration study: family-wise error of subject-level sign-flip inference
# (max-T, TFCE and cluster mass) under the global null, for each alternative.
#
#   FMRIGDS_THREADS=8 Rscript inst/examples/fmrigds_perm_spatial_calibration.R
#
# Environment variables (defaults in brackets):
#   CAL_N_DATA [1000]  simulated datasets
#   CAL_N_PERM [500]   sign flips per dataset
#   CAL_N_SUBJ [20]    subjects
#   CAL_OUT            optional .rds file for the per-dataset minima
#
# Each dataset: 20 subject maps of Gaussian noise on a 30 x 30 x 20 grid,
# smoothed with sigma = 1.5 voxels (FWHM ~ 3.5 voxels), mean zero
# everywhere, restricted to an ellipsoid mask (~9,900 voxels). A dataset
# commits a family-wise error at level alpha when any voxel has corrected
# p <= alpha. Cluster-forming threshold: one-tailed p = 0.01.
suppressPackageStartupMessages(library(fmrigds))

n_data <- as.integer(Sys.getenv("CAL_N_DATA", "1000"))
n_perm <- as.integer(Sys.getenv("CAL_N_PERM", "500"))
n_subj <- as.integer(Sys.getenv("CAL_N_SUBJ", "20"))
out_file <- Sys.getenv("CAL_OUT", "")

smooth3 <- function(a, sigma) {
  r <- ceiling(3 * sigma)
  k <- stats::dnorm(-r:r, sd = sigma)
  k <- k / sum(k)
  d <- dim(a)
  f <- function(x) as.numeric(stats::filter(c(rep(0, r), x, rep(0, r)), k, sides = 2))[(r + 1):(r + length(x))]
  a <- array(apply(a, c(2, 3), f), d)
  a <- aperm(array(apply(aperm(a, c(2, 1, 3)), c(2, 3), f), d[c(2, 1, 3)]), c(2, 1, 3))
  aperm(array(apply(aperm(a, c(3, 1, 2)), c(2, 3), f), d[c(3, 1, 2)]), c(2, 3, 1))
}

dims <- c(30L, 30L, 20L)
g <- expand.grid(i = seq_len(dims[1]), j = seq_len(dims[2]), k = seq_len(dims[3]))
c0 <- (dims + 1) / 2
r <- 0.95 * (dims - 1) / 2
mask_idx <- which(((g$i - c0[1]) / r[1])^2 + ((g$j - c0[2]) / r[2])^2 +
                    ((g$k - c0[3]) / r[3])^2 <= 1)
space <- space_voxel(dims, diag(4), mask_idx = mask_idx, storage = "packed")
cat(sprintf("%d datasets, %d sign flips, %d subjects, %d voxels, threads = %d\n",
            n_data, n_perm, n_subj, length(mask_idx), fmrigds:::.resolve_thread_count()))

alts <- c("greater", "less", "two.sided")
res <- array(NA_real_, c(n_data, length(alts), 3),
             dimnames = list(NULL, alts, c("max_t", "tfce", "cluster_mass")))
set.seed(20261107)
t0 <- proc.time()[["elapsed"]]
for (d in seq_len(n_data)) {
  beta <- vapply(seq_len(n_subj), function(s) {
    a <- smooth3(array(stats::rnorm(prod(dims)), dims), 1.5)
    a[mask_idx]
  }, numeric(length(mask_idx)))
  x <- new_gds(
    list(beta = array(beta, c(length(mask_idx), n_subj, 1L)),
         var = array(1, c(length(mask_idx), n_subj, 1L))),
    space, paste0("s", seq_len(n_subj)), "c1"
  )
  for (alt in alts) {
    out <- compute(reduce(x, method = "perm:onesample", weights = "equal",
                          options = list(n_perm = n_perm, seed = d, alternative = alt,
                                         spatial = c("tfce", "cluster"),
                                         cluster_p = 0.01)))
    res[d, alt, ] <- c(min(assay(out, "p_fwer"), na.rm = TRUE),
                       min(assay(out, "p_tfce_fwer"), na.rm = TRUE),
                       min(assay(out, "p_cluster_fwer"), na.rm = TRUE))
  }
  if (d %% 50 == 0) {
    cat(sprintf("  %d datasets, %.0f s\n", d, proc.time()[["elapsed"]] - t0))
  }
}
if (nzchar(out_file)) saveRDS(res, out_file)

summarise <- function(alpha) {
  fw <- apply(res <= alpha, c(2, 3), mean)
  se <- sqrt(alpha * (1 - alpha) / n_data)
  cat(sprintf("\nFamily-wise error rate at alpha = %.2f (Monte Carlo SE %.4f; 95%% band %.3f-%.3f)\n",
              alpha, se, alpha - 1.96 * se, alpha + 1.96 * se))
  print(round(fw, 4))
}
summarise(0.05)
summarise(0.10)
