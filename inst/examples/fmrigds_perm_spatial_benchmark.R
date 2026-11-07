# Benchmark: whole-brain subject-level sign-flip inference with
# reduce(method = "perm:onesample", options = list(spatial = ...)).
#
#   FMRIGDS_THREADS=8 Rscript inst/examples/fmrigds_perm_spatial_benchmark.R
#
# Environment variables (defaults in brackets):
#   BENCH_N_PERM [5000]  number of sign flips
#   BENCH_N_SUBJ [35]    subjects
#   BENCH_SPATIAL [tfce,cluster]
#   BENCH_STAT [z]       spatial statistic, "z" or "t"
#
# Data: a 91 x 109 x 91 grid (2 mm MNI size) with an ellipsoid mask of
# ~244k voxels; each subject map is Gaussian noise smoothed with
# sigma = 2 voxels (FWHM ~ 9.4 mm).
suppressPackageStartupMessages(library(fmrigds))

n_perm <- as.integer(Sys.getenv("BENCH_N_PERM", "5000"))
n_subj <- as.integer(Sys.getenv("BENCH_N_SUBJ", "35"))
spatial <- strsplit(Sys.getenv("BENCH_SPATIAL", "tfce,cluster"), ",")[[1]]
stat <- Sys.getenv("BENCH_STAT", "z")

smooth3 <- function(a, sigma) {
  # Separable Gaussian smoothing with zero padding (base R only).
  r <- ceiling(3 * sigma)
  k <- stats::dnorm(-r:r, sd = sigma)
  k <- k / sum(k)
  d <- dim(a)
  f <- function(x) as.numeric(stats::filter(c(rep(0, r), x, rep(0, r)), k, sides = 2))[(r + 1):(r + length(x))]
  a <- array(apply(a, c(2, 3), f), d)
  a <- aperm(array(apply(aperm(a, c(2, 1, 3)), c(2, 3), f), d[c(2, 1, 3)]), c(2, 1, 3))
  aperm(array(apply(aperm(a, c(3, 1, 2)), c(2, 3), f), d[c(3, 1, 2)]), c(2, 3, 1))
}

dims <- c(91L, 109L, 91L)
g <- expand.grid(i = seq_len(dims[1]), j = seq_len(dims[2]), k = seq_len(dims[3]))
c0 <- (dims + 1) / 2
r <- 0.81 * (dims - 1) / 2
mask <- ((g$i - c0[1]) / r[1])^2 + ((g$j - c0[2]) / r[2])^2 + ((g$k - c0[3]) / r[3])^2 <= 1
mask_idx <- which(mask)
cat(sprintf("grid %s, %d in-mask voxels, %d subjects, %d sign flips, spatial = %s, stat = %s\n",
            paste(dims, collapse = "x"), length(mask_idx), n_subj, n_perm,
            paste(spatial, collapse = "+"), stat))
cat(sprintf("fmrigds threads: %s\n", fmrigds:::.resolve_thread_count()))

set.seed(1)
t_sim <- system.time({
  beta <- vapply(seq_len(n_subj), function(s) {
    a <- smooth3(array(stats::rnorm(prod(dims)), dims), 2)
    a[mask_idx] / stats::sd(a[mask_idx])
  }, numeric(length(mask_idx)))
})
cat(sprintf("simulated data in %.1f s\n", t_sim[["elapsed"]]))

gds_obj <- new_gds(
  list(beta = array(beta, c(length(mask_idx), n_subj, 1)),
       var = array(1, c(length(mask_idx), n_subj, 1))),
  space_voxel(dims, diag(c(2, 2, 2, 1)), mask_idx = mask_idx, storage = "packed"),
  paste0("s", seq_len(n_subj)), "c1"
)
rm(beta)

timing <- system.time({
  res <- compute(reduce(gds_obj, method = "perm:onesample", weights = "equal",
                        options = list(n_perm = n_perm, seed = 1L,
                                       alternative = "two.sided",
                                       spatial = spatial, spatial_stat = stat)))
})
cat(sprintf("perm:onesample (max-T + %s): %.1f s elapsed, %.1f s CPU\n",
            paste(spatial, collapse = " + "), timing[["elapsed"]],
            timing[["user.self"]] + timing[["sys.self"]]))
att <- metadata(res)$attachments[["reduce/perm:onesample/contrast=c1"]]
if (!is.null(att)) {
  cat(sprintf("null maxima (95th pct): |t| %.2f", stats::quantile(att$null$max_t, 0.95)))
  if (!is.null(att$null$tfce)) cat(sprintf(", TFCE %.0f", stats::quantile(att$null$tfce, 0.95)))
  if (!is.null(att$null$cluster_mass)) cat(sprintf(", cluster mass %.1f", stats::quantile(att$null$cluster_mass, 0.95)))
  cat("\n")
}
