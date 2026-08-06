#!/usr/bin/env Rscript

if (!("fmrigds" %in% loadedNamespaces())) {
  suppressPackageStartupMessages(library(fmrigds))
}

make_knownvar_benchmark_df <- function(n_samples = 50L, n_subjects = 30L) {
  set.seed(20260806)
  subjects <- paste0("sub-", sprintf("%03d", seq_len(n_subjects)))
  contrasts <- c("t0", "t1", "t2")
  time <- c(-1, 0, 1)
  samples <- paste0("ROI_", sprintf("%03d", seq_len(n_samples)))
  df <- expand.grid(
    sample = samples,
    subject = subjects,
    contrast = contrasts,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  sample_idx <- match(df$sample, samples)
  subject_idx <- match(df$subject, subjects)
  time_idx <- match(df$contrast, contrasts)
  b0 <- rnorm(n_subjects, sd = 0.25)
  b1 <- rnorm(n_subjects, sd = 0.10)
  df$time <- time[time_idx]
  df$var <- runif(nrow(df), min = 0.008, max = 0.05)
  df$beta <- 0.8 +
    0.01 * sample_idx +
    b0[subject_idx] +
    (0.45 + b1[subject_idx]) * df$time +
    rnorm(nrow(df), sd = sqrt(df$var + 0.03))
  df
}

run_fmrigds_lmm_knownvar_benchmark <- function(n_samples = 50L, n_subjects = 30L) {
  df <- make_knownvar_benchmark_df(n_samples, n_subjects)
  elapsed <- system.time({
    result <- gds(df, contrast_data_cols = "time") |>
      reduce(
        method = "lmm:ri_slope1_knownvar",
        formula = ~ time,
        options = list(slope = "time", covariance = "diag", fit = "REML")
      ) |>
      compute()
  })[["elapsed"]]

  list(
    elapsed = elapsed,
    samples = n_samples,
    subjects = n_subjects,
    converged = sum(as.numeric(assay(result, "converged")) == 1)
  )
}

if (sys.nframe() == 0L) {
  print(run_fmrigds_lmm_knownvar_benchmark())
}
