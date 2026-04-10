#!/usr/bin/env Rscript

if (!("fmrigds" %in% loadedNamespaces())) {
  suppressPackageStartupMessages({
    library(fmrigds)
  })
}

make_rm_benchmark_df <- function(n_samples = 200L, n_subjects = 40L, time = c(-1, 0, 1)) {
  stopifnot(n_samples >= 1L, n_subjects >= 2L, length(time) >= 2L)

  subjects <- paste0("sub-", sprintf("%03d", seq_len(n_subjects)))
  contrasts <- paste0("t", seq_along(time) - 1L)
  samples <- paste0("ROI_", seq_len(n_samples))

  sample_shift <- seq(0, 0.6, length.out = n_samples)
  sample_slope <- seq(0.3, 0.9, length.out = n_samples)
  sample_theta <- seq(0.05, 0.35, length.out = n_samples)
  subject_intercept <- stats::rnorm(n_subjects, sd = 0.25)
  subject_slope <- stats::rnorm(n_subjects, sd = 0.12)

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
  time_val <- time[time_idx]

  df$time <- time_val
  df$beta <- 0.8 +
    sample_shift[sample_idx] +
    subject_intercept[subject_idx] +
    (sample_slope[sample_idx] + subject_slope[subject_idx]) * time_val +
    stats::rnorm(nrow(df), sd = sample_theta[sample_idx])
  df$var <- 0.04
  df
}

run_fmrigds_lmm_benchmark <- function(n_samples = 200L, n_subjects = 40L, theta_modes = c("pooled", "voxelwise")) {
  df <- make_rm_benchmark_df(n_samples = n_samples, n_subjects = n_subjects)
  timings <- vector("list", length(theta_modes))
  names(timings) <- theta_modes

  for (mode in theta_modes) {
    elapsed <- system.time({
      gout <- gds(df, contrast_data_cols = "time") |>
        reduce(
          method = "lmm:ri_slope1",
          formula = ~ time,
          options = list(
            slope = "time",
            covariance = "diag",
            fit = "REML",
            theta_mode = mode
          )
        ) |>
        compute()
    })[["elapsed"]]

    timings[[mode]] <- list(
      elapsed = elapsed,
      assays = names(assays(gout))
    )
  }

  if (requireNamespace("lme4", quietly = TRUE)) {
    n_ref <- min(10L, n_samples)
    df_ref <- df[df$sample %in% paste0("ROI_", seq_len(n_ref)), , drop = FALSE]
    ref_elapsed <- system.time({
      split(df_ref, df_ref$sample)
      invisible(lapply(split(df_ref, df_ref$sample), function(d) {
        lme4::lmer(beta ~ time + (1 + time | subject), data = d, REML = TRUE)
      }))
    })[["elapsed"]]
    timings$lme4_reference <- list(
      elapsed = ref_elapsed,
      samples = n_ref
    )
  }

  timings
}

if (sys.nframe() == 0L) {
  res <- run_fmrigds_lmm_benchmark()
  print(res)
}
