#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
})

pkgload::load_all(".", compile = FALSE, recompile = FALSE)
source("inst/examples/fmrigds_lmm_benchmark.R")

res <- run_fmrigds_lmm_benchmark(
  n_samples = as.integer(Sys.getenv("FMRIGDS_BENCH_SAMPLES", "64")),
  n_subjects = as.integer(Sys.getenv("FMRIGDS_BENCH_SUBJECTS", "24"))
)

saveRDS(res, "lmm-benchmark.rds")

summary_lines <- c(
  "# Repeated-measures LMM benchmark",
  "",
  sprintf("- pooled elapsed: %.3fs", res$pooled$elapsed),
  sprintf("- voxelwise elapsed: %.3fs", res$voxelwise$elapsed)
)

if (!is.null(res$lme4_reference)) {
  summary_lines <- c(
    summary_lines,
    sprintf("- lme4 reference elapsed (%d samples): %.3fs", res$lme4_reference$samples, res$lme4_reference$elapsed)
  )
}

writeLines(summary_lines, "lmm-benchmark-summary.md")
cat(paste(summary_lines, collapse = "\n"), "\n")
