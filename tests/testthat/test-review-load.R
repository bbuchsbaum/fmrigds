# Regression: compiled reducers must work in a fresh session where only
# fmrigds has been attached (Rcpp must be imported, not merely LinkingTo).
test_that("C++ reducers run in a fresh session without Rcpp preloaded", {
  skip_on_cran()
  skip_if(
    isNamespaceLoaded("pkgload") && pkgload::is_dev_package("fmrigds"),
    "requires the installed package"
  )
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript))
  code <- paste(
    "library(fmrigds)",
    "tmp <- tempfile(fileext = '.csv')",
    "writeLines(c('sample,subject,contrast,beta,var',",
    "  'ROI_1,sub-01,faces,0.5,0.04', 'ROI_1,sub-02,faces,0.7,0.05'), tmp)",
    "for (i in 1:2) r <- compute(reduce(gds(tmp), method = 'fixed'))",
    "cat('FMRIGDS_OK')",
    sep = "\n"
  )
  script <- tempfile(fileext = ".R")
  writeLines(code, script)
  out <- suppressWarnings(system2(rscript, c("--vanilla", script), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("FMRIGDS_OK", out)), info = paste(out, collapse = "\n"))
})
