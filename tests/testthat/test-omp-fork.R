# Forked workers must complete when a single native thread is requested via the
# standard OMP_NUM_THREADS variable (#32). Runs in a fresh, time-bounded process.

test_that("OMP_NUM_THREADS=1 keeps forked permutation workers from hanging", {
  skip_on_cran()
  skip_on_os(c("windows", "mac", "solaris"))
  skip_if(Sys.which("timeout") == "", "timeout(1) not available")
  skip_if_not(requireNamespace("fmrigds", quietly = TRUE), "fmrigds not installed")

  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "library(fmrigds)",
    "set.seed(11)",
    "Y <- matrix(rnorm(16 * 61), 16)",
    "signs <- matrix(sample(c(-1L, 1L), 1000 * 16, TRUE), 1000, 16)",
    "f <- get_reducer('perm:onesample')$fun",
    "run <- function() f(Y, NULL, NULL, NULL, NULL, NULL, NULL, NULL,",
    "  list(signs = signs, n_perm = 1000, alternative = 'two.sided', min_subjects = 16))",
    "invisible(run())",
    "z <- parallel::mclapply(1:2, function(i) run(), mc.cores = 2)",
    "cat('workers complete\\n')"
  ), script)

  out <- suppressWarnings(system2(
    "timeout", c("60", file.path(R.home("bin"), "Rscript"), script),
    stdout = TRUE, stderr = TRUE,
    env = c("OMP_NUM_THREADS=1", "FMRIGDS_THREADS=", "OMP_THREAD_LIMIT=")
  ))
  expect_null(attr(out, "status"))
  expect_true("workers complete" %in% out)
})
