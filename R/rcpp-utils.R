.set_threads_from_option <- function() {
  n <- getOption("fmrigds.threads", default = NULL)
  if (is.null(n)) return(invisible())
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0) return(invisible())
  if (exists("set_omp_threads", mode = "function")) {
    try(set_omp_threads(as.integer(n)), silent = TRUE)
  }
}

