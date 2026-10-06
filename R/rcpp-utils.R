# OpenMP thread control -----------------------------------------------------
#
# By default fmrigds uses at most 2 OpenMP threads so that loading the package
# never grabs every core on a shared machine. The thread count is resolved as:
#   1. `options(fmrigds.threads = n)` if set to a positive number;
#   2. otherwise the `FMRIGDS_THREADS` environment variable if set;
#   3. otherwise the default of 2.
# The result is capped by `OMP_THREAD_LIMIT` when that is set and smaller.

.positive_int_or_null <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NULL)
  x <- suppressWarnings(as.integer(x))
  if (is.na(x) || x <= 0L) return(NULL)
  x
}

.resolve_thread_count <- function(default = 2L) {
  n <- .positive_int_or_null(getOption("fmrigds.threads", default = NULL))
  if (is.null(n)) {
    env <- Sys.getenv("FMRIGDS_THREADS", unset = "")
    if (nzchar(env)) n <- .positive_int_or_null(env)
  }
  if (is.null(n)) n <- as.integer(default)
  limit <- Sys.getenv("OMP_THREAD_LIMIT", unset = "")
  if (nzchar(limit)) {
    limit <- .positive_int_or_null(limit)
    if (!is.null(limit)) n <- min(n, limit)
  }
  n
}

# Compiled kernels look up Rcpp runtime callables (e.g. RNGScope) via
# R_GetCCallable("Rcpp", ...), which fails unless the Rcpp namespace is loaded;
# a failed lookup inside a static initializer then hangs every later call.
# Rcpp is only in LinkingTo, so make sure its namespace is loaded at load time.
.ensure_rcpp_runtime <- function() {
  pkg <- "Rcpp"
  if (!isNamespaceLoaded(pkg)) {
    try(loadNamespace(pkg), silent = TRUE)
  }
  invisible(isNamespaceLoaded(pkg))
}

.set_threads_from_option <- function() {
  .ensure_rcpp_runtime()
  n <- .resolve_thread_count()
  if (exists("set_omp_threads", mode = "function")) {
    try(set_omp_threads(as.integer(n)), silent = TRUE)
  }
  invisible(n)
}
