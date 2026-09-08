# Deterministic soft performance baseline for the experimental diagnostics.
# Run from the package root with:
#   Rscript inst/benchmarks/experimental-cancellation.R

pkgload::load_all(quiet = TRUE)
set.seed(4201)

dims3 <- c(32L, 32L, 16L)
n_voxel <- prod(dims3)
n_subject <- 32L
subjects <- paste0("s", seq_len(n_subject))

beta <- array(rnorm(n_voxel * n_subject, sd = 0.05),
              c(n_voxel, n_subject, 1L))
variance <- array(0.05^2, dim(beta))

signal <- array(0, dims3)
signal[15:18, 15:18, 8:9] <- 0.4
signal[16:17, 16:17, 8:9] <- 0.9
shifts <- rep(c(-1L, 0L, 1L), length.out = n_subject)
shift_x <- function(x, dx) {
  out <- array(0, dim(x))
  target <- arrayInd(seq_len(length(x)), dim(x))
  source <- sweep(target, 2L, c(dx, 0L, 0L), "+")
  ok <- rowSums(source < 1L | source > rep(dim(x), each = nrow(source))) == 0L
  source_idx <- source[ok, 1L] +
    (source[ok, 2L] - 1L) * dim(x)[1L] +
    (source[ok, 3L] - 1L) * dim(x)[1L] * dim(x)[2L]
  out[ok] <- x[source_idx]
  out
}
for (i in seq_len(n_subject)) {
  truth <- as.numeric(shift_x(signal, shifts[i]))
  beta[, i, 1L] <- truth + rnorm(n_voxel, sd = 0.05)
}

subject_maps <- new_gds(
  list(beta = beta, var = variance),
  space_voxel(dims3, diag(4), storage = "dense"),
  subjects,
  "task"
)

timing <- system.time({
  result <- experimental_cancellation(
    subject_maps,
    delta = 0.15,
    equivalence = 0.2,
    shift_radius = 1L,
    patch_radius = 1L
  )
})

cat(
  sprintf(
    "voxels=%d subjects=%d elapsed_seconds=%.3f finite_rescue=%d\n",
    n_voxel,
    n_subject,
    unname(timing[["elapsed"]]),
    sum(is.finite(assay(result, "shift_rescue_fraction")))
  )
)
