# Deterministic soft performance baseline for the exact LDR ROI reference.
# Run from the package root with:
#   Rscript inst/benchmarks/local-displacement-rescue.R

pkgload::load_all(quiet = TRUE)
set.seed(5301)

dims3 <- c(11L, 11L, 11L)
affine <- diag(c(2, 2, 2, 1))
center_coord <- c(6L, 6L, 6L)
center <- center_coord[1L] +
  (center_coord[2L] - 1L) * dims3[1L] +
  (center_coord[3L] - 1L) * dims3[1L] * dims3[2L]
n_subject <- 24L
n_voxel <- prod(dims3)

signal <- array(0, dims3)
signal[5:7, 5:7, 5:7] <- -0.4
signal[6L, 6L, 6L] <- 1.2
shift_x <- function(x, dx) {
  out <- array(0, dim(x))
  target <- arrayInd(seq_len(length(x)), dim(x))
  source <- sweep(target, 2L, c(dx, 0L, 0L), "-")
  valid <- rowSums(source < 1L | source > rep(dim(x), each = nrow(source))) == 0L
  source_idx <- source[valid, 1L] +
    (source[valid, 2L] - 1L) * dim(x)[1L] +
    (source[valid, 3L] - 1L) * dim(x)[1L] * dim(x)[2L]
  out[valid] <- x[source_idx]
  out
}

beta_a <- beta_b <- array(0, c(n_voxel, n_subject, 1L))
shifts <- rep(c(-1L, 0L, 1L), length.out = n_subject)
for (i in seq_len(n_subject)) {
  truth <- as.numeric(shift_x(signal, shifts[i]))
  beta_a[, i, 1L] <- truth + stats::rnorm(n_voxel, sd = 0.1)
  beta_b[, i, 1L] <- truth + stats::rnorm(n_voxel, sd = 0.1)
}
variance <- array(0.1^2, dim(beta_a))
make_split <- function(beta) {
  new_gds(
    list(beta = beta, var = variance),
    space_voxel(dims3, affine, storage = "dense"),
    paste0("s", seq_len(n_subject)),
    "task"
  )
}

timing <- system.time({
  result <- local_displacement_rescue(
    make_split(beta_a),
    make_split(beta_b),
    center = center,
    patch_radius_mm = 4,
    shift_radius_mm = 2,
    tau_grid_mm = c(1, 2),
    folds = 3L,
    min_subjects = 12L,
    iterations = 8L
  )
})

cat(sprintf(
  paste0(
    "patch=3d subjects=%d elapsed_seconds=%.3f ",
    "activation_score=%.3f displacement_score=%.3f mni_loss=%.3f\n"
  ),
  n_subject,
  unname(timing[["elapsed"]]),
  metadata(result)$ldr$activation_score,
  metadata(result)$ldr$displacement_score,
  assay(result, "mni_loss")[center, 1L, 1L]
))
