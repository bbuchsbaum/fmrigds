# Deterministic 3D benchmark for the LDR v0.2 tangent search.

suppressPackageStartupMessages(library(fmrigds))

set.seed(5601L)
dims3 <- c(9L, 9L, 9L)
affine <- diag(c(2, 2, 2, 1))
center_coord <- c(5L, 5L, 5L)
center <- center_coord[1L] +
  (center_coord[2L] - 1L) * dims3[1L] +
  (center_coord[3L] - 1L) * dims3[1L] * dims3[2L]
n_subject <- 28L
subjects <- paste0("s", seq_len(n_subject))
noise_sd <- 0.08

feature <- array(0, dims3)
feature[center] <- 1.2
neighbor <- rbind(
  c(-1L, 0L, 0L), c(1L, 0L, 0L),
  c(0L, -1L, 0L), c(0L, 1L, 0L),
  c(0L, 0L, -1L), c(0L, 0L, 1L)
)
neighbor_coord <- sweep(neighbor, 2L, center_coord, "+")
neighbor_index <- neighbor_coord[, 1L] +
  (neighbor_coord[, 2L] - 1L) * dims3[1L] +
  (neighbor_coord[, 3L] - 1L) * dims3[1L] * dims3[2L]
feature[neighbor_index] <- -0.35
subject_shift <- neighbor[c(1L, 2L, 3L, 4L, 5L, 6L, 1L), , drop = FALSE]
subject_shift[7L, ] <- 0L

shift_array <- function(x, displacement) {
  source <- arrayInd(seq_len(prod(dim(x))), dim(x))
  target <- sweep(source, 2L, displacement, "+")
  valid <- rowSums(target < 1L | target > rep(dim(x), each = nrow(target))) == 0L
  source_index <- which(valid)
  target <- target[valid, , drop = FALSE]
  target_index <- target[, 1L] +
    (target[, 2L] - 1L) * dim(x)[1L] +
    (target[, 3L] - 1L) * dim(x)[1L] * dim(x)[2L]
  out <- array(0, dim(x))
  out[target_index] <- x[source_index]
  out
}

beta_a <- beta_b <- matrix(0, prod(dims3), n_subject)
for (i in seq_len(n_subject)) {
  truth <- shift_array(
    feature,
    subject_shift[(i - 1L) %% nrow(subject_shift) + 1L, ]
  )
  beta_a[, i] <- truth + rnorm(prod(dims3), sd = noise_sd)
  beta_b[, i] <- truth + rnorm(prod(dims3), sd = noise_sd)
}
variance <- matrix(noise_sd^2, prod(dims3), n_subject)
make_split <- function(beta) {
  new_gds(
    assays = list(
      beta = array(beta, c(prod(dims3), n_subject, 1L)),
      var = array(variance, c(prod(dims3), n_subject, 1L))
    ),
    space = space_voxel(dims3, affine),
    subjects = subjects,
    contrasts = "task"
  )
}

search_coord <- as.matrix(expand.grid(
  x = 4:6, y = 4:6, z = 4:6
))
search_centers <- search_coord[, 1L] +
  (search_coord[, 2L] - 1L) * dims3[1L] +
  (search_coord[, 3L] - 1L) * dims3[1L] * dims3[2L]

elapsed <- system.time({
  result <- local_displacement_rescue_map(
    make_split(beta_a),
    make_split(beta_b),
    centers = search_centers,
    patch_radius_mm = 4,
    shift_radius_mm = 2,
    folds = 4L,
    min_subjects = 12L,
    covariance_shrinkage = 0.5,
    n_resamples = 3L,
    seed = 5602L,
    alpha = 0.5,
    iterations = 8L,
    max_refits = 4L
  )
})

cat("subjects:", n_subject, "\n")
cat("searched centers:", length(search_centers), "\n")
cat("resamples:", 3L, "\n")
cat("elapsed seconds:", unname(elapsed[["elapsed"]]), "\n")
cat("maximum tangent activation t:",
    max(assay(result, "tangent_activation_t"), na.rm = TRUE), "\n")
cat("maximum tangent displacement t:",
    max(assay(result, "tangent_displacement_t"), na.rm = TRUE), "\n")
cat("nonzero LDR flags:", sum(assay(result, "ldr_flag") > 0), "\n")
