.ldr_test_fixture <- function(mode = c("shifted", "aligned", "polarity", "null"),
                              seed = 5101L,
                              n_subject = 30L,
                              noise_sd = 0.08) {
  mode <- match.arg(mode)
  set.seed(seed)
  dims3 <- c(11L, 1L, 1L)
  affine <- diag(c(2, 2, 2, 1))
  center <- 6L
  geometry <- .ldr_patch_geometry(
    center, dims3, affine,
    patch_radius_mm = 8,
    shift_radius_mm = 2
  )
  feature <- numeric(prod(dims3))
  feature[5:7] <- c(-0.6, 1.2, -0.6)
  shifts <- rep(c(-1L, 0L, 1L), length.out = n_subject)
  signs <- rep(c(-1, 1), length.out = n_subject)
  shift_vector <- function(x, dx) {
    out <- numeric(length(x))
    target <- seq_along(x)
    source <- target - dx
    valid <- source >= 1L & source <= length(x)
    out[valid] <- x[source[valid]]
    out
  }
  beta_a <- beta_b <- matrix(0, nrow = prod(dims3), ncol = n_subject)
  for (i in seq_len(n_subject)) {
    truth <- switch(
      mode,
      shifted = shift_vector(feature, shifts[i]),
      aligned = feature,
      polarity = signs[i] * feature,
      null = numeric(length(feature))
    )
    beta_a[, i] <- truth + stats::rnorm(length(truth), sd = noise_sd)
    beta_b[, i] <- truth + stats::rnorm(length(truth), sd = noise_sd)
  }
  variance <- matrix(noise_sd^2, nrow = prod(dims3), ncol = n_subject)
  patch <- .ldr_prepare_patch(
    beta_a[geometry$sample_idx, , drop = FALSE],
    beta_b[geometry$sample_idx, , drop = FALSE],
    variance[geometry$sample_idx, , drop = FALSE],
    variance[geometry$sample_idx, , drop = FALSE],
    geometry
  )
  list(
    patch = patch,
    subjects = paste0("s", seq_len(n_subject)),
    geometry = geometry
  )
}

.ldr_test_fit <- function(fixture) {
  .ldr_crossfit(
    fixture$patch,
    fixture$subjects,
    tau_grid_mm = c(1, 2, 3),
    folds = 3L,
    iterations = 8L
  )
}
