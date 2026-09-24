# Local displacement rescue: ROI resampling calibration -------------------

.ldr_with_seed <- function(seed, code) {
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.ldr_rebuild_patch <- function(patch, y_a, y_b) {
  .ldr_prepare_patch(
    y_a,
    y_b,
    patch$var_a,
    patch$var_b,
    patch$geometry
  )
}

.ldr_sign_flip_patch <- function(patch, signs) {
  signs <- as.numeric(signs)
  if (length(signs) != ncol(patch$y_a) || any(!signs %in% c(-1, 1))) {
    stop("LDR sign flips must supply one -1/+1 value per subject.",
         call. = FALSE)
  }
  .ldr_rebuild_patch(
    patch,
    sweep(patch$y_a, 2L, signs, "*"),
    sweep(patch$y_b, 2L, signs, "*")
  )
}

.ldr_simulate_aligned_patch <- function(patch, fit) {
  n_subject <- ncol(patch$y_a)
  n_voxel <- nrow(patch$y_a)
  active <- stats::rbinom(n_subject, size = 1L, prob = fit$prevalence)
  amplitude <- switch(
    fit$amplitude_family,
    positive = abs(stats::rnorm(n_subject, sd = fit$amplitude_scale)),
    signed = stats::rnorm(n_subject, sd = fit$amplitude_scale),
    stop("Unknown LDR amplitude family.", call. = FALSE)
  )
  latent <- outer(fit$template, active * amplitude)
  draw_split <- function(variance) {
    noise <- matrix(
      stats::rnorm(n_voxel * n_subject),
      nrow = n_voxel,
      ncol = n_subject
    ) * sqrt(variance)
    noise[!is.finite(variance) | variance <= 0] <- NA_real_
    latent + noise
  }
  .ldr_rebuild_patch(
    patch,
    draw_split(patch$var_a),
    draw_split(patch$var_b)
  )
}

.ldr_resample_p <- function(observed, null) {
  null <- null[is.finite(null)]
  if (length(observed) != 1L || is.na(observed) || !length(null)) {
    return(NA_real_)
  }
  (1 + sum(null >= observed)) / (length(null) + 1)
}

.ldr_calibrate <- function(patch,
                           subject_ids,
                           observed,
                           tau_grid_mm,
                           folds,
                           iterations,
                           ridge,
                           tolerance,
                           n_resamples,
                           seed) {
  n_resamples <- as.integer(n_resamples)
  if (length(n_resamples) != 1L || is.na(n_resamples) ||
      n_resamples < 1L) {
    stop("`n_resamples` must be one positive integer for calibration.",
         call. = FALSE)
  }
  indices <- seq_along(subject_ids)
  aligned_null_fit <- .ldr_fit_model(
    patch, indices, tau_grid_mm = 0,
    iterations = iterations,
    ridge = ridge,
    tolerance = tolerance,
    amplitude_family = "positive"
  )
  heterogeneity_null_fit <- .ldr_fit_model(
    patch, indices, tau_grid_mm = 0,
    iterations = iterations,
    ridge = ridge,
    tolerance = tolerance,
    amplitude_family = "signed"
  )
  run_crossfit <- function(null_patch) {
    .ldr_crossfit(
      null_patch,
      subject_ids = subject_ids,
      tau_grid_mm = tau_grid_mm,
      folds = folds,
      iterations = iterations,
      ridge = ridge,
      tolerance = tolerance
    )
  }

  null <- list(
    activation = rep(NA_real_, n_resamples),
    displacement_aligned = rep(NA_real_, n_resamples),
    displacement_heterogeneity = rep(NA_real_, n_resamples),
    loss_aligned = rep(NA_real_, n_resamples),
    loss_heterogeneity = rep(NA_real_, n_resamples)
  )
  .ldr_with_seed(seed, {
    for (b in seq_len(n_resamples)) {
      signs <- sample(c(-1, 1), length(subject_ids), replace = TRUE)
      activation <- run_crossfit(.ldr_sign_flip_patch(patch, signs))
      aligned <- run_crossfit(
        .ldr_simulate_aligned_patch(patch, aligned_null_fit)
      )
      heterogeneity <- run_crossfit(
        .ldr_simulate_aligned_patch(patch, heterogeneity_null_fit)
      )
      null$activation[b] <- activation$activation_t
      null$displacement_aligned[b] <-
        aligned$displacement_t_components[["aligned"]]
      null$displacement_heterogeneity[b] <-
        heterogeneity$displacement_t_components[["heterogeneity"]]
      null$loss_aligned[b] <- aligned$mni_loss
      null$loss_heterogeneity[b] <- heterogeneity$mni_loss
    }
  })

  component_p <- c(
    activation = .ldr_resample_p(
      observed$activation_t, null$activation
    ),
    displacement_aligned = .ldr_resample_p(
      observed$displacement_t_components[["aligned"]],
      null$displacement_aligned
    ),
    displacement_heterogeneity = .ldr_resample_p(
      observed$displacement_t_components[["heterogeneity"]],
      null$displacement_heterogeneity
    ),
    loss_aligned = .ldr_resample_p(
      observed$mni_loss, null$loss_aligned
    ),
    loss_heterogeneity = .ldr_resample_p(
      observed$mni_loss, null$loss_heterogeneity
    )
  )
  p_displacement <- max(
    component_p[["displacement_aligned"]],
    component_p[["displacement_heterogeneity"]]
  )
  p_loss <- max(
    component_p[["loss_aligned"]],
    component_p[["loss_heterogeneity"]]
  )
  p_ldr <- max(
    component_p[["activation"]],
    p_displacement,
    p_loss
  )
  list(
    p_ldr = p_ldr,
    p_activation = component_p[["activation"]],
    p_displacement = p_displacement,
    p_loss = p_loss,
    component_p = component_p,
    null = null,
    n_resamples = n_resamples,
    seed = as.integer(seed),
    minimum_p = 1 / (n_resamples + 1)
  )
}
