# Variance-aware repeated-measures Gaussian LMM reducers -------------------

.lmm_knownvar_components <- function(par, kind = c("ri", "ri_slope1"), covariance = c("diag", "full")) {
  kind <- match.arg(kind)
  covariance <- match.arg(covariance)

  if (identical(kind, "ri")) {
    if (length(par) != 2L) stop("known-variance random-intercept model requires two parameters", call. = FALSE)
    sd_intercept <- exp(par[[1L]])
    sd_resid <- exp(par[[2L]])
    factor_mat <- matrix(sd_intercept, nrow = 1L, ncol = 1L)
    return(list(
      factor = factor_mat,
      residual_var = sd_resid^2,
      vc_intercept = sd_intercept^2,
      vc_slope = NULL,
      vc_cov_intercept_slope = NULL,
      corr_intercept_slope = NULL
    ))
  }

  expected <- if (identical(covariance, "diag")) 3L else 4L
  if (length(par) != expected) {
    stop("known-variance random-slope model received the wrong parameter count", call. = FALSE)
  }
  sd_intercept <- exp(par[[1L]])
  if (identical(covariance, "diag")) {
    sd_slope <- exp(par[[2L]])
    corr <- 0
    factor_mat <- diag(c(sd_intercept, sd_slope), nrow = 2L)
    sd_resid <- exp(par[[3L]])
  } else {
    corr <- tanh(par[[2L]])
    sd_slope <- exp(par[[3L]])
    factor_mat <- matrix(
      c(
        sd_intercept, 0,
        corr * sd_slope, sd_slope * sqrt(max(1 - corr^2, 0))
      ),
      nrow = 2L,
      byrow = TRUE
    )
    sd_resid <- exp(par[[4L]])
  }

  list(
    factor = factor_mat,
    residual_var = sd_resid^2,
    vc_intercept = sd_intercept^2,
    vc_slope = sd_slope^2,
    vc_cov_intercept_slope = corr * sd_intercept * sd_slope,
    corr_intercept_slope = corr
  )
}

.lmm_knownvar_start <- function(y, sampling_var, W_block, kind, covariance) {
  observed_var <- stats::var(y)
  sampling_scale <- stats::median(sampling_var)
  total_scale <- max(observed_var, sampling_scale, 1e-8)
  excess_scale <- max(observed_var - sampling_scale, 0.1 * total_scale, 1e-8)
  residual_var <- 0.35 * excess_scale

  if (identical(kind, "ri")) {
    return(c(
      log(sqrt(0.65 * excess_scale)),
      log(sqrt(residual_var))
    ))
  }

  slope_scale <- stats::var(W_block[, 2L])
  if (!is.finite(slope_scale) || slope_scale < 1e-4) slope_scale <- 1
  slope_var <- 0.25 * excess_scale / slope_scale
  intercept_var <- 0.40 * excess_scale
  if (identical(covariance, "diag")) {
    c(log(sqrt(intercept_var)), log(sqrt(slope_var)), log(sqrt(residual_var)))
  } else {
    c(log(sqrt(intercept_var)), 0, log(sqrt(slope_var)), log(sqrt(residual_var)))
  }
}

.lmm_knownvar_bounds <- function(kind, covariance) {
  if (identical(kind, "ri")) {
    return(list(lower = c(-14, -14), upper = c(10, 10)))
  }
  if (identical(covariance, "diag")) {
    return(list(lower = c(-14, -14, -14), upper = c(10, 10, 10)))
  }
  list(lower = c(-14, -5, -14, -14), upper = c(10, 5, 10, 10))
}

.lmm_knownvar_objective <- function(par, y, sampling_var, X, W_block, fit, kind, covariance) {
  if (!all(is.finite(par)) || any(abs(par) > 25)) return(Inf)
  components <- .lmm_knownvar_components(par, kind = kind, covariance = covariance)
  U_block <- unname(W_block %*% components$factor)
  if (!all(is.finite(U_block)) || !is.finite(components$residual_var)) return(Inf)
  lmm_knownvar_objective_cpp(
    y,
    sampling_var,
    X,
    U_block = U_block,
    residual_var = components$residual_var,
    fit = fit
  )
}

.optimize_lmm_knownvar <- function(y, sampling_var, X, W_block, fit, kind, covariance) {
  base <- .lmm_knownvar_start(y, sampling_var, W_block, kind, covariance)
  residual_index <- length(base)
  starts <- list(base)

  random_small <- base
  random_indices <- if (identical(kind, "ri")) {
    1L
  } else if (identical(covariance, "diag")) {
    c(1L, 2L)
  } else {
    c(1L, 3L)
  }
  random_small[random_indices] <- -8
  starts[[length(starts) + 1L]] <- random_small

  residual_small <- base
  residual_small[[residual_index]] <- -8
  starts[[length(starts) + 1L]] <- residual_small

  if (identical(kind, "ri_slope1")) {
    slope_index <- if (identical(covariance, "diag")) 2L else 3L
    slope_small <- base
    slope_small[[slope_index]] <- -8
    starts[[length(starts) + 1L]] <- slope_small
    if (identical(covariance, "full")) {
      corr_positive <- base
      corr_positive[[2L]] <- atanh(0.5)
      corr_negative <- base
      corr_negative[[2L]] <- atanh(-0.5)
      starts <- c(starts, list(corr_positive, corr_negative))
    }
  }

  bounds <- .lmm_knownvar_bounds(kind, covariance)
  fn <- function(par) .lmm_knownvar_objective(
    par,
    y = y,
    sampling_var = sampling_var,
    X = X,
    W_block = W_block,
    fit = fit,
    kind = kind,
    covariance = covariance
  )

  best <- NULL
  best_value <- Inf
  for (start in starts) {
    start <- pmin(pmax(start, bounds$lower), bounds$upper)
    opt <- tryCatch(
      stats::optim(
        start,
        fn,
        method = "L-BFGS-B",
        lower = bounds$lower,
        upper = bounds$upper,
        control = list(maxit = 300, factr = 1e7)
      ),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$value)) next
    if (opt$value < best_value) {
      best <- opt
      best_value <- opt$value
    }
  }

  if (is.null(best)) {
    stop("Failed to optimize known-variance LMM covariance", call. = FALSE)
  }
  best
}

.lmm_fit_knownvar_sample <- function(y, sampling_var, design, fit, kind, covariance) {
  W_block <- if (identical(kind, "ri")) {
    matrix(1, nrow = design$n_repeat, ncol = 1L)
  } else {
    design$W_block
  }
  opt <- .optimize_lmm_knownvar(
    y,
    sampling_var,
    X = design$X,
    W_block = W_block,
    fit = fit,
    kind = kind,
    covariance = covariance
  )
  components <- .lmm_knownvar_components(opt$par, kind = kind, covariance = covariance)
  fit_out <- lmm_knownvar_fit_cpp(
    y,
    sampling_var,
    design$X,
    U_block = unname(W_block %*% components$factor),
    residual_var = components$residual_var,
    fit = fit
  )
  fit_out$t_coef <- fit_out$coef / fit_out$se_coef
  fit_out$p_coef <- 2 * stats::pt(-abs(fit_out$t_coef), df = fit_out$df_res)
  fit_out$components <- components
  fit_out$converged <- as.numeric(identical(opt$convergence, 0L))
  fit_out
}

.fit_lmm_knownvar_reducer <- function(arrays, design, opts, kind = c("ri", "ri_slope1")) {
  kind <- match.arg(kind)
  beta <- arrays$beta
  sampling_var <- arrays$var
  method <- if (identical(kind, "ri")) "lmm:ri_knownvar" else "lmm:ri_slope1_knownvar"
  if (is.null(beta) || is.null(sampling_var)) {
    stop(method, " requires both `beta` and `var` assays", call. = FALSE)
  }
  if (!identical(dim(beta), dim(sampling_var)) || length(dim(beta)) != 3L) {
    stop(method, " expects aligned beta and var arrays [sample x subject x contrast]", call. = FALSE)
  }
  if (any(sampling_var <= 0, na.rm = TRUE)) {
    stop(method, " requires strictly positive sampling variances", call. = FALSE)
  }

  Y <- .stack_joint_beta_cube(beta)
  S <- .stack_joint_beta_cube(sampling_var)
  valid <- colSums(!is.finite(Y) | !is.finite(S)) == 0L
  n_samples <- ncol(Y)
  coef_names <- design$X_colnames %||% paste0("X", seq_len(ncol(design$X)))
  p <- length(coef_names)
  fit_mode <- match.arg(opts$fit %||% "REML", c("REML", "ML"))
  theta_mode <- match.arg(opts$theta_mode %||% "voxelwise", "voxelwise")
  covariance_mode <- if (identical(kind, "ri")) {
    "diag"
  } else {
    match.arg(opts$covariance %||% "diag", c("diag", "full"))
  }

  fit_out <- list(
    coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    se_coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    t_coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    p_coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    sigma2 = rep(NA_real_, n_samples),
    vc_intercept = rep(NA_real_, n_samples),
    vc_resid = rep(NA_real_, n_samples),
    sampling_var_mean = rep(NA_real_, n_samples),
    df_res = rep(NA_real_, n_samples),
    logLik = rep(NA_real_, n_samples),
    converged = rep(0, n_samples)
  )
  if (identical(kind, "ri_slope1")) {
    fit_out$vc_slope <- rep(NA_real_, n_samples)
    fit_out$vc_cov_intercept_slope <- rep(NA_real_, n_samples)
    fit_out$corr_intercept_slope <- rep(NA_real_, n_samples)
  }

  for (j in which(valid)) {
    sample_fit <- .lmm_fit_knownvar_sample(
      Y[, j],
      S[, j],
      design = design,
      fit = fit_mode,
      kind = kind,
      covariance = covariance_mode
    )
    fit_out$coef[, j] <- sample_fit$coef
    fit_out$se_coef[, j] <- sample_fit$se_coef
    fit_out$t_coef[, j] <- sample_fit$t_coef
    fit_out$p_coef[, j] <- sample_fit$p_coef
    fit_out$sigma2[[j]] <- sample_fit$components$residual_var
    fit_out$vc_resid[[j]] <- sample_fit$components$residual_var
    fit_out$vc_intercept[[j]] <- sample_fit$components$vc_intercept
    fit_out$sampling_var_mean[[j]] <- mean(S[, j])
    fit_out$df_res[[j]] <- sample_fit$df_res
    fit_out$logLik[[j]] <- sample_fit$logLik
    fit_out$converged[[j]] <- sample_fit$converged
    if (identical(kind, "ri_slope1")) {
      fit_out$vc_slope[[j]] <- sample_fit$components$vc_slope
      fit_out$vc_cov_intercept_slope[[j]] <- sample_fit$components$vc_cov_intercept_slope
      fit_out$corr_intercept_slope[[j]] <- sample_fit$components$corr_intercept_slope
    }
  }

  design_info <- list(
    method = method,
    formula = design$formula,
    columns = coef_names,
    slope = if (identical(kind, "ri_slope1")) design$slope else NULL,
    covariance = covariance_mode,
    theta_mode = theta_mode,
    sampling_variance = "known_diagonal",
    hash = digest::digest(list(X = design$X, W_block = design$W_block)),
    portable = .portable_design_receipt(design$X, coef_names)
  )

  list(
    arrays = .lmm_build_result_arrays(fit_out, coef_names, n_samples),
    subjects = "meta",
    contrasts = "model",
    contrast_data = data.frame(label = "model", row.names = "model", stringsAsFactors = FALSE),
    design_info = design_info,
    attachments = setNames(list(list(
      type = "lmm_design",
      subjects = design$subjects,
      contrasts = design$contrasts,
      slope = design_info$slope,
      fit = fit_mode,
      theta_mode = theta_mode,
      covariance = covariance_mode,
      sampling_variance = "known_diagonal"
    )), paste0("reduce/", method, "/model"))
  )
}

.fit_lmm_ri_knownvar_reducer <- function(arrays, design, opts) {
  .fit_lmm_knownvar_reducer(arrays, design, opts, kind = "ri")
}

.fit_lmm_ri_slope1_knownvar_reducer <- function(arrays, design, opts) {
  .fit_lmm_knownvar_reducer(arrays, design, opts, kind = "ri_slope1")
}
