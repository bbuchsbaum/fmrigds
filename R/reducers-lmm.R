# Restricted repeated-measures Gaussian LMM reducers -------------------------

.lmm_default_contrast_data <- function(contrasts) {
  data.frame(
    contrast = contrasts,
    row.names = contrasts,
    stringsAsFactors = FALSE
  )
}

.lmm_parse_formula <- function(formula) {
  f <- tryCatch(if (inherits(formula, "formula")) formula else stats::as.formula(formula), error = function(e) NULL)
  if (is.null(f)) stop("`formula` must be a valid formula", call. = FALSE)
  if (length(f) != 2L) {
    stop("LMM reducers require a one-sided fixed-effects formula, e.g. `~ group * condition`", call. = FALSE)
  }
  ftxt <- paste(deparse(f), collapse = " ")
  if (grepl("|", ftxt, fixed = TRUE)) {
    stop("LMM reducers do not support lmer-style random-effects syntax; use `method` and reducer options instead", call. = FALSE)
  }
  f
}

.build_joint_reducer_design <- function(reducer, formula, subjects, contrasts, col_data = NULL, contrast_data = NULL, opts = list()) {
  if (!grepl("^lmm:", reducer$name)) {
    stop("No joint design builder registered for reducer ", reducer$name, call. = FALSE)
  }
  f <- .lmm_parse_formula(formula %||% "~ 1")
  subjects <- as.character(subjects %||% character())
  if (!length(subjects)) stop("Joint reducers require subject identifiers", call. = FALSE)

  contrasts <- as.character(contrasts %||% character())
  if (!length(contrasts)) stop("Joint reducers require contrast identifiers", call. = FALSE)

  subj_df <- if (is.null(col_data)) {
    data.frame(row.names = subjects)
  } else {
    .align_col_data_for_subjects(col_data, subjects, warn_extra = FALSE, context = reducer$name)
  }
  if ("subject" %in% names(subj_df)) subj_df$subject <- NULL

  rep_df <- if (is.null(contrast_data)) .lmm_default_contrast_data(contrasts) else {
    .align_contrast_data_for_contrasts(contrast_data, contrasts, warn_extra = FALSE, context = reducer$name)
  }
  if ("contrast" %in% names(rep_df)) rep_df$contrast <- NULL
  rep_df$contrast <- contrasts

  slope_name <- opts$slope %||% NULL
  slope_reducer <- reducer$name %in% c("lmm:ri_slope1", "lmm:ri_slope1_knownvar")
  if (slope_reducer) {
    if (is.null(slope_name) || !nzchar(slope_name)) {
      stop(
        "`options$slope` is required for method = '", reducer$name, "'",
        call. = FALSE
      )
    }
    if (!slope_name %in% names(rep_df)) {
      stop("Slope variable '", slope_name, "' is not present in contrast_data", call. = FALSE)
    }
    if (!is.numeric(rep_df[[slope_name]])) {
      stop("Slope variable '", slope_name, "' must be numeric", call. = FALSE)
    }
    if (isTRUE(opts$center_slope)) {
      rep_df[[slope_name]] <- rep_df[[slope_name]] - mean(rep_df[[slope_name]])
    }
    if (length(unique(rep_df[[slope_name]])) < 2L) {
      stop("Slope variable '", slope_name, "' must vary across repeated-measure levels", call. = FALSE)
    }
  }

  obs_blocks <- lapply(seq_along(subjects), function(i) {
    block_subj <- subj_df[rep(i, nrow(rep_df)), , drop = FALSE]
    data.frame(
      subject = rep(subjects[[i]], nrow(rep_df)),
      block_subj,
      rep_df,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  obs_data <- do.call(rbind, obs_blocks)
  rownames(obs_data) <- NULL
  obs_data$subject <- factor(obs_data$subject, levels = subjects)
  obs_data$contrast <- factor(obs_data$contrast, levels = contrasts)

  X <- stats::model.matrix(f, data = obs_data)
  if (!is.matrix(X) || !nrow(X) || !ncol(X)) {
    stop("Failed to build fixed-effects design matrix for LMM reducer", call. = FALSE)
  }
  if (!all(is.finite(X))) {
    stop("LMM design matrix contains non-finite values", call. = FALSE)
  }

  list(
    formula = paste(deparse(f), collapse = " "),
    subjects = subjects,
    contrasts = contrasts,
    contrast_data = rep_df[contrasts, , drop = FALSE],
    observation_data = obs_data,
    X = unname(X),
    X_colnames = colnames(X),
    W_block = if (slope_reducer) {
      unname(cbind(Intercept = 1, slope = rep_df[[slope_name]]))
    } else {
      NULL
    },
    n_repeat = length(contrasts),
    slope = slope_name
  )
}

.stack_joint_beta_cube <- function(beta) {
  dims <- dim(beta)
  stopifnot(length(dims) == 3L)
  matrix(
    aperm(beta, c(3L, 2L, 1L)),
    nrow = dims[2L] * dims[3L],
    ncol = dims[1L]
  )
}

.optimize_lmm_ri_theta <- function(Y, X, n_repeat, fit = c("REML", "ML"), eta_hint = NULL, eta_window = 4) {
  fit <- match.arg(fit)
  obj0 <- lmm_ri_objective_cpp(Y, X, n_repeat = n_repeat, lambda = 0, fit = fit)
  fn <- function(eta) lmm_ri_objective_cpp(Y, X, n_repeat = n_repeat, lambda = exp(eta), fit = fit)
  intervals <- list(c(-14, 10))
  if (!is.null(eta_hint) && length(eta_hint) && is.finite(eta_hint)) {
    lower <- max(-14, eta_hint - eta_window)
    upper <- min(10, eta_hint + eta_window)
    intervals <- c(list(c(lower, upper)), intervals)
  }

  best <- NULL
  best_obj <- Inf
  for (interval in intervals) {
    opt <- stats::optimize(fn, interval = interval)
    if (is.finite(opt$objective) && opt$objective < best_obj) {
      best <- opt
      best_obj <- opt$objective
    }
  }

  lambda <- exp(best$minimum)
  if (is.finite(obj0) && obj0 <= best_obj) {
    lambda <- 0
  }
  lambda
}

.lmm_add_inference <- function(fit) {
  fit$t_coef <- fit$coef / fit$se_coef
  dfmat <- matrix(fit$df_res, nrow = nrow(fit$t_coef), ncol = ncol(fit$t_coef), byrow = TRUE)
  fit$p_coef <- 2 * stats::pt(-abs(fit$t_coef), df = dfmat)
  fit
}

.lmm_build_result_arrays <- function(fit, coef_names, n_samples) {
  arr3 <- function(x) array(as.numeric(x), dim = c(n_samples, 1L, 1L))
  out <- list()
  scalar_names <- setdiff(names(fit), c("coef", "se_coef", "t_coef", "p_coef"))
  for (nm in scalar_names) {
    out[[nm]] <- arr3(fit[[nm]])
  }

  for (j in seq_along(coef_names)) {
    nm <- coef_names[[j]]
    out[[paste0("coef:", nm)]] <- arr3(fit$coef[j, ])
    out[[paste0("se_coef:", nm)]] <- arr3(fit$se_coef[j, ])
    out[[paste0("t_coef:", nm)]] <- arr3(fit$t_coef[j, ])
    out[[paste0("p_coef:", nm)]] <- arr3(fit$p_coef[j, ])
  }
  out
}

.lmm_ris1_factor <- function(par, covariance = c("diag", "full")) {
  covariance <- match.arg(covariance)
  if (identical(covariance, "diag")) {
    if (length(par) != 2L) stop("diag covariance requires two parameters", call. = FALSE)
    return(diag(exp(par), nrow = 2L, ncol = 2L))
  }
  if (length(par) != 3L) stop("full covariance requires three parameters", call. = FALSE)
  matrix(
    c(exp(par[1L]), 0, par[2L], exp(par[3L])),
    nrow = 2L,
    byrow = TRUE
  )
}

.lmm_ris1_relative_cov <- function(factor_mat) {
  tcrossprod(factor_mat)
}

.lmm_ris1_relative_components <- function(factor_mat) {
  rel_cov <- .lmm_ris1_relative_cov(factor_mat)
  rel_intercept <- rel_cov[1L, 1L]
  rel_slope <- rel_cov[2L, 2L]
  rel_cov_is <- rel_cov[1L, 2L]
  corr <- if (rel_intercept <= 0 || rel_slope <= 0) {
    0
  } else {
    rel_cov_is / sqrt(rel_intercept * rel_slope)
  }
  list(
    rel_cov = rel_cov,
    lambda_intercept = rel_intercept,
    lambda_slope = rel_slope,
    lambda_cov_intercept_slope = rel_cov_is,
    corr_intercept_slope = corr
  )
}

.lmm_ris1_u_block <- function(W_block, factor_mat) {
  unname(W_block %*% factor_mat)
}

.lmm_ris1_objective <- function(par, Y, X, W_block, fit = c("REML", "ML"), covariance = c("diag", "full")) {
  fit <- match.arg(fit)
  covariance <- match.arg(covariance)
  if (!all(is.finite(par)) || any(abs(par) > 25)) {
    return(Inf)
  }
  factor_mat <- .lmm_ris1_factor(par, covariance = covariance)
  if (!all(is.finite(factor_mat))) {
    return(Inf)
  }
  U_block <- .lmm_ris1_u_block(W_block, factor_mat)
  if (!all(is.finite(U_block))) {
    return(Inf)
  }
  lmm_lowrank_objective_cpp(Y, X, U_block = U_block, fit = fit)
}

.lmm_ris1_bounds <- function(covariance = c("diag", "full")) {
  covariance <- match.arg(covariance)
  if (identical(covariance, "diag")) {
    return(list(lower = c(-14, -14), upper = c(10, 10)))
  }
  list(lower = c(-14, -10, -14), upper = c(10, 10, 10))
}

.optimize_lmm_ris1_theta <- function(Y, X, W_block, fit = c("REML", "ML"), covariance = c("diag", "full"), lambda_ri = NULL, start_par = NULL) {
  fit <- match.arg(fit)
  covariance <- match.arg(covariance)
  lambda_ri <- as.numeric(lambda_ri %||% 0)
  lambda_ri <- if (length(lambda_ri)) max(lambda_ri, 1e-8) else 1e-8
  log_sd_intercept <- 0.5 * log(lambda_ri)

  starts <- if (identical(covariance, "diag")) {
    list(
      c(log_sd_intercept, -6),
      c(log_sd_intercept, -2),
      c(-4, -4)
    )
  } else {
    list(
      c(log_sd_intercept, 0, -6),
      c(log_sd_intercept, 0, -2),
      c(-4, 0, -4)
    )
  }
  if (!is.null(start_par)) {
    starts <- c(list(as.numeric(start_par)), starts)
  }
  starts <- Filter(function(x) all(is.finite(x)), starts)
  bounds <- .lmm_ris1_bounds(covariance)

  fn <- function(par) .lmm_ris1_objective(par, Y = Y, X = X, W_block = W_block, fit = fit, covariance = covariance)
  best <- NULL
  best_val <- Inf

  for (start in starts) {
    start <- pmin(pmax(start, bounds$lower), bounds$upper)
    opt <- tryCatch(
      stats::optim(
        start,
        fn,
        method = "L-BFGS-B",
        lower = bounds$lower,
        upper = bounds$upper,
        control = list(maxit = 200, factr = 1e7)
      ),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$value) || !identical(opt$convergence, 0L)) {
      opt <- tryCatch(
        stats::optim(start, fn, method = "Nelder-Mead", control = list(maxit = 1000, reltol = 1e-9)),
        error = function(e) NULL
      )
    }
    if (!is.null(opt) && is.finite(opt$value) && opt$value < best_val) {
      best <- opt
      best_val <- opt$value
    }
  }

  if (is.null(best)) {
    stop("Failed to optimize pooled theta for lmm:ri_slope1", call. = FALSE)
  }
  best$par
}

.lmm_fit_ri_pooled <- function(Y_valid, design, fit_mode) {
  lambda <- .optimize_lmm_ri_theta(Y_valid, design$X, n_repeat = design$n_repeat, fit = fit_mode)
  fit_valid <- lmm_ri_fit_cpp(Y_valid, design$X, n_repeat = design$n_repeat, lambda = lambda, fit = fit_mode)
  fit_valid <- .lmm_add_inference(fit_valid)
  fit_valid
}

.lmm_fit_ri_voxelwise <- function(Y_valid, design, fit_mode) {
  n_valid <- ncol(Y_valid)
  p <- ncol(design$X)
  fit_valid <- list(
    coef = matrix(NA_real_, nrow = p, ncol = n_valid),
    se_coef = matrix(NA_real_, nrow = p, ncol = n_valid),
    sigma2 = rep(NA_real_, n_valid),
    vc_intercept = rep(NA_real_, n_valid),
    vc_resid = rep(NA_real_, n_valid),
    df_res = rep(NA_real_, n_valid),
    logLik = rep(NA_real_, n_valid),
    converged = rep(0, n_valid),
    lambda = rep(NA_real_, n_valid),
    t_coef = matrix(NA_real_, nrow = p, ncol = n_valid),
    p_coef = matrix(NA_real_, nrow = p, ncol = n_valid)
  )

  lambda_start <- NULL
  if (n_valid > 1L) {
    lambda_start <- tryCatch(
      .optimize_lmm_ri_theta(Y_valid, design$X, n_repeat = design$n_repeat, fit = fit_mode),
      error = function(e) NULL
    )
  }
  eta_start <- if (is.null(lambda_start) || !is.finite(lambda_start) || lambda_start <= 0) NULL else log(lambda_start)

  for (j in seq_len(n_valid)) {
    yj <- Y_valid[, j, drop = FALSE]
    lambda_j <- .optimize_lmm_ri_theta(
      yj,
      design$X,
      n_repeat = design$n_repeat,
      fit = fit_mode,
      eta_hint = eta_start
    )
    fit_j <- lmm_ri_fit_cpp(yj, design$X, n_repeat = design$n_repeat, lambda = lambda_j, fit = fit_mode)
    fit_j <- .lmm_add_inference(fit_j)

    fit_valid$coef[, j] <- fit_j$coef[, 1]
    fit_valid$se_coef[, j] <- fit_j$se_coef[, 1]
    fit_valid$sigma2[j] <- fit_j$sigma2[[1L]]
    fit_valid$vc_intercept[j] <- fit_j$vc_intercept[[1L]]
    fit_valid$vc_resid[j] <- fit_j$vc_resid[[1L]]
    fit_valid$df_res[j] <- fit_j$df_res[[1L]]
    fit_valid$logLik[j] <- fit_j$logLik[[1L]]
    fit_valid$converged[j] <- fit_j$converged[[1L]]
    fit_valid$lambda[j] <- fit_j$lambda[[1L]]
    fit_valid$t_coef[, j] <- fit_j$t_coef[, 1]
    fit_valid$p_coef[, j] <- fit_j$p_coef[, 1]
    eta_start <- if (is.finite(lambda_j) && lambda_j > 0) log(lambda_j) else NULL
  }

  fit_valid
}

.lmm_fit_ris1_pooled <- function(Y_valid, design, fit_mode, covariance_mode) {
  lambda_ri <- .optimize_lmm_ri_theta(Y_valid, design$X, n_repeat = design$n_repeat, fit = fit_mode)
  theta_par <- .optimize_lmm_ris1_theta(
    Y_valid,
    design$X,
    W_block = design$W_block,
    fit = fit_mode,
    covariance = covariance_mode,
    lambda_ri = lambda_ri
  )
  theta_factor <- .lmm_ris1_factor(theta_par, covariance = covariance_mode)
  theta_components <- .lmm_ris1_relative_components(theta_factor)
  fit_valid <- lmm_lowrank_fit_cpp(
    Y_valid,
    design$X,
    U_block = .lmm_ris1_u_block(design$W_block, theta_factor),
    fit = fit_mode
  )
  fit_valid <- .lmm_add_inference(fit_valid)
  list(fit = fit_valid, theta = theta_components, par = theta_par)
}

.lmm_fit_ris1_voxelwise <- function(Y_valid, design, fit_mode, covariance_mode) {
  n_valid <- ncol(Y_valid)
  p <- ncol(design$X)
  fit_valid <- list(
    coef = matrix(NA_real_, nrow = p, ncol = n_valid),
    se_coef = matrix(NA_real_, nrow = p, ncol = n_valid),
    sigma2 = rep(NA_real_, n_valid),
    vc_intercept = rep(NA_real_, n_valid),
    vc_slope = rep(NA_real_, n_valid),
    vc_cov_intercept_slope = rep(NA_real_, n_valid),
    vc_resid = rep(NA_real_, n_valid),
    df_res = rep(NA_real_, n_valid),
    logLik = rep(NA_real_, n_valid),
    converged = rep(0, n_valid),
    lambda_intercept = rep(NA_real_, n_valid),
    lambda_slope = rep(NA_real_, n_valid),
    lambda_cov_intercept_slope = rep(NA_real_, n_valid),
    corr_intercept_slope = rep(NA_real_, n_valid),
    t_coef = matrix(NA_real_, nrow = p, ncol = n_valid),
    p_coef = matrix(NA_real_, nrow = p, ncol = n_valid)
  )

  pooled_start <- NULL
  lambda_start <- NULL
  if (n_valid > 1L) {
    pooled_fit <- tryCatch(
      .lmm_fit_ris1_pooled(Y_valid, design, fit_mode = fit_mode, covariance_mode = covariance_mode),
      error = function(e) NULL
    )
    if (!is.null(pooled_fit)) {
      pooled_start <- pooled_fit$par
      lambda_start <- pooled_fit$theta$lambda_intercept
    }
  }
  eta_start <- if (is.null(lambda_start) || !is.finite(lambda_start) || lambda_start <= 0) NULL else 0.5 * log(lambda_start)
  theta_start <- pooled_start

  for (j in seq_len(n_valid)) {
    yj <- Y_valid[, j, drop = FALSE]
    lambda_ri_j <- .optimize_lmm_ri_theta(
      yj,
      design$X,
      n_repeat = design$n_repeat,
      fit = fit_mode,
      eta_hint = eta_start
    )
    theta_par_j <- .optimize_lmm_ris1_theta(
      yj,
      design$X,
      W_block = design$W_block,
      fit = fit_mode,
      covariance = covariance_mode,
      lambda_ri = lambda_ri_j,
      start_par = theta_start
    )
    theta_factor_j <- .lmm_ris1_factor(theta_par_j, covariance = covariance_mode)
    theta_components_j <- .lmm_ris1_relative_components(theta_factor_j)
    fit_j <- lmm_lowrank_fit_cpp(
      yj,
      design$X,
      U_block = .lmm_ris1_u_block(design$W_block, theta_factor_j),
      fit = fit_mode
    )
    fit_j <- .lmm_add_inference(fit_j)

    sigma2_j <- fit_j$sigma2[[1L]]
    fit_valid$coef[, j] <- fit_j$coef[, 1]
    fit_valid$se_coef[, j] <- fit_j$se_coef[, 1]
    fit_valid$sigma2[j] <- sigma2_j
    fit_valid$vc_intercept[j] <- theta_components_j$lambda_intercept * sigma2_j
    fit_valid$vc_slope[j] <- theta_components_j$lambda_slope * sigma2_j
    fit_valid$vc_cov_intercept_slope[j] <- theta_components_j$lambda_cov_intercept_slope * sigma2_j
    fit_valid$vc_resid[j] <- sigma2_j
    fit_valid$df_res[j] <- fit_j$df_res[[1L]]
    fit_valid$logLik[j] <- fit_j$logLik[[1L]]
    fit_valid$converged[j] <- fit_j$converged[[1L]]
    fit_valid$lambda_intercept[j] <- theta_components_j$lambda_intercept
    fit_valid$lambda_slope[j] <- theta_components_j$lambda_slope
    fit_valid$lambda_cov_intercept_slope[j] <- theta_components_j$lambda_cov_intercept_slope
    fit_valid$corr_intercept_slope[j] <- theta_components_j$corr_intercept_slope
    fit_valid$t_coef[, j] <- fit_j$t_coef[, 1]
    fit_valid$p_coef[, j] <- fit_j$p_coef[, 1]
    eta_start <- if (is.finite(lambda_ri_j) && lambda_ri_j > 0) log(lambda_ri_j) else NULL
    theta_start <- theta_par_j
  }

  fit_valid
}

.fit_lmm_ri_reducer <- function(arrays, design, opts) {
  beta <- arrays$beta
  if (is.null(beta)) stop("lmm:ri requires a `beta` assay", call. = FALSE)
  dims <- dim(beta)
  if (length(dims) != 3L) stop("lmm:ri expects beta with dimensions [sample x subject x contrast]", call. = FALSE)
  if (dims[2L] != length(design$subjects) || dims[3L] != design$n_repeat) {
    stop("LMM design metadata is not aligned with beta dimensions", call. = FALSE)
  }

  Y <- .stack_joint_beta_cube(beta)
  valid <- colSums(!is.finite(Y)) == 0L
  n_samples <- ncol(Y)
  coef_names <- design$X_colnames %||% paste0("X", seq_len(ncol(design$X)))
  p <- length(coef_names)

  fit_mode <- match.arg(opts$fit %||% "REML", c("REML", "ML"))
  theta_mode <- match.arg(opts$theta_mode %||% "pooled", c("pooled", "voxelwise"))

  fit_out <- list(
    coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    se_coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    sigma2 = rep(NA_real_, n_samples),
    vc_intercept = rep(NA_real_, n_samples),
    vc_resid = rep(NA_real_, n_samples),
    df_res = rep(NA_real_, n_samples),
    logLik = rep(NA_real_, n_samples),
    converged = rep(0, n_samples),
    lambda = rep(NA_real_, n_samples)
  )

  if (any(valid)) {
    Y_valid <- Y[, valid, drop = FALSE]
    fit_valid <- if (identical(theta_mode, "pooled")) {
      .lmm_fit_ri_pooled(Y_valid, design, fit_mode = fit_mode)
    } else {
      .lmm_fit_ri_voxelwise(Y_valid, design, fit_mode = fit_mode)
    }

    fit_out$coef[, valid] <- fit_valid$coef
    fit_out$se_coef[, valid] <- fit_valid$se_coef
    fit_out$sigma2[valid] <- fit_valid$sigma2
    fit_out$vc_intercept[valid] <- fit_valid$vc_intercept
    fit_out$vc_resid[valid] <- fit_valid$vc_resid
    fit_out$df_res[valid] <- fit_valid$df_res
    fit_out$logLik[valid] <- fit_valid$logLik
    fit_out$converged[valid] <- fit_valid$converged
    fit_out$lambda[valid] <- fit_valid$lambda
    fit_out$t_coef <- matrix(NA_real_, nrow = p, ncol = n_samples)
    fit_out$p_coef <- matrix(NA_real_, nrow = p, ncol = n_samples)
    fit_out$t_coef[, valid] <- fit_valid$t_coef
    fit_out$p_coef[, valid] <- fit_valid$p_coef
  } else {
    fit_out$t_coef <- matrix(NA_real_, nrow = p, ncol = n_samples)
    fit_out$p_coef <- matrix(NA_real_, nrow = p, ncol = n_samples)
  }

  design_info <- list(
    method = "lmm:ri",
    formula = design$formula,
    columns = coef_names,
    theta_mode = theta_mode,
    hash = digest::digest(design$X),
    portable = .portable_design_receipt(design$X, coef_names)
  )

  list(
    arrays = .lmm_build_result_arrays(fit_out, coef_names, n_samples),
    subjects = "meta",
    contrasts = "model",
    contrast_data = data.frame(label = "model", row.names = "model", stringsAsFactors = FALSE),
    design_info = design_info,
    attachments = list(
      "reduce/lmm:ri/model" = list(
        type = "lmm_design",
        subjects = design$subjects,
        contrasts = design$contrasts,
        fit = fit_mode,
        theta_mode = theta_mode
      )
    )
  )
}

.fit_lmm_ri_slope1_reducer <- function(arrays, design, opts) {
  beta <- arrays$beta
  if (is.null(beta)) stop("lmm:ri_slope1 requires a `beta` assay", call. = FALSE)
  dims <- dim(beta)
  if (length(dims) != 3L) stop("lmm:ri_slope1 expects beta with dimensions [sample x subject x contrast]", call. = FALSE)
  if (dims[2L] != length(design$subjects) || dims[3L] != design$n_repeat) {
    stop("LMM design metadata is not aligned with beta dimensions", call. = FALSE)
  }

  Y <- .stack_joint_beta_cube(beta)
  valid <- colSums(!is.finite(Y)) == 0L
  n_samples <- ncol(Y)
  coef_names <- design$X_colnames %||% paste0("X", seq_len(ncol(design$X)))
  p <- length(coef_names)
  fit_mode <- match.arg(opts$fit %||% "REML", c("REML", "ML"))
  theta_mode <- match.arg(opts$theta_mode %||% "pooled", c("pooled", "voxelwise"))
  covariance_mode <- match.arg(opts$covariance %||% "diag", c("diag", "full"))

  fit_out <- list(
    coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    se_coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    sigma2 = rep(NA_real_, n_samples),
    vc_intercept = rep(NA_real_, n_samples),
    vc_slope = rep(NA_real_, n_samples),
    vc_cov_intercept_slope = rep(NA_real_, n_samples),
    vc_resid = rep(NA_real_, n_samples),
    df_res = rep(NA_real_, n_samples),
    logLik = rep(NA_real_, n_samples),
    converged = rep(0, n_samples),
    lambda_intercept = rep(NA_real_, n_samples),
    lambda_slope = rep(NA_real_, n_samples),
    lambda_cov_intercept_slope = rep(NA_real_, n_samples),
    corr_intercept_slope = rep(NA_real_, n_samples),
    t_coef = matrix(NA_real_, nrow = p, ncol = n_samples),
    p_coef = matrix(NA_real_, nrow = p, ncol = n_samples)
  )

  if (any(valid)) {
    Y_valid <- Y[, valid, drop = FALSE]
    fit_valid <- if (identical(theta_mode, "pooled")) {
      pooled <- .lmm_fit_ris1_pooled(Y_valid, design, fit_mode = fit_mode, covariance_mode = covariance_mode)
      rel_intercept <- pooled$theta$lambda_intercept
      rel_slope <- pooled$theta$lambda_slope
      rel_cov_is <- pooled$theta$lambda_cov_intercept_slope
      corr_is <- pooled$theta$corr_intercept_slope
      sigma2_valid <- pooled$fit$sigma2
      pooled$fit$vc_intercept <- rel_intercept * sigma2_valid
      pooled$fit$vc_slope <- rel_slope * sigma2_valid
      pooled$fit$vc_cov_intercept_slope <- rel_cov_is * sigma2_valid
      pooled$fit$vc_resid <- sigma2_valid
      pooled$fit$lambda_intercept <- rep(rel_intercept, length(sigma2_valid))
      pooled$fit$lambda_slope <- rep(rel_slope, length(sigma2_valid))
      pooled$fit$lambda_cov_intercept_slope <- rep(rel_cov_is, length(sigma2_valid))
      pooled$fit$corr_intercept_slope <- rep(corr_is, length(sigma2_valid))
      pooled$fit
    } else {
      .lmm_fit_ris1_voxelwise(Y_valid, design, fit_mode = fit_mode, covariance_mode = covariance_mode)
    }

    fit_out$coef[, valid] <- fit_valid$coef
    fit_out$se_coef[, valid] <- fit_valid$se_coef
    fit_out$t_coef[, valid] <- fit_valid$t_coef
    fit_out$p_coef[, valid] <- fit_valid$p_coef
    fit_out$sigma2[valid] <- fit_valid$sigma2
    fit_out$vc_intercept[valid] <- fit_valid$vc_intercept
    fit_out$vc_slope[valid] <- fit_valid$vc_slope
    fit_out$vc_cov_intercept_slope[valid] <- fit_valid$vc_cov_intercept_slope
    fit_out$vc_resid[valid] <- fit_valid$vc_resid
    fit_out$df_res[valid] <- fit_valid$df_res
    fit_out$logLik[valid] <- fit_valid$logLik
    fit_out$converged[valid] <- fit_valid$converged
    fit_out$lambda_intercept[valid] <- fit_valid$lambda_intercept
    fit_out$lambda_slope[valid] <- fit_valid$lambda_slope
    fit_out$lambda_cov_intercept_slope[valid] <- fit_valid$lambda_cov_intercept_slope
    fit_out$corr_intercept_slope[valid] <- fit_valid$corr_intercept_slope
  }

  design_info <- list(
    method = "lmm:ri_slope1",
    formula = design$formula,
    columns = coef_names,
    slope = design$slope,
    covariance = covariance_mode,
    theta_mode = theta_mode,
    hash = digest::digest(list(X = design$X, W_block = design$W_block)),
    portable = .portable_design_receipt(design$X, coef_names)
  )

  list(
    arrays = .lmm_build_result_arrays(fit_out, coef_names, n_samples),
    subjects = "meta",
    contrasts = "model",
    contrast_data = data.frame(label = "model", row.names = "model", stringsAsFactors = FALSE),
    design_info = design_info,
    attachments = list(
      "reduce/lmm:ri_slope1/model" = list(
        type = "lmm_design",
        subjects = design$subjects,
        contrasts = design$contrasts,
        slope = design$slope,
        fit = fit_mode,
        theta_mode = theta_mode,
        covariance = covariance_mode
      )
    )
  )
}

register_lmm_reducers <- function() {
  register_reducer(
    name = "lmm:ri",
    fun = .fit_lmm_ri_reducer,
    requires = c("beta"),
    provides = c("coef", "se_coef", "t_coef", "p_coef", "sigma2", "vc_intercept", "vc_resid", "df_res", "logLik", "converged", "lambda"),
    options_schema = list(fit = c("REML", "ML"), theta_mode = c("pooled", "voxelwise")),
    input_shape = "joint_contrast"
  )
  register_reducer(
    name = "lmm:ri_slope1",
    fun = .fit_lmm_ri_slope1_reducer,
    requires = c("beta"),
    provides = c(
      "coef", "se_coef", "t_coef", "p_coef", "sigma2",
      "vc_intercept", "vc_slope", "vc_cov_intercept_slope", "vc_resid",
      "df_res", "logLik", "converged",
      "lambda_intercept", "lambda_slope", "lambda_cov_intercept_slope",
      "corr_intercept_slope"
    ),
    options_schema = list(
      fit = c("REML", "ML"),
      theta_mode = c("pooled", "voxelwise"),
      covariance = c("diag", "full")
    ),
    input_shape = "joint_contrast"
  )
  register_reducer(
    name = "lmm:ri_knownvar",
    fun = .fit_lmm_ri_knownvar_reducer,
    requires = c("beta", "var"),
    provides = c(
      "coef", "se_coef", "t_coef", "p_coef", "sigma2",
      "vc_intercept", "vc_resid", "sampling_var_mean",
      "df_res", "logLik", "converged"
    ),
    options_schema = list(fit = c("REML", "ML"), theta_mode = "voxelwise"),
    input_shape = "joint_contrast",
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "model_specific",
      missingness = "complete_case",
      synthetic_variance = "forbid",
      deletion = "unsupported"
    )
  )
  register_reducer(
    name = "lmm:ri_slope1_knownvar",
    fun = .fit_lmm_ri_slope1_knownvar_reducer,
    requires = c("beta", "var"),
    provides = c(
      "coef", "se_coef", "t_coef", "p_coef", "sigma2",
      "vc_intercept", "vc_slope", "vc_cov_intercept_slope", "vc_resid",
      "corr_intercept_slope", "sampling_var_mean",
      "df_res", "logLik", "converged"
    ),
    options_schema = list(
      fit = c("REML", "ML"),
      theta_mode = "voxelwise",
      covariance = c("diag", "full")
    ),
    input_shape = "joint_contrast",
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "model_specific",
      missingness = "complete_case",
      synthetic_variance = "forbid",
      deletion = "unsupported"
    )
  )
}
