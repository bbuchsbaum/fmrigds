.frame_fit_selection <- function(frame) {
  if (inherits(frame, "fmri_view")) {
    list(
      base = frame$base,
      observations = frame$observation_index,
      features = frame$feature_index
    )
  } else {
    list(
      base = frame,
      observations = seq_len(nrow(frame)),
      features = seq_len(ncol(frame))
    )
  }
}

.frame_random_intercept_only <- function(compiled) {
  formula <- compiled$random_formula
  if (is.null(formula)) {
    stop("Frame group fitting requires a random-intercept formula.", call. = FALSE)
  }
  rhs <- formula[[2L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|"))) {
    stop("Frame group fitting requires random syntax `~ 1 | group`.", call. = FALSE)
  }
  if (!identical(rhs[[2L]], 1)) {
    stop("The walking-skeleton fitter supports a random intercept only.", call. = FALSE)
  }
  grouping <- multidesign::grouping_data(compiled)
  if (ncol(grouping) != 1L || anyNA(grouping[[1L]])) {
    stop("A random-intercept fit requires one complete grouping variable.", call. = FALSE)
  }
  factor(grouping[[1L]], levels = unique(grouping[[1L]]))
}

.frame_ri_gls <- function(y, sampling_var, X, groups, variances, fit) {
  tau2 <- variances[[1L]]
  residual_var <- variances[[2L]]
  Z <- stats::model.matrix(~ 0 + groups)
  V <- diag(sampling_var + residual_var, nrow = length(y)) +
    tau2 * tcrossprod(Z)
  chol_v <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(chol_v)) {
    return(NULL)
  }
  inverse_v <- chol2inv(chol_v)
  information <- crossprod(X, inverse_v %*% X)
  chol_information <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(chol_information)) {
    return(NULL)
  }
  covariance <- chol2inv(chol_information)
  coefficient <- covariance %*% crossprod(X, inverse_v %*% y)
  residual <- as.numeric(y - X %*% coefficient)
  objective <- 2 * sum(log(diag(chol_v))) +
    as.numeric(crossprod(residual, inverse_v %*% residual))
  if (identical(fit, "REML")) {
    objective <- objective + 2 * sum(log(diag(chol_information)))
  }
  list(
    objective = objective,
    coefficient = as.numeric(coefficient),
    covariance = covariance
  )
}

.frame_ri_starts <- function(y, sampling_var) {
  observed <- stats::var(y)
  sampling <- stats::median(sampling_var)
  excess <- max(observed - sampling, observed * 0.1, 1e-8)
  base <- log(c(0.65 * excess, 0.35 * excess))
  list(base, c(-12, base[[2L]]), c(base[[1L]], -12), c(-6, -6))
}

.fit_frame_ri_feature <- function(y, sampling_var, X, groups, fit = "REML") {
  valid <- is.finite(y) & is.finite(sampling_var) & sampling_var > 0 &
    stats::complete.cases(X) & !is.na(groups)
  y <- y[valid]
  sampling_var <- sampling_var[valid]
  X <- X[valid, , drop = FALSE]
  groups <- droplevels(groups[valid])
  p <- ncol(X)
  if (length(y) <= p || nlevels(groups) < 2L || qr(X)$rank < p) {
    return(NULL)
  }

  objective <- function(parameters) {
    if (!all(is.finite(parameters))) {
      return(Inf)
    }
    gls <- .frame_ri_gls(
      y,
      sampling_var,
      X,
      groups,
      variances = exp(parameters),
      fit = fit
    )
    if (is.null(gls) || !is.finite(gls$objective)) Inf else gls$objective
  }
  best <- NULL
  for (start in .frame_ri_starts(y, sampling_var)) {
    candidate <- tryCatch(
      stats::optim(
        start,
        objective,
        method = "L-BFGS-B",
        lower = c(-20, -20),
        upper = c(10, 10),
        control = list(maxit = 300, factr = 1e7)
      ),
      error = function(e) NULL
    )
    if (is.null(candidate) || !is.finite(candidate$value)) next
    if (is.null(best) || candidate$value < best$value) best <- candidate
  }
  if (is.null(best)) {
    return(NULL)
  }

  variances <- exp(best$par)
  gls <- .frame_ri_gls(y, sampling_var, X, groups, variances, fit)
  standard_error <- sqrt(diag(gls$covariance))
  statistic <- gls$coefficient / standard_error
  degrees_freedom <- max(length(y) - p, 1L)
  list(
    coefficient = gls$coefficient,
    standard_error = standard_error,
    statistic = statistic,
    p_value = 2 * stats::pt(-abs(statistic), df = degrees_freedom),
    vc_intercept = variances[[1L]],
    vc_resid = variances[[2L]],
    log_likelihood = -0.5 * gls$objective,
    degrees_freedom = degrees_freedom,
    converged = identical(best$convergence, 0L)
  )
}

.frame_fit_block_size <- function(n_observation, n_feature, memory_budget, block_size) {
  if (!is.numeric(memory_budget) || length(memory_budget) != 1L ||
    !is.finite(memory_budget) || memory_budget <= 0) {
    stop("`memory_budget` must be one positive finite byte count.", call. = FALSE)
  }
  maximum <- max(1L, floor(memory_budget / max(16 * n_observation, 1)))
  if (is.null(block_size)) {
    return(as.integer(min(n_feature, maximum)))
  }
  block_size <- as.integer(block_size)
  if (length(block_size) != 1L || is.na(block_size) || block_size < 1L) {
    stop("`block_size` must be one positive integer.", call. = FALSE)
  }
  if (block_size > maximum) {
    stop("`block_size` exceeds the requested memory_budget.", call. = FALSE)
  }
  min(block_size, n_feature)
}

.frame_result_rows <- function(coefficient_names) {
  data.frame(
    .obs_id = c(coefficient_names, "vc_intercept", "vc_resid"),
    result_type = c(
      rep("fixed_effect", length(coefficient_names)),
      "variance_component",
      "variance_component"
    ),
    term = c(coefficient_names, "random_intercept", "residual"),
    stringsAsFactors = FALSE
  )
}

#' Fit a variance-aware random-intercept model over frame feature blocks
#'
#' `fit_group()` is the first frame-native `fmrigds` execution path. It consumes
#' aligned estimate and sampling-variance assays, compiles a `multidesign`
#' specification when necessary, and reads bounded feature blocks. The model is
#' `diag(sampling_variance) + vc_resid I + vc_intercept Z Z'`.
#'
#' @param frame An `fmri_frame` or synchronized view.
#' @param estimate Name of the effect assay.
#' @param variance Name of the known sampling-variance assay.
#' @param design A `design_spec` or `compiled_design`.
#' @param memory_budget Maximum bytes for the two input feature blocks.
#' @param block_size Optional explicit feature-block width.
#' @param fit Likelihood criterion, `"REML"` or `"ML"`.
#' @return An `fmri_group_fit` whose `result` is an `fmri_frame` sharing the
#'   input feature IDs and spatial domain.
#' @export
fit_group <- function(
  frame,
  estimate = "beta",
  variance = "variance",
  design,
  memory_budget = 256 * 1024^2,
  block_size = NULL,
  fit = c("REML", "ML")
) {
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  fit <- match.arg(fit)
  compiled <- if (inherits(design, "design_spec")) {
    multidesign::compile_design(frame, design)
  } else {
    design
  }
  if (!inherits(compiled, "compiled_design")) {
    stop("`design` must be a design_spec or compiled_design.", call. = FALSE)
  }
  groups <- .frame_random_intercept_only(compiled)
  X <- multidesign::model_matrix(compiled)
  if (!identical(
    compiled$observation_ids,
    fmridataset::observation_ids(frame)
  )) {
    stop("Compiled design observation IDs do not align with the frame.", call. = FALSE)
  }
  if (qr(X)$rank < ncol(X)) {
    stop("The fixed-effects design matrix is rank deficient.", call. = FALSE)
  }

  selection <- .frame_fit_selection(frame)
  estimate_assay <- fmridataset::assay(selection$base, estimate)
  variance_assay <- fmridataset::assay(selection$base, variance)
  feature_count <- length(selection$features)
  width <- .frame_fit_block_size(nrow(X), feature_count, memory_budget, block_size)
  coefficient_names <- colnames(X)
  row_data <- .frame_result_rows(coefficient_names)
  result_count <- nrow(row_data)
  estimate_result <- matrix(NA_real_, nrow = result_count, ncol = feature_count)
  standard_error <- matrix(NA_real_, nrow = result_count, ncol = feature_count)
  statistic <- matrix(NA_real_, nrow = result_count, ncol = feature_count)
  p_value <- matrix(NA_real_, nrow = result_count, ncol = feature_count)
  convergence <- logical(feature_count)
  log_likelihood <- rep(NA_real_, feature_count)

  starts <- seq.int(1L, feature_count, by = width)
  for (start in starts) {
    local_features <- seq.int(start, min(start + width - 1L, feature_count))
    source_features <- selection$features[local_features]
    beta <- fmridataset::source_read(
      estimate_assay$source,
      observations = selection$observations,
      features = source_features
    )
    sampling_variance <- fmridataset::source_read(
      variance_assay$source,
      observations = selection$observations,
      features = source_features
    )
    if (any(!is.finite(sampling_variance)) || any(sampling_variance <= 0)) {
      stop("Known sampling variances must be finite and strictly positive.",
        call. = FALSE
      )
    }
    for (offset in seq_along(local_features)) {
      fitted <- .fit_frame_ri_feature(
        beta[, offset],
        sampling_variance[, offset],
        X,
        groups,
        fit = fit
      )
      if (is.null(fitted)) next
      feature <- local_features[[offset]]
      fixed <- seq_along(coefficient_names)
      estimate_result[fixed, feature] <- fitted$coefficient
      standard_error[fixed, feature] <- fitted$standard_error
      statistic[fixed, feature] <- fitted$statistic
      p_value[fixed, feature] <- fitted$p_value
      estimate_result[length(fixed) + 1L, feature] <- fitted$vc_intercept
      estimate_result[length(fixed) + 2L, feature] <- fitted$vc_resid
      convergence[[feature]] <- fitted$converged
      log_likelihood[[feature]] <- fitted$log_likelihood
    }
  }

  result <- fmridataset::fmri_frame(
    assays = list(
      estimate = estimate_result,
      std_error = standard_error,
      statistic = statistic,
      p_value = p_value
    ),
    observations = row_data,
    features = fmridataset::feature_axis(frame),
    active_assay = "estimate",
    metadata = list(
      method = "variance_aware_random_intercept",
      fit = fit,
      term_data = multidesign::term_data(compiled),
      convergence = convergence,
      log_likelihood = log_likelihood,
      source_observation_ids = compiled$observation_ids
    ),
    provenance = list(
      operation = "fmrigds::fit_group",
      estimate_assay = estimate,
      variance_assay = variance,
      memory_budget = memory_budget,
      block_size = width
    )
  )
  structure(
    list(
      result = result,
      design = compiled,
      diagnostics = list(
        converged = convergence,
        log_likelihood = log_likelihood,
        block_size = width
      )
    ),
    class = "fmri_group_fit"
  )
}

#' @export
print.fmri_group_fit <- function(x, ...) {
  cat("<fmri_group_fit>", ncol(x$result), "features\n")
  cat("  model: variance-aware random intercept\n")
  cat(
    "  converged:", sum(x$diagnostics$converged), "/",
    length(x$diagnostics$converged), "\n"
  )
  invisible(x)
}
