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

.frame_fit_block_size <- function(n_observation, n_feature, n_input_arrays,
                                  memory_budget, block_size) {
  if (!is.numeric(memory_budget) || length(memory_budget) != 1L ||
    !is.finite(memory_budget) || memory_budget <= 0) {
    stop("`memory_budget` must be one positive finite byte count.", call. = FALSE)
  }
  bytes_per_feature <- max(8 * n_observation * n_input_arrays, 1)
  maximum <- floor(memory_budget / bytes_per_feature)
  if (maximum < 1L) {
    stop(
      "`memory_budget` cannot hold one required assay feature block.",
      call. = FALSE
    )
  }
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

.frame_fixed_result_rows <- function(coefficient_names) {
  data.frame(
    .obs_id = coefficient_names,
    result_type = rep("fixed_effect", length(coefficient_names)),
    term = coefficient_names,
    stringsAsFactors = FALSE
  )
}

.fit_frame_ri_knownvar_block <- function(arrays, design, options) {
  beta <- arrays$beta
  sampling_variance <- arrays$var
  X <- design$X
  groups <- design$groups
  fit_mode <- match.arg(options$fit %||% "REML", c("REML", "ML"))
  coefficient_names <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
  rows <- .frame_result_rows(coefficient_names)
  n_result <- nrow(rows)
  n_feature <- ncol(beta)
  estimate <- matrix(NA_real_, nrow = n_result, ncol = n_feature)
  std_error <- matrix(NA_real_, nrow = n_result, ncol = n_feature)
  statistic <- matrix(NA_real_, nrow = n_result, ncol = n_feature)
  p_value <- matrix(NA_real_, nrow = n_result, ncol = n_feature)
  converged <- logical(n_feature)
  log_likelihood <- rep(NA_real_, n_feature)
  n_obs <- integer(n_feature)

  nonpositive <- is.finite(sampling_variance) & sampling_variance <= 0
  if (any(nonpositive)) {
    stop("Known sampling variances must be strictly positive when finite.",
      call. = FALSE
    )
  }
  for (feature in seq_len(n_feature)) {
    valid <- is.finite(beta[, feature]) &
      is.finite(sampling_variance[, feature]) &
      sampling_variance[, feature] > 0 &
      stats::complete.cases(X) & !is.na(groups)
    n_obs[[feature]] <- sum(valid)
    fitted <- .fit_frame_ri_feature(
      beta[, feature],
      sampling_variance[, feature],
      X,
      groups,
      fit = fit_mode
    )
    if (is.null(fitted)) next
    fixed <- seq_along(coefficient_names)
    estimate[fixed, feature] <- fitted$coefficient
    std_error[fixed, feature] <- fitted$standard_error
    statistic[fixed, feature] <- fitted$statistic
    p_value[fixed, feature] <- fitted$p_value
    estimate[length(fixed) + 1L, feature] <- fitted$vc_intercept
    estimate[length(fixed) + 2L, feature] <- fitted$vc_resid
    converged[[feature]] <- fitted$converged
    log_likelihood[[feature]] <- fitted$log_likelihood
  }
  list(
    assays = list(
      estimate = estimate,
      std_error = std_error,
      statistic = statistic,
      p_value = p_value
    ),
    observations = rows,
    diagnostics = list(
      converged = converged,
      log_likelihood = log_likelihood,
      n_obs = n_obs
    )
  )
}

.fit_frame_ols_block <- function(arrays, design, options) {
  beta <- arrays$beta
  X <- design$X
  result <- get_reducer("ols:voxelwise")$fun(
    beta, NULL, X, NULL, NULL, NULL, NULL, NULL, options
  )
  coefficient_names <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
  list(
    assays = list(
      estimate = result$coef,
      std_error = result$se_coef,
      statistic = result$t_coef,
      p_value = result$p_coef
    ),
    observations = .frame_fixed_result_rows(coefficient_names),
    diagnostics = list(
      sigma2 = as.numeric(result$sigma2),
      degrees_freedom = as.numeric(result$df_res),
      n_obs = as.integer(result$n_obs)
    )
  )
}

.frame_assay_names <- function(frame) names(fmridataset::assays(frame))

.frame_reducer_digest <- function(reducer) {
  digest::digest(
    list(
      name = reducer$name,
      requires = reducer$requires,
      provides = reducer$provides,
      options_schema = reducer$options_schema,
      input_shape = reducer$input_shape,
      model_contract = reducer$model_contract,
      frame_formals = formals(reducer$frame_fun),
      frame_body = deparse(body(reducer$frame_fun), width.cutoff = 500L)
    ),
    algo = "xxhash64"
  )
}

.frame_plan_assay_fingerprint <- function(frame, name) {
  if (is.null(name)) return(NULL)
  selection <- .frame_fit_selection(frame)
  descriptor <- fmridataset::assay(selection$base, name)
  fmridataset::source_fingerprint(descriptor$source)
}

.validate_frame_reducer_assays <- function(reducer, estimate, variance, frame) {
  required <- setdiff(reducer$requires %||% character(), "X")
  supported <- list(beta = estimate, var = variance)
  unavailable <- required[
    !required %in% names(supported) |
      vapply(required, function(name) is.null(supported[[name]]), logical(1))
  ]
  if (length(unavailable)) {
    stop(
      "Reducer '", reducer$name, "' requires frame inputs: ",
      paste(unavailable, collapse = ", "), ".",
      call. = FALSE
    )
  }
  chosen <- unname(unlist(supported[intersect(required, names(supported))],
    use.names = FALSE
  ))
  missing <- setdiff(chosen, .frame_assay_names(frame))
  if (length(missing)) {
    stop("Frame assay not found: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Construct a lazy frame-native group-analysis plan
#'
#' `group_plan()` compiles design metadata and resolves stable axis IDs without
#' reading imaging assays. Numerical values are read only by [compute()], in
#' bounded feature blocks, and dispatched through the selected reducer's frame
#' execution contract.
#'
#' @param frame An `fmri_frame` or synchronized view.
#' @param design A `multidesign::design_spec` or compiled design.
#' @param method Registered reducer name with frame execution support.
#' @param estimate Effect-assay name.
#' @param variance Known sampling-variance assay name, or `NULL` for reducers
#'   that do not consume it.
#' @param memory_budget Maximum bytes for input assay blocks.
#' @param block_size Optional explicit feature-block width.
#' @param options Reducer options.
#' @return A serializable `fmri_group_plan`.
#' @export
group_plan <- function(
  frame,
  design,
  method = "lmm:ri_knownvar",
  estimate = "beta",
  variance = "variance",
  memory_budget = 256 * 1024^2,
  block_size = NULL,
  options = list()
) {
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  if (!is.list(options)) stop("`options` must be a list.", call. = FALSE)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)) {
    stop("`method` must be one non-empty reducer name.", call. = FALSE)
  }
  method <- .normalize_reducer_name(method)
  reducer <- get_reducer(method)
  if (is.null(reducer)) stop("Unknown reducer: ", method, ".", call. = FALSE)
  if (!is.function(reducer$frame_fun)) {
    stop("Reducer '", method, "' does not support fmri_frame execution.",
      call. = FALSE
    )
  }
  if (!"var" %in% (reducer$requires %||% character())) variance <- NULL
  .validate_frame_reducer_assays(reducer, estimate, variance, frame)
  compiled <- if (inherits(design, "design_spec")) {
    multidesign::compile_design(frame, design)
  } else {
    design
  }
  if (!inherits(compiled, "compiled_design")) {
    stop("`design` must be a design_spec or compiled_design.", call. = FALSE)
  }
  frame_ids <- fmridataset::observation_ids(frame)
  observation_index <- match(compiled$observation_ids, frame_ids)
  if (anyNA(observation_index) || anyDuplicated(observation_index)) {
    stop("Compiled design observation IDs do not align uniquely with the frame.",
      call. = FALSE
    )
  }
  X <- multidesign::model_matrix(compiled)
  if (qr(X)$rank < ncol(X)) {
    stop("The fixed-effects design matrix is rank deficient.", call. = FALSE)
  }
  groups <- if (identical(method, "lmm:ri_knownvar")) {
    .frame_random_intercept_only(compiled)
  } else {
    NULL
  }
  if (ncol(frame) < 1L) {
    stop("Frame group plans require at least one feature.", call. = FALSE)
  }
  width <- .frame_fit_block_size(
    nrow(X),
    ncol(frame),
    length(intersect(reducer$requires %||% character(), c("beta", "var"))),
    memory_budget,
    block_size
  )
  structure(
    list(
      frame = frame,
      design = compiled,
      method = method,
      estimate = estimate,
      variance = variance,
      observation_index = as.integer(observation_index),
      observation_ids = compiled$observation_ids,
      feature_ids = fmridataset::feature_ids(frame),
      space_digest = fmridataset::space_digest(fmridataset::space(frame)),
      block_size = width,
      memory_budget = memory_budget,
      options = validate_reducer_options(
        reducer$options_schema %||% list(), options
      ),
      execution_design = list(X = X, groups = groups),
      estimate_source_fingerprint = .frame_plan_assay_fingerprint(frame, estimate),
      variance_source_fingerprint = .frame_plan_assay_fingerprint(frame, variance),
      reducer_digest = .frame_reducer_digest(reducer),
      schema_version = 1L
    ),
    class = "fmri_group_plan"
  )
}

.combine_frame_block_results <- function(blocks) {
  if (!length(blocks)) stop("Frame group plans require at least one feature.",
    call. = FALSE
  )
  observations <- blocks[[1L]]$observations
  if (!all(vapply(blocks, function(x) identical(x$observations, observations), logical(1)))) {
    stop("Frame reducer returned inconsistent result observation metadata.",
      call. = FALSE
    )
  }
  assay_names <- names(blocks[[1L]]$assays)
  if (!all(vapply(blocks, function(x) identical(names(x$assays), assay_names), logical(1)))) {
    stop("Frame reducer returned inconsistent assay names across blocks.",
      call. = FALSE
    )
  }
  assays <- setNames(lapply(assay_names, function(name) {
    do.call(cbind, lapply(blocks, function(x) x$assays[[name]]))
  }), assay_names)
  diagnostic_names <- unique(unlist(lapply(blocks, function(x) names(x$diagnostics))))
  diagnostics <- setNames(lapply(diagnostic_names, function(name) {
    unlist(lapply(blocks, function(x) x$diagnostics[[name]]), use.names = FALSE)
  }), diagnostic_names)
  list(assays = assays, observations = observations, diagnostics = diagnostics)
}

.compute_group_plan <- function(plan) {
  if (!inherits(plan, "fmri_group_plan") || !identical(plan$schema_version, 1L)) {
    stop("`plan` must be a valid fmri_group_plan.", call. = FALSE)
  }
  reducer <- get_reducer(plan$method)
  if (is.null(reducer) || !is.function(reducer$frame_fun)) {
    stop("The planned frame reducer is unavailable: ", plan$method, ".",
      call. = FALSE
    )
  }
  current_observation_ids <- fmridataset::observation_ids(plan$frame)
  current_feature_ids <- fmridataset::feature_ids(plan$frame)
  current_space_digest <- fmridataset::space_digest(
    fmridataset::space(plan$frame)
  )
  if (!identical(current_observation_ids[plan$observation_index], plan$observation_ids) ||
      !identical(current_feature_ids, plan$feature_ids) ||
      !identical(current_space_digest, plan$space_digest)) {
    stop("Frame plan axes or spatial identity changed after planning.",
      call. = FALSE
    )
  }
  if (!identical(.frame_reducer_digest(reducer), plan$reducer_digest)) {
    stop("Frame reducer contract changed after planning.", call. = FALSE)
  }
  if (!identical(
    .frame_plan_assay_fingerprint(plan$frame, plan$estimate),
    plan$estimate_source_fingerprint
  ) || !identical(
    .frame_plan_assay_fingerprint(plan$frame, plan$variance),
    plan$variance_source_fingerprint
  )) {
    stop("Frame assay source changed after planning.", call. = FALSE)
  }
  selection <- .frame_fit_selection(plan$frame)
  source_observations <- selection$observations[plan$observation_index]
  estimate_assay <- fmridataset::assay(selection$base, plan$estimate)
  variance_assay <- if (is.null(plan$variance)) NULL else {
    fmridataset::assay(selection$base, plan$variance)
  }
  starts <- seq.int(1L, length(selection$features), by = plan$block_size)
  blocks <- lapply(starts, function(start) {
    local_features <- seq.int(
      start, min(start + plan$block_size - 1L, length(selection$features))
    )
    source_features <- selection$features[local_features]
    arrays <- list(beta = fmridataset::source_read(
      estimate_assay$source,
      observations = source_observations,
      features = source_features
    ))
    if (!is.null(variance_assay)) {
      arrays$var <- fmridataset::source_read(
        variance_assay$source,
        observations = source_observations,
        features = source_features
      )
    }
    reducer$frame_fun(arrays, plan$execution_design, plan$options)
  })
  combined <- .combine_frame_block_results(blocks)
  result <- fmridataset::fmri_frame(
    assays = combined$assays,
    observations = combined$observations,
    features = fmridataset::feature_axis(plan$frame),
    active_assay = "estimate",
    metadata = list(
      method = plan$method,
      term_data = multidesign::term_data(plan$design),
      diagnostics = combined$diagnostics,
      source_observation_ids = plan$observation_ids
    ),
    provenance = list(
      operation = "fmrigds::compute_group_plan",
      method = plan$method,
      estimate_assay = plan$estimate,
      variance_assay = plan$variance,
      memory_budget = plan$memory_budget,
      block_size = plan$block_size
    )
  )
  structure(
    list(
      result = result,
      design = plan$design,
      plan = plan,
      diagnostics = c(
        combined$diagnostics,
        list(block_size = plan$block_size)
      )
    ),
    class = "fmri_group_fit"
  )
}

#' @export
print.fmri_group_plan <- function(x, ...) {
  cat("<fmri_group_plan>", x$method, "\n")
  cat("  observations:", length(x$observation_ids), "\n")
  cat("  features:", length(x$feature_ids), "\n")
  cat("  block size:", x$block_size, "\n")
  invisible(x)
}

#' @export
explain.fmri_group_plan <- function(x, ...) {
  list(
    class = "fmri_group_plan",
    method = x$method,
    observations = length(x$observation_ids),
    features = length(x$feature_ids),
    assays = Filter(Negate(is.null), list(
      estimate = x$estimate,
      variance = x$variance
    )),
    space_digest = x$space_digest,
    block_size = x$block_size,
    memory_budget = x$memory_budget,
    digest = digest_plan(x)
  )
}

#' Fit a registered group reducer over frame feature blocks
#'
#' `fit_group()` is the eager convenience wrapper for [group_plan()] followed
#' by [compute()]. The default known-variance random-intercept model is
#' `diag(sampling_variance) + vc_resid I + vc_intercept Z Z'`; other registered
#' frame kernels may use different input and model contracts.
#'
#' @param frame An `fmri_frame` or synchronized view.
#' @param estimate Name of the effect assay.
#' @param variance Name of the known sampling-variance assay.
#' @param design A `design_spec` or `compiled_design`.
#' @param memory_budget Maximum bytes for the two input feature blocks.
#' @param block_size Optional explicit feature-block width.
#' @param fit Likelihood criterion, `"REML"` or `"ML"`.
#' @param method Registered frame-capable reducer. The default is
#'   `"lmm:ri_knownvar"`; `"ols:voxelwise"` provides unweighted fixed-effects
#'   fitting when no sampling-variance assay is available.
#' @param options Named reducer options.
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
  fit = c("REML", "ML"),
  method = "lmm:ri_knownvar",
  options = list()
) {
  if (!is.list(options)) stop("`options` must be a list.", call. = FALSE)
  fit <- match.arg(fit)
  normalized_method <- if (is.character(method) && length(method) == 1L &&
      !is.na(method)) .normalize_reducer_name(method) else ""
  if (startsWith(normalized_method, "lmm:") && is.null(options$fit)) {
    options$fit <- fit
  }
  compute(group_plan(
    frame = frame,
    design = design,
    method = method,
    estimate = estimate,
    variance = variance,
    memory_budget = memory_budget,
    block_size = block_size,
    options = options
  ))
}

#' @export
print.fmri_group_fit <- function(x, ...) {
  cat("<fmri_group_fit>", ncol(x$result), "features\n")
  cat("  reducer:", x$plan$method, "\n")
  if (!is.null(x$diagnostics$converged)) {
    cat(
      "  converged:", sum(x$diagnostics$converged), "/",
      length(x$diagnostics$converged), "\n"
    )
  }
  invisible(x)
}
