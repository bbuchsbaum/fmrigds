# Deterministic two-pass residual geometry ---------------------------------

.geometry_projection <- function(feature_ids, dimension, seed, nonzero = 3L) {
  n_feature <- length(feature_ids)
  if (!n_feature || dimension < 1L) return(matrix(numeric(), n_feature, 0L))
  nonzero <- min(as.integer(nonzero), dimension)
  omega <- matrix(0, n_feature, dimension)
  seed_code <- utf8ToInt(seed)
  seed_hash <- sum(as.numeric(seed_code) * seq_along(seed_code)) %% 2147483629
  for (i in seq_along(feature_ids)) {
    code <- utf8ToInt(feature_ids[i])
    base <- (sum(as.numeric(code) * seq_along(code)) + seed_hash) %% 2147483629
    used <- integer()
    step <- 1L
    while (length(used) < nonzero) {
      mixed <- (base * 48271 + step * 69621 + step^2 * 1237) %% 2147483629
      candidate <- as.integer(mixed %% dimension) + 1L
      if (!candidate %in% used) used <- c(used, candidate)
      step <- step + 1L
    }
    signs <- ifelse(
      ((base + seq_len(nonzero) * 130363) %% 2) == 0,
      1,
      -1
    )
    omega[i, used] <- signs / sqrt(nonzero)
  }
  omega
}

.geometry_residual_block <- function(diagnostic,
                                     sample_labels,
                                     contrast,
                                     state,
                                     control) {
  residual <- diagnostic$predictive_resid
  eligible <- diagnostic$surprise_eligible & is.finite(residual)
  residual[!eligible] <- 0
  residual <- pmin(residual, control$geometry$cap)
  residual <- pmax(residual, -control$geometry$cap)
  if (isTRUE(control$geometry$balance_contrasts)) {
    residual <- residual / sqrt(max(1, state$n_sample))
  }
  list(
    E = residual,
    feature_ids = paste(contrast, sample_labels, sep = "|")
  )
}

.accumulate_geometry_pass1 <- function(state,
                                       diagnostic,
                                       sample_labels,
                                       contrast,
                                       seed,
                                       control) {
  if (state$geometry_projection_dimension < 1L) return(state)
  block <- .geometry_residual_block(
    diagnostic, sample_labels, contrast, state, control
  )
  omega <- .geometry_projection(
    block$feature_ids,
    state$geometry_projection_dimension,
    seed
  )
  state$geometry_Y <- state$geometry_Y + block$E %*% omega
  state$geometry_pass1_energy <- state$geometry_pass1_energy + sum(block$E^2)
  state
}

.prepare_geometry_basis <- function(state, control) {
  Y <- state$geometry_Y
  if (!length(Y) || !any(is.finite(Y)) || max(abs(Y), na.rm = TRUE) <=
      control$tolerance$degeneracy) {
    return(list(
      Q = matrix(numeric(), nrow(Y), 0L),
      sketch_rank = 0L,
      requested_rank = control$geometry$rank,
      pass1_energy = state$geometry_pass1_energy,
      status = "degenerate_observed"
    ))
  }
  qr_y <- qr(Y, tol = control$tolerance$rank, LAPACK = FALSE)
  rank <- min(qr_y$rank, ncol(Y), nrow(Y))
  if (rank < 1L) {
    Q <- matrix(numeric(), nrow(Y), 0L)
  } else {
    Q <- qr.Q(qr_y, complete = FALSE)[, seq_len(rank), drop = FALSE]
  }
  list(
    Q = Q,
    sketch_rank = as.integer(rank),
    requested_rank = control$geometry$rank,
    pass1_energy = state$geometry_pass1_energy,
    status = if (rank > 0L) "available" else "degenerate_observed"
  )
}

.initialize_geometry_pass2 <- function(basis,
                                       state,
                                       control,
                                       selected_subjects = character(),
                                       model_context = NULL) {
  rank <- ncol(basis$Q)
  n_split <- control$geometry$stability_replicates
  out <- list(
    basis = basis,
    C = matrix(0, rank, rank),
    split_C = array(0, c(rank, rank, n_split)),
    total_energy = 0,
    split_energy = numeric(n_split),
    n_sample = state$n_sample,
    subjects = state$subjects,
    contrasts = state$contrasts,
    estimands = state$estimands
  )
  if (!is.null(model_context)) {
    out <- .initialize_selected_pass(out, selected_subjects, model_context)
  }
  out
}

.update_geometry_pass2_from_arrays <- function(state,
                                               arrays,
                                               block,
                                               model_context,
                                               control) {
  reducer <- model_context$reducer
  arrays <- .ensure_required_arrays(arrays, reducer$requires %||% character())
  opts <- validate_reducer_options(
    reducer$options_schema %||% list(),
    model_context$options %||% list()
  )
  split <- .examination_feature_split(
    block$sample_label,
    state$contrasts,
    model_context$source_plan_digest,
    control$geometry$stability_replicates
  )
  for (k in seq_along(state$contrasts)) {
    beta <- if (!is.null(arrays$beta)) .slice_subjects_samples(arrays$beta, k) else NULL
    var <- if (!is.null(arrays$var)) .slice_subjects_samples(arrays$var, k) else NULL
    z <- if (!is.null(arrays$z)) .slice_subjects_samples(arrays$z, k) else NULL
    p <- if (!is.null(arrays$p)) .slice_subjects_samples(arrays$p, k) else NULL
    fit <- reducer$fun(
      beta, var, model_context$X, z, p,
      arrays$df %||% NULL, arrays$df1 %||% NULL, arrays$df2 %||% NULL,
      opts
    )
    diagnostic <- model_context$diagnostics$fun(
      fit = fit,
      beta = beta,
      var = var,
      X = model_context$X,
      estimands = model_context$estimand_matrix,
      opts = opts,
      tolerance = control$tolerance
    )
    state <- .accumulate_geometry_pass2(
      state,
      diagnostic,
      block$sample_label,
      state$contrasts[k],
      if (length(split)) split[, k] else integer(),
      control
    )
    state <- .accumulate_selected_subject_maps(
      state,
      beta,
      var,
      diagnostic,
      fit,
      contrast_index = k,
      ordinal = block$ordinal,
      model_context = model_context,
      control = control
    )
  }
  state
}

.accumulate_geometry_pass2 <- function(state,
                                       diagnostic,
                                       sample_labels,
                                       contrast,
                                       split,
                                       control) {
  if (!ncol(state$basis$Q)) return(state)
  block <- .geometry_residual_block(
    diagnostic, sample_labels, contrast, state, control
  )
  projected <- crossprod(state$basis$Q, block$E)
  state$C <- state$C + tcrossprod(projected)
  state$total_energy <- state$total_energy + sum(block$E^2)
  for (r in seq_len(dim(state$split_C)[3L])) {
    use <- split == r
    if (!any(use)) next
    projected_split <- projected[, use, drop = FALSE]
    state$split_C[, , r] <- state$split_C[, , r] + tcrossprod(projected_split)
    state$split_energy[r] <- state$split_energy[r] + sum(block$E[, use, drop = FALSE]^2)
  }
  state
}

.finalize_residual_geometry <- function(state, control) {
  Q <- state$basis$Q
  if (!ncol(Q) || !is.finite(state$total_energy) ||
      state$total_energy <= control$tolerance$degeneracy) {
    return(list(
      coordinates = data.frame(subject = state$subjects, stringsAsFactors = FALSE),
      explained_energy = numeric(),
      captured_energy = 0,
      sketch_rank = state$basis$sketch_rank,
      requested_rank = state$basis$requested_rank,
      residual_cap = control$geometry$cap,
      contrast_balanced = control$geometry$balance_contrasts,
      status = "degenerate_observed",
      stability_method = "deterministic_split_features"
    ))
  }
  decomposition <- eigen((state$C + t(state$C)) / 2, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  keep <- which(values > control$tolerance$degeneracy * max(1, values[1L]))
  keep <- head(keep, control$geometry$rank)
  if (!length(keep)) {
    return(list(
      coordinates = data.frame(subject = state$subjects, stringsAsFactors = FALSE),
      explained_energy = numeric(),
      captured_energy = 0,
      sketch_rank = state$basis$sketch_rank,
      requested_rank = state$basis$requested_rank,
      residual_cap = control$geometry$cap,
      contrast_balanced = control$geometry$balance_contrasts,
      status = "degenerate_observed",
      stability_method = "deterministic_split_features"
    ))
  }
  values <- values[keep]
  vectors <- decomposition$vectors[, keep, drop = FALSE]
  coordinates <- Q %*% vectors %*% diag(sqrt(values), nrow = length(values))
  coordinates <- .canonicalize_geometry_signs(coordinates)
  stability <- .geometry_split_stability(
    coordinates,
    Q,
    state$split_C,
    length(values),
    control
  )
  coordinate_data <- data.frame(
    subject = state$subjects,
    coordinates,
    stability = stability,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(coordinate_data)[seq_len(ncol(coordinates)) + 1L] <-
    paste0("dimension", seq_len(ncol(coordinates)))
  explained_energy <- values / state$total_energy
  captured_energy <- sum(explained_energy)
  if (captured_energy > 1) {
    # Q is orthonormal, so projected energy cannot exceed total residual
    # energy mathematically. Normalize boundary roundoff so the public
    # fidelity remains a fraction and its components retain the same sum.
    explained_energy <- explained_energy / captured_energy
    captured_energy <- 1
  }
  list(
    coordinates = coordinate_data,
    explained_energy = explained_energy,
    captured_energy = captured_energy,
    total_residual_energy = state$total_energy,
    sketch_rank = state$basis$sketch_rank,
    requested_rank = state$basis$requested_rank,
    residual_cap = control$geometry$cap,
    contrast_balanced = control$geometry$balance_contrasts,
    status = "available",
    stability_method = "deterministic_split_features"
  )
}

.canonicalize_geometry_signs <- function(coordinates) {
  if (!ncol(coordinates)) return(coordinates)
  for (column in seq_len(ncol(coordinates))) {
    anchor <- which.max(abs(coordinates[, column]))
    if (coordinates[anchor, column] < 0) {
      coordinates[, column] <- -coordinates[, column]
    }
  }
  coordinates
}

.geometry_split_stability <- function(full, Q, split_C, rank, control) {
  n_split <- dim(split_C)[3L]
  if (!n_split || !ncol(full)) return(rep(NA_real_, nrow(full)))
  scores <- matrix(NA_real_, nrow(full), n_split)
  for (r in seq_len(n_split)) {
    C <- split_C[, , r, drop = TRUE]
    if (!length(C) || max(abs(C), na.rm = TRUE) <= control$tolerance$degeneracy) next
    decomposition <- eigen((C + t(C)) / 2, symmetric = TRUE)
    values <- pmax(decomposition$values, 0)
    use <- seq_len(min(rank, sum(values > control$tolerance$degeneracy), ncol(Q)))
    if (!length(use)) next
    split_coordinates <- Q %*% decomposition$vectors[, use, drop = FALSE] %*%
      diag(sqrt(values[use]), nrow = length(use))
    full_use <- full[, seq_len(length(use)), drop = FALSE]
    alignment <- svd(crossprod(split_coordinates, full_use))
    rotation <- alignment$u %*% t(alignment$v)
    aligned <- split_coordinates %*% rotation
    numerator <- sqrt(rowSums((aligned - full_use)^2))
    denominator <- sqrt(rowSums(aligned^2)) + sqrt(rowSums(full_use^2))
    scores[, r] <- pmax(0, 1 - numerator / pmax(denominator, control$tolerance$degeneracy))
  }
  stability <- rowMeans(scores, na.rm = TRUE)
  stability[is.nan(stability)] <- NA_real_
  stability
}
