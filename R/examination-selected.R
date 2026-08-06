# Selected-subject localization and random-effects refits ------------------

.initialize_selected_pass <- function(state, selected_subjects, model_context) {
  selected_subjects <- as.character(selected_subjects %||% character())
  selected_index <- match(selected_subjects, state$subjects)
  if (anyNA(selected_index)) {
    stop("Selected localization subjects do not match the model context.", call. = FALSE)
  }
  state$selected_subjects <- selected_subjects
  state$selected_index <- as.integer(selected_index)
  state$selected_maps <- list()
  state$selected_map_modes <- character()
  n_selected <- length(selected_subjects)
  n_contrast <- length(state$contrasts)
  n_estimand <- nrow(model_context$estimand_matrix)
  state$exact_influence_sum_sq <- array(0, c(n_selected, n_contrast, n_estimand))
  state$exact_influence_count <- array(0, c(n_selected, n_contrast, n_estimand))
  state$exact_influence_max <- array(NA_real_, c(n_selected, n_contrast, n_estimand))
  state
}

.accumulate_selected_subject_maps <- function(state,
                                              beta,
                                              var,
                                              diagnostic,
                                              fit,
                                              contrast_index,
                                              ordinal,
                                              model_context,
                                              control) {
  if (!length(state$selected_index)) return(state)
  mode <- diagnostic$mode
  state <- .selected_map_assign(
    state, "observed", t(beta[state$selected_index, , drop = FALSE]),
    ordinal, contrast_index, "observed"
  )
  state <- .selected_map_assign(
    state, "expected", t(diagnostic$expected[state$selected_index, , drop = FALSE]),
    ordinal, contrast_index, mode
  )
  state <- .selected_map_assign(
    state, "predictive_residual",
    t(diagnostic$predictive_resid[state$selected_index, , drop = FALSE]),
    ordinal, contrast_index, mode
  )
  for (e in seq_along(state$estimands)) {
    name <- state$estimands[e]
    contribution <- matrix(
      diagnostic$fit_contribution[state$selected_index, e, ],
      nrow = length(state$selected_index),
      ncol = ncol(beta)
    )
    delta_effect <- matrix(
      diagnostic$delta_effect[state$selected_index, e, ],
      nrow = length(state$selected_index),
      ncol = ncol(beta)
    )
    delta_stat <- matrix(
      diagnostic$delta_stat[state$selected_index, e, ],
      nrow = length(state$selected_index),
      ncol = ncol(beta)
    )
    state <- .selected_map_assign(
      state,
      paste0("fit_contribution:", name),
      t(contribution),
      ordinal, contrast_index, mode
    )
    state <- .selected_map_assign(
      state,
      paste0("delta_effect:", name),
      t(delta_effect),
      ordinal, contrast_index, mode
    )
    state <- .selected_map_assign(
      state,
      paste0("delta_stat:", name),
      t(delta_stat),
      ordinal, contrast_index, mode
    )
  }

  if (model_context$method %in% c("meta:re", "meta:re_reg")) {
    exact <- .exact_random_deletion_block(
      beta,
      var,
      model_context$X,
      model_context$estimand_matrix,
      fit,
      state$selected_index,
      model_context$reducer,
      model_context$options,
      control$tolerance
    )
    state <- .selected_map_assign(
      state, "expected_exact", t(exact$expected),
      ordinal, contrast_index, "tau2_refit_exact"
    )
    state <- .selected_map_assign(
      state, "predictive_residual_exact", t(exact$predictive_resid),
      ordinal, contrast_index, "tau2_refit_exact"
    )
    state <- .selected_map_assign(
      state, "tau2_deleted_exact", t(exact$tau2_deleted),
      ordinal, contrast_index, "tau2_refit_exact"
    )
    for (e in seq_along(state$estimands)) {
      name <- state$estimands[e]
      exact_delta_effect <- matrix(
        exact$delta_effect[, e, ],
        nrow = length(state$selected_index),
        ncol = ncol(beta)
      )
      exact_delta_stat <- matrix(
        exact$delta_stat[, e, ],
        nrow = length(state$selected_index),
        ncol = ncol(beta)
      )
      state <- .selected_map_assign(
        state,
        paste0("delta_effect_exact:", name),
        t(exact_delta_effect),
        ordinal, contrast_index, "tau2_refit_exact"
      )
      state <- .selected_map_assign(
        state,
        paste0("delta_stat_exact:", name),
        t(exact_delta_stat),
        ordinal, contrast_index, "tau2_refit_exact"
      )
      for (j in seq_along(state$selected_index)) {
        values <- exact$delta_stat[j, e, ]
        values <- values[is.finite(values)]
        if (!length(values)) next
        state$exact_influence_sum_sq[j, contrast_index, e] <-
          state$exact_influence_sum_sq[j, contrast_index, e] +
          sum(pmin(values^2, control$geometry$cap^2))
        state$exact_influence_count[j, contrast_index, e] <-
          state$exact_influence_count[j, contrast_index, e] + length(values)
        state$exact_influence_max[j, contrast_index, e] <- max(
          c(state$exact_influence_max[j, contrast_index, e], abs(values)),
          na.rm = TRUE
        )
      }
    }
  }
  state
}

.selected_map_assign <- function(state,
                                 name,
                                 values,
                                 ordinal,
                                 contrast_index,
                                 mode) {
  n_selected <- length(state$selected_subjects)
  if (!is.matrix(values)) {
    values <- matrix(values, nrow = length(ordinal), ncol = n_selected)
  }
  if (!identical(dim(values), c(length(ordinal), n_selected))) {
    stop("Selected-subject map '", name, "' has invalid block dimensions.", call. = FALSE)
  }
  if (is.null(state$selected_maps[[name]])) {
    state$selected_maps[[name]] <- array(
      NA_real_,
      c(state$n_sample, n_selected, length(state$contrasts)),
      dimnames = list(NULL, state$selected_subjects, state$contrasts)
    )
    state$selected_map_modes[[name]] <- mode
  }
  state$selected_maps[[name]][ordinal, , contrast_index] <- values
  state
}

.exact_random_deletion_block <- function(beta,
                                         var,
                                         X,
                                         estimands,
                                         full_fit,
                                         selected_index,
                                         reducer,
                                         options,
                                         tolerance) {
  n_selected <- length(selected_index)
  n_sample <- ncol(beta)
  n_estimand <- nrow(estimands)
  out <- list(
    expected = matrix(NA_real_, n_selected, n_sample),
    predictive_resid = matrix(NA_real_, n_selected, n_sample),
    tau2_deleted = matrix(NA_real_, n_selected, n_sample),
    delta_effect = array(NA_real_, c(n_selected, n_estimand, n_sample)),
    delta_stat = array(NA_real_, c(n_selected, n_estimand, n_sample))
  )
  options <- validate_reducer_options(
    reducer$options_schema %||% list(),
    options %||% list()
  )
  for (j in seq_along(selected_index)) {
    i <- selected_index[j]
    keep <- seq_len(nrow(beta)) != i
    X_minus <- if (is.null(X)) NULL else X[keep, , drop = FALSE]
    deleted <- reducer$fun(
      beta[keep, , drop = FALSE],
      var[keep, , drop = FALSE],
      X_minus,
      NULL, NULL, NULL, NULL, NULL,
      options
    )
    if (identical(reducer$name, "meta:re")) {
      full_effect <- as.numeric(full_fit$beta_g)
      full_stat <- as.numeric(full_fit$z_g)
      deleted_effect <- as.numeric(deleted$beta_g)
      deleted_stat <- as.numeric(deleted$z_g)
      out$expected[j, ] <- deleted_effect
      predictive_variance <- var[i, ] + deleted$tau2 + deleted$var_g
      eligible <- is.finite(beta[i, ]) & is.finite(predictive_variance) &
        predictive_variance > tolerance$degeneracy
      out$predictive_resid[j, eligible] <-
        (beta[i, eligible] - deleted_effect[eligible]) /
        sqrt(predictive_variance[eligible])
      out$delta_effect[j, 1L, ] <- full_effect - deleted_effect
      out$delta_stat[j, 1L, ] <- full_stat - deleted_stat
      out$tau2_deleted[j, ] <- deleted$tau2
    } else {
      full_effect <- estimands %*% full_fit$coef
      full_se <- sqrt(vapply(seq_len(n_sample), function(b) {
        tau2 <- full_fit$tau2[b]
        valid <- is.finite(beta[, b]) & is.finite(var[, b]) & var[, b] > 0
        if (!is.finite(tau2) || sum(valid) < ncol(X)) return(rep(NA_real_, n_estimand))
        w <- 1 / (var[valid, b] + tau2)
        A <- .diagnostic_inverse(crossprod(X[valid, , drop = FALSE] * sqrt(w)))
        if (is.null(A)) return(rep(NA_real_, n_estimand))
        rowSums((estimands %*% A) * estimands)
      }, numeric(n_estimand)))
      if (n_estimand == 1L) full_se <- matrix(full_se, 1L, n_sample)
      full_stat <- full_effect / full_se
      deleted_effect <- estimands %*% deleted$coef
      deleted_se <- matrix(NA_real_, n_estimand, n_sample)
      for (b in seq_len(n_sample)) {
        tau2 <- deleted$tau2[b]
        valid <- is.finite(beta[keep, b]) & is.finite(var[keep, b]) & var[keep, b] > 0
        if (!is.finite(tau2) || sum(valid) < ncol(X_minus)) next
        w <- 1 / (var[keep, b][valid] + tau2)
        Xm <- X_minus[valid, , drop = FALSE]
        A <- .diagnostic_inverse(crossprod(Xm * sqrt(w)))
        if (is.null(A)) next
        deleted_se[, b] <- sqrt(rowSums((estimands %*% A) * estimands))
        if (is.finite(beta[i, b]) && is.finite(var[i, b]) && var[i, b] > 0) {
          theta <- deleted$coef[, b]
          expected <- sum(X[i, ] * theta)
          prediction_variance <- var[i, b] + tau2 +
            drop(X[i, , drop = FALSE] %*% A %*% X[i, ])
          if (is.finite(prediction_variance) && prediction_variance > tolerance$degeneracy) {
            out$expected[j, b] <- expected
            out$predictive_resid[j, b] <-
              (beta[i, b] - expected) / sqrt(prediction_variance)
          }
        }
      }
      deleted_stat <- deleted_effect / deleted_se
      out$delta_effect[j, , ] <- full_effect - deleted_effect
      out$delta_stat[j, , ] <- full_stat - deleted_stat
      out$tau2_deleted[j, ] <- deleted$tau2
    }
  }
  out
}

.finalize_selected_pass <- function(state, compiled, model_context, control) {
  subject_maps <- if (!length(state$selected_subjects)) {
    NULL
  } else {
    components <- .examination_output_components(compiled)
    new_gds(
      assays = state$selected_maps,
      space = components$space,
      subjects = state$selected_subjects,
      contrasts = state$contrasts,
      row_data = components$row_data,
      metadata = list(
        examination = list(
          assay_modes = state$selected_map_modes,
          selected_subjects = state$selected_subjects,
          selection_frozen_before_exact_refit = TRUE,
          model_context_digest = model_context$digest
        )
      )
    )
  }
  exact <- list()
  index <- 1L
  if (model_context$method %in% c("meta:re", "meta:re_reg")) {
    for (k in seq_along(state$contrasts)) {
      for (e in seq_along(state$estimands)) {
        for (j in seq_along(state$selected_subjects)) {
          count <- state$exact_influence_count[j, k, e]
          exact[[index]] <- data.frame(
            subject = state$selected_subjects[j],
            contrast = state$contrasts[k],
            estimand = state$estimands[e],
            metric = "delta_stat",
            mode = "tau2_refit_exact",
            influence_energy = if (count > 0) {
              sqrt(state$exact_influence_sum_sq[j, k, e] / count)
            } else {
              NA_real_
            },
            influence_cap = control$geometry$cap,
            max_abs_delta_stat = state$exact_influence_max[j, k, e],
            abs_delta_q90_approx = NA_real_,
            abs_delta_q95_approx = NA_real_,
            abs_delta_q99_approx = NA_real_,
            eligible_n = count,
            status = if (count > 0) "available" else "nonestimable",
            stability = NA_real_,
            stable = NA,
            ranking_stage = "selected_refit",
            stringsAsFactors = FALSE
          )
          index <- index + 1L
        }
      }
    }
  }
  list(
    subject_maps = subject_maps,
    exact_estimand_data = if (length(exact)) do.call(rbind, exact) else NULL,
    retained_map_count = length(state$selected_maps)
  )
}

.examination_output_components <- function(compiled) {
  sample_index <- compiled$axis_selection$sample
  source_space <- compiled$plan$source$probe$space
  full_n <- compiled$plan$source$probe$dims[["sample"]]
  output_space <- if (length(sample_index) == full_n &&
                      identical(sample_index, seq_len(full_n))) {
    source_space
  } else if (inherits(source_space, "space_voxel")) {
    space_subset(source_space, sample_index, pack = TRUE)
  } else {
    space_subset(source_space, sample_index)
  }
  source_row_data <- row_data(compiled$plan)
  output_row_data <- if (is.null(source_row_data)) NULL else
    source_row_data[sample_index, , drop = FALSE]
  list(space = output_space, row_data = output_row_data)
}
