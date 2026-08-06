# Group-examination accumulation and result object -------------------------

.initialize_examination_accumulator <- function(scan_context,
                                                model_context,
                                                control,
                                                quality_data,
                                                retain,
                                                preflight) {
  n_subject <- length(model_context$subjects)
  n_contrast <- length(scan_context$contrasts)
  n_estimand <- nrow(model_context$estimand_matrix)
  n_split <- control$geometry$stability_replicates
  histogram_upper <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8, 16, Inf)
  matrix_sc <- function(value = 0) {
    matrix(value, n_subject, n_contrast)
  }
  array_sce <- function(value = 0) {
    array(value, c(n_subject, n_contrast, n_estimand))
  }

  list(
    subjects = model_context$subjects,
    contrasts = scan_context$contrasts,
    estimands = rownames(model_context$estimand_matrix),
    n_sample = length(scan_context$source_samples),
    maps = list(),
    quality_data = quality_data,
    retain = retain %||% character(),
    preflight = preflight,
    screening_mode = if (model_context$method %in% c("meta:re", "meta:re_reg")) {
      "tau2_fixed_full"
    } else {
      "exact"
    },
    beta_finite = matrix_sc(),
    se_sum = matrix_sc(),
    se_count = matrix_sc(),
    contrast_feature_count = numeric(n_contrast),
    surprise_sum_sq = matrix_sc(),
    surprise_count = matrix_sc(),
    surprise_tail = matrix_sc(),
    surprise_split_sum_sq = array(0, c(n_subject, n_contrast, n_split)),
    surprise_split_count = array(0, c(n_subject, n_contrast, n_split)),
    surprise_split_tail = array(0, c(n_subject, n_contrast, n_split)),
    agreement_sw = matrix_sc(),
    agreement_swx = matrix_sc(),
    agreement_swy = matrix_sc(),
    agreement_swxx = matrix_sc(),
    agreement_swyy = matrix_sc(),
    agreement_swxy = matrix_sc(),
    agreement_sign = matrix_sc(),
    influence_sum_sq = array_sce(),
    influence_count = array_sce(),
    influence_max = array_sce(NA_real_),
    influence_split_sum_sq = array(0, c(n_subject, n_contrast, n_estimand, n_split)),
    influence_split_count = array(0, c(n_subject, n_contrast, n_estimand, n_split)),
    influence_hist = array(
      0,
      c(n_subject, n_contrast, n_estimand, length(histogram_upper))
    ),
    histogram_upper = histogram_upper,
    geometry_projection_dimension = min(
      n_subject,
      control$geometry$rank + control$geometry$oversample
    ),
    geometry_Y = matrix(
      0,
      n_subject,
      min(n_subject, control$geometry$rank + control$geometry$oversample)
    ),
    geometry_pass1_energy = 0
  )
}

.update_examination_accumulator <- function(state,
                                            arrays,
                                            block,
                                            scan_context,
                                            model_context,
                                            control) {
  reducer <- model_context$reducer
  arrays <- .ensure_required_arrays(arrays, reducer$requires %||% character())
  n_contrast <- length(state$contrasts)
  n_block <- dim(arrays[[1L]])[1L]
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

  for (k in seq_len(n_contrast)) {
    beta <- if (!is.null(arrays$beta)) .slice_subjects_samples(arrays$beta, k) else NULL
    var <- if (!is.null(arrays$var)) .slice_subjects_samples(arrays$var, k) else NULL
    z <- if (!is.null(arrays$z)) .slice_subjects_samples(arrays$z, k) else NULL
    p <- if (!is.null(arrays$p)) .slice_subjects_samples(arrays$p, k) else NULL
    fit <- reducer$fun(
      beta, var, model_context$X, z, p,
      arrays$df %||% NULL,
      arrays$df1 %||% NULL,
      arrays$df2 %||% NULL,
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
    state <- .accumulate_examination_maps(
      state,
      diagnostic$maps,
      block$ordinal,
      k,
      n_contrast
    )
    state <- .accumulate_examination_subjects(
      state,
      beta,
      var,
      diagnostic,
      contrast_index = k,
      split = if (length(split)) split[, k] else integer(),
      control = control
    )
    state <- .accumulate_geometry_pass1(
      state,
      diagnostic,
      block$sample_label,
      state$contrasts[k],
      model_context$digest,
      control
    )
  }
  state
}

.accumulate_examination_maps <- function(state,
                                         maps,
                                         ordinal,
                                         contrast_index,
                                         n_contrast) {
  for (name in names(maps)) {
    values <- as.vector(maps[[name]])
    if (length(values) != length(ordinal)) {
      stop("Diagnostic map '", name, "' has the wrong block length.", call. = FALSE)
    }
    if (is.null(state$maps[[name]])) {
      missing_value <- if (is.integer(maps[[name]])) NA_integer_ else NA_real_
      state$maps[[name]] <- array(
        missing_value,
        c(state$n_sample, 1L, n_contrast),
        dimnames = list(NULL, "meta", state$contrasts)
      )
    }
    state$maps[[name]][ordinal, 1L, contrast_index] <- values
  }
  state
}

.accumulate_examination_subjects <- function(state,
                                             beta,
                                             var,
                                             diagnostic,
                                             contrast_index,
                                             split,
                                             control) {
  n_subject <- nrow(beta)
  n_estimand <- length(state$estimands)
  cap <- control$geometry$cap
  residual_threshold <- control$review$surprise$residual_threshold
  state$contrast_feature_count[contrast_index] <-
    state$contrast_feature_count[contrast_index] + ncol(beta)

  for (i in seq_len(n_subject)) {
    finite_beta <- is.finite(beta[i, ])
    state$beta_finite[i, contrast_index] <-
      state$beta_finite[i, contrast_index] + sum(finite_beta)
    if (!is.null(var)) {
      valid_se <- is.finite(var[i, ]) & var[i, ] > 0
      state$se_sum[i, contrast_index] <- state$se_sum[i, contrast_index] +
        sum(sqrt(var[i, valid_se]))
      state$se_count[i, contrast_index] <- state$se_count[i, contrast_index] +
        sum(valid_se)
    }

    surprise_ok <- diagnostic$surprise_eligible[i, ] &
      is.finite(diagnostic$predictive_resid[i, ])
    if (any(surprise_ok)) {
      residual <- diagnostic$predictive_resid[i, surprise_ok]
      capped_sq <- pmin(residual^2, cap^2)
      state$surprise_sum_sq[i, contrast_index] <-
        state$surprise_sum_sq[i, contrast_index] + sum(capped_sq)
      state$surprise_count[i, contrast_index] <-
        state$surprise_count[i, contrast_index] + length(residual)
      state$surprise_tail[i, contrast_index] <-
        state$surprise_tail[i, contrast_index] +
        sum(abs(residual) >= residual_threshold)

      for (r in seq_len(dim(state$surprise_split_sum_sq)[3L])) {
        in_split <- surprise_ok & split == r
        if (!any(in_split)) next
        split_residual <- diagnostic$predictive_resid[i, in_split]
        state$surprise_split_sum_sq[i, contrast_index, r] <-
          state$surprise_split_sum_sq[i, contrast_index, r] +
          sum(pmin(split_residual^2, cap^2))
        state$surprise_split_count[i, contrast_index, r] <-
          state$surprise_split_count[i, contrast_index, r] + length(split_residual)
        state$surprise_split_tail[i, contrast_index, r] <-
          state$surprise_split_tail[i, contrast_index, r] +
          sum(abs(split_residual) >= residual_threshold)
      }
    }

    agreement_ok <- surprise_ok &
      is.finite(diagnostic$expected[i, ]) &
      is.finite(beta[i, ]) &
      is.finite(diagnostic$predictive_weight[i, ]) &
      diagnostic$predictive_weight[i, ] > 0
    if (any(agreement_ok)) {
      expected <- diagnostic$expected[i, agreement_ok]
      observed <- beta[i, agreement_ok]
      weight <- diagnostic$predictive_weight[i, agreement_ok]
      state$agreement_sw[i, contrast_index] <-
        state$agreement_sw[i, contrast_index] + sum(weight)
      state$agreement_swx[i, contrast_index] <-
        state$agreement_swx[i, contrast_index] + sum(weight * expected)
      state$agreement_swy[i, contrast_index] <-
        state$agreement_swy[i, contrast_index] + sum(weight * observed)
      state$agreement_swxx[i, contrast_index] <-
        state$agreement_swxx[i, contrast_index] + sum(weight * expected^2)
      state$agreement_swyy[i, contrast_index] <-
        state$agreement_swyy[i, contrast_index] + sum(weight * observed^2)
      state$agreement_swxy[i, contrast_index] <-
        state$agreement_swxy[i, contrast_index] + sum(weight * expected * observed)
      state$agreement_sign[i, contrast_index] <-
        state$agreement_sign[i, contrast_index] +
        sum(weight * (sign(expected) == sign(observed)))
    }

    for (e in seq_len(n_estimand)) {
      influence_ok <- diagnostic$influence_eligible[i, e, ] &
        is.finite(diagnostic$delta_stat[i, e, ])
      if (!any(influence_ok)) next
      values <- abs(diagnostic$delta_stat[i, e, influence_ok])
      state$influence_sum_sq[i, contrast_index, e] <-
        state$influence_sum_sq[i, contrast_index, e] +
        sum(pmin(values^2, cap^2))
      state$influence_count[i, contrast_index, e] <-
        state$influence_count[i, contrast_index, e] + length(values)
      current_max <- state$influence_max[i, contrast_index, e]
      state$influence_max[i, contrast_index, e] <-
        max(c(current_max, values), na.rm = TRUE)

      bins <- vapply(values, function(value) {
        which(value <= state$histogram_upper)[1L]
      }, integer(1))
      for (bin in bins) {
        state$influence_hist[i, contrast_index, e, bin] <-
          state$influence_hist[i, contrast_index, e, bin] + 1
      }
      for (r in seq_len(dim(state$influence_split_sum_sq)[4L])) {
        in_split <- influence_ok & split == r
        if (!any(in_split)) next
        split_values <- diagnostic$delta_stat[i, e, in_split]
        state$influence_split_sum_sq[i, contrast_index, e, r] <-
          state$influence_split_sum_sq[i, contrast_index, e, r] +
          sum(pmin(split_values^2, cap^2))
        state$influence_split_count[i, contrast_index, e, r] <-
          state$influence_split_count[i, contrast_index, e, r] +
          length(split_values)
      }
    }
  }
  state
}

.examination_feature_split <- function(sample_labels,
                                       contrasts,
                                       digest,
                                       n_split) {
  if (n_split < 1L) return(matrix(integer(), nrow = length(sample_labels), ncol = 0L))
  out <- matrix(NA_integer_, length(sample_labels), length(contrasts))
  for (k in seq_along(contrasts)) {
    out[, k] <- vapply(sample_labels, function(label) {
      code <- utf8ToInt(paste(label, contrasts[k], digest, sep = "|"))
      as.integer(sum(as.numeric(code %% 104729L) * seq_along(code)) %% n_split) + 1L
    }, integer(1))
  }
  out
}

.finalize_examination_accumulator <- function(state,
                                              receipt,
                                              scan_context,
                                              compiled,
                                              model_context,
                                              control,
                                              preflight) {
  contrast_data <- .finalize_contrast_data(state, control)
  estimand_data <- .finalize_estimand_data(state, control)
  subject_data <- .finalize_subject_data(
    state,
    contrast_data,
    estimand_data,
    control,
    model_context
  )
  availability <- .finalize_examination_availability(
    state,
    contrast_data,
    estimand_data,
    preflight,
    model_context
  )
  group_maps <- .finalize_examination_group_maps(
    state,
    compiled,
    model_context
  )
  posthoc_availability <- .posthoc_case_deletion_preflight(compiled$conclusion_tail)
  scientific_control <- control
  scientific_control$block_size <- NULL
  scientific_control$staging <- NULL
  provenance <- list(
    source_plan_digest = compiled$source_plan_digest,
    examination_digest = digest::digest(
      list(
        model = model_context$digest,
        control = scientific_control,
        subjects = state$subjects,
        contrasts = state$contrasts,
        estimands = state$estimands
      ),
      algo = "xxhash64"
    ),
    model_context_digest = model_context$digest,
    reducer_node_id = model_context$reducer_node_id,
    node_ids = compiled$origin_node_ids %||% compiled$node_ids,
    discarded_write_node_ids = vapply(
      compiled$origin_discarded_writes %||% compiled$discarded_writes,
      function(node) node$node_id %||% NA_character_,
      character(1)
    ),
    excluded_subjects = model_context$excluded_subjects,
    scan_receipt = receipt,
    software = list(
      package = "fmrigds",
      version = .pkg_version(),
      R_version = as.character(getRversion())
    )
  )
  cohort <- list(
    n_selected = length(model_context$selected_subject_index),
    n_included = length(state$subjects),
    n_excluded_covariates = nrow(model_context$excluded_subjects),
    n_samples = state$n_sample,
    n_contrasts = length(state$contrasts),
    model = model_context$method,
    formula = if (is.null(model_context$formula)) NULL else
      paste(deparse(model_context$formula), collapse = " "),
    estimands = state$estimands,
    variance_mode = model_context$variance_mode,
    diagnostic_modes = model_context$diagnostics$modes,
    review_n = sum(subject_data$review_status == "review"),
    insufficient_n = sum(subject_data$review_status == "insufficient")
  )
  sensitivity <- estimand_data[, c(
    "subject", "contrast", "estimand", "mode", "ranking_stage", "influence_energy",
    "max_abs_delta_stat", "abs_delta_q95_approx", "eligible_n", "stable"
  ), drop = FALSE]

  structure(
    list(
      subject_data = subject_data,
      contrast_data = contrast_data,
      estimand_data = estimand_data,
      availability = availability,
      cohort = cohort,
      embedding = NULL,
      group_maps = group_maps,
      subject_maps = NULL,
      sensitivity = sensitivity,
      conclusion = list(availability = posthoc_availability, results = NULL),
      pairwise = NULL,
      config = list(
        method = model_context$method,
        formula = cohort$formula,
        estimand_matrix = model_context$estimand_matrix,
        weight_contract = model_context$weight_contract,
        variance_mode = model_context$variance_mode,
        control = control,
        retained_subjects = subject_data$subject[subject_data$retained]
      ),
      provenance = provenance
    ),
    class = "gds_examination"
  )
}

.finalize_contrast_data <- function(state, control) {
  rows <- vector("list", length(state$subjects) * length(state$contrasts))
  index <- 1L
  tolerance <- control$tolerance$degeneracy
  for (k in seq_along(state$contrasts)) {
    for (i in seq_along(state$subjects)) {
      count <- state$surprise_count[i, k]
      surprise_energy <- if (count > 0) {
        sqrt(state$surprise_sum_sq[i, k] / count)
      } else {
        NA_real_
      }
      tail_extent <- if (count > 0) state$surprise_tail[i, k] / count else NA_real_
      split_energy <- .split_energy(
        state$surprise_split_sum_sq[i, k, ],
        state$surprise_split_count[i, k, ]
      )
      split_tail <- ifelse(
        state$surprise_split_count[i, k, ] > 0,
        state$surprise_split_tail[i, k, ] /
          state$surprise_split_count[i, k, ],
        NA_real_
      )
      stable <- .gate_stability(
        split_energy >= control$review$surprise$energy_threshold &
          split_tail >= control$review$surprise$tail_threshold
      )

      sw <- state$agreement_sw[i, k]
      correlation <- cosine <- gain <- sign_concordance <- NA_real_
      correlation_status <- cosine_status <- gain_status <- "insufficient_samples"
      if (is.finite(sw) && sw > 0) {
        sx <- state$agreement_swx[i, k]
        sy <- state$agreement_swy[i, k]
        sxx <- state$agreement_swxx[i, k]
        syy <- state$agreement_swyy[i, k]
        sxy <- state$agreement_swxy[i, k]
        expected_var <- sxx - sx^2 / sw
        observed_var <- syy - sy^2 / sw
        if (expected_var <= tolerance) {
          correlation_status <- "degenerate_expected"
        } else if (observed_var <= tolerance) {
          correlation_status <- "degenerate_observed"
        } else {
          correlation <- (sxy - sx * sy / sw) /
            sqrt(expected_var * observed_var)
          correlation_status <- "available"
        }
        if (sxx <= tolerance) {
          cosine_status <- gain_status <- "degenerate_expected"
        } else {
          gain <- sxy / sxx
          gain_status <- "available"
          if (syy <= tolerance) {
            cosine_status <- "degenerate_observed"
          } else {
            cosine <- sxy / sqrt(sxx * syy)
            cosine_status <- "available"
          }
        }
        sign_concordance <- state$agreement_sign[i, k] / sw
      }
      total <- state$contrast_feature_count[k]
      rows[[index]] <- data.frame(
        subject = state$subjects[i],
        contrast = state$contrasts[k],
        coverage_fraction = if (total > 0) state$beta_finite[i, k] / total else NA_real_,
        mean_first_level_se = if (state$se_count[i, k] > 0) {
          state$se_sum[i, k] / state$se_count[i, k]
        } else {
          NA_real_
        },
        surprise_energy = surprise_energy,
        tail_extent = tail_extent,
        surprise_eligible_n = count,
        surprise_status = if (count > 0) "available" else "insufficient_samples",
        surprise_stability = stable,
        weighted_correlation = correlation,
        correlation_status = correlation_status,
        weighted_cosine = cosine,
        cosine_status = cosine_status,
        zero_intercept_gain = gain,
        gain_status = gain_status,
        sign_concordance = sign_concordance,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  do.call(rbind, rows)
}

.finalize_estimand_data <- function(state, control) {
  rows <- vector(
    "list",
    length(state$subjects) * length(state$contrasts) * length(state$estimands)
  )
  index <- 1L
  for (k in seq_along(state$contrasts)) {
    for (e in seq_along(state$estimands)) {
      for (i in seq_along(state$subjects)) {
        count <- state$influence_count[i, k, e]
        energy <- if (count > 0) {
          sqrt(state$influence_sum_sq[i, k, e] / count)
        } else {
          NA_real_
        }
        split_energy <- .split_energy(
          state$influence_split_sum_sq[i, k, e, ],
          state$influence_split_count[i, k, e, ]
        )
        stable <- .gate_stability(
          split_energy >= control$review$influence$energy_threshold
        )
        histogram <- state$influence_hist[i, k, e, ]
        rows[[index]] <- data.frame(
          subject = state$subjects[i],
          contrast = state$contrasts[k],
          estimand = state$estimands[e],
          metric = "delta_stat",
          mode = state$screening_mode,
          influence_energy = energy,
          influence_cap = control$geometry$cap,
          max_abs_delta_stat = state$influence_max[i, k, e],
          abs_delta_q90_approx = .histogram_quantile(
            histogram, state$histogram_upper, 0.90,
            state$influence_max[i, k, e]
          ),
          abs_delta_q95_approx = .histogram_quantile(
            histogram, state$histogram_upper, 0.95,
            state$influence_max[i, k, e]
          ),
          abs_delta_q99_approx = .histogram_quantile(
            histogram, state$histogram_upper, 0.99,
            state$influence_max[i, k, e]
          ),
          eligible_n = count,
          status = if (count > 0) "available" else "nonestimable",
          stability = stable,
          stable = is.finite(stable) && stable >= control$review$min_stability,
          ranking_stage = "screening",
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
  }
  do.call(rbind, rows)
}

.split_energy <- function(sum_sq, count) {
  if (!length(count)) return(numeric())
  ifelse(count > 0, sqrt(sum_sq / count), NA_real_)
}

.gate_stability <- function(trigger) {
  available <- !is.na(trigger)
  if (!any(available)) return(NA_real_)
  mean(trigger[available])
}

.histogram_quantile <- function(counts, upper, probability, maximum) {
  total <- sum(counts)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  bin <- which(cumsum(counts) >= probability * total)[1L]
  value <- upper[bin]
  if (is.infinite(value)) maximum else value
}

.finalize_subject_data <- function(state,
                                   contrast_data,
                                   estimand_data,
                                   control,
                                   model_context = NULL) {
  subjects <- state$subjects
  out <- data.frame(subject = subjects, stringsAsFactors = FALSE)
  subject_metadata <- model_context$subject_metadata %||%
    model_context$model_frame %||% NULL
  if (!is.null(subject_metadata)) {
    subject_metadata <- subject_metadata[subjects, , drop = FALSE]
    for (name in setdiff(names(subject_metadata), "subject")) {
      value <- subject_metadata[[name]]
      if (is.atomic(value) && length(value) == length(subjects)) {
        out[[name]] <- value
      }
    }
  }
  for (name in names(state$quality_data)) {
    out[[name]] <- state$quality_data[subjects, name]
  }
  out$coverage_fraction <- vapply(subjects, function(subject) {
    min(contrast_data$coverage_fraction[contrast_data$subject == subject], na.rm = TRUE)
  }, numeric(1))
  out$surprise_score <- vapply(subjects, function(subject) {
    .safe_max(contrast_data$surprise_energy[contrast_data$subject == subject])
  }, numeric(1))
  out$influence_score <- vapply(subjects, function(subject) {
    .safe_max(estimand_data$influence_energy[estimand_data$subject == subject])
  }, numeric(1))
  out$quality_score <- 1 - out$coverage_fraction
  quality_risks <- list(out$quality_score)
  for (name in names(control$review$quality)) {
    if (identical(name, "coverage_fraction") || !name %in% names(out)) next
    spec <- control$review$quality[[name]]
    value <- suppressWarnings(as.numeric(out[[name]]))
    quality_risks[[length(quality_risks) + 1L]] <-
      if (identical(spec$direction, "high")) value else -value
  }
  if (length(quality_risks) > 1L) {
    risk_percentiles <- do.call(cbind, lapply(quality_risks, .percentile_score))
    out$quality_percentile <- apply(risk_percentiles, 1L, .safe_max)
  } else {
    out$quality_percentile <- .percentile_score(out$quality_score)
  }
  out$surprise_percentile <- .percentile_score(out$surprise_score)
  out$influence_percentile <- .percentile_score(out$influence_score)
  out$review_priority <- apply(
    cbind(out$quality_percentile, out$surprise_percentile, out$influence_percentile),
    1L,
    .safe_max
  )

  out$review_status <- "none"
  out$review_source <- NA_character_
  out$review_reason <- "No absolute review criterion met."
  for (i in seq_along(subjects)) {
    subject <- subjects[i]
    contrast_rows <- contrast_data[contrast_data$subject == subject, , drop = FALSE]
    estimand_rows <- estimand_data[estimand_data$subject == subject, , drop = FALSE]
    surprise_trigger <- with(
      contrast_rows,
      is.finite(surprise_energy) &
        surprise_energy >= control$review$surprise$energy_threshold &
        tail_extent >= control$review$surprise$tail_threshold &
        surprise_stability >= control$review$min_stability
    )
    influence_trigger <- with(
      estimand_rows,
      is.finite(influence_energy) &
        influence_energy >= control$review$influence$energy_threshold &
        max_abs_delta_stat >= control$review$influence$max_abs_threshold &
        stability >= control$review$min_stability
    )
    quality_trigger <- FALSE
    quality_reason <- NULL
    for (name in names(control$review$quality)) {
      spec <- control$review$quality[[name]]
      value <- if (identical(name, "coverage_fraction")) {
        out$coverage_fraction[i]
      } else if (name %in% names(out)) {
        suppressWarnings(as.numeric(out[[name]][i]))
      } else {
        NA_real_
      }
      triggered <- is.finite(value) && if (identical(spec$direction, "high")) {
        value >= spec$threshold
      } else {
        value <= spec$threshold
      }
      if (triggered) {
        quality_trigger <- TRUE
        quality_reason <- paste0("Data-validity criterion met for ", name, ".")
        break
      }
    }

    availability <- c(
      any(contrast_rows$surprise_status == "available"),
      any(estimand_rows$status == "available"),
      is.finite(out$coverage_fraction[i])
    )
    if (!any(availability)) {
      out$review_status[i] <- "insufficient"
      out$review_reason[i] <- "Insufficient eligible data for examination."
      next
    }
    triggers <- c(
      quality = quality_trigger,
      surprise = any(surprise_trigger),
      influence = any(influence_trigger)
    )
    if (!any(triggers)) next
    source_scores <- c(
      quality = out$quality_percentile[i],
      surprise = out$surprise_percentile[i],
      influence = out$influence_percentile[i]
    )
    source_scores[!triggers] <- -Inf
    source <- names(source_scores)[which.max(source_scores)]
    out$review_status[i] <- "review"
    out$review_source[i] <- source
    if (identical(source, "quality")) {
      out$review_reason[i] <- quality_reason %||% "Data-validity criterion met."
    } else if (identical(source, "surprise")) {
      row <- contrast_rows[which(surprise_trigger)[1L], , drop = FALSE]
      if (is.finite(row$zero_intercept_gain) && row$zero_intercept_gain < 0) {
        out$review_reason[i] <- paste0(
          "High unexpectedness with negative map gain in contrast ",
          row$contrast, "."
        )
      } else {
        out$review_reason[i] <- paste0(
          "High model surprise in contrast ", row$contrast, "."
        )
      }
    } else {
      row <- estimand_rows[which(influence_trigger)[1L], , drop = FALSE]
      out$review_reason[i] <- paste0(
        "High group-statistic influence for estimand ", row$estimand,
        " in contrast ", row$contrast, "."
      )
    }
  }

  ordering <- order(-out$review_priority, out$subject, na.last = TRUE)
  automatic <- head(out$subject[ordering], control$retain_n)
  out$retained <- out$subject %in% unique(c(
    state$retain,
    out$subject[out$review_status == "review"],
    automatic
  ))
  out
}

.safe_max <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) NA_real_ else max(value)
}

.percentile_score <- function(value) {
  out <- rep(NA_real_, length(value))
  ok <- is.finite(value)
  if (!any(ok)) return(out)
  if (sum(ok) == 1L) {
    out[ok] <- 1
  } else {
    out[ok] <- (rank(value[ok], ties.method = "average") - 1) / (sum(ok) - 1)
  }
  out
}

.finalize_examination_availability <- function(state,
                                               contrast_data,
                                               estimand_data,
                                               preflight,
                                               model_context) {
  rows <- list()
  index <- 1L
  for (diagnostic in unique(preflight$availability$diagnostic)) {
    uses_estimand <- diagnostic %in% c("coefficient_deletion", "statistic_deletion")
    targets <- if (uses_estimand) state$estimands else NA_character_
    for (contrast in state$contrasts) {
      for (estimand in targets) {
        status <- preflight$availability$status[
          match(diagnostic, preflight$availability$diagnostic)
        ]
        if (identical(status, "available")) {
          if (uses_estimand) {
            subset <- estimand_data$contrast == contrast &
              estimand_data$estimand == estimand
            if (!any(estimand_data$status[subset] == "available")) status <- "nonestimable"
          } else if (diagnostic %in% c("prediction", "surprise")) {
            subset <- contrast_data$contrast == contrast
            if (!any(contrast_data$surprise_status[subset] == "available")) {
              status <- "insufficient_samples"
            }
          }
        }
        rows[[index]] <- data.frame(
          diagnostic = diagnostic,
          contrast = contrast,
          estimand = estimand,
          mode = paste(model_context$diagnostics$modes, collapse = ","),
          status = status,
          reason = if (identical(status, "available")) NA_character_ else
            paste0("Diagnostic status: ", status, "."),
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
  }
  agreement_specs <- list(
    weighted_correlation = "correlation_status",
    weighted_cosine = "cosine_status",
    zero_intercept_gain = "gain_status"
  )
  for (diagnostic in names(agreement_specs)) {
    status_column <- agreement_specs[[diagnostic]]
    for (contrast in state$contrasts) {
      statuses <- contrast_data[[status_column]][contrast_data$contrast == contrast]
      status <- if (any(statuses == "available")) {
        "available"
      } else if (any(statuses == "degenerate_expected")) {
        "degenerate_expected"
      } else if (any(statuses == "degenerate_observed")) {
        "degenerate_observed"
      } else {
        "insufficient_samples"
      }
      rows[[index]] <- data.frame(
        diagnostic = diagnostic,
        contrast = contrast,
        estimand = NA_character_,
        mode = "exact",
        status = status,
        reason = if (identical(status, "available")) NA_character_ else
          paste0("Diagnostic status: ", status, "."),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  rows[[index]] <- data.frame(
    diagnostic = "coverage_fraction",
    contrast = NA_character_,
    estimand = NA_character_,
    mode = "descriptive",
    status = if (any(is.finite(contrast_data$coverage_fraction))) {
      "available"
    } else {
      "insufficient_samples"
    },
    reason = NA_character_,
    stringsAsFactors = FALSE
  )
  index <- index + 1L
  for (quality_name in names(state$quality_data)) {
    value <- state$quality_data[[quality_name]]
    status <- if (any(is.finite(value))) "available" else "insufficient_samples"
    rows[[index]] <- data.frame(
      diagnostic = paste0("quality:", quality_name),
      contrast = NA_character_,
      estimand = NA_character_,
      mode = "descriptive",
      status = status,
      reason = if (identical(status, "available")) NA_character_ else
        "Quality metric has no finite values.",
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
  do.call(rbind, rows)
}

.finalize_examination_group_maps <- function(state, compiled, model_context) {
  components <- .examination_output_components(compiled)
  categorical <- grep(
    "^(argmax_delta_stat|tie_count_delta_stat):",
    names(state$maps),
    value = TRUE
  )
  new_gds(
    assays = state$maps,
    space = components$space,
    subjects = "meta",
    contrasts = state$contrasts,
    row_data = components$row_data,
    metadata = list(
      examination = list(
        categorical_assays = categorical,
        non_interpolable_assays = categorical,
        responsible_subject_lookup = data.frame(
          index = seq_along(state$subjects),
          subject = state$subjects,
          stringsAsFactors = FALSE
        ),
        model_context_digest = model_context$digest
      )
    )
  )
}

#' Summarize a group examination
#'
#' @param object A `gds_examination`.
#' @param ... Unused.
#'
#' @return A compact `summary.gds_examination` object.
#' @export
summary.gds_examination <- function(object, ...) {
  structure(
    list(
      cohort = object$cohort,
      review_queue = object$subject_data[
        object$subject_data$review_status == "review",
        c("subject", "review_priority", "review_source", "review_reason"),
        drop = FALSE
      ],
      availability = object$availability,
      retained_subjects = object$config$retained_subjects
    ),
    class = "summary.gds_examination"
  )
}

#' @export
print.summary.gds_examination <- function(x, ...) {
  cat("Group examination\n")
  cat("  model: ", x$cohort$model, "\n", sep = "")
  if (!is.null(x$cohort$formula)) {
    cat("  formula: ", x$cohort$formula, "\n", sep = "")
  }
  cat(
    "  subjects: ", x$cohort$n_included,
    " included, ", x$cohort$n_excluded_covariates,
    " excluded for covariates\n",
    sep = ""
  )
  cat("  variance mode: ", x$cohort$variance_mode, "\n", sep = "")
  cat("  review cases: ", nrow(x$review_queue), "\n", sep = "")
  if (nrow(x$review_queue)) print(x$review_queue, row.names = FALSE)
  invisible(x)
}

#' @export
print.gds_examination <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
