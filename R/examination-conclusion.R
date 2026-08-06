# Selected case-deletion conclusions --------------------------------------

.examination_statistic_arrays <- function(examination,
                                          estimand,
                                          model_context,
                                          deleted_subject = NULL) {
  statistic_name <- paste0("stat:", estimand)
  if (!statistic_name %in% names(examination$group_maps$assays)) {
    stop("Full statistic map is unavailable for estimand ", estimand, ".", call. = FALSE)
  }
  statistic <- examination$group_maps$assays[[statistic_name]]
  mode <- "full"
  if (!is.null(deleted_subject)) {
    subject_index <- match(deleted_subject, examination$subject_maps$subjects)
    if (is.na(subject_index)) {
      stop("Deleted conclusion requested for a subject without retained maps.", call. = FALSE)
    }
    exact_name <- paste0("delta_stat_exact:", estimand)
    analytic_name <- paste0("delta_stat:", estimand)
    delta_name <- if (exact_name %in% names(examination$subject_maps$assays)) {
      exact_name
    } else if (analytic_name %in% names(examination$subject_maps$assays)) {
      analytic_name
    } else {
      stop("Deletion statistic map is unavailable for estimand ", estimand, ".", call. = FALSE)
    }
    delta <- examination$subject_maps$assays[[delta_name]][
      , subject_index, , drop = FALSE
    ]
    statistic <- statistic - delta
    mode <- examination$subject_maps$metadata$examination$assay_modes[[delta_name]] %||%
      if (grepl("_exact:", delta_name, fixed = TRUE)) "tau2_refit_exact" else "exact"
  }
  arrays <- list()
  alternative <- model_context$options$alternative %||% "two.sided"
  if (identical(model_context$method, "ols:voxelwise")) {
    arrays$t <- statistic
    df <- examination$group_maps$assays$df_res
    if (!is.null(deleted_subject)) df <- df - 1
    arrays$df <- df
    arrays$p <- .p_from_t(statistic, df, alternative)
  } else {
    arrays$z <- statistic
    arrays$p <- .p_from_z(statistic, alternative)
  }
  arrays[[paste0("p_coef:", estimand)]] <- arrays$p
  list(arrays = arrays, mode = mode, statistic = statistic)
}

.p_from_t <- function(statistic, df, alternative = "two.sided") {
  if (identical(alternative, "greater")) {
    stats::pt(statistic, df = df, lower.tail = FALSE)
  } else if (identical(alternative, "less")) {
    stats::pt(statistic, df = df, lower.tail = TRUE)
  } else {
    2 * stats::pt(-abs(statistic), df = df)
  }
}

.conclusion_source_is_applicable <- function(node, estimand) {
  source <- node$options$source %||% NULL
  is.null(source) || identical(source, "p") ||
    identical(source, paste0("p_coef:", estimand))
}

.apply_examination_posthoc <- function(node, arrays, examination) {
  apply_posthoc(
    node,
    arrays,
    context = list(
      space = space(examination$group_maps),
      subjects = "meta",
      contrasts = contrasts(examination$group_maps),
      row_data = row_data(examination$group_maps)
    )
  )$arrays
}

.extract_examination_conclusion <- function(arrays, node, estimand, alpha) {
  q_name <- if (!is.null(arrays[[paste0("q_coef:", estimand)]])) {
    paste0("q_coef:", estimand)
  } else if (!is.null(arrays$q)) {
    "q"
  } else {
    NA_character_
  }
  adjusted <- if (!is.na(q_name)) arrays[[q_name]] else NULL
  significant <- if (!is.null(arrays$sig_mask)) {
    out <- array(NA, dim(arrays$sig_mask))
    ok <- is.finite(arrays$sig_mask)
    out[ok] <- arrays$sig_mask[ok] != 0
    out
  } else if (!is.null(adjusted)) {
    out <- array(NA, dim(adjusted))
    ok <- is.finite(adjusted)
    out[ok] <- adjusted[ok] <= alpha
    out
  } else {
    NULL
  }
  source_assay <- node$options$source %||%
    if (!is.null(arrays$p)) "p" else paste0("p_coef:", estimand)
  list(
    adjusted = adjusted,
    significant = significant,
    source_assay = source_assay,
    result_assay = q_name
  )
}

.conclusion_map_name <- function(kind, method, estimand) {
  method <- gsub("[^A-Za-z0-9_.-]+", "_", method)
  paste(kind, method, estimand, sep = ":")
}

.assign_conclusion_map <- function(store,
                                   name,
                                   values,
                                   subject_index,
                                   dimensions,
                                   subjects,
                                   contrasts) {
  if (is.null(values)) return(store)
  if (is.null(store[[name]])) {
    store[[name]] <- array(
      NA_real_,
      dimensions,
      dimnames = list(NULL, subjects, contrasts)
    )
  }
  store[[name]][, subject_index, ] <- values[, 1L, ]
  store
}

.conclusion_summary_rows <- function(full,
                                     deleted,
                                     subject,
                                     contrasts,
                                     estimand,
                                     method,
                                     alpha,
                                     mode,
                                     seed,
                                     status = "available",
                                     reason = NA_character_) {
  if (!identical(status, "available") || is.null(full$significant) ||
      is.null(deleted$significant)) {
    return(data.frame(
      subject = subject,
      contrast = NA_character_,
      estimand = estimand,
      method = method,
      source_assay = full$source_assay %||% NA_character_,
      result_assay = full$result_assay %||% NA_character_,
      alpha = alpha,
      full_significant_n = NA_integer_,
      deleted_significant_n = NA_integer_,
      transition_count = NA_integer_,
      gained_n = NA_integer_,
      lost_n = NA_integer_,
      mode = mode,
      seed = seed,
      status = status,
      reason = reason,
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(seq_along(contrasts), function(k) {
    full_mask <- full$significant[, 1L, k]
    deleted_mask <- deleted$significant[, 1L, k]
    eligible <- !is.na(full_mask) & !is.na(deleted_mask)
    data.frame(
      subject = subject,
      contrast = contrasts[k],
      estimand = estimand,
      method = method,
      source_assay = full$source_assay %||% deleted$source_assay %||% NA_character_,
      result_assay = full$result_assay %||% deleted$result_assay %||% NA_character_,
      alpha = alpha,
      full_significant_n = sum(full_mask, na.rm = TRUE),
      deleted_significant_n = sum(deleted_mask, na.rm = TRUE),
      transition_count = sum(full_mask[eligible] != deleted_mask[eligible]),
      gained_n = sum(!full_mask[eligible] & deleted_mask[eligible]),
      lost_n = sum(full_mask[eligible] & !deleted_mask[eligible]),
      mode = mode,
      seed = seed,
      status = status,
      reason = reason,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.unavailable_conclusion_rows <- function(examination,
                                         node,
                                         contract,
                                         status,
                                         reason) {
  selected <- examination$config$retained_subjects %||% NA_character_
  estimands <- examination$cohort$estimands %||% NA_character_
  alpha <- as.numeric(node$options$alpha %||% 0.05)
  rows <- lapply(selected, function(subject) {
    lapply(estimands, function(estimand) {
      .conclusion_summary_rows(
        full = list(), deleted = list(), subject = subject,
        contrasts = examination$group_maps$contrasts,
        estimand = estimand, method = node$method, alpha = alpha,
        mode = contract$mode, seed = node$options$seed %||% NA_real_,
        status = status, reason = reason
      )
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

.compute_examination_conclusions <- function(examination,
                                             compiled,
                                             model_context) {
  nodes <- compiled$conclusion_tail
  availability <- examination$conclusion$availability
  if (!length(nodes)) {
    return(list(
      availability = availability,
      results = NULL,
      full_maps = NULL,
      deleted_maps = NULL
    ))
  }
  selected <- examination$config$retained_subjects %||% character()
  n_sample <- dim(examination$group_maps$assays[[1L]])[1L]
  n_contrast <- length(examination$group_maps$contrasts)
  full_maps <- list()
  deleted_maps <- list()
  summary_rows <- list()
  row_index <- 1L
  full_dimensions <- c(n_sample, 1L, n_contrast)
  deleted_dimensions <- c(n_sample, length(selected), n_contrast)

  for (node in nodes) {
    posthoc <- get_posthoc(node$method)
    contract <- posthoc$case_deletion %||% .normalize_posthoc_case_deletion()
    if (!isTRUE(contract$supported)) {
      summary_rows[[row_index]] <- .unavailable_conclusion_rows(
        examination, node, contract, "unsupported_method",
        "The post-hoc method does not declare selected case-deletion support."
      )
      row_index <- row_index + 1L
      next
    }
    seed <- node$options$seed %||% NA_real_
    if (!isTRUE(contract$deterministic) && !is.finite(seed)) {
      summary_rows[[row_index]] <- .unavailable_conclusion_rows(
        examination, node, contract, "seed_required",
        paste0("Selected recomputation requires an explicit ", contract$seed_contract, ".")
      )
      row_index <- row_index + 1L
      next
    }
    if (!length(selected) || is.null(examination$subject_maps)) {
      summary_rows[[row_index]] <- .unavailable_conclusion_rows(
        examination, node, contract, "no_selected_subjects",
        "No retained subject maps are available for selected recomputation."
      )
      row_index <- row_index + 1L
      next
    }

    for (estimand in examination$cohort$estimands) {
      if (!.conclusion_source_is_applicable(node, estimand)) next
      alpha <- as.numeric(node$options$alpha %||% 0.05)
      full_input <- tryCatch(
        .examination_statistic_arrays(
          examination, estimand, model_context, deleted_subject = NULL
        ),
        error = identity
      )
      if (inherits(full_input, "error")) {
        summary_rows[[row_index]] <- .unavailable_conclusion_rows(
          examination, node, contract, "nonestimable",
          conditionMessage(full_input)
        )
        row_index <- row_index + 1L
        next
      }
      full_result <- tryCatch(
        .apply_examination_posthoc(node, full_input$arrays, examination),
        error = identity
      )
      if (inherits(full_result, "error")) {
        summary_rows[[row_index]] <- .unavailable_conclusion_rows(
          examination, node, contract, "recompute_error",
          conditionMessage(full_result)
        )
        row_index <- row_index + 1L
        next
      }
      full <- .extract_examination_conclusion(full_result, node, estimand, alpha)
      adjusted_name <- .conclusion_map_name("adjusted_p", node$method, estimand)
      significant_name <- .conclusion_map_name("significant", node$method, estimand)
      full_maps <- .assign_conclusion_map(
        full_maps, adjusted_name, full$adjusted, 1L,
        full_dimensions, "meta", examination$group_maps$contrasts
      )
      full_maps <- .assign_conclusion_map(
        full_maps, significant_name,
        if (is.null(full$significant)) NULL else 1 * full$significant,
        1L, full_dimensions, "meta", examination$group_maps$contrasts
      )

      for (j in seq_along(selected)) {
        deleted_input <- tryCatch(
          .examination_statistic_arrays(
            examination, estimand, model_context,
            deleted_subject = selected[j]
          ),
          error = identity
        )
        if (inherits(deleted_input, "error")) {
          summary_rows[[row_index]] <- .conclusion_summary_rows(
            full, list(), selected[j], examination$group_maps$contrasts,
            estimand, node$method, alpha, "unavailable", seed,
            status = "nonestimable", reason = conditionMessage(deleted_input)
          )
          row_index <- row_index + 1L
          next
        }
        deleted_result <- tryCatch(
          .apply_examination_posthoc(node, deleted_input$arrays, examination),
          error = identity
        )
        if (inherits(deleted_result, "error")) {
          summary_rows[[row_index]] <- .conclusion_summary_rows(
            full, list(), selected[j], examination$group_maps$contrasts,
            estimand, node$method, alpha, deleted_input$mode, seed,
            status = "recompute_error", reason = conditionMessage(deleted_result)
          )
          row_index <- row_index + 1L
          next
        }
        deleted <- .extract_examination_conclusion(
          deleted_result, node, estimand, alpha
        )
        deleted_maps <- .assign_conclusion_map(
          deleted_maps, adjusted_name, deleted$adjusted, j,
          deleted_dimensions, selected, examination$group_maps$contrasts
        )
        deleted_maps <- .assign_conclusion_map(
          deleted_maps, significant_name,
          if (is.null(deleted$significant)) NULL else 1 * deleted$significant,
          j, deleted_dimensions, selected, examination$group_maps$contrasts
        )
        summary_rows[[row_index]] <- .conclusion_summary_rows(
          full, deleted, selected[j], examination$group_maps$contrasts,
          estimand, node$method, alpha, deleted_input$mode, seed
        )
        row_index <- row_index + 1L
      }
    }
  }

  components <- list(
    space = space(examination$group_maps),
    row_data = row_data(examination$group_maps)
  )
  categorical_full <- grep("^significant:", names(full_maps), value = TRUE)
  categorical_deleted <- grep("^significant:", names(deleted_maps), value = TRUE)
  list(
    availability = availability,
    results = if (length(summary_rows)) do.call(rbind, summary_rows) else NULL,
    full_maps = if (length(full_maps)) {
      new_gds(
        full_maps, components$space, "meta", examination$group_maps$contrasts,
        row_data = components$row_data,
        metadata = list(examination = list(
          categorical_assays = categorical_full,
          non_interpolable_assays = categorical_full,
          scope = "full_group_conclusion"
        ))
      )
    } else NULL,
    deleted_maps = if (length(deleted_maps)) {
      new_gds(
        deleted_maps, components$space, selected, examination$group_maps$contrasts,
        row_data = components$row_data,
        metadata = list(examination = list(
          categorical_assays = categorical_deleted,
          non_interpolable_assays = categorical_deleted,
          scope = "selected_case_deletion_conclusion"
        ))
      )
    } else NULL
  )
}
