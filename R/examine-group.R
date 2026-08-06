# Terminal group examination -----------------------------------------------

#' Examine subjects under an intended group model
#'
#' `examine_group()` branches from subject-level data at a group reducer. It
#' reports model-conditioned predictive surprise separately from each
#' subject's exact influence on requested group estimands. It never deletes a
#' subject or rewrites the source analysis.
#'
#' When `x` already contains one reducer, the examination inherits that
#' reducer's method, formula, options, and weighting semantics. Supplying model
#' overrides in that case is an error. When no reducer is present, `method` is
#' required and formula-driven reducers default to `~ 1`.
#'
#' @section Result contract:
#' The returned object keeps spatial maps in ordinary GDS objects and keeps
#' subject, contrast, estimand, availability, sensitivity, and provenance data
#' in tables or named lists. `subject_data$review_priority` orders inspection;
#' `subject_data$review_status` changes to `"review"` only when an absolute,
#' stability-aware criterion is met. Random-effects screening is labeled
#' `"tau2_fixed_full"`; exact retained-subject refits are labeled
#' `"tau2_refit_exact"` and never alter the original review queue.
#'
#' @examples
#' subject_ids <- paste0("sub-", 1:8)
#' participant_data <- data.frame(
#'   group = factor(rep(c("control", "patient"), each = 4)),
#'   row.names = subject_ids
#' )
#' beta <- array(rnorm(12 * 8, sd = 0.2), c(12, 8, 1))
#' beta[, participant_data$group == "patient", 1] <-
#'   beta[, participant_data$group == "patient", 1] + 0.4
#' variance <- array(0.1, dim(beta))
#' subject_gds <- new_gds(
#'   assays = list(beta = beta, var = variance),
#'   space = space_sample_labels(paste0("feature-", 1:12)),
#'   subjects = subject_ids,
#'   contrasts = "task",
#'   col_data = participant_data
#' )
#' analysis <- as_plan(subject_gds) |>
#'   reduce(method = "meta:re_reg", formula = ~ group)
#' exam <- examine_group(
#'   analysis,
#'   estimands = "grouppatient",
#'   control = examination_control(retain_n = 2L)
#' )
#' summary(exam)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(exam)
#'
#' @param x A subject-level [`gds_plan`], [`gds_source`], or [`gds`], optionally
#'   containing one group reducer and a post-hoc conclusion tail.
#' @param method Reducer method when `x` has no reducer.
#' @param formula Model formula when `x` has no reducer.
#' @param options Reducer options when `x` has no reducer.
#' @param quality Optional character vector naming first-level quality columns
#'   in `col_data(x)`. Metrics are descriptive unless an explicit direction and
#'   absolute threshold are declared in `control$review$quality`.
#' @param estimands Named linear estimands or design-column names. See
#'   [examination_control()] for execution controls.
#' @param retain Optional subject IDs retained for later localization regardless
#'   of review rank.
#' @param control A validated object from [examination_control()].
#' @param na_action Missing-covariate policy. The default fails; `"exclude"`
#'   deterministically records and excludes affected subjects from this
#'   examination only.
#'
#' @return A `gds_examination` object.
#' @export
examine_group <- function(x,
                          method = NULL,
                          formula = NULL,
                          options = list(),
                          quality = NULL,
                          estimands = NULL,
                          retain = NULL,
                          control = examination_control(),
                          na_action = c("fail", "exclude")) {
  na_action <- match.arg(na_action)
  if (!inherits(control, "gds_examination_control")) {
    stop("control must be created by examination_control().", call. = FALSE)
  }
  if (!is.list(options)) stop("options must be a list.", call. = FALSE)
  if (!is.null(quality) &&
      (!is.character(quality) || anyNA(quality) || any(!nzchar(quality)))) {
    stop("quality must be NULL or a character vector of col_data columns.", call. = FALSE)
  }
  if (!is.null(retain) &&
      (!is.character(retain) || anyNA(retain) || any(!nzchar(retain)))) {
    stop("retain must be NULL or a character vector of subject IDs.", call. = FALSE)
  }
  configured_quality <- setdiff(
    names(control$review$quality),
    "coverage_fraction"
  )
  undeclared_quality <- setdiff(configured_quality, quality %||% character())
  if (length(undeclared_quality)) {
    stop(
      "Quality review criteria require the corresponding columns in `quality`: ",
      paste(undeclared_quality, collapse = ", "), ".",
      call. = FALSE
    )
  }

  compiled <- compile_examination_plan(
    x,
    method = method,
    formula = formula,
    options = options
  )
  staging <- NULL
  stage_path <- NULL
  if (!identical(compiled$scan_strategy, "direct")) {
    staging <- .stage_examination_plan(compiled, control)
    stage_path <- staging$path
    on.exit({
      if (!is.null(stage_path) && file.exists(stage_path)) {
        .cleanup_examination_stage(stage_path)
      }
    }, add = TRUE)
    compiled <- staging$compiled
  }
  model_context <- .build_reducer_model_context(
    compiled,
    estimands = estimands,
    na_action = na_action
  )
  preflight <- .preflight_reducer_diagnostics(model_context)
  diagnostics <- model_context$diagnostics
  if (is.null(diagnostics) || !is.function(diagnostics$fun)) {
    stop(
      "Reducer '", model_context$method,
      "' has no implemented group-examination diagnostic kernel. ",
      "Its model surprise and influence outputs are unavailable.",
      call. = FALSE
    )
  }
  if (isTRUE(model_context$synthetic_variance) &&
      identical(model_context$model_contract$synthetic_variance, "forbid")) {
    stop(
      "Reducer '", model_context$method,
      "' requires genuine first-level variance; synthetic unit variance cannot ",
      "be interpreted as measured precision.",
      call. = FALSE
    )
  }
  # Missing-covariate exclusion is resolved before any adapter read. The model
  # context retains both selected and included indices for provenance.
  compiled$axis_selection$subject <-
    compiled$axis_selection$subject[model_context$included_subject_index]

  all_subjects <- as.character(compiled$plan$source$probe$subjects)
  missing_retain <- setdiff(retain %||% character(), model_context$subjects)
  if (length(missing_retain)) {
    stop("Unknown retained subjects: ", paste(missing_retain, collapse = ", "), call. = FALSE)
  }
  quality_data <- .examination_quality_data(
    compiled$plan,
    model_context$subjects,
    quality
  )
  required_assays <- .examination_required_assays(
    compiled$plan,
    model_context$reducer
  )

  first_scan <- .scan_compiled_plan(
    compiled,
    assays = required_assays,
    block_size = control$block_size,
    initialize = function(scan_context) {
      .initialize_examination_accumulator(
        scan_context,
        model_context,
        control,
        quality_data,
        retain,
        preflight
      )
    },
    update = function(state, arrays, block, scan_context) {
      .update_examination_accumulator(
        state,
        arrays,
        block,
        scan_context,
        model_context,
        control
      )
    },
    finalize = function(state, receipt, scan_context) {
      list(state = state, receipt = receipt, scan_context = scan_context)
    }
  )
  first <- first_scan$value
  examination <- .finalize_examination_accumulator(
    first$state,
    first$receipt,
    first$scan_context,
    compiled,
    model_context,
    control,
    preflight
  )
  geometry_basis <- .prepare_geometry_basis(first$state, control)
  selected_subjects <- examination$subject_data$subject[
    examination$subject_data$retained
  ]
  second_scan <- .scan_compiled_plan(
    compiled,
    assays = required_assays,
    block_size = control$block_size,
    initialize = function(scan_context) {
      .initialize_geometry_pass2(
        geometry_basis,
        first$state,
        control,
        selected_subjects = selected_subjects,
        model_context = model_context
      )
    },
    update = function(state, arrays, block, scan_context) {
      .update_geometry_pass2_from_arrays(
        state, arrays, block, model_context, control
      )
    },
    finalize = function(state, receipt, scan_context) {
      list(
        embedding = .finalize_residual_geometry(state, control),
        selected = .finalize_selected_pass(
          state, compiled, model_context, control
        ),
        receipt = receipt
      )
    }
  )
  examination$embedding <- second_scan$value$embedding
  examination$subject_maps <- second_scan$value$selected$subject_maps
  exact_data <- second_scan$value$selected$exact_estimand_data
  if (!is.null(exact_data)) {
    examination$estimand_data <- rbind(examination$estimand_data, exact_data)
    examination$sensitivity <- examination$estimand_data[, c(
      "subject", "contrast", "estimand", "mode", "ranking_stage",
      "influence_energy", "max_abs_delta_stat", "abs_delta_q95_approx",
      "eligible_n", "stable"
    ), drop = FALSE]
  }
  examination$availability <- rbind(
    examination$availability,
    data.frame(
      diagnostic = "residual_geometry",
      contrast = NA_character_,
      estimand = NA_character_,
      mode = "deterministic_randomized_svd",
      status = examination$embedding$status,
      reason = if (identical(examination$embedding$status, "available")) {
        NA_character_
      } else {
        paste0("Diagnostic status: ", examination$embedding$status, ".")
      },
      stringsAsFactors = FALSE
    )
  )
  examination$conclusion <- .compute_examination_conclusions(
    examination,
    compiled,
    model_context
  )
  examination$provenance$scan_receipt <- .combine_examination_receipts(
    first$receipt,
    second_scan$value$receipt,
    staging = staging$record %||% NULL,
    retained_map_count = second_scan$value$selected$retained_map_count
  )
  examination$provenance$source_fingerprint <- compiled$source_fingerprint
  if (!is.null(staging)) {
    staging$record$cleanup_succeeded <- .cleanup_examination_stage(stage_path)
    stage_path <- NULL
    examination$provenance$staging <- staging$record
    examination$provenance$scan_receipt$staging <- staging$record
  }
  examination
}

.combine_examination_receipts <- function(pass1,
                                          pass2,
                                          staging = NULL,
                                          retained_map_count = 0L) {
  rss <- c(pass1$peak_rss_bytes, pass2$peak_rss_bytes)
  peak_rss <- if (any(is.finite(rss))) max(rss[is.finite(rss)]) else NA_real_
  list(
    pass1 = pass1,
    pass2 = pass2,
    staging = staging,
    adapter_reads = pass1$adapter_reads + pass2$adapter_reads +
      (staging$adapter_reads %||% 0L),
    bytes_read = pass1$bytes_read + pass2$bytes_read +
      (staging$bytes_read %||% 0),
    elapsed_seconds = pass1$elapsed_seconds + pass2$elapsed_seconds +
      (staging$elapsed_seconds %||% 0),
    peak_rss_bytes = peak_rss,
    stage_size_bytes = staging$stage_size_bytes %||% 0,
    retained_map_count = as.integer(retained_map_count)
  )
}

.examination_required_assays <- function(plan, reducer) {
  available <- as.character(plan$source$probe$assays %||% character())
  required <- setdiff(reducer$requires %||% character(), "X")
  requested <- intersect(required, available)
  if ("var" %in% required && !"var" %in% available && "se" %in% available) {
    requested <- c(requested, "se")
  }
  if ("se" %in% required && !"se" %in% available && "var" %in% available) {
    requested <- c(requested, "var")
  }
  derivable <- c(
    if ("var" %in% required && "se" %in% available) "var" else character(),
    if ("se" %in% required && "var" %in% available) "se" else character()
  )
  unavailable <- setdiff(required, c(requested, derivable))
  if (length(unavailable)) {
    stop(
      "Examination source cannot provide required assays for reducer '",
      reducer$name, "': ", paste(unavailable, collapse = ", "), ".",
      call. = FALSE
    )
  }
  unique(requested)
}

.examination_quality_data <- function(plan, subjects, quality) {
  data <- col_data(plan)
  if (is.null(data)) {
    if (length(quality)) {
      stop("quality columns were requested but col_data is absent.", call. = FALSE)
    }
    return(data.frame(row.names = subjects))
  }
  data <- .align_col_data_for_subjects(
    data,
    subjects,
    warn_extra = FALSE,
    context = "group examination"
  )
  missing <- setdiff(quality %||% character(), names(data))
  if (length(missing)) {
    stop("Unknown quality columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  nonnumeric <- (quality %||% character())[
    !vapply(data[, quality %||% character(), drop = FALSE], is.numeric, logical(1))
  ]
  if (length(nonnumeric)) {
    stop("Quality columns must be numeric: ", paste(nonnumeric, collapse = ", "), call. = FALSE)
  }
  data[, quality %||% character(), drop = FALSE]
}
