# Authoritative reducer model context --------------------------------------

.examination_availability_statuses <- c(
  "available",
  "insufficient_samples",
  "degenerate_expected",
  "degenerate_observed",
  "nonestimable",
  "synthetic_variance",
  "unsupported_reducer"
)

.formula_from_spec <- function(formula) {
  if (is.null(formula)) return(NULL)
  if (inherits(formula, "formula")) return(formula)
  if (is.character(formula) && length(formula)) {
    return(tryCatch(
      stats::as.formula(paste(formula, collapse = " ")),
      error = function(e) {
        stop("Invalid reducer formula: ", conditionMessage(e), call. = FALSE)
      }
    ))
  }
  stop("Reducer formula must be NULL, a formula, or a character representation.", call. = FALSE)
}

.build_execution_design <- function(formula,
                                    subjects,
                                    col_data = NULL,
                                    na_action = c("fail", "exclude"),
                                    context = "reducer") {
  na_action <- match.arg(na_action)
  subjects <- as.character(subjects)
  formula <- .formula_from_spec(formula)
  if (is.null(formula)) {
    return(list(
      formula = NULL,
      subjects = subjects,
      included_subject_index = seq_along(subjects),
      excluded_subjects = .empty_excluded_subjects(),
      model_frame = NULL,
      X = NULL,
      terms = NULL,
      factor_levels = list(),
      factor_contrasts = NULL
    ))
  }

  if (is.null(col_data)) {
    data <- data.frame(subject = subjects, row.names = subjects, check.names = FALSE)
  } else {
    data <- .align_col_data_for_subjects(
      col_data,
      subjects,
      warn_extra = FALSE,
      context = context
    )
    if (!"subject" %in% names(data)) data$subject <- subjects
  }
  rownames(data) <- subjects

  model_frame <- tryCatch(
    stats::model.frame(
      formula,
      data = data,
      na.action = stats::na.pass,
      drop.unused.levels = FALSE
    ),
    error = function(e) {
      stop(
        "Could not build the design model frame from formula ",
        paste(deparse(formula), collapse = " "), ": ", conditionMessage(e),
        if (is.null(col_data)) {
          " (no col_data supplied; attach subject covariates with with_col_data())"
        } else {
          ""
        },
        call. = FALSE
      )
    }
  )
  if (nrow(model_frame) != length(subjects)) {
    stop("Design model frame did not preserve the realized subject axis.", call. = FALSE)
  }
  rownames(model_frame) <- subjects
  complete <- if (ncol(model_frame)) {
    stats::complete.cases(model_frame)
  } else {
    rep(TRUE, length(subjects))
  }
  missing_subjects <- subjects[!complete]
  if (length(missing_subjects) && identical(na_action, "fail")) {
    stop(
      "Missing model covariates for subjects: ",
      paste(missing_subjects, collapse = ", "),
      ". Set na_action = 'exclude' to exclude and record them explicitly.",
      call. = FALSE
    )
  }
  included <- which(complete)
  excluded <- if (length(missing_subjects)) {
    data.frame(
      subject = missing_subjects,
      reason = rep("missing_covariate", length(missing_subjects)),
      subject_index = which(!complete),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    .empty_excluded_subjects()
  }
  model_frame <- model_frame[included, , drop = FALSE]
  included_subjects <- subjects[included]
  rownames(model_frame) <- included_subjects

  terms <- stats::terms(model_frame)
  X <- tryCatch(
    stats::model.matrix(terms, data = model_frame),
    error = function(e) {
      stop("Could not build the design matrix: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.matrix(X) || nrow(X) != length(included_subjects)) {
    stop("Design matrix did not preserve the included subject axis.", call. = FALSE)
  }
  rownames(X) <- included_subjects
  factor_columns <- vapply(model_frame, is.factor, logical(1))
  factor_levels <- lapply(model_frame[factor_columns], levels)

  list(
    formula = formula,
    subjects = included_subjects,
    included_subject_index = as.integer(included),
    excluded_subjects = excluded,
    model_frame = model_frame,
    X = X,
    terms = terms,
    factor_levels = factor_levels,
    factor_contrasts = attr(X, "contrasts") %||% NULL
  )
}

.empty_excluded_subjects <- function() {
  data.frame(
    subject = character(),
    reason = character(),
    subject_index = integer(),
    stringsAsFactors = FALSE
  )
}

.build_reducer_model_context <- function(compiled,
                                         estimands = NULL,
                                         na_action = c("fail", "exclude")) {
  if (!inherits(compiled, "gds_examination_plan")) {
    stop("compiled must be a gds_examination_plan.", call. = FALSE)
  }
  na_action <- match.arg(na_action)
  plan <- compiled$plan
  reducer_name <- .normalize_reducer_name(compiled$reducer$method)
  reducer <- get_reducer(reducer_name)
  if (is.null(reducer)) {
    stop("Unknown reducer in examination plan: ", reducer_name, call. = FALSE)
  }
  contract <- reducer$model_contract %||% NULL
  source_subjects <- as.character(plan$source$probe$subjects)
  source_subject_index <- as.integer(compiled$axis_selection$subject)
  selected_subjects <- source_subjects[source_subject_index]
  formula <- .formula_from_spec(compiled$reducer$formula %||% NULL)

  if (!is.null(contract) && !isTRUE(contract$uses_X) && !is.null(formula)) {
    stop(
      "Reducer '", reducer_name,
      "' does not consume a design matrix; a covariate formula would be ignored.",
      call. = FALSE
    )
  }
  if (is.null(contract) && !is.null(formula)) {
    stop(
      "Reducer '", reducer_name,
      "' lacks a model contract, so its formula cannot be interpreted for group examination.",
      call. = FALSE
    )
  }

  design <- if (!is.null(contract) && isTRUE(contract$uses_X)) {
    .build_execution_design(
      formula = formula %||% stats::as.formula("~ 1"),
      subjects = selected_subjects,
      col_data = col_data(plan),
      na_action = na_action,
      context = "group examination"
    )
  } else {
    list(
      formula = formula,
      subjects = selected_subjects,
      included_subject_index = seq_along(selected_subjects),
      excluded_subjects = .empty_excluded_subjects(),
      model_frame = if (is.null(col_data(plan))) NULL else
        .align_col_data_for_subjects(
          col_data(plan), selected_subjects,
          warn_extra = FALSE, context = "group examination"
        ),
      X = NULL,
      terms = NULL,
      factor_levels = list(),
      factor_contrasts = NULL
    )
  }

  estimand_matrix <- .resolve_reducer_estimands(
    estimands,
    contract = contract,
    X = design$X
  )
  if (!is.null(design$X) && nrow(estimand_matrix)) {
    .validate_estimability(
      design$X,
      estimand_matrix,
      tolerance = sqrt(.Machine$double.eps)
    )
  }

  synthetic_variance <- .plan_has_synthetic_variance(plan)
  variance_mode <- if (synthetic_variance) {
    if (!is.null(contract) && identical(contract$synthetic_variance, "allow_effect_only")) {
      "effect_only_synthetic"
    } else {
      "synthetic"
    }
  } else if ("var" %in% (reducer$requires %||% character())) {
    "measured_first_level"
  } else {
    "effect_only"
  }

  included_source_index <- source_subject_index[design$included_subject_index]
  subject_metadata <- if (is.null(col_data(plan))) {
    data.frame(row.names = design$subjects)
  } else {
    selected_metadata <- .align_col_data_for_subjects(
      col_data(plan),
      selected_subjects,
      warn_extra = FALSE,
      context = "group examination"
    )
    selected_metadata[design$included_subject_index, , drop = FALSE]
  }
  rownames(subject_metadata) <- design$subjects
  context <- list(
    method = reducer_name,
    reducer = reducer,
    formula = design$formula,
    subjects = design$subjects,
    source_subject_index = included_source_index,
    selected_subject_index = source_subject_index,
    included_subject_index = design$included_subject_index,
    excluded_subjects = design$excluded_subjects,
    model_frame = design$model_frame,
    subject_metadata = subject_metadata,
    X = design$X,
    terms = design$terms,
    contrasts = design$factor_contrasts,
    factor_levels = design$factor_levels,
    factor_contrasts = design$factor_contrasts,
    estimand_matrix = estimand_matrix,
    variance_mode = variance_mode,
    weight_contract = contract$weight_mode %||% NULL,
    missing_covariate_policy = na_action,
    synthetic_variance = synthetic_variance,
    model_contract = contract,
    diagnostics = reducer$diagnostics %||% NULL,
    options = compiled$reducer$options %||% list(),
    reducer_node_id = compiled$reducer$node_id %||% NULL,
    source_plan_digest = compiled$source_plan_digest
  )
  context$digest <- digest::digest(
    list(
      source_plan_digest = compiled$source_plan_digest,
      source_fingerprint_digest = compiled$source_fingerprint_digest,
      source_hash = compiled$origin_source_hash %||% plan$source$hash,
      selected_subject_index = source_subject_index,
      subjects = context$subjects,
      excluded_subjects = context$excluded_subjects,
      formula = if (is.null(context$formula)) NULL else paste(deparse(context$formula), collapse = " "),
      model_frame = context$model_frame,
      subject_metadata = context$subject_metadata,
      X = context$X,
      estimand_matrix = context$estimand_matrix,
      variance_mode = context$variance_mode,
      reducer = reducer_name,
      options = context$options
    ),
    algo = "xxhash64"
  )
  class(context) <- "reducer_model_context"
  context
}

.plan_has_synthetic_variance <- function(plan) {
  probe_flag <- isTRUE(plan$source$probe$metadata$synthetic_var)
  source_flag <- isTRUE(tryCatch(
    plan$source$source$metadata$synthetic_var,
    error = function(e) FALSE
  ))
  probe_flag || source_flag
}

.resolve_reducer_estimands <- function(estimands, contract, X = NULL) {
  kind <- contract$estimands %||% "none"
  if (identical(kind, "none")) {
    if (!is.null(estimands)) {
      stop("This reducer does not define linear estimands.", call. = FALSE)
    }
    return(matrix(numeric(), nrow = 0L, ncol = 0L))
  }
  if (identical(kind, "intercept")) {
    if (is.null(estimands)) {
      return(matrix(
        1,
        nrow = 1L,
        dimnames = list("pooled_effect", "pooled_effect")
      ))
    }
    allowed <- c("pooled_effect", "(Intercept)")
    if (is.character(estimands) && length(estimands) == 1L && estimands %in% allowed) {
      return(matrix(
        1,
        nrow = 1L,
        dimnames = list(as.character(estimands), "pooled_effect")
      ))
    }
    stop("Intercept-only reducers support only the pooled_effect estimand.", call. = FALSE)
  }
  if (is.null(X)) {
    stop("Linear estimands require an execution-time design matrix.", call. = FALSE)
  }
  columns <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
  if (is.null(estimands)) {
    out <- diag(length(columns))
    dimnames(out) <- list(columns, columns)
    return(out)
  }
  if (is.character(estimands)) {
    unknown <- setdiff(estimands, columns)
    if (length(unknown)) {
      stop("Unknown estimand or design column: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- matrix(0, nrow = length(estimands), ncol = length(columns))
    rownames(out) <- estimands
    colnames(out) <- columns
    out[cbind(seq_along(estimands), match(estimands, columns))] <- 1
    return(out)
  }
  if (is.list(estimands) && !is.data.frame(estimands)) {
    if (is.null(names(estimands)) || any(!nzchar(names(estimands))) || anyDuplicated(names(estimands))) {
      stop("Estimand lists must have unique non-empty names.", call. = FALSE)
    }
    out <- matrix(0, nrow = length(estimands), ncol = length(columns))
    dimnames(out) <- list(names(estimands), columns)
    for (i in seq_along(estimands)) {
      value <- estimands[[i]]
      if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
        stop("Estimand '", names(estimands)[i], "' must be a finite numeric vector.", call. = FALSE)
      }
      if (is.null(names(value))) {
        if (length(value) != length(columns)) {
          stop("Unnamed estimands must have one value per design column.", call. = FALSE)
        }
        out[i, ] <- value
      } else {
        unknown <- setdiff(names(value), columns)
        if (length(unknown)) {
          stop("Unknown estimand coefficient: ", paste(unknown, collapse = ", "), call. = FALSE)
        }
        out[i, match(names(value), columns)] <- value
      }
    }
    return(out)
  }
  if (is.numeric(estimands)) {
    estimands <- matrix(estimands, nrow = 1L)
  }
  if (is.matrix(estimands) && is.numeric(estimands)) {
    if (ncol(estimands) != length(columns)) {
      stop("Estimand matrix must have one column per design column.", call. = FALSE)
    }
    if (anyNA(estimands) || any(!is.finite(estimands))) {
      stop("Estimand matrix must contain finite numeric values.", call. = FALSE)
    }
    colnames(estimands) <- columns
    if (is.null(rownames(estimands))) {
      rownames(estimands) <- paste0("estimand", seq_len(nrow(estimands)))
    }
    return(estimands)
  }
  stop("Unsupported estimands specification.", call. = FALSE)
}

.validate_estimability <- function(X, estimands, tolerance) {
  decomp <- svd(X, nu = 0L, nv = ncol(X))
  scale <- if (length(decomp$d)) max(decomp$d) else 0
  rank <- if (scale > 0) sum(decomp$d > tolerance * scale) else 0L
  if (rank < ncol(X)) {
    null <- decomp$v[, seq.int(rank + 1L, ncol(X)), drop = FALSE]
    residual <- abs(estimands %*% null)
    bad <- apply(residual, 1L, max) > tolerance * pmax(1, apply(abs(estimands), 1L, max))
    if (any(bad)) {
      stop(
        "Linear estimand is not estimable under the realized design: ",
        paste(rownames(estimands)[bad], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.preflight_reducer_diagnostics <- function(context) {
  if (!inherits(context, "reducer_model_context")) {
    stop("context must be a reducer_model_context.", call. = FALSE)
  }
  diagnostics <- context$diagnostics %||% NULL
  capabilities <- diagnostics$capabilities %||% "model_diagnostics"
  if (!length(capabilities)) capabilities <- "model_diagnostics"
  modes <- diagnostics$modes %||% NA_character_
  if (!length(modes)) modes <- NA_character_
  rows <- expand.grid(
    diagnostic = capabilities,
    mode = modes,
    stringsAsFactors = FALSE
  )
  implemented <- !is.null(diagnostics) && is.function(diagnostics$fun)
  status <- if (!implemented) {
    "unsupported_reducer"
  } else if (isTRUE(context$synthetic_variance) &&
             identical(context$model_contract$synthetic_variance, "forbid")) {
    "synthetic_variance"
  } else {
    "available"
  }
  rows$status <- status
  rows$reason <- switch(
    status,
    unsupported_reducer = "Reducer has no implemented diagnostic kernel.",
    synthetic_variance = "Reducer requires genuine first-level variance.",
    available = NA_character_
  )
  rows$reducer <- context$method
  list(
    model_contract = context$model_contract,
    diagnostics = diagnostics,
    availability = rows[, c("diagnostic", "mode", "reducer", "status", "reason"), drop = FALSE]
  )
}
