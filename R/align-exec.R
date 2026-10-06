# Alignment execution -----------------------------------------------------

apply_align <- function(family, arrays, subjects, space) {
  operators <- family$by_subject
  uncertainty <- family$uncertainty %||% UncertaintyRule("independent")

  if (!all(c("beta", "var") %in% names(arrays))) {
    stop("align() requires beta and var assays", call. = FALSE)
  }

  subject_names <- subjects
  if (is.null(names(operators))) names(operators) <- subject_names

  target_space <- family$to
  first_mat <- .align_matrix(operators[[subject_names[1]]])
  n_target <- nrow(first_mat)
  dims <- dim(arrays$beta)
  n_subjects <- dims[2]
  n_contrasts <- dims[3]

  had_stats <- any(c("z", "p") %in% names(arrays))
  satterthwaite <- identical(uncertainty$df_rule, "satterthwaite")

  new_beta <- array(NA_real_, dim = c(n_target, n_subjects, n_contrasts))
  new_var <- array(NA_real_, dim = c(n_target, n_subjects, n_contrasts))
  new_df <- if (!is.null(arrays$df) && satterthwaite) {
    array(NA_real_, dim = c(n_target, n_subjects, n_contrasts))
  } else {
    NULL
  }

  for (j in seq_len(n_subjects)) {
    subj_id <- subject_names[j]
    mat <- .align_matrix(operators[[subj_id]])
    if (nrow(mat) != n_target) {
      stop("align(): all subject operators must map to the same number of target samples", call. = FALSE)
    }
    .check_map_columns(mat, arrays, context = sprintf("align() [subject '%s']", subj_id))
    beta_j <- arrays$beta[, j, , drop = FALSE]
    var_j <- arrays$var[, j, , drop = FALSE]

    res <- if (uncertainty$mode == "cov_provider") {
      propagate_variance_covariance(mat, beta_j, var_j, uncertainty$cov_provider)
    } else {
      propagate_variance_independent(mat, beta_j, var_j)
    }

    new_beta[, j, ] <- res$beta[, 1, ]
    new_var[, j, ] <- res$var[, 1, ]

    if (!is.null(new_df)) {
      df_res <- aggregate_df_satterthwaite(mat, var_j, arrays$df[, j, , drop = FALSE])
      new_df[, j, ] <- df_res[, 1, ]
    }
  }

  # Without Satterthwaite aggregation, source-space df can only be carried
  # over when it is constant per subject/contrast; otherwise it is dropped.
  if (is.null(new_df) && !is.null(arrays$df)) {
    new_df <- .broadcast_constant_df(arrays$df, n_target)
  }

  out <- list(beta = new_beta, var = new_var, df = new_df)
  arrays <- .reset_after_transform(out, c("beta", "var", "df"), had_stats = had_stats)

  list(arrays = arrays, space = target_space)
}

.align_matrix <- function(entry) {
  if (is.matrix(entry)) return(entry)
  if (inherits(entry, "Matrix")) return(as.matrix(entry))
  if (is.list(entry) && !is.null(entry$matrix)) return(as.matrix(entry$matrix))
  stop("Subject map must be matrix or list(matrix = ..)", call. = FALSE)
}
