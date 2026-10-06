# Reducer execution -------------------------------------------------------

apply_reduce <- function(node, arrays, weights, subjects, col_data = NULL, contrast_data = NULL, contrasts = NULL) {
  name <- .normalize_reducer_name(node$method)
  reducer <- get_reducer(name)
  if (is.null(reducer)) {
    stop(.unknown_reducer_message(node$method), call. = FALSE)
  }

  # Ensure required inputs exist (derive z/p if necessary)
  req <- reducer$requires %||% character(0)
  # Consumption-site backstop for direct apply_reduce() calls that pass a tagged
  # `var` array (the array attribute does not survive pipeline slicing/masking,
  # so realised plans are instead protected by the metadata-flag guard in the
  # reduce() verb). Keyed on the attribute, not on values, so a legitimate real
  # variance that happens to equal 1 is not falsely tripped.
  if ("var" %in% req && isTRUE(attr(arrays$var, "synthetic_unit_variance"))) {
    stop(sprintf(
      paste0("Reducer '%s' is variance-weighted, but the `var` assay is a synthetic ",
             "unit-variance placeholder (beta/stat maps ingested without standard errors), ",
             "so the group standard errors would be meaningless. Use an unweighted reducer ",
             "such as method = \"ols:voxelwise\", or supply real standard errors."),
      node$method), call. = FALSE)
  }
  arrays <- .ensure_required_arrays(arrays, req)

  if (identical(reducer$input_shape %||% "contrastwise", "joint_contrast")) {
    return(.apply_joint_reducer(
      node,
      reducer,
      arrays,
      subjects,
      contrasts = contrasts,
      col_data = col_data,
      contrast_data = contrast_data
    ))
  }

  beta <- arrays$beta
  var  <- arrays$var
  z    <- arrays$z
  p    <- arrays$p
  dims <- dim(beta %||% var %||% z %||% p)
  n_samples <- dims[1]; n_subject <- dims[2]; n_contrast <- dims[3]

  # Prepare outputs
  out_fields <- reducer$provides
  # Prepare base fields (scalar/vector per sample) and allow dynamic param fields later
  reg_fields <- intersect(
    out_fields,
    c("coef", "se_coef", "z_coef", "t_coef", "p_coef")
  )
  base_fields <- setdiff(out_fields, reg_fields)
  out_arrays <- lapply(base_fields, function(.) array(NA_real_, dim = c(n_samples, 1, n_contrast)))
  names(out_arrays) <- base_fields

  is_reg <- grepl("^meta:(fe|re)_reg$", reducer$name) || grepl("^ols:voxelwise$", reducer$name)
  design_info <- NULL
  attachments <- list()

  param_arrays <- list()
  # Validate and prepare reducer options once
  opts_root <- validate_reducer_options(reducer$options_schema %||% list(), node$options %||% list())
  # A formula is authoritative and is always resolved against the realized
  # subject axis. Never reuse a plan-construction-time X after subject
  # subsetting or reordering.
  if (!is.null(node$formula)) {
    contract <- reducer$model_contract %||% NULL
    if (!is.null(contract) && !isTRUE(contract$uses_X)) {
      stop(
        "Reducer '", reducer$name,
        "' does not consume a design matrix; a formula would be ignored.",
        call. = FALSE
      )
    }
    design <- .build_execution_design(
      formula = node$formula,
      subjects = subjects,
      col_data = col_data,
      na_action = "fail",
      context = reducer$name
    )
    opts_root$X <- design$X
  }
  if (is_reg && !is.null(opts_root$X)) {
    cols <- colnames(opts_root$X)
    if (is.null(cols)) cols <- paste0("X", seq_len(ncol(opts_root$X)))
    design_info <- list(
      method = reducer$name,
      formula = node$formula %||% NULL,
      columns = cols,
      hash = digest::digest(opts_root$X)
    )
  }
  for (k in seq_len(n_contrast)) {
    # slice and transpose to [subjects x samples]
    beta_mat <- if (!is.null(beta)) .slice_subjects_samples(beta, k) else NULL
    var_mat  <- if (!is.null(var))  .slice_subjects_samples(var, k)  else NULL
    z_mat    <- if (!is.null(z))    .slice_subjects_samples(z, k)    else NULL
    p_mat    <- if (!is.null(p))    .slice_subjects_samples(p, k)    else NULL
    # Lancaster dfw auto-derive if missing. Degrees of freedom may differ per
    # sample (e.g. after map_to()/align() with Satterthwaite df), so derive a
    # [samples x subjects] df matrix and run the kernel per distinct df row.
    opts_local <- opts_root
    lancaster_dfmat <- NULL
    if (identical(reducer$name, "combine:lancaster") && is.null(opts_local$dfw)) {
      df_src <- arrays$df %||% arrays$df1
      if (!is.null(df_src)) {
        lancaster_dfmat <- matrix(
          as.integer(round(as.numeric(df_src[, , k]))),
          nrow = n_samples,
          ncol = n_subject
        )
        opts_local$dfw <- lancaster_dfmat[1L, ]
      } else {
        opts_local$dfw <- rep.int(1L, n_subject)
      }
    }
    X <- opts_local$X %||% NULL
    if (is_reg && is.null(design_info) && !is.null(X)) {
      cols <- colnames(X)
      if (is.null(cols)) cols <- paste0("X", seq_len(ncol(X)))
      design_info <- list(
        method = reducer$name,
        formula = node$formula %||% NULL,
        columns = cols,
        hash = digest::digest(X)
      )
    }

    res <- if (!is.null(lancaster_dfmat)) {
      .lancaster_by_df_rows(reducer, p_mat, lancaster_dfmat, opts_local, arrays)
    } else {
      reducer$fun(beta_mat, var_mat, X, z_mat, p_mat, arrays$df, arrays$df1, arrays$df2, opts_local)
    }
    .warn_on_reduced_effective_n(reducer$name, res$n_eff, n_subject)
    #
    # Handle regression/parametric results (matrix outputs): expand into param-suffixed assays
    if ((is.matrix(res$coef) && nrow(res$coef) >= 1) || (is.matrix(res$se_coef) && nrow(res$se_coef) >= 1)) {
      par_names <- colnames(X) %||% paste0("X", seq_len(nrow(res$coef %||% res$se_coef)))
      add_matrix <- function(M, prefix) {
        if (is.null(M)) return()
        stopifnot(is.matrix(M), nrow(M) == length(par_names))
        pa <- param_arrays
        for (pi in seq_len(nrow(M))) {
          nm <- paste0(prefix, par_names[pi])
          if (is.null(pa[[nm]])) pa[[nm]] <- array(NA_real_, dim = c(n_samples, 1, n_contrast))
          pa[[nm]][, 1, k] <- as.numeric(M[pi, ])
        }
        param_arrays <<- pa
      }
      add_matrix(res$coef,    "coef:")
      add_matrix(res$se_coef, "se_coef:")
      add_matrix(res$z_coef %||% NULL,  "z_coef:")
      add_matrix(res$t_coef %||% NULL,  "t_coef:")
      add_matrix(res$p_coef %||% NULL,  "p_coef:")
      # debug: uncomment if needed
      # warning(sprintf("emit params for %s: %s", reducer$name, paste(names(out_arrays), collapse=",")))
      if (is.null(design_info) && !is.null(X)) {
        cols <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
        design_info <- list(
          method = reducer$name,
          formula = node$formula %||% NULL,
          columns = cols,
          hash = digest::digest(X)
        )
      }
    }
    # Scalars/vectors like Q, df_res, tau2 are still handled as usual
    for (nm in intersect(names(res), base_fields)) {
      vec <- as.numeric(res[[nm]])
      if (length(vec) != n_samples) vec <- rep_len(vec, n_samples)
      out_arrays[[nm]][, 1, k] <- vec
    }
    # Optional covariance output: emit assays per parameter pair and also attach packed triangles
    if (!is.null(res$cov_tri)) {
      c_name <- as.character(k)
      key <- paste0("reduce/", reducer$name, "/contrast=", c_name)
      par_names <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
      # Promote to assays: cov:<term_i>:<term_j> with [samples x 1 x contrasts]
      t_idx <- 1L
      p_len <- length(par_names)
      for (pi in seq_len(p_len)) {
        for (pj in pi:p_len) {
          nm <- paste0("cov:", par_names[pi], ":", par_names[pj])
          if (is.null(param_arrays[[nm]])) param_arrays[[nm]] <- array(NA_real_, dim = c(n_samples, 1, n_contrast))
          # cov_tri rows are packed upper triangles per sample; each row maps to a (pi,pj) pair
          param_arrays[[nm]][, 1, k] <- as.numeric(res$cov_tri[t_idx, ])
          t_idx <- t_idx + 1L
        }
      }
      # Keep attachments for backwards compatibility and richer metadata
      attachments[[key]] <- list(type = "cov_tri", terms = par_names, pack = "upper", cov_tri = res$cov_tri)
    }
  }

  # Merge param arrays collected across contrasts
  out_arrays <- c(out_arrays, param_arrays)
  # Merge back into arrays; if param-suffixed fields exist, replace arrays entirely
  #
  has_param <- length(param_arrays) > 0
  if (has_param) {
    arrays <- out_arrays
  } else {
    arrays[names(out_arrays)] <- out_arrays
  }
  # Maintain common convenience derivations and legacy field mapping
  if (all(c("beta_g", "var_g") %in% names(out_arrays))) {
    arrays$se_g <- arrays$se_g %||% sqrt(arrays$var_g)
    arrays$z_g <- arrays$z_g %||% (arrays$beta_g / arrays$se_g)
    # Only compute two-sided p if reducer didn't provide p_g
    if (is.null(arrays$p_g)) {
      arrays$p_g <- 2 * stats::pnorm(-abs(arrays$z_g))
    }
    # Legacy: overwrite beta/var for group-level result
    arrays$beta <- arrays$beta_g
    arrays$var  <- arrays$var_g
  }
  if (all(c("z_g", "p_g") %in% names(out_arrays))) {
    arrays$z <- arrays$z_g
    arrays$p <- arrays$p_g
  }
  # After reduction, collapse to group-level: drop any leftover assays with subject dimension > 1
  arrays <- Filter(function(a) {
    is.array(a) && length(dim(a)) == 3L && dim(a)[2L] == 1L
  }, arrays)
  list(arrays = arrays, subjects = "meta", design_info = design_info, attachments = attachments)
}

.warn_on_reduced_effective_n <- function(reducer_name, n_eff, n_subject) {
  if (!(reducer_name %in% c("meta:fe", "meta:re")) || is.null(n_eff)) {
    return(invisible(NULL))
  }

  n_eff <- as.numeric(n_eff)
  reduced <- is.finite(n_eff) & n_eff < n_subject
  if (any(reduced)) {
    warning(sprintf(
      "%s used fewer than %d finite subjects for %d of %d samples; inspect `n_eff` for effective subject counts.",
      reducer_name, n_subject, sum(reduced), length(n_eff)
    ), call. = FALSE)
  }
  invisible(NULL)
}

.apply_joint_reducer <- function(node, reducer, arrays, subjects, contrasts = NULL, col_data = NULL, contrast_data = NULL) {
  opts_root <- validate_reducer_options(reducer$options_schema %||% list(), node$options %||% list())
  res <- reducer$fun(
    arrays = arrays,
    design = .build_joint_reducer_design(
      reducer = reducer,
      formula = node$formula %||% opts_root$formula %||% "~ 1",
      subjects = subjects,
      contrasts = contrasts %||% dimnames(arrays$beta)[[3L]] %||% paste0("contrast", seq_len(dim(arrays$beta)[3L])),
      col_data = col_data,
      contrast_data = contrast_data,
      opts = opts_root
    ),
    opts = opts_root
  )
  if (!is.list(res) || is.null(res$arrays)) {
    stop("Joint reducer must return a list with an `arrays` element", call. = FALSE)
  }
  res$subjects <- res$subjects %||% "meta"
  res$contrasts <- res$contrasts %||% "model"
  res
}

.ensure_required_arrays <- function(arrays, requires) {
  # Support se<->var auto-derivation when reducers declare needs
  need_var <- ("var" %in% requires) && is.null(arrays$var) && !is.null(arrays$se)
  if (need_var) arrays$var <- derive_var(arrays)
  need_se <- ("se" %in% requires) && is.null(arrays$se) && !is.null(arrays$var)
  if (need_se) arrays$se <- derive_se(arrays)
  need_z <- ("z" %in% requires) && is.null(arrays$z)
  # A per-term `p_coef:<term>` family (regression/LMM output) satisfies the "p"
  # requirement: the FDR post-hoc resolves it itself. Only derive a bare `p`
  # when neither an exact `p` nor any `p_coef:<term>` is available, so that
  # multi-term models do not error trying (and failing) to derive one `p`.
  has_p <- !is.null(arrays[["p"]]) || any(grepl("^p_coef:", names(arrays)))
  need_p <- ("p" %in% requires) && !has_p
  if (need_z) arrays$z <- derive_z(arrays)
  if (need_p) arrays$p <- derive_p(arrays)
  arrays
}

.slice_subjects_samples <- function(x, k) {
  d <- dim(x)
  mat <- x[, , k]
  # force into matrix [samples x subjects]
  mat <- matrix(mat, nrow = d[1], ncol = d[2])
  t(mat)
}

# Run the Lancaster kernel once per group of samples sharing the same
# per-subject df vector, so each sample is combined with its own df.
.lancaster_by_df_rows <- function(reducer, p_mat, dfmat, opts, arrays) {
  n_samples <- nrow(dfmat)
  keys <- apply(dfmat, 1L, paste, collapse = ",")
  groups <- split(seq_len(n_samples), factor(keys, levels = unique(keys)))
  out <- list()
  for (idx in groups) {
    opts_g <- opts
    opts_g$dfw <- dfmat[idx[1L], ]
    res_g <- reducer$fun(NULL, NULL, NULL, NULL, p_mat[, idx, drop = FALSE],
                         arrays$df, arrays$df1, arrays$df2, opts_g)
    for (nm in names(res_g)) {
      if (is.null(out[[nm]])) out[[nm]] <- rep(NA_real_, n_samples)
      out[[nm]][idx] <- as.numeric(res_g[[nm]])
    }
  }
  out
}
