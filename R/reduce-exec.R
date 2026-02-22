# Reducer execution -------------------------------------------------------

apply_reduce <- function(node, arrays, weights, subjects, col_data = NULL) {
  name <- .normalize_reducer_name(node$method)
  reducer <- get_reducer(name)
  if (is.null(reducer)) {
    # fallback to legacy paths
    method <- node$method
    if (method %in% c("fixed", "random")) {
      arrays <- .reduce_effects(arrays, method, node$weights, node$options)
    } else {
      arrays <- .reduce_evidence(arrays, method)
    }
    return(list(arrays = arrays, subjects = "meta"))
  }

  # Ensure required inputs exist (derive z/p if necessary)
  req <- reducer$requires %||% character(0)
  arrays <- .ensure_required_arrays(arrays, req)

  beta <- arrays$beta
  var  <- arrays$var
  z    <- arrays$z
  p    <- arrays$p
  dims <- dim(beta %||% var %||% z %||% p)
  n_samples <- dims[1]; n_subject <- dims[2]; n_contrast <- dims[3]

  # Prepare outputs
  out_fields <- reducer$provides
  # Prepare base fields (scalar/vector per sample) and allow dynamic param fields later
  reg_fields <- intersect(out_fields, c("coef", "se_coef", "t_coef", "p_coef"))
  base_fields <- setdiff(out_fields, reg_fields)
  out_arrays <- lapply(base_fields, function(.) array(NA_real_, dim = c(n_samples, 1, n_contrast)))
  names(out_arrays) <- base_fields

  is_reg <- grepl("^meta:(fe|re)_reg$", reducer$name) || grepl("^ols:voxelwise$", reducer$name)
  design_info <- NULL
  attachments <- list()

  param_arrays <- list()
  # Validate and prepare reducer options once
  opts_root <- validate_reducer_options(reducer$options_schema %||% list(), node$options %||% list())
  # Build X once per reducer node if requested via options or formula
  if (is.null(opts_root$X) && !is.null(node$formula) && !is.null(col_data)) {
    f <- tryCatch(if (inherits(node$formula, "formula")) node$formula else stats::as.formula(node$formula), error = function(e) NULL)
    if (!is.null(f)) {
      # Match subject order; require rownames(col_data) or a 'subject' column
      subj_col <- NULL
      if (!is.null(rownames(col_data)) && length(rownames(col_data))) {
        subj_col <- data.frame(subject = rownames(col_data), col_data, check.names = FALSE)
      } else if ("subject" %in% names(col_data)) {
        subj_col <- col_data
      }
      if (!is.null(subj_col)) {
        idx <- match(subjects, subj_col$subject)
        Xdat <- subj_col[idx, , drop = FALSE]
        opts_root$X <- stats::model.matrix(f, data = Xdat)
      }
    }
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
    # Lancaster dfw auto-derive if missing
    opts_local <- opts_root
    if (identical(reducer$name, "combine:lancaster") && is.null(opts_local$dfw)) {
      if (!is.null(arrays$df)) {
        opts_local$dfw <- as.integer(round(as.numeric(arrays$df[1, , k])))
      } else if (!is.null(arrays$df1)) {
        opts_local$dfw <- as.integer(round(as.numeric(arrays$df1[1, , k])))
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

    res <- reducer$fun(beta_mat, var_mat, X, z_mat, p_mat, arrays$df, arrays$df1, arrays$df2, opts_local)
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

.ensure_required_arrays <- function(arrays, requires) {
  # Support se<->var auto-derivation when reducers declare needs
  need_var <- ("var" %in% requires) && is.null(arrays$var) && !is.null(arrays$se)
  if (need_var) arrays$var <- derive_var(arrays)
  need_se <- ("se" %in% requires) && is.null(arrays$se) && !is.null(arrays$var)
  if (need_se) arrays$se <- derive_se(arrays)
  need_z <- ("z" %in% requires) && is.null(arrays$z)
  need_p <- ("p" %in% requires) && is.null(arrays$p)
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

.reduce_effects <- function(arrays, method, weight_scheme, opts) {
  # Auto-derive var from se if available for legacy paths
  if (is.null(arrays$var) && !is.null(arrays$se)) {
    arrays$var <- derive_var(arrays)
  }
  if (!all(c("beta", "var") %in% names(arrays))) {
    stop("beta and var required for fixed/random reduce", call. = FALSE)
  }
  weights <- .compute_weights(weight_scheme, arrays, opts)
  beta <- arrays$beta
  var <- arrays$var
  dims <- dim(beta)

  tau2 <- NULL
  if (identical(method, "random")) {
    tau2 <- .estimate_tau2(beta, var, weights)
    weights <- 1 / (1 / weights + array(tau2, dim = dim(weights)))
  }

  beta_out <- array(NA_real_, dim = c(dims[1], 1, dims[3]))
  var_out <- array(NA_real_, dim = c(dims[1], 1, dims[3]))
  df_out <- if ("df" %in% names(arrays)) array(NA_real_, dim = c(dims[1], 1, dims[3])) else NULL

  for (i in seq_len(dims[1])) {
    for (k in seq_len(dims[3])) {
      w <- weights[i, , k]
      b <- beta[i, , k]
      v <- var[i, , k]
      w <- replace(w, !is.finite(w), 0)
      denom <- sum(w, na.rm = TRUE)
      if (denom == 0) {
        beta_out[i, 1, k] <- NA
        var_out[i, 1, k] <- NA
        if (!is.null(df_out)) df_out[i, 1, k] <- NA
        next
      }
      beta_out[i, 1, k] <- sum(w * b, na.rm = TRUE) / denom
      var_out[i, 1, k] <- 1 / denom
      if (!is.null(df_out)) {
        df_weight <- .weights_for_df(w)
        df_out[i, 1, k] <- aggregate_df_satterthwaite(diag(df_weight), array(v, c(length(v), 1, 1)), array(arrays$df[i, , k], c(length(v), 1, 1)))[1, 1, 1]
      }
    }
  }

  arrays$beta <- beta_out
  arrays$var <- var_out
  if (!is.null(df_out)) arrays$df <- df_out
  if ("se" %in% names(arrays)) arrays$se <- sqrt(var_out)
  arrays$t <- beta_out / sqrt(var_out)
  if (("df" %in% names(arrays)) || ("p" %in% names(arrays) && "beta" %in% names(arrays))) {
    arrays$z <- derive_z(arrays)
    arrays$p <- derive_p(arrays)
  }
  arrays
}

.reduce_evidence <- function(arrays, method) {
  z <- if ("z" %in% names(arrays)) arrays$z else NULL
  if (is.null(z) && "p" %in% names(arrays)) {
    z <- stats::qnorm(1 - arrays$p / 2) * sign(arrays$beta %||% 1)
  }
  if (is.null(z)) stop("Evidence reducer requires z or p", call. = FALSE)

  dims <- dim(z)
  z_out <- array(NA_real_, dim = c(dims[1], 1, dims[3]))
  p_out <- array(NA_real_, dim = c(dims[1], 1, dims[3]))
  chi_out <- array(NA_real_, dim = c(dims[1], 1, dims[3]))
  df_out <- array(NA_real_, dim = c(dims[1], 1, dims[3]))

  for (i in seq_len(dims[1])) {
    for (k in seq_len(dims[3])) {
      z_vals <- z[i, , k]
      z_vals <- z_vals[is.finite(z_vals)]
      if (!length(z_vals)) next

      if (identical(method, "stouffer")) {
        z_comb <- sum(z_vals) / sqrt(length(z_vals))
        p_comb <- 2 * stats::pnorm(-abs(z_comb))
        z_out[i, 1, k] <- z_comb
        p_out[i, 1, k] <- p_comb
      } else {
        chi <- -2 * sum(log(2 * stats::pnorm(-abs(z_vals))))
        dfv <- 2 * length(z_vals)
        chi_out[i, 1, k] <- chi
        df_out[i, 1, k] <- dfv
        p_out[i, 1, k] <- stats::pchisq(chi, dfv, lower.tail = FALSE)
        z_out[i, 1, k] <- stats::qnorm(1 - p_out[i, 1, k] / 2)
      }
    }
  }

  arrays$z <- z_out
  arrays$p <- p_out
  arrays$chi2 <- chi_out
  arrays$df <- df_out
  arrays
}

.compute_weights <- function(scheme, arrays, opts) {
  dims <- dim(arrays$beta)
  weights <- array(1, dim = dims)
  if (scheme == "equal") return(weights)
  if (scheme == "1/var") {
    weights <- 1 / arrays$var
  } else if (scheme == "n_eff" && "n_eff" %in% names(arrays)) {
    weights <- array(arrays$n_eff, dim = dims)
  } else if (scheme == "custom") {
    if (is.null(opts$custom_weights)) stop("custom weights require options$custom_weights", call. = FALSE)
    weights <- opts$custom_weights
  }
  weights
}

.estimate_tau2 <- function(beta, var, weights) {
  dims <- dim(beta)
  tau2 <- array(0, dim = c(dims[1], 1, dims[3]))
  for (i in seq_len(dims[1])) {
    for (k in seq_len(dims[3])) {
      w <- weights[i, , k]
      b <- beta[i, , k]
      v <- var[i, , k]
      w <- replace(w, !is.finite(w), 0)
      W <- sum(w)
      if (W <= 0) next
      b_bar <- sum(w * b) / W
      Q <- sum(w * (b - b_bar)^2)
      df <- length(b) - 1
      if (df <= 0) next
      C <- W - sum(w^2) / W
      tau <- max((Q - df) / C, 0)
      tau2[i, 1, k] <- tau
    }
  }
  tau2
}

.weights_for_df <- function(w) {
  w / sum(w)
}
