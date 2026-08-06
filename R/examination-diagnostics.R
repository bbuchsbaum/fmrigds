# Exact group-examination diagnostic kernels -------------------------------

.diagnostic_empty <- function(n_subject, n_sample, n_estimand) {
  list(
    expected = matrix(NA_real_, n_subject, n_sample),
    predictive_resid = matrix(NA_real_, n_subject, n_sample),
    predictive_weight = matrix(NA_real_, n_subject, n_sample),
    leverage = matrix(NA_real_, n_subject, n_sample),
    fit_contribution = array(NA_real_, c(n_subject, n_estimand, n_sample)),
    delta_effect = array(NA_real_, c(n_subject, n_estimand, n_sample)),
    delta_stat = array(NA_real_, c(n_subject, n_estimand, n_sample)),
    deleted_stat = array(NA_real_, c(n_subject, n_estimand, n_sample)),
    surprise_eligible = matrix(FALSE, n_subject, n_sample),
    influence_eligible = array(FALSE, c(n_subject, n_estimand, n_sample))
  )
}

.diagnose_meta_fe_block <- function(fit,
                                    beta,
                                    var,
                                    X,
                                    estimands,
                                    opts,
                                    tolerance) {
  n_subject <- nrow(beta)
  n_sample <- ncol(beta)
  out <- .diagnostic_empty(n_subject, n_sample, 1L)
  full_effect <- rep(NA_real_, n_sample)
  full_se <- rep(NA_real_, n_sample)
  full_stat <- matrix(NA_real_, 1L, n_sample)
  coverage <- integer(n_sample)
  n_weight <- max_weight_fraction <- sign_agreement <- rep(NA_real_, n_sample)
  df_res <- rep(NA_real_, n_sample)

  for (b in seq_len(n_sample)) {
    y <- beta[, b]
    v <- var[, b]
    valid <- is.finite(y) & is.finite(v) & v > 0
    coverage[b] <- sum(valid)
    if (coverage[b] < 2L) next
    w <- rep(NA_real_, n_subject)
    w[valid] <- 1 / v[valid]
    W <- sum(w[valid])
    sum_wy <- sum(w[valid] * y[valid])
    if (!is.finite(W) || W <= 0) next
    mu <- sum_wy / W
    se <- sqrt(1 / W)
    z <- mu / se
    full_effect[b] <- mu
    full_se[b] <- se
    full_stat[1L, b] <- z
    df_res[b] <- coverage[b] - 1L
    n_weight[b] <- W^2 / sum(w[valid]^2)
    max_weight_fraction[b] <- max(w[valid]) / W
    direction <- sign(mu)
    if (direction != 0) {
      sign_agreement[b] <- sum(w[valid] * (sign(y[valid]) == direction)) / W
    }

    for (i in which(valid)) {
      W_minus <- W - w[i]
      if (!is.finite(W_minus) || W_minus <= tolerance$degeneracy) next
      mu_minus <- (sum_wy - w[i] * y[i]) / W_minus
      se_minus <- sqrt(1 / W_minus)
      z_minus <- mu_minus / se_minus
      pred_var <- v[i] + 1 / W_minus
      if (!is.finite(pred_var) || pred_var <= tolerance$degeneracy) next
      out$expected[i, b] <- mu_minus
      out$predictive_resid[i, b] <- (y[i] - mu_minus) / sqrt(pred_var)
      out$predictive_weight[i, b] <- 1 / pred_var
      out$leverage[i, b] <- w[i] / W
      out$fit_contribution[i, 1L, b] <- (w[i] / W) * y[i]
      out$delta_effect[i, 1L, b] <- mu - mu_minus
      out$deleted_stat[i, 1L, b] <- z_minus
      out$delta_stat[i, 1L, b] <- z - z_minus
      out$surprise_eligible[i, b] <- TRUE
      out$influence_eligible[i, 1L, b] <- TRUE
    }
  }

  maps <- list(
    coverage_n = as.numeric(coverage),
    df_res = df_res,
    model_rank = as.numeric(coverage > 0),
    rank_failure = as.numeric(coverage < 2L),
    finite_eligible = as.numeric(coverage >= 2L),
    n_weight = n_weight,
    max_weight_fraction = max_weight_fraction,
    sign_agreement = sign_agreement,
    "effect:pooled_effect" = full_effect,
    "se:pooled_effect" = full_se,
    "stat:pooled_effect" = drop(full_stat)
  )
  for (field in intersect(c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "n_eff"), names(fit))) {
    maps[[field]] <- as.numeric(fit[[field]])
  }
  out$full_effect <- matrix(full_effect, 1L, n_sample)
  out$full_se <- matrix(full_se, 1L, n_sample)
  out$full_stat <- full_stat
  out$maps <- .add_influence_group_maps(
    maps,
    out$delta_stat,
    out$influence_eligible,
    "pooled_effect",
    tolerance$degeneracy
  )
  out$statistic <- "z"
  out$mode <- "exact"
  out
}

.diagnose_meta_fe_reg_block <- function(fit,
                                        beta,
                                        var,
                                        X,
                                        estimands,
                                        opts,
                                        tolerance) {
  .diagnose_linear_block(
    fit = fit,
    beta = beta,
    var = var,
    X = X,
    estimands = estimands,
    opts = opts,
    tolerance = tolerance,
    method = "meta:fe_reg"
  )
}

.diagnose_meta_re_block <- function(fit,
                                    beta,
                                    var,
                                    X,
                                    estimands,
                                    opts,
                                    tolerance) {
  n_subject <- nrow(beta)
  n_sample <- ncol(beta)
  out <- .diagnostic_empty(n_subject, n_sample, 1L)
  full_effect <- full_se <- rep(NA_real_, n_sample)
  full_stat <- matrix(NA_real_, 1L, n_sample)
  coverage <- integer(n_sample)
  n_weight <- max_weight_fraction <- sign_agreement <- rep(NA_real_, n_sample)
  df_res <- rep(NA_real_, n_sample)

  for (b in seq_len(n_sample)) {
    y <- beta[, b]
    v <- var[, b]
    tau2 <- fit$tau2[b]
    valid <- is.finite(y) & is.finite(v) & v > 0 & is.finite(tau2) & tau2 >= 0
    coverage[b] <- sum(valid)
    if (coverage[b] < 2L) next
    w <- rep(NA_real_, n_subject)
    w[valid] <- 1 / (v[valid] + tau2)
    W <- sum(w[valid])
    sum_wy <- sum(w[valid] * y[valid])
    if (!is.finite(W) || W <= 0) next
    mu <- sum_wy / W
    se <- sqrt(1 / W)
    z <- mu / se
    full_effect[b] <- mu
    full_se[b] <- se
    full_stat[1L, b] <- z
    df_res[b] <- coverage[b] - 1L
    n_weight[b] <- W^2 / sum(w[valid]^2)
    max_weight_fraction[b] <- max(w[valid]) / W
    direction <- sign(mu)
    if (direction != 0) {
      sign_agreement[b] <- sum(w[valid] * (sign(y[valid]) == direction)) / W
    }
    for (i in which(valid)) {
      W_minus <- W - w[i]
      if (!is.finite(W_minus) || W_minus <= tolerance$degeneracy) next
      mu_minus <- (sum_wy - w[i] * y[i]) / W_minus
      z_minus <- mu_minus * sqrt(W_minus)
      pred_var <- v[i] + tau2 + 1 / W_minus
      if (!is.finite(pred_var) || pred_var <= tolerance$degeneracy) next
      out$expected[i, b] <- mu_minus
      out$predictive_resid[i, b] <- (y[i] - mu_minus) / sqrt(pred_var)
      out$predictive_weight[i, b] <- 1 / pred_var
      out$leverage[i, b] <- w[i] / W
      out$fit_contribution[i, 1L, b] <- (w[i] / W) * y[i]
      out$delta_effect[i, 1L, b] <- mu - mu_minus
      out$deleted_stat[i, 1L, b] <- z_minus
      out$delta_stat[i, 1L, b] <- z - z_minus
      out$surprise_eligible[i, b] <- TRUE
      out$influence_eligible[i, 1L, b] <- TRUE
    }
  }
  maps <- list(
    coverage_n = as.numeric(coverage),
    df_res = df_res,
    model_rank = as.numeric(coverage > 0),
    rank_failure = as.numeric(coverage < 2L),
    finite_eligible = as.numeric(coverage >= 2L),
    n_weight = n_weight,
    max_weight_fraction = max_weight_fraction,
    sign_agreement = sign_agreement,
    "effect:pooled_effect" = full_effect,
    "se:pooled_effect" = full_se,
    "stat:pooled_effect" = drop(full_stat)
  )
  for (field in intersect(
    c("beta_g", "var_g", "se_g", "z_g", "p_g", "tau2", "Q", "I2", "n_eff"),
    names(fit)
  )) {
    maps[[field]] <- as.numeric(fit[[field]])
  }
  out$full_effect <- matrix(full_effect, 1L, n_sample)
  out$full_se <- matrix(full_se, 1L, n_sample)
  out$full_stat <- full_stat
  out$maps <- .add_influence_group_maps(
    maps, out$delta_stat, out$influence_eligible,
    "pooled_effect", tolerance$degeneracy
  )
  out$statistic <- "z"
  out$mode <- "tau2_fixed_full"
  out
}

.diagnose_meta_re_reg_block <- function(fit,
                                        beta,
                                        var,
                                        X,
                                        estimands,
                                        opts,
                                        tolerance) {
  .diagnose_linear_block(
    fit = fit,
    beta = beta,
    var = var,
    X = X,
    estimands = estimands,
    opts = opts,
    tolerance = tolerance,
    method = "meta:re_reg"
  )
}

.diagnose_ols_block <- function(fit,
                                beta,
                                var,
                                X,
                                estimands,
                                opts,
                                tolerance) {
  .diagnose_linear_block(
    fit = fit,
    beta = beta,
    var = NULL,
    X = X,
    estimands = estimands,
    opts = opts,
    tolerance = tolerance,
    method = "ols:voxelwise"
  )
}

.diagnose_linear_block <- function(fit,
                                   beta,
                                   var,
                                   X,
                                   estimands,
                                   opts,
                                   tolerance,
                                   method) {
  n_subject <- nrow(beta)
  n_sample <- ncol(beta)
  n_coef <- ncol(X)
  n_estimand <- nrow(estimands)
  out <- .diagnostic_empty(n_subject, n_sample, n_estimand)
  full_effect <- matrix(NA_real_, n_estimand, n_sample)
  full_se <- matrix(NA_real_, n_estimand, n_sample)
  full_stat <- matrix(NA_real_, n_estimand, n_sample)
  coverage <- integer(n_sample)
  df_res <- model_rank <- rep(NA_real_, n_sample)
  rank_failure <- rep(1, n_sample)
  max_leverage <- rep(NA_real_, n_sample)
  finite_design <- rowSums(!is.finite(X)) == 0L
  is_ols <- identical(method, "ols:voxelwise")
  is_random <- identical(method, "meta:re_reg")
  min_obs <- as.integer(opts$min_subjects %||% (n_coef + 1L))
  min_obs <- max(min_obs, n_coef + 1L)

  for (b in seq_len(n_sample)) {
    y <- beta[, b]
    valid <- is.finite(y) & finite_design
    if (!is_ols) valid <- valid & is.finite(var[, b]) & var[, b] > 0
    coverage[b] <- sum(valid)
    if (coverage[b] < min_obs) next
    Xv <- X[valid, , drop = FALSE]
    rank <- qr(Xv, tol = tolerance$rank)$rank
    model_rank[b] <- rank
    if (rank < n_coef) next
    rank_failure[b] <- 0
    yv <- y[valid]
    tau2 <- if (is_random) fit$tau2[b] else 0
    if (is_random && (!is.finite(tau2) || tau2 < 0)) next
    wv <- if (is_ols) rep(1, length(yv)) else 1 / (var[valid, b] + tau2)
    G <- crossprod(Xv * sqrt(wv))
    A <- .diagnostic_inverse(G)
    if (is.null(A)) {
      rank_failure[b] <- 1
      next
    }
    theta <- drop(A %*% crossprod(Xv, wv * yv))
    residual <- yv - drop(Xv %*% theta)
    leverage <- wv * rowSums((Xv %*% A) * Xv)
    max_leverage[b] <- max(leverage)
    df <- length(yv) - n_coef
    df_res[b] <- df
    sse <- sum(wv * residual^2)
    scale_full <- if (is_ols) sse / df else 1
    if (!is.finite(scale_full) || scale_full < 0) next
    estimand_variance <- drop(rowSums((estimands %*% A) * estimands)) * scale_full
    valid_estimand <- is.finite(estimand_variance) & estimand_variance > tolerance$degeneracy
    psi <- drop(estimands %*% theta)
    full_effect[, b] <- psi
    full_se[valid_estimand, b] <- sqrt(estimand_variance[valid_estimand])
    full_stat[valid_estimand, b] <- psi[valid_estimand] / full_se[valid_estimand, b]

    valid_index <- which(valid)
    for (j in seq_along(valid_index)) {
      i <- valid_index[j]
      h <- leverage[j]
      one_minus_h <- 1 - h
      if (!is.finite(one_minus_h) || one_minus_h <= tolerance$leverage) next
      keep <- seq_along(valid_index) != j
      if (sum(keep) < min_obs ||
          qr(Xv[keep, , drop = FALSE], tol = tolerance$rank)$rank < n_coef) next

      ax <- drop(A %*% Xv[j, ])
      delta_theta <- ax * wv[j] * residual[j] / one_minus_h
      theta_minus <- theta - delta_theta
      A_minus <- A + (wv[j] / one_minus_h) * tcrossprod(ax)
      expected <- sum(Xv[j, ] * theta_minus)

      if (is_ols) {
        df_minus <- df - 1L
        if (df_minus < 1L) next
        sse_minus <- sse - residual[j]^2 / one_minus_h
        if (sse_minus < 0 && abs(sse_minus) <= tolerance$degeneracy * max(1, sse)) {
          sse_minus <- 0
        }
        if (!is.finite(sse_minus) || sse_minus <= tolerance$degeneracy) next
        scale_minus <- sse_minus / df_minus
        predictive_resid <- residual[j] / sqrt(scale_minus * one_minus_h)
        predictive_weight <- one_minus_h / scale_minus
      } else {
        predictive_variance <- var[i, b] + tau2 + drop(
          Xv[j, , drop = FALSE] %*% A_minus %*% Xv[j, ]
        )
        if (!is.finite(predictive_variance) ||
            predictive_variance <= tolerance$degeneracy) next
        scale_minus <- 1
        predictive_resid <- (y[i] - expected) / sqrt(predictive_variance)
        predictive_weight <- 1 / predictive_variance
      }

      psi_minus <- drop(estimands %*% theta_minus)
      variance_minus <- drop(rowSums((estimands %*% A_minus) * estimands)) * scale_minus
      estimand_ok <- is.finite(variance_minus) & variance_minus > tolerance$degeneracy &
        is.finite(full_stat[, b])
      stat_minus <- rep(NA_real_, n_estimand)
      stat_minus[estimand_ok] <- psi_minus[estimand_ok] /
        sqrt(variance_minus[estimand_ok])

      out$expected[i, b] <- expected
      out$predictive_resid[i, b] <- predictive_resid
      out$predictive_weight[i, b] <- predictive_weight
      out$leverage[i, b] <- h
      out$fit_contribution[i, , b] <- drop(
        estimands %*% (ax * wv[j] * y[i])
      )
      out$delta_effect[i, , b] <- drop(estimands %*% delta_theta)
      out$deleted_stat[i, estimand_ok, b] <- stat_minus[estimand_ok]
      out$delta_stat[i, estimand_ok, b] <- full_stat[estimand_ok, b] -
        stat_minus[estimand_ok]
      out$surprise_eligible[i, b] <- TRUE
      out$influence_eligible[i, estimand_ok, b] <- TRUE
    }
  }

  estimand_names <- rownames(estimands)
  maps <- list(
    coverage_n = as.numeric(coverage),
    df_res = df_res,
    model_rank = model_rank,
    rank_failure = rank_failure,
    finite_eligible = as.numeric(rank_failure == 0),
    max_leverage = max_leverage
  )
  for (e in seq_len(n_estimand)) {
    name <- estimand_names[e]
    maps[[paste0("effect:", name)]] <- full_effect[e, ]
    maps[[paste0("se:", name)]] <- full_se[e, ]
    maps[[paste0("stat:", name)]] <- full_stat[e, ]
  }
  if (!is.null(fit$Q)) maps$Q <- as.numeric(fit$Q)
  if (!is.null(fit$tau2)) maps$tau2 <- as.numeric(fit$tau2)
  if (!is.null(fit$df_res)) maps$reducer_df_res <- as.numeric(fit$df_res)
  if (!is.null(fit$sigma2)) maps$sigma2 <- as.numeric(fit$sigma2)
  if (!is.null(fit$n_obs)) maps$n_obs <- as.numeric(fit$n_obs)
  out$full_effect <- full_effect
  out$full_se <- full_se
  out$full_stat <- full_stat
  out$maps <- .add_influence_group_maps(
    maps,
    out$delta_stat,
    out$influence_eligible,
    estimand_names,
    tolerance$degeneracy
  )
  out$statistic <- if (is_ols) "t" else "z"
  out$mode <- if (is_random) "tau2_fixed_full" else "exact"
  out
}

.diagnostic_inverse <- function(G) {
  tryCatch(chol2inv(chol(G)), error = function(e) NULL)
}

.add_influence_group_maps <- function(maps,
                                      delta_stat,
                                      eligible,
                                      estimand_names,
                                      tolerance) {
  n_subject <- dim(delta_stat)[1L]
  n_sample <- dim(delta_stat)[3L]
  for (e in seq_along(estimand_names)) {
    maximum <- rep(NA_real_, n_sample)
    argmax <- rep(NA_integer_, n_sample)
    tie_count <- rep(NA_integer_, n_sample)
    eligible_n <- integer(n_sample)
    for (b in seq_len(n_sample)) {
      ok <- eligible[, e, b] & is.finite(delta_stat[, e, b])
      eligible_n[b] <- sum(ok)
      if (!any(ok)) next
      values <- abs(delta_stat[, e, b])
      maximum[b] <- max(values[ok])
      tied <- which(ok & abs(values - maximum[b]) <= tolerance * max(1, maximum[b]))
      argmax[b] <- tied[1L]
      tie_count[b] <- length(tied)
    }
    suffix <- estimand_names[e]
    maps[[paste0("max_abs_delta_stat:", suffix)]] <- maximum
    maps[[paste0("argmax_delta_stat:", suffix)]] <- argmax
    maps[[paste0("tie_count_delta_stat:", suffix)]] <- tie_count
    maps[[paste0("eligible_delta_stat:", suffix)]] <- as.numeric(eligible_n)
  }
  maps
}
