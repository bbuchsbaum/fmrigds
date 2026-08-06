# Core reducer kernels (R prototypes) ----------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.p_from_z <- function(z, alternative = "two.sided") {
  if (identical(alternative, "greater")) {
    stats::pnorm(z, lower.tail = FALSE)
  } else if (identical(alternative, "less")) {
    stats::pnorm(z, lower.tail = TRUE)
  } else {
    2 * stats::pnorm(-abs(z))
  }
}

.add_meta_reg_z_p <- function(result, alternative = "two.sided") {
  if (is.null(result$coef) || is.null(result$se_coef)) return(result)
  z <- result$coef / result$se_coef
  z[!is.finite(result$se_coef) | result$se_coef <= 0] <- NA_real_
  p <- .p_from_z(as.numeric(z), alternative = alternative)
  dim(p) <- dim(z)
  dimnames(p) <- dimnames(z)
  result$z_coef <- z
  result$p_coef <- p
  result
}

.meta_effective_n <- function(beta, var) {
  colSums(is.finite(beta) & is.finite(var) & var > 0)
}

.colwise_fe <- function(beta, var, eps = 1e-12, alternative = "two.sided", min_subj = 1L) {
  ok <- is.finite(beta) & is.finite(var) & var > 0
  beta_ok <- beta
  beta_ok[!ok] <- NA_real_
  w <- 1 / pmax(var, eps)
  w[!ok] <- NA_real_
  sw <- colSums(w, na.rm = TRUE)
  wy <- colSums(w * beta_ok, na.rm = TRUE)
  mu <- wy / sw
  var_mu <- 1 / sw
  se_mu <- sqrt(var_mu)
  # Q and I2
  resid <- beta_ok - rep(mu, each = nrow(beta))
  Q <- colSums(w * resid * resid, na.rm = TRUE)
  k <- colSums(ok)
  I2 <- pmax(0, (Q - (k - 1)) / pmax(Q, .Machine$double.eps))
  z_g <- mu / se_mu
  p_g <- .p_from_z(z_g, alternative = alternative)
  out <- list(beta_g = mu, var_g = var_mu, se_g = se_mu, z_g = z_g, p_g = p_g, Q = Q, I2 = I2)
  bad <- k < as.integer(min_subj) | !is.finite(sw) | sw <= 0
  if (any(bad)) {
    out <- lapply(out, function(x) {
      x[bad] <- NA_real_
      x
    })
  }
  out
}

.colwise_tau2_dl <- function(beta, var, eps = 1e-12) {
  ok <- is.finite(beta) & is.finite(var) & var > 0
  beta_ok <- beta
  beta_ok[!ok] <- NA_real_
  w <- 1 / pmax(var, eps)
  w[!ok] <- NA_real_
  sw <- colSums(w, na.rm = TRUE)
  wy <- colSums(w * beta_ok, na.rm = TRUE)
  mu_fe <- wy / sw
  Q <- colSums(w * (beta_ok - rep(mu_fe, each = nrow(beta)))^2, na.rm = TRUE)
  k <- colSums(ok)
  C <- sw - colSums(w^2, na.rm = TRUE) / sw
  pmax(0, (Q - (k - 1)) / pmax(C, .Machine$double.eps))
}

core_meta_fe_kernel <- function(beta, var, X = NULL, df = NULL, opts = list()) {
  eps <- opts$eps %||% 1e-12
  alt <- opts$alternative %||% "two.sided"
  min_subj <- opts$min_subjects %||% 2L
  tail <- if (identical(alt, "less")) 1L else if (identical(alt, "greater")) 2L else 0L
  n_eff <- .meta_effective_n(beta, var)
  if (exists("meta_fe_cpp", mode = "function")) { # nocov start
    out <- meta_fe_cpp(beta, var, min_subj = min_subj, eps = eps, tail = tail)
  } else {
    out <- .colwise_fe(beta, var, eps = eps, alternative = alt, min_subj = min_subj)
  } # nocov end
  out$n_eff <- as.numeric(n_eff)
  out
}

core_meta_re_dl_kernel <- function(beta, var, X = NULL, df = NULL, opts = list()) {
  eps <- opts$eps %||% 1e-12
  alt <- opts$alternative %||% "two.sided"
  min_subj <- opts$min_subjects %||% 2L
  tail <- if (identical(alt, "less")) 1L else if (identical(alt, "greater")) 2L else 0L
  n_eff <- .meta_effective_n(beta, var)
  if (exists("meta_re_dl_cpp", mode = "function")) { # nocov start
    out <- meta_re_dl_cpp(beta, var, min_subj = min_subj, eps = eps, tail = tail)
  } else {
    tau2 <- .colwise_tau2_dl(beta, var, eps = eps)
    ok <- is.finite(beta) & is.finite(var) & var > 0
    beta_ok <- beta
    beta_ok[!ok] <- NA_real_
    wstar <- 1 / (pmax(var, eps) + rep(tau2, each = nrow(var)))
    wstar[!ok] <- NA_real_
    sws <- colSums(wstar, na.rm = TRUE)
    wys <- colSums(wstar * beta_ok, na.rm = TRUE)
    mu <- wys / sws
    var_mu <- 1 / sws
    se_mu <- sqrt(var_mu)
    z_g <- mu / se_mu
    p_g <- .p_from_z(z_g, alternative = alt)
    w <- 1 / pmax(var, eps)
    w[!ok] <- NA_real_
    mu_fe <- colSums(w * beta_ok, na.rm = TRUE) / colSums(w, na.rm = TRUE)
    Q <- colSums(w * (beta_ok - rep(mu_fe, each = nrow(beta)))^2, na.rm = TRUE)
    k <- colSums(ok)
    I2 <- pmax(0, (Q - (k - 1)) / pmax(Q, .Machine$double.eps))
    out <- list(beta_g = mu, var_g = var_mu, se_g = se_mu, z_g = z_g, p_g = p_g, tau2 = tau2, Q = Q, I2 = I2)
    bad <- k < as.integer(min_subj) | !is.finite(sws) | sws <= 0
    if (any(bad)) {
      out <- lapply(out, function(x) {
        x[bad] <- NA_real_
        x
      })
    }
  } # nocov end
  out$n_eff <- as.numeric(n_eff)
  out
}

.validate_meta_reg_inputs <- function(beta, var, X) {
  if (!is.matrix(beta) || !is.matrix(var) || !is.matrix(X)) {
    stop("`beta`, `var`, and `X` must be matrices.", call. = FALSE)
  }
  if (!identical(dim(beta), dim(var))) {
    stop("`beta` and `var` must have identical dimensions.", call. = FALSE)
  }
  if (nrow(X) != nrow(beta)) {
    stop("`X` must have one row per subject in `beta`.", call. = FALSE)
  }
  if (ncol(X) < 1L) {
    stop("`X` must contain at least one design column.", call. = FALSE)
  }
  invisible(NULL)
}

.meta_re_reg_dl_r <- function(beta,
                              var,
                              X,
                              min_subj = 2L,
                              eps = 1e-12) {
  .validate_meta_reg_inputs(beta, var, X)

  n_sample <- ncol(beta)
  n_coef <- ncol(X)
  required <- max(as.integer(min_subj), n_coef + 1L)
  finite_design <- rowSums(!is.finite(X)) == 0L

  coef <- matrix(NA_real_, nrow = n_coef, ncol = n_sample)
  se <- matrix(NA_real_, nrow = n_coef, ncol = n_sample)
  tau2 <- rep(NA_real_, n_sample)
  Q <- rep(NA_real_, n_sample)
  df_res <- rep(NA_real_, n_sample)

  for (b in seq_len(n_sample)) {
    y <- beta[, b]
    v <- var[, b]
    ok <- is.finite(y) & is.finite(v) & v > 0 & finite_design
    k <- sum(ok)
    if (k < required) next

    Xok <- X[ok, , drop = FALSE]
    yok <- y[ok]
    vok <- pmax(v[ok], eps)
    wok <- 1 / vok

    if (qr(Xok)$rank < n_coef) next

    G <- crossprod(Xok * sqrt(wok))
    A <- tryCatch(chol2inv(chol(G)), error = function(e) NULL)
    if (is.null(A)) next

    beta_fe <- A %*% crossprod(Xok, wok * yok)
    residual <- yok - drop(Xok %*% beta_fe)
    Qb <- sum(wok * residual^2)
    dfb <- k - n_coef

    # P = W - W X (X' W X)^-1 X' W.  The second W matters:
    # tr(P) = sum(w_i) - sum(w_i^2 x_i' A x_i).
    x_A_x <- rowSums((Xok %*% A) * Xok)
    trace_P <- sum(wok) - sum(wok^2 * x_A_x)
    if (!is.finite(trace_P) || trace_P <= 0) next

    t2 <- max(0, (Qb - dfb) / trace_P)
    wstar <- 1 / (vok + t2)
    Gstar <- crossprod(Xok * sqrt(wstar))
    Astar <- tryCatch(chol2inv(chol(Gstar)), error = function(e) NULL)
    if (is.null(Astar)) next

    coef[, b] <- drop(Astar %*% crossprod(Xok, wstar * yok))
    se[, b] <- sqrt(pmax(diag(Astar), 0))
    tau2[b] <- t2
    Q[b] <- Qb
    df_res[b] <- dfb
  }

  list(coef = coef, se_coef = se, tau2 = tau2, Q = Q, df_res = df_res)
}

core_meta_re_reg_dl_kernel <- function(beta,
                                       var,
                                       X,
                                       df = NULL,
                                       opts = list()) {
  if (is.null(X)) {
    stop("meta:re_reg requires X (subjects x p) in options$X", call. = FALSE)
  }
  .validate_meta_reg_inputs(beta, var, X)
  min_subj <- opts$min_subjects %||% 2L
  eps <- opts$eps %||% 1e-12
  if (exists("meta_re_reg_dl_cpp", mode = "function")) { # nocov start
    meta_re_reg_dl_cpp(beta, var, X, min_subj = min_subj, eps = eps)
  } else {
    .meta_re_reg_dl_r(beta, var, X, min_subj = min_subj, eps = eps)
  } # nocov end
}

.stouffer_fallback <- function(z, weights = NULL, min_subj = 1L) {
  if (!is.null(weights)) {
    if (length(weights) == 1L) weights <- rep(weights, nrow(z))
    if (length(weights) != nrow(z)) stop("weights length must equal number of subjects", call. = FALSE)
    finite_z <- is.finite(z)
    num <- colSums(z * weights, na.rm = TRUE)
    den <- sqrt(colSums((weights^2) * finite_z, na.rm = TRUE))
    k <- colSums(finite_z)
  } else {
    k <- colSums(is.finite(z))
    num <- colSums(z, na.rm = TRUE)
    den <- sqrt(k)
  }
  Z <- num / den
  bad <- k < as.integer(min_subj) | !is.finite(den) | den <= 0
  Z[bad] <- NA_real_
  list(z_g = Z, p_g = .p_from_z(Z))
}

core_stouffer_kernel <- function(beta = NULL, var = NULL, X = NULL, df = NULL, opts = list(), z) {
  w <- opts$weights %||% NULL
  if (exists("stouffer_combine_cpp", mode = "function")) { # nocov start
    return(stouffer_combine_cpp(z, weights = w, min_subj = opts$min_subjects %||% 1L))
  } # nocov end
  .stouffer_fallback(z, weights = w, min_subj = opts$min_subjects %||% 1L)
}

core_fisher_kernel <- function(beta = NULL, var = NULL, X = NULL, df = NULL, opts = list(), p) {
  if (exists("fisher_combine_cpp", mode = "function")) { # nocov start
    return(fisher_combine_cpp(p, min_subj = opts$min_subjects %||% 1L))
  } # nocov end
  k <- colSums(is.finite(p))
  X2 <- -2 * colSums(log(p), na.rm = TRUE)
  dfc <- 2 * k
  bad <- k < as.integer(opts$min_subjects %||% 1L)
  X2[bad] <- NA_real_
  dfc[bad] <- NA_real_
  list(chi2 = X2, df = dfc, p_g = stats::pchisq(X2, dfc, lower.tail = FALSE))
}

.lancaster_fallback <- function(p, dfw, min_subj = 1L) {
  if (length(dfw) != nrow(p)) stop("dfw length must equal number of subjects", call. = FALSE)
  n_cols <- ncol(p)
  X2 <- rep(NA_real_, n_cols)
  dfc <- rep(NA_real_, n_cols)
  k <- colSums(is.finite(p))
  for (j in seq_len(n_cols)) {
    ok <- is.finite(p[, j])
    if (!any(ok)) next
    X2[j] <- sum(stats::qchisq(1 - p[ok, j], df = 2 * dfw[ok]))
    dfc[j] <- 2 * sum(dfw[ok])
  }
  bad <- k < as.integer(min_subj)
  X2[bad] <- NA_real_
  dfc[bad] <- NA_real_
  list(chi2 = X2, df = dfc, p_g = stats::pchisq(X2, dfc, lower.tail = FALSE))
}

core_lancaster_kernel <- function(beta = NULL, var = NULL, X = NULL, df = NULL, opts = list(), p, dfw) {
  w <- as.integer(dfw)
  if (exists("lancaster_combine_cpp", mode = "function")) { # nocov start
    return(lancaster_combine_cpp(p, dfw = w, min_subj = opts$min_subjects %||% 1L))
  } # nocov end
  .lancaster_fallback(p, dfw = w, min_subj = opts$min_subjects %||% 1L)
}

.perm_tail_code <- function(alternative) {
  alt <- alternative %||% "two.sided"
  if (identical(alt, "less")) 1L else if (identical(alt, "greater")) 2L else 0L
}

.perm_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

.make_sign_matrix <- function(n_subject, n_perm, seed = NULL, include_observed = TRUE) {
  n_perm <- as.integer(n_perm)
  if (!is.finite(n_perm) || n_perm < 1L) stop("options$n_perm must be a positive integer", call. = FALSE)
  out <- .perm_with_seed(seed, {
    matrix(sample.int(2L, n_perm * n_subject, replace = TRUE), nrow = n_perm, ncol = n_subject)
  })
  out[out == 1L] <- -1L
  out[out == 2L] <- 1L
  if (isTRUE(include_observed)) out[1L, ] <- 1L
  storage.mode(out) <- "integer"
  out
}

.make_group_matrix <- function(group, n_perm, seed = NULL, include_observed = TRUE) {
  g <- as.integer(group)
  if (!all(g %in% c(0L, 1L))) stop("Two-sample permutation group must encode two groups as 0/1", call. = FALSE)
  n_perm <- as.integer(n_perm)
  if (!is.finite(n_perm) || n_perm < 1L) stop("options$n_perm must be a positive integer", call. = FALSE)
  out <- .perm_with_seed(seed, {
    t(replicate(n_perm, sample(g), simplify = TRUE))
  })
  if (isTRUE(include_observed)) out[1L, ] <- g
  storage.mode(out) <- "integer"
  out
}

.infer_two_sample_group <- function(X, opts = list()) {
  group <- opts$group %||% NULL
  if (!is.null(group)) {
    f <- as.factor(group)
    if (nlevels(f) != 2L) stop("options$group must contain exactly two groups", call. = FALSE)
    return(as.integer(f) - 1L)
  }
  if (is.null(X)) stop("perm:twosample requires X or options$group", call. = FALSE)
  term <- opts$term %||% NULL
  j <- if (is.null(term)) {
    ncol(X)
  } else if (is.character(term)) {
    match(term, colnames(X))
  } else {
    as.integer(term)
  }
  if (length(j) != 1L || is.na(j) || j < 1L || j > ncol(X)) {
    stop("options$term must identify one column of the model matrix", call. = FALSE)
  }
  x <- X[, j]
  ux <- sort(unique(x[is.finite(x)]))
  if (length(ux) != 2L) {
    stop("perm:twosample expects the tested design column to have exactly two finite values", call. = FALSE)
  }
  as.integer(x == ux[2L])
}

.add_parametric_p_from_t <- function(res, alternative = "two.sided") {
  tval <- as.numeric(res$t_g)
  df <- as.numeric(res$df)
  alt <- alternative %||% "two.sided"
  p <- if (identical(alt, "greater")) {
    stats::pt(tval, df = df, lower.tail = FALSE)
  } else if (identical(alt, "less")) {
    stats::pt(tval, df = df, lower.tail = TRUE)
  } else {
    2 * stats::pt(-abs(tval), df = df)
  }
  res$p_g <- p
  res
}

core_perm_onesample_kernel <- function(beta, var = NULL, X = NULL, df = NULL, opts = list()) {
  if (is.null(beta)) stop("perm:onesample requires beta", call. = FALSE)
  n_perm <- opts$n_perm %||% 5000L
  signs <- opts$signs %||% .make_sign_matrix(
    n_subject = nrow(beta),
    n_perm = n_perm,
    seed = opts$seed %||% NULL,
    include_observed = opts$include_observed %||% FALSE
  )
  tail <- .perm_tail_code(opts$alternative %||% "two.sided")
  if (exists("perm_onesample_t_cpp", mode = "function")) { # nocov start
    res <- perm_onesample_t_cpp(
      beta,
      signs,
      tail = tail,
      min_subj = opts$min_subjects %||% 2L
    )
  } else {
    stop("perm:onesample requires compiled C++ support", call. = FALSE)
  } # nocov end
  .add_parametric_p_from_t(res, opts$alternative %||% "two.sided")
}

core_perm_twosample_kernel <- function(beta, var = NULL, X = NULL, df = NULL, opts = list()) {
  if (is.null(beta)) stop("perm:twosample requires beta", call. = FALSE)
  group <- .infer_two_sample_group(X, opts)
  n_perm <- opts$n_perm %||% 5000L
  group_mat <- opts$group_mat %||% .make_group_matrix(
    group = group,
    n_perm = n_perm,
    seed = opts$seed %||% NULL,
    include_observed = opts$include_observed %||% TRUE
  )
  tail <- .perm_tail_code(opts$alternative %||% "two.sided")
  if (exists("perm_twosample_t_cpp", mode = "function")) { # nocov start
    res <- perm_twosample_t_cpp(
      beta,
      group_mat,
      tail = tail,
      welch = identical((opts$variance %||% "welch"), "welch"),
      min_group = opts$min_group %||% 2L
    )
  } else {
    stop("perm:twosample requires compiled C++ support", call. = FALSE)
  } # nocov end
  .add_parametric_p_from_t(res, opts$alternative %||% "two.sided")
}

register_core_reducers <- function() {
  register_reducer(
    name = "meta:fe",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_meta_fe_kernel(beta, var, X, df, opts),
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "n_eff"),
    model_contract = list(
      uses_X = FALSE,
      estimands = "intercept",
      weight_mode = "inverse_variance",
      missingness = "samplewise",
      synthetic_variance = "forbid",
      deletion = "closed_form"
    ),
    diagnostics = list(
      fun = .diagnose_meta_fe_block,
      capabilities = c("prediction", "surprise", "coefficient_deletion", "statistic_deletion"),
      modes = "exact"
    )
  )
  register_reducer(
    name = "meta:re",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_meta_re_dl_kernel(beta, var, X, df, opts),
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "se_g", "z_g", "p_g", "tau2", "Q", "I2", "n_eff"),
    options_schema = list(tau2 = c("DL")),
    model_contract = list(
      uses_X = FALSE,
      estimands = "intercept",
      weight_mode = "inverse_variance",
      missingness = "samplewise",
      synthetic_variance = "forbid",
      deletion = "tau2_fixed_full"
    ),
    diagnostics = list(
      fun = .diagnose_meta_re_block,
      capabilities = c("prediction", "surprise", "coefficient_deletion", "statistic_deletion"),
      modes = c("tau2_fixed_full", "tau2_refit_exact")
    )
  )
  register_reducer(
    name = "meta:fe_reg",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      if (is.null(X)) stop("meta:fe_reg requires X (subjects x p) in options$X", call. = FALSE)
      eps <- opts$eps %||% 1e-12
      S <- nrow(beta); B <- ncol(beta); pcols <- ncol(X)
      coef <- matrix(NA_real_, pcols, B)
      se   <- matrix(NA_real_, pcols, B)
      Q    <- rep(NA_real_, B)
      df_res <- rep(NA_real_, B)
      for (b in seq_len(B)) {
        y <- beta[, b]
        ok <- is.finite(y) & is.finite(var[, b]) & var[, b] > 0 &
          rowSums(!is.finite(X)) == 0L
        if (sum(ok) < (pcols + 1)) next
        Xok <- X[ok, , drop = FALSE]
        wok <- 1 / pmax(var[ok, b], eps)
        yok <- y[ok]
        G <- crossprod(Xok * sqrt(wok), Xok * sqrt(wok))
        A <- tryCatch(solve(G), error = function(e) NULL)
        if (is.null(A)) next
        bh <- A %*% crossprod(Xok * wok, yok)
        # residuals and Q
        r <- yok - as.vector(Xok %*% bh)
        Q[b] <- sum(wok * r * r)
        df_res[b] <- sum(ok) - pcols
        coef[, b] <- as.vector(bh)
        se[, b] <- sqrt(pmax(diag(A), 0))
      }
      .add_meta_reg_z_p(
        list(coef = coef, se_coef = se, Q = Q, df_res = df_res)
      )
    },
    requires = c("beta", "var", "X"),
    provides = c("coef", "se_coef", "z_coef", "p_coef", "Q", "df_res"),
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "inverse_variance",
      missingness = "samplewise",
      synthetic_variance = "forbid",
      deletion = "hat_matrix"
    ),
    diagnostics = list(
      fun = .diagnose_meta_fe_reg_block,
      capabilities = c("prediction", "surprise", "leverage", "coefficient_deletion", "statistic_deletion"),
      modes = "exact"
    )
  )
  register_reducer(
    name = "meta:re_reg",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      .add_meta_reg_z_p(core_meta_re_reg_dl_kernel(beta, var, X, df, opts))
    },
    requires = c("beta", "var", "X"),
    provides = c("coef", "se_coef", "z_coef", "p_coef", "tau2", "Q", "df_res"),
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "inverse_variance",
      missingness = "samplewise",
      synthetic_variance = "forbid",
      deletion = "tau2_fixed_full"
    ),
    diagnostics = list(
      fun = .diagnose_meta_re_reg_block,
      capabilities = c("prediction", "surprise", "leverage", "coefficient_deletion", "statistic_deletion"),
      modes = c("tau2_fixed_full", "tau2_refit_exact")
    )
  )
  register_reducer(
    name = "combine:stouffer",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_stouffer_kernel(beta, var, X, df, opts, z = z),
    requires = c("z"),
    provides = c("z_g", "p_g"),
    model_contract = list(
      uses_X = FALSE,
      estimands = "none",
      weight_mode = "evidence",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "unsupported"
    )
  )
  register_reducer(
    name = "combine:fisher",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_fisher_kernel(beta, var, X, df, opts, p = p),
    requires = c("p"),
    provides = c("p_g", "chi2", "df"),
    model_contract = list(
      uses_X = FALSE,
      estimands = "none",
      weight_mode = "evidence",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "unsupported"
    )
  )
  register_reducer(
    name = "combine:lancaster",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_lancaster_kernel(beta, var, X, df, opts, p = p, dfw = opts$dfw),
    requires = c("p"),
    provides = c("p_g", "chi2", "df"),
    model_contract = list(
      uses_X = FALSE,
      estimands = "none",
      weight_mode = "evidence",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "unsupported"
    )
  )

  register_reducer(
    name = "perm:onesample",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_perm_onesample_kernel(beta, var, X, df, opts),
    requires = c("beta"),
    provides = c("beta_g", "se_g", "t_g", "df", "p_g", "p_perm", "p_fwer"),
    options_schema = list(alternative = c("two.sided", "less", "greater")),
    model_contract = list(
      uses_X = FALSE,
      estimands = "intercept",
      weight_mode = "unweighted",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "unsupported"
    )
  )

  register_reducer(
    name = "perm:twosample",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_perm_twosample_kernel(beta, var, X, df, opts),
    requires = c("beta", "X"),
    provides = c("beta_g", "se_g", "t_g", "df", "p_g", "p_perm", "p_fwer"),
    options_schema = list(
      alternative = c("two.sided", "less", "greater"),
      variance = c("welch", "pooled")
    ),
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "unweighted",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "unsupported"
    )
  )

  # OLS across subjects per sample (voxelwise)
  register_reducer(
    name = "ols:voxelwise",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      if (is.null(X)) stop("ols:voxelwise requires X (subjects x p) in options$X or via formula", call. = FALSE)
      ret <- opts$return_cov %||% "none"
      if (!ret %in% c("none", "tri")) stop("options$return_cov must be 'none' or 'tri'", call. = FALSE)
      res <- ols_voxelwise_cpp(
        beta, X,
        return_cov_tri = identical(ret, "tri"),
        min_obs = opts$min_subjects %||% NULL
      )
      # add t and p for coefficients
      tmat <- res$coef / res$se_coef
      dfv <- matrix(res$df_res, nrow = nrow(tmat), ncol = ncol(tmat), byrow = TRUE)
      pmat <- 2 * stats::pt(-abs(tmat), df = dfv)
      res$t_coef <- tmat
      res$p_coef <- pmat
      res
    },
    requires = c("beta", "X"),
    provides = c("coef", "se_coef", "t_coef", "p_coef", "sigma2", "df_res", "n_obs", "cov_tri"),
    options_schema = list(return_cov = c("none", "tri")),
    model_contract = list(
      uses_X = TRUE,
      estimands = "linear",
      weight_mode = "unweighted",
      missingness = "samplewise",
      synthetic_variance = "allow_effect_only",
      deletion = "hat_matrix"
    ),
    diagnostics = list(
      fun = .diagnose_ols_block,
      capabilities = c("prediction", "surprise", "leverage", "coefficient_deletion", "statistic_deletion"),
      modes = "exact"
    )
  )
}

# R implementation of the voxelwise OLS kernel (no compiled override exists).
#
# Performs per-sample (voxel) listwise deletion: subjects whose effect is
# non-finite at a given sample are dropped for that sample only, and the model
# is refit on the remaining rows. Samples with fewer than `min_obs` finite
# observations (default p + 1, the minimum for a residual df) cannot be
# estimated and are returned as NA. The per-sample finite count is returned in
# `n_obs` so callers can audit effective sample size, and a single summary
# warning is emitted when any sample had non-finite subjects.
ols_voxelwise_cpp <- function(beta, X, return_cov_tri = FALSE, min_obs = NULL) {
  # beta: subjects x samples; X: subjects x p
  N <- nrow(beta); B <- ncol(beta); p <- ncol(X)
  min_obs <- if (is.null(min_obs)) p + 1L else max(as.integer(min_obs), p)
  Xt <- t(X)
  A_full <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  coef <- matrix(NA_real_, nrow = p, ncol = B)
  se   <- matrix(NA_real_, nrow = p, ncol = B)
  sigma2 <- rep(NA_real_, B)
  df_res <- rep(NA_real_, B)
  n_obs  <- rep(0L, B)
  if (return_cov_tri) {
    L <- p * (p + 1) / 2
    cov_tri <- matrix(NA_real_, nrow = L, ncol = B)
  } else {
    cov_tri <- NULL
  }
  pack_tri <- function(S, s2) {
    out <- numeric(length = p * (p + 1) / 2)
    t <- 1L
    for (i in seq_len(p)) {
      for (j in i:p) {
        out[t] <- s2 * S[i, j]
        t <- t + 1L
      }
    }
    out
  }
  for (b in seq_len(B)) {
    y <- beta[, b]
    ok <- is.finite(y)
    nb <- sum(ok)
    n_obs[b] <- nb
    if (nb < min_obs) next
    if (nb == N) {
      A <- A_full; Xb <- X; Xtb <- Xt; yb <- y
    } else {
      Xb <- X[ok, , drop = FALSE]; yb <- y[ok]; Xtb <- t(Xb)
      A <- tryCatch(solve(crossprod(Xb)), error = function(e) NULL)
    }
    if (is.null(A)) next
    bh <- A %*% (Xtb %*% yb)
    r <- as.numeric(yb - Xb %*% bh)
    dff <- max(nb - p, 1)
    s2 <- sum(r * r) / dff
    coef[, b] <- as.numeric(bh)
    se[, b] <- sqrt(pmax(diag(A) * s2, 0))
    sigma2[b] <- s2
    df_res[b] <- dff
    if (!is.null(cov_tri)) cov_tri[, b] <- pack_tri(A, s2)
  }
  n_reduced <- sum(n_obs < N)
  if (n_reduced > 0L) {
    n_dropped <- sum(n_obs < min_obs)
    warning(sprintf(
      "ols:voxelwise: %d of %d samples had non-finite subjects (per-voxel listwise deletion); %d had fewer than %d finite observations and were set to NA.",
      n_reduced, B, n_dropped, min_obs
    ), call. = FALSE)
  }
  out <- list(coef = coef, se_coef = se, sigma2 = sigma2, df_res = df_res, n_obs = as.numeric(n_obs))
  if (!is.null(cov_tri)) out$cov_tri <- cov_tri
  out
}
