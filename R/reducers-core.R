# Core reducer kernels (R prototypes) ----------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.colwise_fe <- function(beta, var, eps = 1e-12) {
  w <- 1 / pmax(var, eps)
  sw <- colSums(w, na.rm = TRUE)
  wy <- colSums(w * beta, na.rm = TRUE)
  mu <- wy / sw
  var_mu <- 1 / sw
  se_mu <- sqrt(var_mu)
  # Q and I2
  resid <- beta - rep(mu, each = nrow(beta))
  Q <- colSums(w * resid * resid, na.rm = TRUE)
  k <- colSums(is.finite(beta) & is.finite(var))
  I2 <- pmax(0, (Q - (k - 1)) / pmax(Q, .Machine$double.eps))
  z_g <- mu / se_mu
  p_g <- 2 * pnorm(-abs(z_g))
  list(beta_g = mu, var_g = var_mu, se_g = se_mu, z_g = z_g, p_g = p_g, Q = Q, I2 = I2)
}

.colwise_tau2_dl <- function(beta, var, eps = 1e-12) {
  w <- 1 / pmax(var, eps)
  sw <- colSums(w, na.rm = TRUE)
  wy <- colSums(w * beta, na.rm = TRUE)
  mu_fe <- wy / sw
  Q <- colSums(w * (beta - rep(mu_fe, each = nrow(beta)))^2, na.rm = TRUE)
  k <- colSums(is.finite(beta) & is.finite(var))
  C <- sw - colSums(w^2, na.rm = TRUE) / sw
  pmax(0, (Q - (k - 1)) / pmax(C, .Machine$double.eps))
}

core_meta_fe_kernel <- function(beta, var, X = NULL, df = NULL, opts = list()) {
  eps <- opts$eps %||% 1e-12
  alt <- opts$alternative %||% "two.sided"
  tail <- if (identical(alt, "less")) 1L else if (identical(alt, "greater")) 2L else 0L
  if (exists("meta_fe_cpp", mode = "function")) {
    return(meta_fe_cpp(beta, var, min_subj = opts$min_subjects %||% 2L, eps = eps, tail = tail))
  }
  .colwise_fe(beta, var, eps = eps)
}

core_meta_re_dl_kernel <- function(beta, var, X = NULL, df = NULL, opts = list()) {
  eps <- opts$eps %||% 1e-12
  alt <- opts$alternative %||% "two.sided"
  tail <- if (identical(alt, "less")) 1L else if (identical(alt, "greater")) 2L else 0L
  if (exists("meta_re_dl_cpp", mode = "function")) {
    return(meta_re_dl_cpp(beta, var, min_subj = opts$min_subjects %||% 2L, eps = eps, tail = tail))
  }
  tau2 <- .colwise_tau2_dl(beta, var, eps = eps)
  wstar <- 1 / (pmax(var, eps) + rep(tau2, each = nrow(var)))
  sws <- colSums(wstar, na.rm = TRUE)
  wys <- colSums(wstar * beta, na.rm = TRUE)
  mu <- wys / sws
  var_mu <- 1 / sws
  se_mu <- sqrt(var_mu)
  z_g <- mu / se_mu
  p_g <- 2 * pnorm(-abs(z_g))
  w <- 1 / pmax(var, eps)
  mu_fe <- colSums(w * beta, na.rm = TRUE) / colSums(w, na.rm = TRUE)
  Q <- colSums(w * (beta - rep(mu_fe, each = nrow(beta)))^2, na.rm = TRUE)
  k <- colSums(is.finite(beta) & is.finite(var))
  I2 <- pmax(0, (Q - (k - 1)) / pmax(Q, .Machine$double.eps))
  list(beta_g = mu, var_g = var_mu, se_g = se_mu, z_g = z_g, p_g = p_g, tau2 = tau2, Q = Q, I2 = I2)
}

core_stouffer_kernel <- function(beta = NULL, var = NULL, X = NULL, df = NULL, opts = list(), z) {
  w <- opts$weights %||% NULL
  if (exists("stouffer_combine_cpp", mode = "function")) {
    return(stouffer_combine_cpp(z, weights = w, min_subj = opts$min_subjects %||% 1L))
  }
  if (!is.null(w)) {
    if (length(w) == 1L) w <- rep(w, nrow(z))
    num <- colSums(z * w, na.rm = TRUE)
    den <- sqrt(sum(w^2, na.rm = TRUE))
  } else {
    num <- colSums(z, na.rm = TRUE)
    den <- sqrt(colSums(is.finite(z)))
  }
  Z <- num / den
  list(z_g = Z, p_g = 2 * pnorm(-abs(Z)))
}

core_fisher_kernel <- function(beta = NULL, var = NULL, X = NULL, df = NULL, opts = list(), p) {
  if (exists("fisher_combine_cpp", mode = "function")) {
    return(fisher_combine_cpp(p, min_subj = opts$min_subjects %||% 1L))
  }
  X2 <- -2 * colSums(log(p), na.rm = TRUE)
  dfc <- 2 * colSums(is.finite(p))
  list(chi2 = X2, df = dfc, p_g = pchisq(X2, dfc, lower.tail = FALSE))
}

core_lancaster_kernel <- function(beta = NULL, var = NULL, X = NULL, df = NULL, opts = list(), p, dfw) {
  w <- as.integer(dfw)
  if (exists("lancaster_combine_cpp", mode = "function")) {
    return(lancaster_combine_cpp(p, dfw = w, min_subj = opts$min_subjects %||% 1L))
  }
  if (length(w) != nrow(p)) stop("dfw length must equal number of subjects", call. = FALSE)
  chi_parts <- apply(p, 2, function(col) mapply(function(pi, wi) stats::qchisq(1 - pi, df = 2 * wi), col, w))
  X2 <- colSums(chi_parts, na.rm = TRUE)
  dfc <- 2 * sum(w)
  list(chi2 = X2, df = rep(dfc, length(X2)), p_g = pchisq(X2, dfc, lower.tail = FALSE))
}

register_core_reducers <- function() {
  register_reducer(
    name = "meta:fe",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_meta_fe_kernel(beta, var, X, df, opts),
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2")
  )
  register_reducer(
    name = "meta:re",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_meta_re_dl_kernel(beta, var, X, df, opts),
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "se_g", "z_g", "p_g", "tau2", "Q", "I2"),
    options_schema = list(tau2 = c("DL"))
  )
  register_reducer(
    name = "meta:fe_reg",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      if (is.null(X)) stop("meta:fe_reg requires X (subjects x p) in options$X", call. = FALSE)
      eps <- opts$eps %||% 1e-12
      S <- nrow(beta); B <- ncol(beta); pcols <- ncol(X)
      coef <- matrix(NA_real_, pcols, B)
      se   <- matrix(NA_real_, pcols, B)
      Q    <- numeric(B)
      df_res <- numeric(B)
      for (b in seq_len(B)) {
        w <- 1 / pmax(var[, b], eps); y <- beta[, b]
        ok <- is.finite(y) & is.finite(w)
        if (sum(ok) < (pcols + 1)) next
        Xok <- X[ok, , drop = FALSE]; wok <- w[ok]; yok <- y[ok]
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
      list(coef = coef, se_coef = se, Q = Q, df_res = df_res)
    },
    requires = c("beta", "var", "X"),
    provides = c("coef", "se_coef", "Q", "df_res")
  )
  register_reducer(
    name = "meta:re_reg",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      if (is.null(X)) stop("meta:re_reg requires X (subjects x p) in options$X", call. = FALSE)
      eps <- opts$eps %||% 1e-12
      S <- nrow(beta); B <- ncol(beta); pcols <- ncol(X)
      coef <- matrix(NA_real_, pcols, B)
      se   <- matrix(NA_real_, pcols, B)
      tau2 <- numeric(B); Q <- numeric(B); df_res <- numeric(B)
      for (b in seq_len(B)) {
        w <- 1 / pmax(var[, b], eps); y <- beta[, b]
        ok <- is.finite(y) & is.finite(w)
        if (sum(ok) < (pcols + 1)) next
        Xok <- X[ok, , drop = FALSE]; wok <- w[ok]; yok <- y[ok]
        # FE step
        G <- crossprod(Xok * sqrt(wok), Xok * sqrt(wok))
        A <- tryCatch(solve(G), error = function(e) NULL)
        if (is.null(A)) next
        bh_fe <- A %*% crossprod(Xok * wok, yok)
        r <- yok - as.vector(Xok %*% bh_fe)
        Qb <- sum(wok * r * r)
        # DL tau2 (regression): C = sum w - tr(H), df = n - p
        trH <- sum(wok * rowSums((Xok %*% A) * Xok))
        C <- sum(wok) - trH
        dfb <- sum(ok) - pcols
        t2 <- max((Qb - dfb) / max(C, .Machine$double.eps), 0)
        # RE step
        wstar <- 1 / (1 / wok + t2)
        Gs <- crossprod(Xok * sqrt(wstar), Xok * sqrt(wstar))
        As <- tryCatch(solve(Gs), error = function(e) NULL)
        if (is.null(As)) next
        bh <- As %*% crossprod(Xok * wstar, yok)
        coef[, b] <- as.vector(bh)
        se[, b] <- sqrt(pmax(diag(As), 0))
        tau2[b] <- t2; Q[b] <- Qb; df_res[b] <- dfb
      }
      list(coef = coef, se_coef = se, tau2 = tau2, Q = Q, df_res = df_res)
    },
    requires = c("beta", "var", "X"),
    provides = c("coef", "se_coef", "tau2", "Q", "df_res")
  )
  register_reducer(
    name = "combine:stouffer",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_stouffer_kernel(beta, var, X, df, opts, z = z),
    requires = c("z"),
    provides = c("z_g", "p_g")
  )
  register_reducer(
    name = "combine:fisher",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_fisher_kernel(beta, var, X, df, opts, p = p),
    requires = c("p"),
    provides = c("p_g", "chi2", "df")
  )
  register_reducer(
    name = "combine:lancaster",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) core_lancaster_kernel(beta, var, X, df, opts, p = p, dfw = opts$dfw),
    requires = c("p"),
    provides = c("p_g", "chi2", "df")
  )

  # OLS across subjects per sample (voxelwise)
  register_reducer(
    name = "ols:voxelwise",
    fun = function(beta, var, X, z, p, df, df1, df2, opts) {
      if (is.null(X)) stop("ols:voxelwise requires X (subjects x p) in options$X or via formula", call. = FALSE)
      ret <- opts$return_cov %||% "none"
      if (!ret %in% c("none", "tri")) stop("options$return_cov must be 'none' or 'tri'", call. = FALSE)
      res <- ols_voxelwise_cpp(beta, X, return_cov_tri = identical(ret, "tri"))
      # add t and p for coefficients
      tmat <- res$coef / res$se_coef
      dfv <- matrix(res$df_res, nrow = nrow(tmat), ncol = ncol(tmat), byrow = TRUE)
      pmat <- 2 * stats::pt(-abs(tmat), df = dfv)
      res$t_coef <- tmat
      res$p_coef <- pmat
      res
    },
    requires = c("beta", "X"),
    provides = c("coef", "se_coef", "t_coef", "p_coef", "sigma2", "df_res", "cov_tri"),
    options_schema = list(return_cov = c("none", "tri"))
  )
}

# R fallback implementation for ols_voxelwise_cpp (overridden by Rcpp if present)
ols_voxelwise_cpp <- function(beta, X, return_cov_tri = FALSE) {
  # beta: subjects x samples; X: subjects x p
  N <- nrow(beta); B <- ncol(beta); p <- ncol(X)
  XtX <- crossprod(X)
  A <- tryCatch(solve(XtX), error = function(e) NULL)
  if (is.null(A)) stop("X'X is not SPD; OLS failed", call. = FALSE)
  Xt <- t(X)
  coef <- matrix(NA_real_, nrow = p, ncol = B)
  se   <- matrix(NA_real_, nrow = p, ncol = B)
  sigma2 <- numeric(B)
  df_res <- numeric(B)
  dA <- diag(A)
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
    bh <- A %*% (Xt %*% y)
    r <- as.numeric(y - X %*% bh)
    dff <- max(N - p, 1)
    s2 <- sum(r * r) / dff
    coef[, b] <- as.numeric(bh)
    se[, b] <- sqrt(pmax(dA * s2, 0))
    sigma2[b] <- s2
    df_res[b] <- dff
    if (!is.null(cov_tri)) cov_tri[, b] <- pack_tri(A, s2)
  }
  out <- list(coef = coef, se_coef = se, sigma2 = sigma2, df_res = df_res)
  if (!is.null(cov_tri)) out$cov_tri <- cov_tri
  out
}
