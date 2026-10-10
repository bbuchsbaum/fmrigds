# Variance propagation -----------------------------------------------------

#' Propagate variance assuming independent samples
#'
#' @param M Mapping matrix \code{[n_target x n_source]}
#' @param beta Effect array \code{[n_source x n_subject x n_contrast]}
#' @param var Variance array (same dims as beta)
#'
#' @return List with `beta` and `var` arrays in target space
#' @keywords internal
propagate_variance_independent <- function(M, beta, var) {
  stopifnot(is.array(beta), is.array(var), identical(dim(beta), dim(var)))
  if (!is.matrix(M) && !inherits(M, "Matrix")) {
    stop("M must be matrix or Matrix", call. = FALSE)
  }
  M <- as.matrix(M)

  dims <- dim(beta)
  n_target <- nrow(M)

  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  M2 <- M^2

  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      beta_out[, j, k] <- M %*% beta[, j, k]
      var_out[, j, k] <- M2 %*% var[, j, k]
    }
  }

  list(beta = beta_out, var = var_out)
}

#' Propagate variance using supplied covariance blocks
#'
#' `cov_provider(idx)` describes the *spatial correlation structure* among the
#' source samples `idx` that feed one target sample. It may return either a
#' correlation matrix or a covariance matrix; in both cases only its correlation
#' structure is used (a covariance is converted as with [stats::cov2cor()]). The
#' variances come from `var`, so each subject/contrast keeps its own
#' uncertainty: for target row `w = M[i, idx]` and
#' `D = diag(sqrt(var[idx, j, k]))`,
#' `Var_target[i, j, k] = w' D R D w`, where `R` is the correlation implied by
#' `cov_provider(idx)`. When `cov_provider` returns a covariance whose diagonal
#' equals `var[idx, j, k]`, this reproduces `w' Sigma w` exactly. Source indices
#' whose provider diagonal is zero or non-finite are treated as uncorrelated
#' with the others.
#'
#' @param M Mapping matrix \code{[n_target x n_source]}
#' @param beta Effect array in source space
#' @param var Variance array (diagonal entries) in source space
#' @param cov_provider Function of source indices returning a square
#'   correlation (or covariance) matrix for those indices
#'
#' @return List with mapped `beta` and `var`
#' @keywords internal
propagate_variance_covariance <- function(M, beta, var, cov_provider) {
  if (!is.function(cov_provider)) {
    stop("cov_provider must be a function", call. = FALSE)
  }
  stopifnot(is.array(beta), is.array(var), identical(dim(beta), dim(var)))
  if (!is.matrix(M) && !inherits(M, "Matrix")) {
    stop("M must be matrix or Matrix", call. = FALSE)
  }
  M <- as.matrix(M)

  dims <- dim(beta)
  n_target <- nrow(M)
  n_source <- dims[1]
  n_jk <- dims[2] * dims[3]

  # Flatten subject x contrast so each target needs one matrix product.
  beta_flat <- matrix(beta, nrow = n_source, ncol = n_jk)
  var_flat <- matrix(var, nrow = n_source, ncol = n_jk)
  sd_flat <- sqrt(pmax(var_flat, 0))
  sd_flat[!is.finite(var_flat)] <- NA_real_

  beta_out <- matrix(NA_real_, nrow = n_target, ncol = n_jk)
  var_out <- matrix(NA_real_, nrow = n_target, ncol = n_jk)

  for (i in seq_len(n_target)) {
    w_full <- M[i, ]
    idx <- which(w_full != 0)
    if (!length(idx)) next
    w <- w_full[idx]
    Sigma <- cov_provider(idx)
    if (!is.matrix(Sigma) || nrow(Sigma) != length(idx) || ncol(Sigma) != length(idx)) {
      stop("cov_provider must return square matrix matching index length", call. = FALSE)
    }
    R <- .cov_to_correlation(Sigma)

    beta_out[i, ] <- drop(crossprod(w, beta_flat[idx, , drop = FALSE]))
    # S[, jk] = w * sd[idx, jk];  Var = diag(S' R S) = colSums(S * (R %*% S))
    S <- w * sd_flat[idx, , drop = FALSE]
    var_out[i, ] <- colSums(S * (R %*% S))
  }

  list(
    beta = array(beta_out, dim = c(n_target, dims[2], dims[3])),
    var = array(var_out, dim = c(n_target, dims[2], dims[3]))
  )
}

.cov_to_correlation <- function(Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  d <- diag(Sigma)
  ok <- is.finite(d) & d > 0
  s <- numeric(length(d))
  s[ok] <- 1 / sqrt(d[ok])
  R <- Sigma * tcrossprod(s)
  R[!ok, ] <- 0
  R[, !ok] <- 0
  diag(R) <- 1
  R
}
