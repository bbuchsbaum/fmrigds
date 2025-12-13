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
#' @param M Mapping matrix \code{[n_target x n_source]}
#' @param beta Effect array in source space
#' @param var Variance array (diagonal entries) in source space
#' @param cov_provider Function returning covariance matrix for source indices
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

  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (i in seq_len(n_target)) {
    w_full <- M[i, ]
    idx <- which(w_full != 0)
    if (!length(idx)) next
    w <- w_full[idx]
    Sigma <- cov_provider(idx)
    if (!is.matrix(Sigma) || nrow(Sigma) != length(idx) || ncol(Sigma) != length(idx)) {
      stop("cov_provider must return square matrix matching index length", call. = FALSE)
    }

    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        beta_out[i, j, k] <- sum(w * beta[idx, j, k])
        var_out[i, j, k] <- as.numeric(t(w) %*% Sigma %*% w)
      }
    }
  }

  list(beta = beta_out, var = var_out)
}
