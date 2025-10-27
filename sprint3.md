# Sprint 3 Plan — Statistical Operations & Derivations (of 8)

**Sprint:** 3 of 8
**Duration:** 3 weeks
**Prerequisites:** Sprints 1-2 complete (core foundations + adapters + compute scaffold)
**Target:** ≥95% test coverage for statistical operations modules

---

## Objectives

1. **Implement complete derivation engine** with all statistical transformations
2. **Build variance propagation machinery** for statistically correct spatial mappings
3. **Enable `map_to()` with uncertainty propagation** for space transformations
4. **Ensure statistical correctness** throughout the pipeline (never average SE, refuse invalid mappings)

---

## Scope & Tasks

### 1. Complete Derivation Engine

**File:** `R/derive-stats.R`

#### Task 1.1: Implement All Derivation Rules

Per TECHNICAL_SPECIFICATION.md §3.3 and §6.3:

**Derivation Rules to Implement:**

```r
# var ↔ se  (already sketched in Sprint 1, finalise here)
var = se^2
se = sqrt(var)

# t from beta + var (+ df separately)
t = beta / sqrt(var)

# z from t + df (two-sided, sign-preserving)
p_two = 2 * pt(-abs(t), df)
z = qnorm(1 - p_two/2) * sign(t)

# t from z + df (inverse of above)
t = sign(z) * qt(pnorm(abs(z)), df)

# p from t + df (two-sided)
p = 2 * pt(-abs(t), df)

# p from z (two-sided)
p = 2 * pnorm(-abs(z))

# p from F + df1 + df2
p = pf(F, df1, df2, lower.tail = FALSE)

# p from chi2 + df
p = pchisq(chi2, df, lower.tail = FALSE)

# z from p (optional when sign is available)
if (sign_effect available) z = qnorm(1 - p/2) * sign_effect
```

**Implementation Functions:**

```r
#' Execute derive operation on arrays (run inside compute())
#'
#' @param arrays Named list of assays [sample × subject × contrast]
#' @param what Character vector of assays to derive
#' @param options List of derivation options (e.g., overwrite = TRUE)
#'
#' @return Modified arrays with derived assays
execute_derive <- function(arrays, what, options = list()) {
  overwrite <- isTRUE(options$overwrite)

  for (target in what) {
    if (target %in% names(arrays) && !overwrite) next

    arrays[[target]] <- switch(target,
      var = derive_var(arrays),
      se  = derive_se(arrays),
      t   = derive_t(arrays),
      z   = derive_z(arrays),
      p   = derive_p(arrays),
      stop("Unknown derivation target: ", target, call. = FALSE)
    )
  }

  arrays
}

#' Derive variance from standard error
#' @keywords internal
derive_var <- function(block_data) {
  if (!"se" %in% names(block_data)) {
    stop("Cannot derive var: se not available", call. = FALSE)
  }
  block_data$se^2
}

#' Derive standard error from variance
#' @keywords internal
derive_se <- function(block_data) {
  if (!"var" %in% names(block_data)) {
    stop("Cannot derive se: var not available", call. = FALSE)
  }
  sqrt(block_data$var)
}

#' Derive t-statistic from beta and variance
#' @keywords internal
derive_t <- function(block_data) {
  if (!all(c("beta", "var") %in% names(block_data))) {
    stop("Cannot derive t: need beta and var", call. = FALSE)
  }

  # Note: df must be provided separately (not derived)
  block_data$beta / sqrt(block_data$var)
}

#' Derive z-score from t-statistic/p-values
#' @keywords internal
derive_z <- function(arrays) {
  if ("t" %in% names(arrays) && "df" %in% names(arrays)) {
    t <- arrays$t
    df <- arrays$df
    if (length(dim(df)) == 2) {
      df <- array(df, dim = dim(t))
    }
    p_two <- 2 * pt(-abs(t), df)
    return(qnorm(1 - p_two/2) * sign(t))
  }

  if ("p" %in% names(arrays)) {
    if (!"beta" %in% names(arrays)) {
      stop("Cannot derive signed z from p without sign information", call. = FALSE)
    }
    p <- arrays$p
    sign_effect <- sign(arrays$beta)
    return(qnorm(1 - p/2) * sign_effect)
  }

  stop("Cannot derive z: need {t, df} or {p, sign}", call. = FALSE)
}

#' Derive p-value from test statistic
#' @keywords internal
derive_p <- function(arrays) {
  if ("t" %in% names(arrays) && "df" %in% names(arrays)) {
    return(2 * pt(-abs(arrays$t), arrays$df))
  }
  if ("z" %in% names(arrays)) {
    return(2 * pnorm(-abs(arrays$z)))
  }
  if ("F" %in% names(arrays) && all(c("df1", "df2") %in% names(arrays))) {
    return(pf(arrays$F, arrays$df1, arrays$df2, lower.tail = FALSE))
  }
  if ("chi2" %in% names(arrays) && "df" %in% names(arrays)) {
    return(pchisq(arrays$chi2, arrays$df, lower.tail = FALSE))
  }
  stop("Cannot derive p: need {t,df}, z, {F,df1,df2}, or {chi2,df}", call. = FALSE)
}
```

#### Task 1.2: Integration with Executor

Update `R/executor.R` to use complete `execute_derive()`:

```r
#' Apply derive operation (complete implementation)
#' @keywords internal
apply_derive <- function(node, block_data) {
  execute_derive(block_data, node$what, node$options)
}
```

#### Task 1.3: Validation & Constraints

Add validation for statistical requirements:

```r
#' Validate assay combinations for derivations
#' @keywords internal
validate_derivation_requirements <- function(block_data, target) {
  switch(target,
    var = if (!"se" %in% names(block_data)) {
      stop("Cannot derive var without se", call. = FALSE)
    },
    se = if (!"var" %in% names(block_data)) {
      stop("Cannot derive se without var", call. = FALSE)
    },
    t = if (!all(c("beta", "var") %in% names(block_data))) {
      stop("Cannot derive t without beta and var", call. = FALSE)
    },
    z = if (!("t" %in% names(block_data) && "df" %in% names(block_data)) &&
            !"p" %in% names(block_data)) {
      stop("Cannot derive z without (t + df) or p", call. = FALSE)
    }
  )
  invisible(TRUE)
}
```

**Tests:** `tests/testthat/test-derive-stats.R`

```r
test_that("derive var from se", {
  se <- array(runif(10 * 5 * 2, 0.1, 2), dim = c(10, 5, 2))
  block_data <- list(se = se)

  result <- execute_derive(block_data, "var")

  expect_equal(result$var, se^2)
})

test_that("derive t from beta and var", {
  beta <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
  var <- array(runif(10 * 5 * 2, 0.1, 2), dim = c(10, 5, 2))
  block_data <- list(beta = beta, var = var)

  result <- execute_derive(block_data, "t")

  expected_t <- beta / sqrt(var)
  expect_equal(result$t, expected_t)
})

test_that("derive z from t and df", {
  t <- array(rnorm(10 * 5 * 2, mean = 2), dim = c(10, 5, 2))
  df <- array(30, dim = c(5, 2))  # Subject × contrast df
  block_data <- list(t = t, df = df)

  result <- execute_derive(block_data, "z")

  # Check dimensions
  expect_equal(dim(result$z), dim(t))

  # Check first element manually
  p_two <- 2 * pt(-abs(t[1, 1, 1]), df[1, 1])
  expected_z <- qnorm(1 - p_two/2) * sign(t[1, 1, 1])
  expect_equal(result$z[1, 1, 1], expected_z)
})

test_that("derive refuses when requirements not met", {
  block_data <- list(beta = array(1, c(10, 5, 2)))

  expect_error(
    execute_derive(block_data, "t"),
    "need beta and var"
  )
})

test_that("overwrite control works", {
  block_data <- list(
    beta = array(1, c(10, 5, 2)),
    var = array(0.25, c(10, 5, 2)),
    t = array(999, c(10, 5, 2))  # Existing t
  )

  # Default: don't overwrite
  result1 <- execute_derive(block_data, "t", options = list())
  expect_equal(result1$t[1, 1, 1], 999)

  # With overwrite: recompute
  result2 <- execute_derive(block_data, "t", options = list(overwrite = TRUE))
  expected_t <- 1 / sqrt(0.25)  # = 2
  expect_equal(result2$t[1, 1, 1], expected_t)
})
```

---

### 2. Variance Propagation

**File:** `R/propagate-variance.R`

#### Task 2.1: Independent Variance Propagation

Per TECHNICAL_SPECIFICATION.md §6.1.1:

```r
#' Propagate variance through linear mapping (independence assumption)
#'
#' @param M Mapping matrix [n_target × n_source], sparse or dense
#' @param beta Effect array [n_source × n_subjects × n_contrasts]
#' @param var Variance array (same dims as beta)
#'
#' @return List with beta_out and var_out
#' @export
propagate_variance_independent <- function(M, beta, var) {
  stopifnot(
    is.matrix(M) || inherits(M, "Matrix"),
    is.array(beta),
    is.array(var),
    identical(dim(beta), dim(var)),
    ncol(M) == nrow(beta)
  )

  dims <- dim(beta)
  n_target <- nrow(M)

  # Allocate output
  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  # Element-wise square of M for variance propagation
  M2 <- M^2

  # Process each subject × contrast
  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      beta_out[, j, k] <- as.vector(M %*% beta[, j, k])
      var_out[, j, k] <- as.vector(M2 %*% var[, j, k])
    }
  }

  # Set dimnames if present
  if (!is.null(dimnames(beta))) {
    dimnames(beta_out) <- list(NULL, dimnames(beta)[[2]], dimnames(beta)[[3]])
    dimnames(var_out) <- list(NULL, dimnames(var)[[2]], dimnames(var)[[3]])
  }

  list(beta = beta_out, var = var_out)
}
```

#### Task 2.2: Covariance-Based Propagation

Per TECHNICAL_SPECIFICATION.md §6.1.2:

```r
#' Propagate variance with full covariance structure
#'
#' @param M Mapping matrix [n_target × n_source]
#' @param beta Effect array [n_source × n_subjects × n_contrasts]
#' @param var Variance array (diagonal of covariance)
#' @param cov_provider Function(indices) -> Sigma_block
#'
#' @return List with beta_out and var_out
#' @export
propagate_variance_covariance <- function(M, beta, var, cov_provider) {
  stopifnot(
    is.matrix(M) || inherits(M, "Matrix"),
    is.function(cov_provider)
  )

  dims <- dim(beta)
  n_target <- nrow(M)

  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  # Beta mapping (same as independent)
  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      beta_out[, j, k] <- as.vector(M %*% beta[, j, k])
    }
  }

  # Variance propagation per target row (slower, accurate)
  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      for (i in seq_len(n_target)) {
        # Find non-zero entries in row i
        if (inherits(M, "Matrix")) {
          idx <- which(M[i, ] != 0)
        } else {
          idx <- which(M[i, ] != 0)
        }

        if (length(idx) == 0) {
          var_out[i, j, k] <- 0
          next
        }

        w <- M[i, idx]

        # Get covariance block for these indices
        Sigma <- cov_provider(idx, j, k)  # May vary by subject/contrast

        # Var(y[i]) = w^T Sigma w
        var_out[i, j, k] <- drop(t(w) %*% Sigma %*% w)
      }
    }
  }

  list(beta = beta_out, var = var_out)
}
```

#### Task 2.3: Satterthwaite DF Aggregation

Per TECHNICAL_SPECIFICATION.md §6.2.1:

```r
#' Satterthwaite approximation for combined degrees of freedom
#'
#' @param weights Coefficients [n] (e.g., row of mapping matrix)
#' @param variances Unbiased variances [n]
#' @param dfs Degrees of freedom [n]
#'
#' @return Approximated combined df (scalar)
#' @export
satterthwaite_df <- function(weights, variances, dfs) {
  stopifnot(
    length(weights) == length(variances),
    length(weights) == length(dfs),
    all(variances >= 0),
    all(dfs > 0)
  )

  # Remove zero-weight entries
  keep <- weights != 0
  w <- weights[keep]
  v <- variances[keep]
  df <- dfs[keep]

  if (length(w) == 0) return(NA_real_)

  # Numerator: (Σ w_i^2 v_i)^2
  numerator <- sum(w^2 * v)^2

  # Denominator: Σ (w_i^4 v_i^2 / df_i)
  denominator <- sum(w^4 * v^2 / df)

  if (denominator == 0) return(Inf)

  numerator / denominator
}

#' Aggregate df through linear mapping
#'
#' @param M Mapping matrix
#' @param var Variance array
#' @param df Degrees of freedom array (may be [subject × contrast])
#'
#' @return df_out array
#' @keywords internal
aggregate_df_satterthwaite <- function(M, var, df) {
  dims <- dim(var)
  n_target <- nrow(M)

  # Determine df dimensions
  df_is_broadcast <- length(dim(df)) == 2

  if (df_is_broadcast) {
    # df is [subject × contrast], constant across samples
    df_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        for (i in seq_len(n_target)) {
          idx <- which(M[i, ] != 0)
          if (length(idx) > 0) {
            w <- M[i, idx]
            v <- var[idx, j, k]
            df_vec <- rep(df[j, k], length(idx))
            df_out[i, j, k] <- satterthwaite_df(w, v, df_vec)
          }
        }
      }
    }
  } else {
    # df is [sample × subject × contrast]
    df_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        for (i in seq_len(n_target)) {
          idx <- which(M[i, ] != 0)
          if (length(idx) > 0) {
            w <- M[i, idx]
            v <- var[idx, j, k]
            df_vec <- df[idx, j, k]
            df_out[i, j, k] <- satterthwaite_df(w, v, df_vec)
          }
        }
      }
    }
  }

  df_out
}
```

**Tests:** `tests/testthat/test-propagate-variance.R`

```r
test_that("independent propagation is correct", {
  # Simple 2→1 mapping
  M <- matrix(c(0.6, 0.8), nrow = 1)
  beta <- array(c(10, 20), dim = c(2, 1, 1))
  var <- array(c(1, 4), dim = c(2, 1, 1))

  result <- propagate_variance_independent(M, beta, var)

  # Expected beta: 0.6*10 + 0.8*20 = 22
  expect_equal(result$beta[1, 1, 1], 22)

  # Expected var: 0.6^2*1 + 0.8^2*4 = 0.36 + 2.56 = 2.92
  expect_equal(result$var[1, 1, 1], 2.92)
})

test_that("variance is always non-negative", {
  set.seed(42)
  M <- matrix(runif(10 * 20, -1, 1), nrow = 10)
  beta <- array(rnorm(20 * 5 * 2), dim = c(20, 5, 2))
  var <- array(runif(20 * 5 * 2, 0.1, 2), dim = c(20, 5, 2))

  result <- propagate_variance_independent(M, beta, var)

  expect_true(all(result$var >= 0, na.rm = TRUE))
})

test_that("covariance propagation with identity Sigma equals independent", {
  M <- matrix(runif(5 * 10), nrow = 5)
  beta <- array(rnorm(10 * 3 * 2), dim = c(10, 3, 2))
  var <- array(runif(10 * 3 * 2, 0.1, 1), dim = c(10, 3, 2))

  # Identity covariance (diagonal)
  cov_provider <- function(idx, j, k) diag(var[idx, j, k])

  result_ind <- propagate_variance_independent(M, beta, var)
  result_cov <- propagate_variance_covariance(M, beta, var, cov_provider)

  expect_equal(result_cov$beta, result_ind$beta)
  expect_equal(result_cov$var, result_ind$var, tolerance = 1e-10)
})

test_that("Satterthwaite df aggregation is correct", {
  # Known example: two groups with different variances and df
  weights <- c(0.5, 0.5)
  variances <- c(1, 4)
  dfs <- c(10, 20)

  df_combined <- satterthwaite_df(weights, variances, dfs)

  # Manual calculation
  numerator <- (0.5^2 * 1 + 0.5^2 * 4)^2  # = 1.5625
  denominator <- (0.5^4 * 1^2 / 10) + (0.5^4 * 4^2 / 20)  # = 0.00625 + 0.02 = 0.02625
  expected <- numerator / denominator  # ≈ 59.5

  expect_equal(df_combined, expected, tolerance = 1e-6)
})
```

---

### 3. map_to() Implementation

**File:** `R/verb-map-to.R` (enhance) and `R/executor.R`

#### Task 3.1: Assay-Aware Mapping Logic

Guiding rules before implementing code:
1. Only assays with roles `location`, `variance`, or `stdev` may be mapped linearly. Others (`t`, `z`, `p`, `evidence`, `log_evidence`, `posterior`, etc.) must remain untouched or be recomputed explicitly after mapping.
2. When beta/var are mapped, recompute any derived statistics (se, t, z) using the mapped values and aggregated df.
3. If beta/var are absent, refuse generic mapping; only allow explicit evidence combiners (Stouffer/Fisher) for test-stat-only workflows.

```r
#' Apply map_to operation to block
#'
#' @param node Operation node from plan
#' @param block_data Named list of arrays
#' @param executor Executor environment
#'
#' @return Transformed block_data
#' @keywords internal
apply_map_to <- function(node, block_data, executor) {
  map <- node$map
  uncertainty <- node$uncertainty
  target_space <- node$target_space
  combine <- node$combine

  # Get mapping operator
  M <- if (inherits(map, "gds_map")) {
    if (!is.null(map$by_subject)) {
      stop("Subject-specific maps require align(); use align() before map_to()",
           call. = FALSE)
    }
    map$operator
  } else if (is.matrix(map) || inherits(map, "Matrix")) {
    map
  } else {
    stop("Invalid map: must be gds_map or matrix", call. = FALSE)
  }

  # Check what assays are available
  has_beta_var <- all(c("beta", "var") %in% names(block_data))
  has_t_df <- all(c("t", "df") %in% names(block_data))
  has_z <- "z" %in% names(block_data)

  # Only location/uncertainty roles may be mapped linearly
  linear_assays <- names(block_data)[vapply(names(block_data), can_map_linear, logical(1))]
  non_linear <- setdiff(names(block_data), c(linear_assays, "t", "z", "df", "p", "chi2"))
  if (length(non_linear)) {
    message("Skipping linear mapping for assays of role: ", paste(non_linear, collapse = ", "))
  }

  if (has_beta_var) {
    # LINEAR MAPPING PATH: Map beta and var, recompute test stats

    # Propagate beta and var
    if (uncertainty$mode == "independent") {
      propagated <- propagate_variance_independent(M, block_data$beta, block_data$var)
    } else if (uncertainty$mode == "cov_provider") {
      propagated <- propagate_variance_covariance(M, block_data$beta, block_data$var,
                                                  uncertainty$cov_provider)
    } else {
      stop("Unsupported uncertainty mode: ", uncertainty$mode, call. = FALSE)
    }

    block_data$beta <- propagated$beta
    block_data$var <- propagated$var

    # Derive se if it existed
    if ("se" %in% names(block_data)) {
      block_data$se <- sqrt(block_data$var)
    }

    # Recompute t if it existed
    if ("t" %in% names(block_data)) {
      block_data$t <- block_data$beta / sqrt(block_data$var)

      # Aggregate df if present
      if ("df" %in% names(block_data) && uncertainty$df_rule == "satterthwaite") {
        block_data$df <- aggregate_df_satterthwaite(M, block_data$var, block_data$df)
      }
    }

    # Recompute z if it existed
    if ("z" %in% names(block_data) && "t" %in% names(block_data) &&
        "df" %in% names(block_data)) {
      block_data$z <- derive_z(block_data)
    }

  } else if (has_z || has_t_df) {
    # TEST STATISTIC PATH: Use explicit combiners

    if (is.null(combine)) {
      stop("Cannot linearly map test statistics without effect scale.\n",
           "  Either provide beta+var, or specify combine='stouffer' or 'fisher'",
           call. = FALSE)
    }

    if (combine == "stouffer" && has_z) {
      # Stouffer combination for z-scores
      block_data$z <- combine_stouffer_spatial(M, block_data$z)
    } else if (combine == "fisher" && has_z) {
      # Fisher combination
      block_data$z <- combine_fisher_spatial(M, block_data$z)
    } else {
      stop("Unsupported combiner '", combine, "' for available assays", call. = FALSE)
    }

    # Remove incompatible assays (t, df cannot be linearly mapped)
    if ("t" %in% names(block_data)) block_data$t <- NULL
    if ("df" %in% names(block_data)) block_data$df <- NULL

  } else {
    stop("No mappable assays found. Need (beta + var) or z", call. = FALSE)
  }

  block_data
}

#' Stouffer combination for spatial mapping
#'
#' @param M Mapping matrix (interpreted as weights)
#' @param z Z-score array
#'
#' @return Combined z-scores
#' @keywords internal
combine_stouffer_spatial <- function(M, z) {
  dims <- dim(z)
  n_target <- nrow(M)
  z_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      for (i in seq_len(n_target)) {
        idx <- which(M[i, ] != 0)
        if (length(idx) > 0) {
          w <- M[i, idx]
          z_vals <- z[idx, j, k]

          # Weighted Stouffer: (Σ w_i z_i) / sqrt(Σ w_i^2)
          z_out[i, j, k] <- sum(w * z_vals) / sqrt(sum(w^2))
        }
      }
    }
  }

  z_out
}

#' Fisher combination for spatial mapping
#'
#' @keywords internal
combine_fisher_spatial <- function(M, z) {
  dims <- dim(z)
  n_target <- nrow(M)
  z_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      for (i in seq_len(n_target)) {
        idx <- which(M[i, ] != 0)
        if (length(idx) > 0) {
          z_vals <- z[idx, j, k]

          # Convert to p-values (two-tailed)
          p_vals <- 2 * pnorm(-abs(z_vals))

          # Fisher's method: chi2 = -2 Σ log(p)
          chi2 <- -2 * sum(log(p_vals))
          df <- 2 * length(idx)

          # Combined p-value
          p_combined <- pchisq(chi2, df, lower.tail = FALSE)

          # Convert back to z (preserve sign)
          sign_combined <- sign(mean(z_vals))
          z_out[i, j, k] <- qnorm(p_combined / 2, lower.tail = FALSE) * sign_combined
        }
      }
    }
  }

  z_out
}
```

#### Task 3.2: Integration with Executor

Update `R/executor.R`:

```r
apply_operation <- function(node, block_data, executor) {
  switch(node$op,
    subset_axis = apply_subset(node, block_data),
    derive = apply_derive(node, block_data),
    map = apply_map_to(node, block_data, executor),  # NEW
    align_to_group = block_data,  # TODO: Sprint 4
    mask_policy = block_data,     # TODO: Sprint 4
    reduce = block_data,           # TODO: Sprint 5
    write = block_data,            # TODO: Sprint 7
    block_data
  )
}
```

**Tests:** `tests/testthat/test-map-to.R`

```r
test_that("map_to propagates beta and var correctly", {
  # 2→1 mapping
  M <- matrix(c(0.6, 0.8), nrow = 1)

  beta <- array(c(10, 20), dim = c(2, 1, 1))
  var <- array(c(1, 4), dim = c(2, 1, 1))

  node <- list(
    op = "map",
    map = M,
    uncertainty = UncertaintyRule("independent"),
    combine = NULL
  )

  block_data <- list(beta = beta, var = var)

  result <- apply_map_to(node, block_data, NULL)

  expect_equal(result$beta[1, 1, 1], 22)
  expect_equal(result$var[1, 1, 1], 2.92)
})

test_that("map_to recomputes t after mapping", {
  M <- matrix(c(0.5, 0.5), nrow = 1)

  beta <- array(c(2, 4), dim = c(2, 1, 1))
  var <- array(c(1, 1), dim = c(2, 1, 1))
  t_orig <- array(c(2, 4), dim = c(2, 1, 1))

  node <- list(
    op = "map",
    map = M,
    uncertainty = UncertaintyRule("independent"),
    combine = NULL
  )

  block_data <- list(beta = beta, var = var, t = t_orig)

  result <- apply_map_to(node, block_data, NULL)

  # Expected: beta'=3, var'=0.5, t'=3/sqrt(0.5)≈4.24
  expect_equal(result$beta[1, 1, 1], 3)
  expect_equal(result$var[1, 1, 1], 0.5)
  expect_equal(result$t[1, 1, 1], 3 / sqrt(0.5), tolerance = 1e-6)
})

test_that("map_to refuses to map t without beta/var", {
  M <- matrix(c(0.5, 0.5), nrow = 1)

  t <- array(c(2, 3), dim = c(2, 1, 1))
  df <- array(30, dim = c(1, 1))

  node <- list(
    op = "map",
    map = M,
    uncertainty = UncertaintyRule("independent"),
    combine = NULL
  )

  block_data <- list(t = t, df = df)

  expect_error(
    apply_map_to(node, block_data, NULL),
    "Cannot linearly map test statistics"
  )
})

test_that("Stouffer combiner works for z-only mapping", {
  M <- matrix(c(0.5, 0.5), nrow = 1)

  z <- array(c(2, 3), dim = c(2, 1, 1))

  node <- list(
    op = "map",
    map = M,
    uncertainty = UncertaintyRule("independent"),
    combine = "stouffer"
  )

  block_data <- list(z = z)

  result <- apply_map_to(node, block_data, NULL)

  # Weighted Stouffer: (0.5*2 + 0.5*3) / sqrt(0.5^2 + 0.5^2)
  # = 2.5 / sqrt(0.5) ≈ 3.536
  expected <- 2.5 / sqrt(0.5)
  expect_equal(result$z[1, 1, 1], expected, tolerance = 1e-6)
})
```

---

### 4. Testing & Quality Assurance

#### Task 4.1: Unit Test Coverage

**Test Files:**
- `test-derive-stats.R` (all derivation rules)
- `test-propagate-variance.R` (independent, covariance, Satterthwaite)
- `test-map-to.R` (mapping with uncertainty, combiners)

**Coverage Target:** ≥95% for:
- `R/derive-stats.R`
- `R/propagate-variance.R`
- `R/verb-map-to.R` (enhanced)

#### Task 4.2: Integration Tests

**File:** `tests/testthat/test-integration-sprint3.R`

```r
test_that("End-to-end: derive + map_to with variance propagation", {
  skip_if_not_installed("data.table")

  # Create ROI data
  df <- expand.grid(
    subject = paste0("sub-", 1:5),
    roi = paste0("ROI_", 1:10),
    contrast = "A"
  )
  df$beta <- rnorm(nrow(df), mean = 1)
  df$se <- runif(nrow(df), 0.1, 0.5)

  tmp <- tempfile(fileext = ".csv")
  write.csv(df, tmp, row.names = FALSE)
  on.exit(unlink(tmp))

  # Create mapping (10 ROIs → 3 clusters)
  M <- Matrix::sparseMatrix(
    i = c(rep(1, 3), rep(2, 3), rep(3, 4)),
    j = c(1:3, 4:6, 7:10),
    x = rep(1/c(3, 3, 4), c(3, 3, 4)),
    dims = c(3, 10)
  )

  cluster_space <- space_parcels(labels = paste0("Cluster_", 1:3))

  # Pipeline
  plan <- gds(tmp) %>%
    derive(c("var", "t")) %>%
    map_to(target_space = cluster_space,
           map = M,
           uncertainty = UncertaintyRule("independent"))

  result <- compute(plan, sink = "memory")

  # Validations
  expect_equal(dim(assay(result, "beta")), c(3, 5, 1))
  expect_true("var" %in% names(result$assays))
  expect_true("t" %in% names(result$assays))
  expect_equal(result$space$type, "parcels")
  expect_equal(result$space$labels, paste0("Cluster_", 1:3))
})

test_that("Statistical correctness: variance never negative", {
  skip_if_not_installed("data.table")

  # Large random dataset
  set.seed(42)
  df <- expand.grid(
    subject = paste0("sub-", 1:20),
    roi = paste0("ROI_", 1:50),
    contrast = c("A", "B")
  )
  df$beta <- rnorm(nrow(df))
  df$var <- runif(nrow(df), 0.01, 2)

  tmp <- tempfile(fileext = ".csv")
  data.table::fwrite(df, tmp)
  on.exit(unlink(tmp))

  # Random mapping
  M <- matrix(runif(25 * 50), nrow = 25)
  M <- M / rowSums(M)  # Row-normalize

  plan <- gds(tmp) %>%
    map_to(target_space = space_parcels(labels = paste0("S_", 1:25)),
           map = M,
           uncertainty = UncertaintyRule("independent"))

  result <- compute(plan)

  # All variances must be non-negative
  expect_true(all(assay(result, "var") >= 0, na.rm = TRUE))
})
```

#### Task 4.3: Edge Case Tests

```r
test_that("Handles NA values in derivations", {
  beta <- array(c(1, NA, 3), dim = c(3, 1, 1))
  var <- array(c(0.25, 0.36, NA), dim = c(3, 1, 1))

  block_data <- list(beta = beta, var = var)
  result <- execute_derive(block_data, "t")

  expect_true(is.na(result$t[2, 1, 1]))
  expect_true(is.na(result$t[3, 1, 1]))
  expect_equal(result$t[1, 1, 1], 1 / sqrt(0.25))
})

test_that("Zero variance handling", {
  beta <- array(c(1, 2), dim = c(2, 1, 1))
  var <- array(c(0, 0.25), dim = c(2, 1, 1))

  block_data <- list(beta = beta, var = var)
  result <- execute_derive(block_data, "t")

  expect_true(is.infinite(result$t[1, 1, 1]))
  expect_equal(result$t[2, 1, 1], 2 / sqrt(0.25))
})
```

---

## Deliverables

### Core Functionality
- ✅ Complete derivation engine (var, se, t, z, p)
- ✅ Independent variance propagation
- ✅ Covariance-based variance propagation
- ✅ Satterthwaite df aggregation
- ✅ map_to() with assay-aware logic
- ✅ Statistical correctness enforcement (refuse invalid mappings)
- ✅ Stouffer and Fisher combiners for z-scores

### Testing & Quality
- ✅ ≥95% test coverage for statistical modules
- ✅ Unit tests for all derivation rules
- ✅ Unit tests for variance propagation modes
- ✅ Integration tests with realistic workflows
- ✅ Edge case tests (NA, zero variance, dimension mismatches)
- ✅ CI passing cleanly

### Documentation
- ✅ Roxygen2 documentation for all exported functions
- ✅ Mathematical formulas in documentation
- ✅ Examples demonstrating statistical correctness
- ✅ NEWS.md entry for Sprint 3

---

## Technical Debt & Deferred Work

**Deferred to Sprint 4:**
- Kernel-based covariance estimation (uncertainty mode = "kernel")
- Optimization of covariance propagation for large blocks
- `align()` execution (subject-specific transforms)

**Deferred to Sprint 5:**
- `reduce()` meta-analysis (uses variance propagation built here)

**Deferred to Later:**
- Random-effects variance propagation
- Spatial correlation models for covariance

---

## Dependencies & Notes

**Technical Specification References:**
- Derivation rules: §3.3, §6.3
- Variance propagation: §6.1.1, §6.1.2
- Satterthwaite: §6.2.1
- map_to() execution: §3.6
- Statistical invariants: §6.5

**Package Dependencies:**
- **Required imports:** Matrix (sparse matrices)
- **No new dependencies**

**Implementation Notes:**
- Variance propagation is the **critical path** for correctness
- Test with known analytical solutions (validate formulas)
- Performance optimization deferred to Sprint 8
- Focus on **statistical correctness** over speed

---

## Acceptance Criteria

**Sprint 3 is complete when:**
1. ✅ All derivation rules work correctly (var, se, t, z, p)
2. ✅ Independent variance propagation is mathematically correct
3. ✅ Covariance propagation works with user-supplied covariance
4. ✅ Satterthwaite df aggregation matches known formulas
5. ✅ map_to() refuses to map t/z without effect scale
6. ✅ map_to() recomputes test statistics after mapping beta/var
7. ✅ Stouffer combiner works for z-only mappings
8. ✅ Fisher combiner works for z-only mappings
9. ✅ All unit tests pass (≥95% coverage)
10. ✅ Integration tests demonstrate realistic workflows
11. ✅ Edge cases handled correctly (NA, zero variance)
12. ✅ CI pipeline passes
13. ✅ Documentation complete and accurate

---

## Example Usage (End of Sprint 3)

```r
library(gdsfmri)

# Example 1: Complete derivation pipeline
plan <- gds("roi_stats.csv") %>%
  derive(c("var", "t", "z"))

result <- compute(plan)
print(result$assays)
# List of 4: beta, var, t, z

# Example 2: Spatial mapping with variance propagation
# 100 voxels → 10 parcels
M <- create_parcel_map(voxel_membership)

plan <- gds("voxel_betas.csv") %>%
  derive("var") %>%
  map_to(target_space = parcel_space,
         map = M,
         uncertainty = UncertaintyRule("independent",
                                      df_rule = "satterthwaite"))

result <- compute(plan)
# Variance correctly propagated through mapping
# t-statistics recomputed with aggregated df

# Example 3: z-only mapping with Stouffer
plan <- gds("z_scores.csv") %>%
  map_to(target_space = cluster_space,
         map = M,
         uncertainty = UncertaintyRule("independent"),
         combine = "stouffer")

result <- compute(plan)
# z-scores combined via weighted Stouffer method
```

---

**Sprint 3 Version:** 1.0
**Last Updated:** 2025-01-XX
**Status:** Ready for implementation
