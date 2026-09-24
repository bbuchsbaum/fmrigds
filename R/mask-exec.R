# Mask execution ----------------------------------------------------------

apply_mask_policy <- function(node, arrays, space) {
  policies <- node$policy
  if (inherits(policies, "gds_mask_policy")) {
    policies <- list(policies)
  }

  dims <- dim(arrays[[1]])
  keep <- matrix(TRUE, nrow = dims[1], ncol = dims[3])
  for (policy in policies) {
    if (isTRUE(policy$zero_is_missing)) {
      reference <- arrays$beta %||% arrays[[1]]
      background <- is.finite(reference) & reference == 0
      if (any(background)) {
        arrays <- lapply(arrays, function(a) {
          a[background] <- NA_real_
          a
        })
      }
    }
    keep <- keep & .mask_compute(policy, arrays)
    if (any(!keep)) {
      arrays <- lapply(arrays, function(a) {
        for (k in seq_len(dims[3])) a[!keep[, k], , k] <- NA_real_
        a
      })
    }
  }

  idx <- which(rowSums(keep) > 0L)
  if (!length(idx)) stop("Mask removed all samples", call. = FALSE)

  arrays <- lapply(arrays, function(a) a[idx, , , drop = FALSE])
  space <- .update_space_after_mask(space, idx)
  arrays <- .sync_derived(arrays)

  list(arrays = arrays, space = space, subset = list(samples = idx))
}

.mask_compute <- function(policy, arrays) {
  dims <- dim(arrays[[1]])
  samples <- dims[1]

  finite_mat <- .finite_mask(arrays)

  if (policy$rule == "custom") {
    if (!is.function(policy$custom)) {
      stop("Custom mask policy requires a function", call. = FALSE)
    }
    res <- policy$custom(arrays)
    if (!is.logical(res) || length(res) != samples || anyNA(res)) {
      stop("Custom mask must return a non-missing logical vector of length sample", call. = FALSE)
    }
    return(matrix(res, nrow = samples, ncol = dims[3]))
  }

  keep <- switch(policy$scope,
    group = .mask_group(policy$rule, finite_mat, policy$threshold),
    subject = matrix(
      .mask_subject(policy$rule, finite_mat, policy$threshold),
      nrow = samples, ncol = dims[3]
    )
  )
  keep
}

.finite_mask <- function(arrays) {
  if ("beta" %in% names(arrays)) {
    finite_vals <- is.finite(arrays$beta)
  } else {
    finite_vals <- is.finite(arrays[[1]])
  }
  finite_vals
}

.mask_group <- function(rule, finite_mat, threshold) {
  dim_s <- dim(finite_mat)
  prop <- matrix(0, nrow = dim_s[1], ncol = dim_s[3])
  for (k in seq_len(dim_s[3])) {
    prop[, k] <- rowMeans(matrix(finite_mat[, , k], nrow = dim_s[1]))
  }
  switch(rule,
    intersection = prop == 1,
    union = prop > 0,
    threshold = prop >= threshold
  )
}

.mask_subject <- function(rule, finite_mat, threshold) {
  # Compute per-sample indicator aggregated over subjects
  keep <- rep(TRUE, dim(finite_mat)[1])
  for (j in seq_len(dim(finite_mat)[2])) {
    subj_mat <- finite_mat[, j, , drop = FALSE]
    reshaped <- matrix(subj_mat, nrow = dim(subj_mat)[1])
    prop <- rowMeans(reshaped)
    rule_keep <- switch(rule,
      intersection = prop == 1,
      union = prop > 0,
      threshold = prop >= threshold
    )
    keep <- keep & rule_keep
  }
  keep
}

.update_space_after_mask <- function(space, idx) {
  if (inherits(space, "space_voxel")) {
    if (!is.null(space$mask_idx)) {
      space$mask_idx <- space$mask_idx[idx]
    } else {
      space$mask_idx <- idx
    }
    space$storage <- "packed"
    space$mask_bitmap <- NULL
  }
  if (inherits(space, c("space_parcels", "space_sample_labels"))) {
    space$labels <- space$labels[idx]
  }
  space
}
