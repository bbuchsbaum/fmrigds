# Mask execution ----------------------------------------------------------

apply_mask_policy <- function(node, arrays, space) {
  policies <- node$policy
  if (inherits(policies, "gds_mask_policy")) {
    policies <- list(policies)
  }

  keep <- rep(TRUE, dim(arrays[[1]])[1])
  for (policy in policies) {
    keep <- keep & .mask_compute(policy, arrays)
  }

  idx <- which(keep)
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
    if (length(res) != samples) stop("Custom mask must return vector of length sample", call. = FALSE)
    return(as.logical(res))
  }

  keep <- switch(policy$scope,
    group = .mask_group(policy$rule, finite_mat, policy$threshold),
    subject = .mask_subject(policy$rule, finite_mat, policy$threshold)
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
  reshaped <- array(finite_mat, dim = c(dim_s[1], dim_s[2] * dim_s[3]))
  prop <- rowMeans(reshaped)
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
  if (inherits(space, "space_parcels")) {
    space$labels <- space$labels[idx]
  }
  space
}
