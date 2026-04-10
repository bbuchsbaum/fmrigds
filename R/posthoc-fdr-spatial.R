# Spatial FDR (group-wise) ---------------------------------------------------

# Implements a simple spatial FDR by aggregating p-values within spatial
# groups (e.g., parcels/clusters) using Simes' method, then applying a
# weighted BH procedure across groups with weights proportional to group size
# (normalized to mean 1). The resulting group-level q-values are broadcast
# back to samples.

.posthoc_spatial_fdr <- function() {
  function(arrays, opts) {
    p <- arrays$p
    if (is.null(p)) {
      if (!is.null(arrays$z)) {
        p <- derive_p(list(z = arrays$z))
      } else if (!is.null(arrays$t) && !is.null(arrays$df)) {
        p <- derive_p(list(t = arrays$t, df = arrays$df))
      } else {
        stop("spatial FDR requires p (or z/t+df)", call. = FALSE)
      }
    }
    dims <- dim(p)

    group <- opts$group %||% .resolve_spatial_fdr_group(opts$row_data %||% NULL)
    if (is.null(group)) {
      stop("spatial FDR requires options$group or row_data with a spatial grouping column", call. = FALSE)
    }
    if (length(group) != dims[1L]) stop("length(options$group) must equal sample dimension", call. = FALSE)
    # Normalize group IDs to consecutive integers; preserve NA groups
    gfac <- suppressWarnings(as.factor(group))
    glevels <- levels(gfac)
    gid <- as.integer(gfac)
    n_groups <- length(glevels)
    # Group sizes (exclude NA)
    valid_gid <- !is.na(gid)
    sizes <- tabulate(gid[valid_gid], nbins = n_groups)
    # Normalize weights to mean 1 (sum = n_groups)
    w <- if (n_groups > 0) sizes / mean(sizes[sizes > 0]) else numeric(0)

    out <- array(NA_real_, dim = dims)

    # Helper: Simes aggregation for a vector of p-values in a group
    simes <- function(pv) {
      pv <- pv[is.finite(pv)]
      m <- length(pv)
      if (m == 0) return(NA_real_)
      o <- order(pv)
      ranks <- seq_len(m)
      min(pv[o] * m / ranks)
    }

    # Weighted BH adjusted p-values for group-level p's
    wbh_adjust <- function(p_g, w_g) {
      ok <- is.finite(p_g) & is.finite(w_g) & (w_g > 0)
      if (!any(ok)) return(rep(NA_real_, length(p_g)))
      idx <- which(ok)
      m <- length(idx)
      ord <- order(p_g[idx] / w_g[idx])
      p_sorted <- p_g[idx][ord]
      w_sorted <- w_g[idx][ord]
      # compute adjusted values on the sorted scale
      adj_sorted <- (n_groups * (p_sorted / w_sorted)) / seq_len(m)
      # enforce monotonicity and cap at 1
      adj_sorted <- rev(cummin(rev(adj_sorted)))
      adj_sorted <- pmin(adj_sorted, 1)
      outv <- rep(NA_real_, length(p_g))
      outv[idx[ord]] <- adj_sorted
      outv
    }

    for (j in seq_len(dims[2L])) {
      for (k in seq_len(dims[3L])) {
        pv <- p[, j, k]
        # Aggregate to group-level p-values with Simes
        p_g <- rep(NA_real_, n_groups)
        if (n_groups > 0) {
          for (g in seq_len(n_groups)) {
            members <- which(valid_gid & gid == g)
            if (!length(members)) next
            p_g[g] <- simes(pv[members])
          }
        }
        # Weighted BH across groups
        q_g <- wbh_adjust(p_g, w)
        # Broadcast back to samples within each group
        qv <- rep(NA_real_, length(pv))
        for (g in seq_len(n_groups)) {
          members <- which(valid_gid & gid == g)
          if (!length(members)) next
          qv[members] <- q_g[g]
        }
        out[, j, k] <- qv
      }
    }

    list(q = out)
  }
}

.resolve_spatial_fdr_group <- function(row_data) {
  .sample_groups_from_row_data(row_data)
}
