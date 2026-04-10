# neurothresh post-hoc bridge ------------------------------------------------- # nocov start

.posthoc_resolve_voxel_spec <- function(opts, n_samples, method_name = "nt:tfce_fwer") {
  sp <- opts$space %||% NULL
  if (is.null(sp)) {
    stop(
      method_name,
      " requires voxel-space metadata. Pass options$space (space_voxel) ",
      "or run via compute() so context is supplied automatically.",
      call. = FALSE
    )
  }
  if (!inherits(sp, "space_voxel")) {
    stop(method_name, " requires voxel space; got ", paste(class(sp), collapse = "/"), call. = FALSE)
  }
  if (is.null(sp$dim) || length(sp$dim) != 3L) {
    stop(method_name, " requires 3D voxel dimensions in space$dim", call. = FALSE)
  }

  vdim <- as.integer(sp$dim)
  n_full <- prod(vdim)
  mask_idx <- sp$mask_idx

  if (is.null(mask_idx)) {
    if (n_samples != n_full) {
      stop(
        method_name,
        " sample dimension (", n_samples, ") does not match dense voxel count (", n_full, ")",
        call. = FALSE
      )
    }
    mask_idx <- seq_len(n_full)
  } else {
    mask_idx <- as.integer(mask_idx)
    if (length(mask_idx) != n_samples) {
      stop(
        method_name,
        " sample dimension (", n_samples, ") must match length(space$mask_idx) (", length(mask_idx), ")",
        call. = FALSE
      )
    }
    if (any(mask_idx < 1L | mask_idx > n_full)) {
      stop(method_name, " space$mask_idx has indices outside 1:prod(space$dim)", call. = FALSE)
    }
  }

  mask_full <- rep.int(FALSE, n_full)
  mask_full[mask_idx] <- TRUE

  list(
    dim = vdim,
    n_full = n_full,
    mask_idx = mask_idx,
    mask = array(mask_full, dim = vdim)
  )
}

.posthoc_unpack_to_full <- function(sample_vec, spec, fill = NA_real_) {
  if (length(sample_vec) != length(spec$mask_idx)) {
    stop("sample vector length does not match voxel mask length", call. = FALSE)
  }
  full <- rep.int(fill, spec$n_full)
  full[spec$mask_idx] <- as.numeric(sample_vec)
  full
}

.posthoc_pack_from_full <- function(full_vec, spec) {
  if (length(full_vec) != spec$n_full) {
    stop("full volume vector length does not match voxel dimensions", call. = FALSE)
  }
  as.numeric(full_vec[spec$mask_idx])
}

.posthoc_pick_null_fun <- function(null_spec, contrast_idx, n_contrast) {
  if (is.null(null_spec)) return(NULL)
  if (is.function(null_spec)) return(null_spec)
  if (!is.list(null_spec)) {
    stop("options$null_fun must be a function or a list of per-contrast functions", call. = FALSE)
  }

  if (!is.null(names(null_spec))) {
    key <- as.character(contrast_idx)
    if (key %in% names(null_spec)) {
      fn <- null_spec[[key]]
      if (!is.function(fn)) {
        stop("options$null_fun[['", key, "']] must be a function", call. = FALSE)
      }
      return(fn)
    }
  }

  if (length(null_spec) == n_contrast) {
    fn <- null_spec[[contrast_idx]]
    if (!is.function(fn)) {
      stop("options$null_fun[[", contrast_idx, "]] must be a function", call. = FALSE)
    }
    return(fn)
  }

  stop(
    "When options$null_fun is a list, provide one function per contrast ",
    "(by position or named with contrast index strings).",
    call. = FALSE
  )
}

.posthoc_tfce_corrected_p <- function(tfce_res, mask_idx, n_full) {
  p_corr <- tfce_res$p_corr %||% NULL
  if (!is.null(p_corr)) {
    p_vec <- as.numeric(p_corr)
    if (length(p_vec) != n_full) {
      stop("neurothresh::tfce_fwer returned p_corr with unexpected length", call. = FALSE)
    }
    return(p_vec)
  }

  tfce_map <- tfce_res$tfce_2 %||% tfce_res$tfce %||% NULL
  max_null <- as.numeric(tfce_res$max_null %||% numeric(0))
  if (is.null(tfce_map) || !length(max_null)) {
    stop("Unable to construct corrected p-values from neurothresh output", call. = FALSE)
  }
  tfce_vec <- as.numeric(tfce_map)
  if (length(tfce_vec) != n_full) {
    stop("neurothresh TFCE map length does not match expected voxel count", call. = FALSE)
  }

  p_vec <- rep.int(1, n_full)
  denom <- length(max_null) + 1
  for (idx in mask_idx) {
    p_vec[idx] <- (1 + sum(max_null >= tfce_vec[idx])) / denom
  }
  p_vec
}

.posthoc_neurothresh_tfce_fwer <- function() {
  function(arrays, opts) {
    if (!requireNamespace("neurothresh", quietly = TRUE)) {
      stop("posthoc method 'nt:tfce_fwer' requires optional package 'neurothresh'", call. = FALSE)
    }

    z <- arrays$z
    if (is.null(z)) {
      stop("nt:tfce_fwer requires z scores", call. = FALSE)
    }

    dims <- dim(z)
    if (length(dims) != 3L) {
      stop("z must be a 3D assay [sample x subject x contrast]", call. = FALSE)
    }
    if (dims[2L] != 1L) {
      stop(
        "nt:tfce_fwer currently expects group-level arrays with subject dimension = 1. ",
        "Run reduce() before this post-hoc method.",
        call. = FALSE
      )
    }

    spec <- .posthoc_resolve_voxel_spec(opts, n_samples = dims[1L], method_name = "nt:tfce_fwer")

    n_perm <- as.integer(opts$n_perm %||% 5000L)
    if (!is.finite(n_perm) || n_perm < 1L) stop("options$n_perm must be a positive integer", call. = FALSE)
    alpha <- as.numeric(opts$alpha %||% 0.05)
    if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("options$alpha must be in (0, 1)", call. = FALSE)
    H <- as.numeric(opts$H %||% 2.0)
    E <- as.numeric(opts$E %||% 0.5)
    dh <- as.numeric(opts$dh %||% 0.1)
    if (!is.finite(H) || H <= 0) stop("options$H must be > 0", call. = FALSE)
    if (!is.finite(E) || E <= 0) stop("options$E must be > 0", call. = FALSE)
    if (!is.finite(dh) || dh <= 0) stop("options$dh must be > 0", call. = FALSE)
    tail <- match.arg(opts$tail %||% "two", c("pos", "neg", "two"))
    allow_voxel_signflip <- isTRUE(opts$allow_voxel_signflip)
    null_spec <- opts$null_fun %||% NULL

    q_out <- array(NA_real_, dim = dims)
    sig_out <- array(0, dim = dims)
    tfce_out <- array(NA_real_, dim = dims)
    attachments <- list()

    for (k in seq_len(dims[3L])) {
      z_vec <- as.numeric(z[, 1L, k])
      if (any(!is.finite(z_vec))) {
        stop("nt:tfce_fwer requires finite z-values in active voxels (contrast ", k, ")", call. = FALSE)
      }

      stat_full <- .posthoc_unpack_to_full(z_vec, spec, fill = NA_real_)
      stat_vol <- array(stat_full, dim = spec$dim)

      null_fun_k <- .posthoc_pick_null_fun(null_spec, contrast_idx = k, n_contrast = dims[3L])
      if (is.null(null_fun_k) && !allow_voxel_signflip) {
        stop(
          "nt:tfce_fwer requires options$null_fun unless allow_voxel_signflip = TRUE.",
          call. = FALSE
        )
      }

      seed_k <- NULL
      if (!is.null(opts$seed)) {
        seed_raw <- opts$seed
        if (length(seed_raw) == 1L) {
          seed_k <- as.integer(seed_raw) + (k - 1L)
        } else if (length(seed_raw) == dims[3L]) {
          seed_k <- as.integer(seed_raw[[k]])
        } else {
          stop("options$seed must be scalar or length equal to number of contrasts", call. = FALSE)
        }
      }

      tfce_res <- neurothresh::tfce_fwer(
        stat_vol = stat_vol,
        mask = spec$mask,
        n_perm = n_perm,
        H = H,
        E = E,
        dh = dh,
        alpha = alpha,
        seed = seed_k,
        null_fun = null_fun_k,
        allow_voxel_signflip = allow_voxel_signflip,
        tail = tail
      )

      q_full <- .posthoc_tfce_corrected_p(tfce_res, mask_idx = spec$mask_idx, n_full = spec$n_full)
      q_out[, 1L, k] <- .posthoc_pack_from_full(q_full, spec)

      sig_full <- as.numeric(as.logical(tfce_res$sig_mask))
      if (length(sig_full) != spec$n_full) {
        stop("neurothresh::tfce_fwer returned sig_mask with unexpected length", call. = FALSE)
      }
      sig_out[, 1L, k] <- .posthoc_pack_from_full(sig_full, spec)

      tfce_map <- tfce_res$tfce %||% tfce_res$tfce_2 %||% NULL
      if (is.null(tfce_map)) {
        stop("neurothresh::tfce_fwer returned no TFCE map", call. = FALSE)
      }
      tfce_full <- as.numeric(tfce_map)
      if (length(tfce_full) != spec$n_full) {
        stop("neurothresh TFCE map length does not match expected voxel count", call. = FALSE)
      }
      tfce_out[, 1L, k] <- .posthoc_pack_from_full(tfce_full, spec)

      att_key <- sprintf("posthoc/nt:tfce_fwer/subject=1/contrast=%d", k)
      att_params <- tfce_res$params
      if (!is.list(att_params)) att_params <- list()
      attachments[[att_key]] <- list(
        type = "neurothresh_tfce_fwer",
        method = as.character(tfce_res$method %||% "tfce_fwer"),
        contrast_index = as.integer(k),
        threshold = as.numeric(tfce_res$threshold %||% NA_real_),
        n_sig = as.integer(sum(sig_out[, 1L, k] > 0.5, na.rm = TRUE)),
        null_model = if (!is.null(null_fun_k)) "explicit_null_fun" else "voxel_signflip_fallback",
        tail = tail,
        params = c(att_params, list(alpha = alpha, n_perm = n_perm, H = H, E = E, dh = dh))
      )
    }

    list(
      arrays = list(
        q = q_out,
        sig_mask = sig_out,
        tfce = tfce_out
      ),
      attachments = attachments
    )
  }
}

.posthoc_cluster_q_from_labels <- function(cluster_labels, cluster_ids, cluster_q, n_full) {
  if (length(cluster_labels) != n_full) {
    stop("cluster label vector length does not match voxel dimensions", call. = FALSE)
  }
  q_full <- rep.int(1, n_full)
  if (!length(cluster_ids)) return(q_full)
  for (i in seq_along(cluster_ids)) {
    vox <- which(cluster_labels == cluster_ids[i])
    if (!length(vox)) next
    q_full[vox] <- pmin(q_full[vox], as.numeric(cluster_q[i]))
  }
  q_full
}

.posthoc_neurothresh_cluster_fdr_perm <- function() {
  function(arrays, opts) {
    if (!requireNamespace("neurothresh", quietly = TRUE)) {
      stop("posthoc method 'nt:cluster_fdr_perm' requires optional package 'neurothresh'", call. = FALSE)
    }

    z <- arrays$z
    if (is.null(z)) stop("nt:cluster_fdr_perm requires z scores", call. = FALSE)
    dims <- dim(z)
    if (length(dims) != 3L) {
      stop("z must be a 3D assay [sample x subject x contrast]", call. = FALSE)
    }
    if (dims[2L] != 1L) {
      stop(
        "nt:cluster_fdr_perm currently expects group-level arrays with subject dimension = 1. ",
        "Run reduce() before this post-hoc method.",
        call. = FALSE
      )
    }

    spec <- .posthoc_resolve_voxel_spec(opts, n_samples = dims[1L], method_name = "nt:cluster_fdr_perm")

    q_level <- as.numeric(opts$q %||% opts$alpha %||% 0.05)
    if (!is.finite(q_level) || q_level <= 0 || q_level >= 1) stop("options$q must be in (0, 1)", call. = FALSE)
    cluster_thresh <- as.numeric(opts$cluster_thresh %||% 3.0)
    if (!is.finite(cluster_thresh) || cluster_thresh <= 0) {
      stop("options$cluster_thresh must be > 0", call. = FALSE)
    }
    n_perm <- as.integer(opts$n_perm %||% 5000L)
    if (!is.finite(n_perm) || n_perm < 1L) stop("options$n_perm must be a positive integer", call. = FALSE)
    tail <- match.arg(opts$tail %||% "pos", c("pos", "neg", "two"))
    two_sided_policy <- match.arg(opts$two_sided_policy %||% "BH_all", c("BH_all", "split_q"))
    allow_voxel_signflip <- isTRUE(opts$allow_voxel_signflip)
    null_spec <- opts$null_fun %||% NULL

    q_out <- array(NA_real_, dim = dims)
    sig_out <- array(0, dim = dims)
    attachments <- list()

    for (k in seq_len(dims[3L])) {
      z_vec <- as.numeric(z[, 1L, k])
      if (any(!is.finite(z_vec))) {
        stop("nt:cluster_fdr_perm requires finite z-values in active voxels (contrast ", k, ")", call. = FALSE)
      }
      stat_full <- .posthoc_unpack_to_full(z_vec, spec, fill = NA_real_)
      stat_vol <- array(stat_full, dim = spec$dim)

      null_fun_k <- .posthoc_pick_null_fun(null_spec, contrast_idx = k, n_contrast = dims[3L])
      if (is.null(null_fun_k) && !allow_voxel_signflip) {
        stop(
          "nt:cluster_fdr_perm requires options$null_fun unless allow_voxel_signflip = TRUE.",
          call. = FALSE
        )
      }

      seed_k <- NULL
      if (!is.null(opts$seed)) {
        seed_raw <- opts$seed
        if (length(seed_raw) == 1L) {
          seed_k <- as.integer(seed_raw) + (k - 1L)
        } else if (length(seed_raw) == dims[3L]) {
          seed_k <- as.integer(seed_raw[[k]])
        } else {
          stop("options$seed must be scalar or length equal to number of contrasts", call. = FALSE)
        }
      }

      if (tail %in% c("pos", "neg")) {
        c_res <- neurothresh::cluster_fdr(
          stat_vol = stat_vol,
          mask = spec$mask,
          cluster_thresh = cluster_thresh,
          q = q_level,
          p_method = "perm",
          n_perm = n_perm,
          null_fun = null_fun_k,
          allow_voxel_signflip = allow_voxel_signflip,
          tail = tail,
          seed = seed_k
        )

        c_tab <- c_res$table
        if (!is.data.frame(c_tab) || nrow(c_tab) == 0L) {
          q_full <- rep.int(1, spec$n_full)
          c_tab <- data.frame()
        } else {
          c_tab$q_cluster <- stats::p.adjust(c_tab$p_unc, method = "BH")
          q_full <- .posthoc_cluster_q_from_labels(
            cluster_labels = as.integer(c_res$cluster_labels),
            cluster_ids = as.integer(c_tab$cluster_id),
            cluster_q = as.numeric(c_tab$q_cluster),
            n_full = spec$n_full
          )
        }
        sig_full <- as.numeric(as.logical(c_res$sig_mask))
        if (length(sig_full) != spec$n_full) {
          stop("neurothresh::cluster_fdr returned sig_mask with unexpected length", call. = FALSE)
        }
        result_table <- c_tab
      } else {
        q_tail <- if (identical(two_sided_policy, "split_q")) q_level / 2 else q_level
        pos <- neurothresh::cluster_fdr(
          stat_vol = stat_vol,
          mask = spec$mask,
          cluster_thresh = cluster_thresh,
          q = q_tail,
          p_method = "perm",
          n_perm = n_perm,
          null_fun = null_fun_k,
          allow_voxel_signflip = allow_voxel_signflip,
          tail = "pos",
          seed = seed_k
        )
        neg <- neurothresh::cluster_fdr(
          stat_vol = stat_vol,
          mask = spec$mask,
          cluster_thresh = cluster_thresh,
          q = q_tail,
          p_method = "perm",
          n_perm = n_perm,
          null_fun = null_fun_k,
          allow_voxel_signflip = allow_voxel_signflip,
          tail = "neg",
          seed = if (is.null(seed_k)) NULL else seed_k + 100000L
        )

        pos_tab <- pos$table
        neg_tab <- neg$table
        if (!is.data.frame(pos_tab)) pos_tab <- data.frame()
        if (!is.data.frame(neg_tab)) neg_tab <- data.frame()

        if (nrow(pos_tab) + nrow(neg_tab) > 0L) {
          if (identical(two_sided_policy, "split_q")) {
            if (nrow(pos_tab) > 0L) pos_tab$q_cluster <- stats::p.adjust(pos_tab$p_unc, method = "BH")
            if (nrow(neg_tab) > 0L) neg_tab$q_cluster <- stats::p.adjust(neg_tab$p_unc, method = "BH")
          } else {
            p_all <- c(pos_tab$p_unc %||% numeric(0), neg_tab$p_unc %||% numeric(0))
            q_all <- stats::p.adjust(p_all, method = "BH")
            if (nrow(pos_tab) > 0L) pos_tab$q_cluster <- q_all[seq_len(nrow(pos_tab))]
            if (nrow(neg_tab) > 0L) {
              neg_tab$q_cluster <- q_all[seq.int(nrow(pos_tab) + 1L, length(q_all))]
            }
          }
        } else {
          pos_tab$q_cluster <- numeric(0)
          neg_tab$q_cluster <- numeric(0)
        }

        q_full <- rep.int(1, spec$n_full)
        if (nrow(pos_tab) > 0L) {
          q_pos <- .posthoc_cluster_q_from_labels(
            cluster_labels = as.integer(pos$cluster_labels),
            cluster_ids = as.integer(pos_tab$cluster_id),
            cluster_q = as.numeric(pos_tab$q_cluster),
            n_full = spec$n_full
          )
          q_full <- pmin(q_full, q_pos)
        }
        if (nrow(neg_tab) > 0L) {
          q_neg <- .posthoc_cluster_q_from_labels(
            cluster_labels = as.integer(neg$cluster_labels),
            cluster_ids = as.integer(neg_tab$cluster_id),
            cluster_q = as.numeric(neg_tab$q_cluster),
            n_full = spec$n_full
          )
          q_full <- pmin(q_full, q_neg)
        }

        sig_full <- as.numeric(as.logical(pos$sig_mask) | as.logical(neg$sig_mask))
        if (nrow(pos_tab) == 0L && nrow(neg_tab) == 0L) {
          result_table <- data.frame()
        } else if (nrow(pos_tab) == 0L) {
          result_table <- neg_tab
        } else if (nrow(neg_tab) == 0L) {
          result_table <- pos_tab
        } else {
          result_table <- rbind(pos_tab, neg_tab)
        }
      }

      q_out[, 1L, k] <- .posthoc_pack_from_full(q_full, spec)
      sig_out[, 1L, k] <- .posthoc_pack_from_full(sig_full, spec)

      att_key <- sprintf("posthoc/nt:cluster_fdr_perm/subject=1/contrast=%d", k)
      attachments[[att_key]] <- list(
        type = "neurothresh_cluster_fdr_perm",
        contrast_index = as.integer(k),
        n_sig_vox = as.integer(sum(sig_out[, 1L, k] > 0.5, na.rm = TRUE)),
        n_sig_clusters = if (!is.data.frame(result_table) || !("q_cluster" %in% names(result_table))) {
          0L
        } else {
          as.integer(sum(result_table$q_cluster <= q_level, na.rm = TRUE))
        },
        null_model = if (!is.null(null_fun_k)) "explicit_null_fun" else "voxel_signflip_fallback",
        tail = tail,
        two_sided_policy = if (tail == "two") two_sided_policy else NA_character_,
        params = list(q = q_level, cluster_thresh = cluster_thresh, n_perm = n_perm, tail = tail, p_method = "perm"),
        table = result_table
      )
    }

    list(
      arrays = list(
        q = q_out,
        sig_mask = sig_out
      ),
      attachments = attachments
    )
  }
}

# nocov end
