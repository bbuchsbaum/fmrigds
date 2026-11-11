# Spatial (TFCE / cluster) sign-flip inference for perm:onesample ------------

.perm_spatial_choices <- c("none", "tfce", "cluster")

# Normalise and validate the spatial options of perm:onesample. Returns NULL
# when no spatial statistic was requested.
.perm_spatial_options <- function(opts, n_subject) {
  spatial <- unique(as.character(opts$spatial %||% "none"))
  bad <- setdiff(spatial, .perm_spatial_choices)
  if (length(bad)) {
    stop(sprintf(
      "perm:onesample: unknown options$spatial value(s) %s; use \"tfce\" and/or \"cluster\".",
      paste(sprintf("'%s'", bad), collapse = ", ")
    ), call. = FALSE)
  }
  spatial <- setdiff(spatial, "none")
  if (!length(spatial)) return(NULL)

  num1 <- function(x, name, positive = TRUE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || (positive && x <= 0)) {
      stop(sprintf("perm:onesample: options$%s must be one %snumber.",
                   name, if (positive) "positive " else "finite "), call. = FALSE)
    }
    as.numeric(x)
  }
  stat <- opts$spatial_stat %||% "z"
  if (!(length(stat) == 1L && stat %in% c("z", "t"))) {
    stop("perm:onesample: options$spatial_stat must be \"z\" or \"t\".", call. = FALSE)
  }
  connectivity <- suppressWarnings(as.integer(opts$connectivity %||% 6L))
  if (length(connectivity) != 1L || is.na(connectivity) ||
      !connectivity %in% c(6L, 18L, 26L)) {
    stop("perm:onesample: options$connectivity must be 6, 18 or 26.", call. = FALSE)
  }
  cluster_stat <- opts$cluster_stat %||% "mass"
  if (!(length(cluster_stat) == 1L && cluster_stat %in% c("mass", "extent"))) {
    stop("perm:onesample: options$cluster_stat must be \"mass\" or \"extent\".", call. = FALSE)
  }
  threshold <- NA_real_
  cluster_p <- NA_real_
  if ("cluster" %in% spatial) {
    if (!is.null(opts$cluster_threshold) && !is.null(opts$cluster_p)) {
      stop("perm:onesample: give options$cluster_threshold or options$cluster_p, not both.",
           call. = FALSE)
    }
    if (!is.null(opts$cluster_threshold)) {
      threshold <- num1(opts$cluster_threshold, "cluster_threshold")
    } else {
      cluster_p <- num1(opts$cluster_p %||% 0.001, "cluster_p")
      if (cluster_p >= 0.5) {
        stop("perm:onesample: options$cluster_p must be below 0.5.", call. = FALSE)
      }
      threshold <- if (identical(stat, "z")) {
        stats::qnorm(cluster_p, lower.tail = FALSE)
      } else {
        stats::qt(cluster_p, df = n_subject - 1, lower.tail = FALSE)
      }
    }
  }
  threads <- opts$threads %||% 0L
  threads <- suppressWarnings(as.integer(threads))
  if (length(threads) != 1L || is.na(threads) || threads < 0L) {
    stop("perm:onesample: options$threads must be a non-negative integer.", call. = FALSE)
  }
  list(
    tfce = "tfce" %in% spatial,
    cluster = "cluster" %in% spatial,
    stat = stat,
    H = num1(opts$tfce_H %||% 2, "tfce_H"),
    E = num1(opts$tfce_E %||% 0.5, "tfce_E"),
    dh = num1(opts$tfce_dh %||% 0.1, "tfce_dh"),
    connectivity = connectivity,
    cluster_stat = cluster_stat,
    cluster_threshold = threshold,
    cluster_p = cluster_p,
    threads = threads
  )
}

# Sample axis -> voxel grid for a voxel space.
.perm_spatial_grid <- function(space, n_samples) {
  if (is.null(space)) {
    stop("perm:onesample with options$spatial needs the voxel space of the data; ",
         "run it through compute() on a voxel-space GDS.", call. = FALSE)
  }
  if (!inherits(space, "space_voxel")) {
    stop("perm:onesample with options$spatial requires a voxel space; got ",
         paste(class(space), collapse = "/"), ".", call. = FALSE)
  }
  .posthoc_resolve_voxel_spec(list(space = space), n_samples,
                              method_name = "perm:onesample (spatial)")
}

# Max-statistic FWER p-values, (1 + #{null >= obs}) / (1 + n_null).
.perm_maxnull_p <- function(obs, null) {
  B <- length(null)
  below <- findInterval(obs, sort(null), left.open = TRUE)
  p <- (1 + B - below) / (B + 1)
  p[is.na(obs)] <- NA_real_
  p
}

# Cluster table plus per-sample cluster id and cluster FWER p.
.perm_spatial_clusters <- function(res, sp, grid, null_mass, null_extent, space, n_samples) {
  parts <- list()
  label_full <- integer(n_samples)
  for (side in c("pos", "neg")) {
    cl <- res[[paste0("clusters_", side)]]
    if (!length(cl) || !length(cl$sizes)) next
    sgn <- if (side == "pos") 1L else -1L
    parts[[side]] <- list(
      sign = sgn, labels = cl$labels, sizes = cl$sizes, mass = cl$mass, peak = cl$peak
    )
  }
  tab <- data.frame(
    cluster_id = integer(0), sign = integer(0), size = integer(0),
    mass = numeric(0), p_fwer_mass = numeric(0), p_fwer_extent = numeric(0),
    peak_stat = numeric(0), peak_t = numeric(0), peak_sample = integer(0),
    peak_i = integer(0), peak_j = integer(0), peak_k = integer(0),
    peak_x = numeric(0), peak_y = numeric(0), peak_z = numeric(0)
  )
  p_vox <- rep(1, n_samples)
  if (!length(parts)) {
    return(list(table = tab, cluster_id = label_full, p = p_vox))
  }
  rows <- lapply(parts, function(pt) {
    data.frame(
      sign = rep(pt$sign, length(pt$sizes)),
      local = seq_along(pt$sizes),
      size = pt$sizes,
      mass = pt$mass,
      peak_sample = pt$peak
    )
  })
  tab0 <- do.call(rbind, rows)
  side_of <- ifelse(tab0$sign > 0, "pos", "neg")
  # Largest mass first; ties by extent then by position.
  ord <- order(-tab0$mass, -tab0$size, tab0$peak_sample)
  tab0 <- tab0[ord, , drop = FALSE]
  side_of <- side_of[ord]
  tab0$cluster_id <- seq_len(nrow(tab0))
  tab0$p_fwer_mass <- .perm_maxnull_p(tab0$mass, null_mass)
  tab0$p_fwer_extent <- .perm_maxnull_p(tab0$size, null_extent)
  p_use <- if (identical(sp$cluster_stat, "mass")) tab0$p_fwer_mass else tab0$p_fwer_extent
  for (r in seq_len(nrow(tab0))) {
    lab <- parts[[side_of[r]]]$labels
    hit <- which(lab == tab0$local[r])
    label_full[hit] <- tab0$cluster_id[r]
    p_vox[hit] <- p_use[r]
  }
  ijk <- arrayInd(grid$mask_idx[tab0$peak_sample], grid$dim)
  aff <- space$affine
  xyz <- if (is.matrix(aff) && all(dim(aff) == c(4L, 4L))) {
    t(aff %*% rbind(t(ijk - 1), 1))[, 1:3, drop = FALSE]
  } else {
    matrix(NA_real_, nrow(ijk), 3)
  }
  tab <- data.frame(
    cluster_id = tab0$cluster_id,
    sign = as.integer(tab0$sign),
    size = as.integer(tab0$size),
    mass = tab0$mass,
    p_fwer_mass = tab0$p_fwer_mass,
    p_fwer_extent = tab0$p_fwer_extent,
    peak_stat = res$stat[tab0$peak_sample],
    peak_t = res$t_g[tab0$peak_sample],
    peak_sample = as.integer(tab0$peak_sample),
    peak_i = as.integer(ijk[, 1]),
    peak_j = as.integer(ijk[, 2]),
    peak_k = as.integer(ijk[, 3]),
    peak_x = xyz[, 1],
    peak_y = xyz[, 2],
    peak_z = xyz[, 3]
  )
  rownames(tab) <- NULL
  list(table = tab, cluster_id = label_full, p = p_vox)
}

# Subject-level sign-flip inference with TFCE and/or cluster FWER.
.perm_onesample_spatial <- function(beta, signs, skip_first, opts, sp) {
  n_samples <- ncol(beta)
  space <- opts$space
  grid <- .perm_spatial_grid(space, n_samples)
  tail <- .perm_tail_code(opts$alternative %||% "two.sided")
  min_subj <- as.integer(opts$min_subjects %||% 2L)
  res <- perm_onesample_spatial_cpp(
    beta,
    signs,
    dims = as.integer(grid$dim),
    voxels = as.integer(grid$mask_idx),
    connectivity = sp$connectivity,
    tail = tail,
    min_subj = min_subj,
    skip_first = skip_first,
    do_tfce = sp$tfce,
    do_cluster = sp$cluster,
    H = sp$H, E = sp$E, dh = sp$dh,
    cluster_threshold = if (sp$cluster) sp$cluster_threshold else 0,
    stat_type = if (identical(sp$stat, "z")) 1L else 0L,
    threads = sp$threads
  )
  keep <- if (skip_first) -1L else seq_len(nrow(signs))
  null_t <- res$max_t_null[keep]
  excluded <- !is.finite(res$t_g)

  out <- res[c("beta_g", "se_g", "t_g", "df", "p_perm", "p_fwer")]
  out$z_g <- if (identical(sp$stat, "z")) res$stat else .z_from_t(res$t_g, res$df)
  out$n_eff <- as.numeric(colSums(is.finite(beta)))

  null <- list(max_t = null_t)
  if (sp$tfce) {
    null_tfce <- res$max_tfce_null[keep]
    score <- pmax(res$tfce_pos, res$tfce_neg)
    tfce <- res$tfce_pos - res$tfce_neg
    tfce[excluded] <- NA_real_
    score[excluded] <- NA_real_
    out$tfce <- tfce
    out$p_tfce_fwer <- .perm_maxnull_p(score, null_tfce)
    null$tfce <- null_tfce
  }
  clusters <- NULL
  if (sp$cluster) {
    null_mass <- res$max_mass_null[keep]
    null_extent <- res$max_extent_null[keep]
    cl <- .perm_spatial_clusters(res, sp, grid, null_mass, null_extent, space, n_samples)
    cid <- as.numeric(cl$cluster_id)
    pc <- cl$p
    cid[excluded] <- NA_real_
    pc[excluded] <- NA_real_
    out$cluster_id <- cid
    out$p_cluster_fwer <- pc
    null$cluster_mass <- null_mass
    null$cluster_extent <- null_extent
    clusters <- cl$table
  }
  out$attachment <- list(
    type = "perm_spatial",
    reducer = "perm:onesample",
    alternative = opts$alternative %||% "two.sided",
    n_perm = nrow(signs),
    n_null = length(null_t),
    statistic = sp$stat,
    connectivity = sp$connectivity,
    tfce = if (sp$tfce) list(H = sp$H, E = sp$E, dh = sp$dh) else NULL,
    cluster = if (sp$cluster) {
      list(threshold = sp$cluster_threshold, cluster_p = sp$cluster_p,
           statistic = sp$cluster_stat)
    } else {
      NULL
    },
    min_subjects = min_subj,
    n_excluded = sum(excluded),
    seed = opts$seed %||% NULL,
    null = null,
    clusters = clusters
  )
  out
}
