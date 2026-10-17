# SABHA-style spatial FDR (fdr:sabha) -----------------------------------------
#
# Voxel-level structure-adaptive weighted Benjamini-Hochberg with groupwise
# Storey pi0 weights and an Efron central-matching empirical null. Ported from
# fmrireg::spatial_fdr() (review 2026-10 contract) and numerically identical to
# it. This is distinct from `fdr:spatial`, which aggregates p-values within
# parcels by Simes' method and returns parcel-level q-values.

#' Structure-adaptive spatial FDR (SABHA)
#'
#' Voxel-level spatially aware FDR control by a structure-adaptive weighted
#' Benjamini-Hochberg procedure. Each feature belongs to a spatial group
#' (parcel, block, cluster); the proportion of nulls \eqn{\pi_0} is estimated
#' per group (Storey), optionally smoothed over neighbouring groups, and the
#' weighted BH procedure of Genovese, Roeder & Wasserman (2006) is applied with
#' weights \eqn{1/\pi_0} normalised to sum to the number of finite p-values:
#' reject \eqn{p_i \le w_i \alpha k / n}. Unlike [posthoc()] method
#' `"fdr:spatial"` (Simes aggregation, one q-value per parcel), every feature
#' keeps its own q-value. Results are identical to `fmrireg::spatial_fdr()`.
#'
#' **Group numbering.** Groups are numbered 1..G by sorted unique value, or by
#' factor level (unused levels included) when `group` is a factor. This
#' numbering indexes `neighbors`, `pi0_raw`, `pi0_smooth` and `m_g`. Features
#' with `NA` group are excluded from \eqn{\pi_0} estimation and receive weight
#' 1 before normalisation.
#'
#' **Empirical null.** When `z` is supplied with `empirical_null = TRUE`,
#' \eqn{\mu_0} and \eqn{\sigma_0} are estimated by Efron's (2004) central
#' matching: a quadratic is fitted by Poisson regression to the log histogram
#' counts of z within \eqn{\pm 1.5} robust SDs of the median. With fewer than
#' 500 finite z-values (or if the fit fails) the median and normal-consistent
#' MAD of all finite z-values are used instead; with fewer than 10 the
#' theoretical null is used. p-values are then two-sided normal p-values of
#' \eqn{(z - \mu_0)/\sigma_0}.
#'
#' @param z Numeric vector (or features x maps matrix) of z-values. Supply
#'   either `z` or `p`.
#' @param p Numeric vector (or matrix) of two-sided p-values.
#' @param group Integer, character or factor vector of group ids, one per
#'   feature (row).
#' @param alpha FDR level in (0, 1); sets `reject`, `threshold` and `k`
#'   (q-values do not depend on it).
#' @param tau Storey threshold in (0, 1) for \eqn{\pi_0} estimation.
#' @param lambda Non-negative smoothing strength of \eqn{\pi_0} across
#'   neighbours (one step of \eqn{(\pi_g + \lambda \bar\pi_{N(g)})/(1+\lambda)}).
#' @param neighbors Optional list of length G; element g holds the 1-based ids
#'   (in the numbering above) of the neighbours of group g.
#' @param min_pi0 Lower bound in \[0, 1\] for the smoothed \eqn{\pi_0}.
#' @param empirical_null Logical; estimate \eqn{\mu_0, \sigma_0} from `z`.
#'
#' @return For vector input, a list with `reject`, `q`, `p` (the p-values
#'   tested), `weights` (normalised), `pi0_raw`, `pi0_smooth`, `threshold`,
#'   `k`, `mu0`, `sigma0`, `groups` (1..G ids), `G`, `m_g` (features per
#'   group), `alpha`, `tau`, `lambda`. For matrix input every column is
#'   processed independently: `reject`, `q`, `p`, `weights` are matrices,
#'   `pi0_raw`/`pi0_smooth` are G x columns, and `threshold`, `k`, `mu0`,
#'   `sigma0` are per-column vectors.
#' @references
#' Benjamini & Hochberg (1995); Storey (2002); Genovese, Roeder & Wasserman
#' (2006) False discovery control with p-value weighting; Efron (2004)
#' Large-scale simultaneous hypothesis testing: the choice of a null
#' hypothesis.
#' @seealso [posthoc-methods] (`"fdr:sabha"`), [posthoc()]
#' @export
#' @examples
#' set.seed(1)
#' grp <- rep(1:10, each = 100)
#' z <- rnorm(1000)
#' z[grp <= 2] <- z[grp <= 2] + 3
#' res <- spatial_fdr_sabha(z = z, group = grp, alpha = 0.05)
#' sum(res$reject)
#' res$pi0_smooth
spatial_fdr_sabha <- function(z = NULL,
                              p = NULL,
                              group,
                              alpha = 0.05,
                              tau = 0.5,
                              lambda = 1.0,
                              neighbors = NULL,
                              min_pi0 = 0.05,
                              empirical_null = TRUE) {
  if (is.null(z) && is.null(p)) stop("Must provide either z or p values", call. = FALSE)
  if (!is.null(z) && !is.null(p)) stop("Provide either z or p, not both", call. = FALSE)
  x <- z %||% p
  if (is.matrix(x)) {
    cols <- lapply(seq_len(ncol(x)), function(j) {
      args <- if (is.null(z)) list(p = x[, j]) else list(z = x[, j])
      do.call(spatial_fdr_sabha, c(args, list(
        group = group, alpha = alpha, tau = tau, lambda = lambda,
        neighbors = neighbors, min_pi0 = min_pi0, empirical_null = empirical_null
      )))
    })
    if (!length(cols)) stop("`z`/`p` matrix has no columns", call. = FALSE)
    bind <- function(f) {
      out <- do.call(cbind, lapply(cols, `[[`, f))
      colnames(out) <- colnames(x)
      out
    }
    vec <- function(f) stats::setNames(vapply(cols, `[[`, numeric(1), f), colnames(x))
    first <- cols[[1L]]
    return(list(
      reject = bind("reject"), q = bind("q"), p = bind("p"), weights = bind("weights"),
      pi0_raw = bind("pi0_raw"), pi0_smooth = bind("pi0_smooth"),
      threshold = vec("threshold"), k = vec("k"), mu0 = vec("mu0"), sigma0 = vec("sigma0"),
      groups = first$groups, G = first$G, m_g = first$m_g,
      alpha = alpha, tau = tau, lambda = first$lambda
    ))
  }

  m <- length(x)
  if (length(group) != m) stop("group length must match z/p length", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)", call. = FALSE)
  if (tau <= 0 || tau >= 1) stop("tau must be in (0, 1)", call. = FALSE)
  if (lambda < 0) {
    warning("lambda < 0 not recommended; setting to 0", call. = FALSE)
    lambda <- 0
  }
  if (min_pi0 < 0 || min_pi0 > 1) stop("min_pi0 must be in [0, 1]", call. = FALSE)

  # Authoritative group numbering 1..G: factor levels (unused levels kept, so a
  # neighbors list indexed by level stays aligned), otherwise sorted unique
  # values. pi0, smoothing, weights and `groups` all use it.
  gfac <- if (is.factor(group)) group else factor(group)
  G <- nlevels(gfac)
  gid <- as.integer(gfac)
  if (G == 0L) stop("No valid groups found", call. = FALSE)

  mu0 <- 0
  sigma0 <- 1
  if (is.null(p)) {
    zz <- as.numeric(z)
    if (isTRUE(empirical_null)) {
      en <- .sabha_empirical_null(zz)
      mu0 <- en$mu0
      sigma0 <- en$sigma0
    }
    pv <- 2 * stats::pnorm(abs((zz - mu0) / sigma0), lower.tail = FALSE)
  } else {
    pv <- pmin(pmax(as.numeric(p), 0), 1)
  }

  gstats <- sabha_group_pi0_counts_cpp(pv, gid, tau, G)
  m_g <- as.integer(gstats$m_g)
  pi0_raw <- as.numeric(gstats$pi0_raw)
  pi0_raw[!is.finite(pi0_raw)] <- 1

  if (!is.null(neighbors)) {
    if (length(neighbors) != G) {
      stop("neighbors must be a list of length G (", G, ")", call. = FALSE)
    }
    pi0_sm <- sabha_pi0_smooth_cpp(pi0_raw, neighbors, lambda = lambda, iters = 1L)
  } else {
    pi0_sm <- pi0_raw
  }
  pi0_sm <- pmin(pmax(pi0_sm, min_pi0), 1.0)

  w <- (1 / pi0_sm)[gid]
  w[is.na(w)] <- 1
  wb <- sabha_weighted_bh_cpp(pv, w, alpha = alpha)
  q <- sabha_bh_qvalues_scaled_cpp(pv / wb$w_norm)

  list(
    reject = wb$reject,
    q = q,
    p = pv,
    weights = wb$w_norm,
    pi0_raw = pi0_raw,
    pi0_smooth = pi0_sm,
    threshold = wb$threshold,
    k = wb$k,
    mu0 = mu0,
    sigma0 = sigma0,
    groups = gid,
    G = G,
    m_g = m_g,
    alpha = alpha,
    tau = tau,
    lambda = lambda
  )
}

# Central matching (Efron 2004): fit log f(z) = b0 + b1 z + b2 z^2 by Poisson
# regression on histogram counts within +/- half_width robust SDs of the
# median, so sigma0 = sqrt(-1 / (2 b2)) and mu0 = b1 sigma0^2. Falls back to
# median / normal-consistent MAD of all finite z with fewer than `min_cm`
# values or an unusable fit, and to the theoretical null with fewer than 10.
.sabha_empirical_null <- function(z, half_width = 1.5, nbins = 40L, min_cm = 500L) {
  z <- z[is.finite(z)]
  n <- length(z)
  if (n < 10L) return(list(mu0 = 0, sigma0 = 1, method = "theoretical"))

  mu_r <- stats::median(z)
  s_r <- stats::mad(z, center = mu_r)
  if (!is.finite(s_r) || s_r <= 0) s_r <- stats::sd(z)
  if (!is.finite(s_r) || s_r <= 0) {
    warning("Could not estimate sigma0; using 1", call. = FALSE)
    return(list(mu0 = mu_r, sigma0 = 1, method = "median"))
  }

  if (n >= min_cm) {
    lo <- mu_r - half_width * s_r
    hi <- mu_r + half_width * s_r
    br <- seq(lo, hi, length.out = nbins + 1L)
    mids <- (br[-1L] + br[-length(br)]) / 2
    inwin <- z >= lo & z < hi
    y <- tabulate(findInterval(z[inwin], br), nbins)
    fit <- tryCatch(
      suppressWarnings(stats::glm.fit(cbind(1, mids, mids^2), y, family = stats::poisson())),
      error = function(e) NULL
    )
    if (!is.null(fit) && isTRUE(fit$converged)) {
      b <- unname(stats::coef(fit))
      if (all(is.finite(b)) && b[3] < 0) {
        s_cm <- sqrt(-1 / (2 * b[3]))
        mu_cm <- b[2] * s_cm^2
        if (is.finite(s_cm) && s_cm > 0.5 * s_r && s_cm < 2 * s_r &&
            abs(mu_cm - mu_r) < half_width * s_r) {
          return(list(mu0 = mu_cm, sigma0 = s_cm, method = "central matching"))
        }
      }
    }
  }
  list(mu0 = mu_r, sigma0 = s_r, method = "median/MAD")
}

# Resolve the statistic the post-hoc adjusts: named list (output q-name ->
# list(kind = "z"|"p", values = array)).
.sabha_sources <- function(arrays, opts) {
  statistic <- opts$statistic %||% "auto"
  if (!statistic %in% c("auto", "z", "p")) {
    stop("fdr:sabha: options$statistic must be \"auto\", \"z\" or \"p\".", call. = FALSE)
  }
  qname <- function(src) if (grepl("^[pz]_[A-Za-z0-9]+:", src)) sub("^[pz]_", "q_", src) else "q"
  wrap <- function(lst, kind) lapply(lst, function(a) list(kind = kind, values = a))
  src <- opts$source %||% NULL
  if (!is.null(src)) {
    if (is.null(arrays[[src]])) {
      stop("fdr:sabha: options$source assay not found: ", src, call. = FALSE)
    }
    kind <- if (grepl("^z($|_)", src)) "z" else "p"
    return(stats::setNames(wrap(list(arrays[[src]]), kind), qname(src)))
  }
  if (!identical(statistic, "p")) {
    if (!is.null(arrays[["z"]])) return(list(q = list(kind = "z", values = arrays[["z"]])))
    zc <- grep("^z_coef:", names(arrays), value = TRUE)
    if (length(zc)) {
      return(stats::setNames(wrap(lapply(zc, function(n) arrays[[n]]), "z"), qname(zc)))
    }
    if (identical(statistic, "z")) {
      stop("fdr:sabha: statistic = \"z\" but no z or z_coef:<term> assay is present.",
           call. = FALSE)
    }
  }
  wrap(.fdr_p_sources(arrays, opts), "p")
}

.posthoc_sabha_fdr <- function() {
  function(arrays, opts) {
    sources <- .sabha_sources(arrays, opts)
    dims <- dim(sources[[1L]]$values)
    group <- opts$group %||% .resolve_spatial_fdr_group(opts$row_data %||% NULL)
    if (is.null(group)) {
      stop(
        "fdr:sabha needs a spatial grouping: one parcel/block id per sample.\n",
        sprintf("  Supply options = list(group = <length-%d vector>) or attach a row_data ", dims[1L]),
        "column (feature_group, spatial_group, group, parcel, parcel_id).",
        call. = FALSE
      )
    }
    if (length(group) != dims[1L]) {
      stop(sprintf("length(options$group) (%d) must equal the sample dimension (%d).",
                   length(group), dims[1L]), call. = FALSE)
    }
    args <- list(
      group = group,
      alpha = opts$alpha %||% 0.05,
      tau = opts$tau %||% 0.5,
      lambda = opts$lambda %||% 1.0,
      neighbors = opts$neighbors %||% NULL,
      min_pi0 = opts$min_pi0 %||% 0.05,
      empirical_null = opts$empirical_null %||% TRUE
    )
    out <- list()
    attachments <- list()
    for (out_name in names(sources)) {
      src <- sources[[out_name]]
      q <- array(NA_real_, dim = dims)
      info <- list(
        type = "sabha_fdr", statistic = src$kind, alpha = args$alpha,
        mu0 = matrix(NA_real_, dims[2L], dims[3L]),
        sigma0 = matrix(NA_real_, dims[2L], dims[3L]),
        threshold = matrix(NA_real_, dims[2L], dims[3L]),
        k = matrix(NA_real_, dims[2L], dims[3L]),
        pi0_smooth = NULL
      )
      for (j in seq_len(dims[2L])) {
        for (k in seq_len(dims[3L])) {
          val <- list(src$values[, j, k])
          names(val) <- src$kind
          res <- do.call(spatial_fdr_sabha, c(val, args))
          q[, j, k] <- res$q
          info$mu0[j, k] <- res$mu0
          info$sigma0[j, k] <- res$sigma0
          info$threshold[j, k] <- res$threshold
          info$k[j, k] <- res$k
          if (is.null(info$pi0_smooth)) {
            info$pi0_smooth <- array(NA_real_, c(res$G, dims[2L], dims[3L]))
            info$groups <- if (is.factor(group)) levels(group) else sort(unique(group))
          }
          info$pi0_smooth[, j, k] <- res$pi0_smooth
        }
      }
      out[[out_name]] <- q
      attachments[[paste0("posthoc/fdr:sabha/", out_name)]] <- info
    }
    list(arrays = out, attachments = attachments)
  }
}
