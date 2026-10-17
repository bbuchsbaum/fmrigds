# fdr:sabha / spatial_fdr_sabha(): port of fmrireg::spatial_fdr() (SABHA-style
# weighted BH, Efron central-matching empirical null, sorted-label group
# numbering). Reference values below were produced by fmrireg's
# spatial_fdr() (branch claude/simplify-retire-cruft, 27d1543) on exactly
# these inputs; the port is identical to it field by field.

sabha_case_inputs <- function() {
  set.seed(20261017)
  # A: central matching path (n >= 500), chain neighbours, smoothing
  gA <- rep(1:20, each = 100)
  zA <- rnorm(2000)
  zA[gA <= 3] <- zA[gA <= 3] + 3
  nbA <- lapply(1:20, function(g) setdiff(c(g - 1, g + 1), c(0, 21)))
  # B: p input, NA p, NA group, factor with an unused level
  pB <- c(rbeta(60, 0.2, 5), runif(140))
  pB[c(3, 50, 177)] <- NA
  gB <- factor(rep(c("a", "c", "d", "e"), each = 50), levels = c("a", "b", "c", "d", "e"))
  gB[c(10, 120)] <- NA
  nbB <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
  # C: small n (median/MAD fallback), character labels, shifted null
  zC <- rnorm(200, 0.3, 1.3)
  zC[1:20] <- zC[1:20] + 7
  gC <- sprintf("p%02d", rep(c(7, 2, 11, 5), each = 50))
  list(zA = zA, gA = gA, nbA = nbA, pB = pB, gB = gB, nbB = nbB, zC = zC, gC = gC)
}

test_that("spatial_fdr_sabha reproduces fmrireg::spatial_fdr() reference values", {
  I <- sabha_case_inputs()
  check <- function(g, mu0, sigma0, k, thr, sumq, sumrej, pi0, idx, q) {
    expect_equal(g$mu0, mu0, tolerance = 1e-10)
    expect_equal(g$sigma0, sigma0, tolerance = 1e-10)
    expect_identical(as.integer(g$k), as.integer(k))
    expect_equal(g$threshold, thr, tolerance = 1e-10)
    expect_equal(sum(g$q, na.rm = TRUE), sumq, tolerance = 1e-10)
    expect_identical(sum(which(g$reject)), as.integer(sumrej))
    expect_equal(g$pi0_smooth, pi0, tolerance = 1e-10)
    expect_equal(g$q[idx], q, tolerance = 1e-10)
  }
  A <- spatial_fdr_sabha(z = I$zA, group = I$gA, neighbors = I$nbA, lambda = 1, alpha = 0.05)
  check(A, 0.026116920444394138, 1.0882101848346413, 173, 0.0043249999999999999,
        1704.8425973417411, 21942,
        c(0.05, 0.05, 0.275, 0.76, 1, 1, 1, 1, 1, 1, 0.995, 0.99, 0.995, 1, 1, 1,
          0.98, 0.96, 0.98, 1),
        c(1, 287, 572, 858, 1143, 1429, 1714, 2000),
        c(0.028879194857822286, 0.058150449184670841, 1, 1, 1, 1, 1, 1))
  B <- spatial_fdr_sabha(p = I$pB, group = I$gB, neighbors = I$nbB, lambda = 0.5,
                         tau = 0.4, min_pi0 = 0.1, alpha = 0.1, empirical_null = FALSE)
  # Case B has missing p-values: Storey group sizes count finite tests only
  # (fmrireg's spatial_fdr() counted missing ones until the same fix).
  check(B, 0, 1, 43, 0.021827411167512689, 119.43344461373988, 1316,
        c(0.33333333333333331, 0.80555555555555547, 0.85827664399092962,
          0.84750566893424040, 0.92970521541950113),
        c(1, 29, 58, 86, 115, 143, 172, 200),
        c(0.059278982187982152, 0.454274904599392737, 0.424392815888882724,
          0.915841462519309535, 0.333922384510125847, 1, 0.822430146986976029,
          0.278076269382573849))
  expect_identical(B$m_g, c(47L, 0L, 50L, 49L, 49L))
  C <- spatial_fdr_sabha(z = I$zC, group = I$gC, alpha = 0.05)
  check(C, 0.68492872593075682, 1.6960777500164816, 20, 0.005000000000000001,
        173.44465975698444, 210, c(1, 1, 0.56, 1),
        c(1, 29, 58, 86, 115, 143, 172, 200),
        c(0.007701367236807059, 0.96053211274840455, 1, 0.90971783918949967, 1, 1, 1, 1))
  # character labels are numbered by sorted label: p02 < p05 < p07 < p11
  expect_identical(C$groups, as.integer(factor(I$gC)))
})

# ---------------------------------------------------------------------------
# Ported from fmrireg tests (test-review-A-spatial-fdr.R, test_spatial_fdr.R)
# ---------------------------------------------------------------------------

wbh_reference <- function(p, w, alpha) {
  valid <- is.finite(p)
  n <- sum(valid)
  wn <- numeric(length(p))
  wn[valid] <- w[valid] * n / sum(w[valid])
  ps <- rep(NA_real_, length(p))
  ps[valid] <- p[valid] / wn[valid]
  srt <- sort(ps[valid])
  ok <- which(srt <= alpha * seq_len(n) / n)
  k <- if (length(ok)) max(ok) else 0L
  thr <- if (k > 0) alpha * k / n else 0
  reject <- valid & !is.na(ps) & ps <= thr & k > 0
  q <- rep(NA_real_, length(p))
  q[valid] <- stats::p.adjust(pmin(ps[valid], 1), method = "BH")
  list(reject = reject, q = q, w_norm = wn, k = k, threshold = thr)
}

pi0_reference <- function(p, group, tau) {
  f <- if (is.factor(group)) group else factor(group)
  tail <- tapply(is.finite(p) & p > tau, f, sum)
  m <- tapply(is.finite(p), f, sum)
  est <- as.numeric(tail) / ((1 - tau) * as.numeric(m))
  est[is.na(est)] <- 1
  pmin(pmax(est, 0), 1)
}

test_that("empirical null is consistent under N(0,1) and with 10% signal", {
  set.seed(20261011)
  m <- 1e5
  grp <- rep(seq_len(100), each = m / 100)
  res <- spatial_fdr_sabha(z = rnorm(m), group = grp)
  expect_gte(res$sigma0, 0.95)
  expect_lte(res$sigma0, 1.05)
  expect_lt(abs(res$mu0), 0.05)
  expect_lte(sum(res$reject), 5)

  set.seed(42)
  is_sig <- grp <= 10
  z <- rnorm(m)
  z[is_sig] <- rnorm(sum(is_sig), mean = 4)
  res <- spatial_fdr_sabha(z = z, group = grp)
  oracle <- spatial_fdr_sabha(z = z, group = grp, empirical_null = FALSE)
  expect_gte(res$sigma0, 0.95)
  expect_lte(res$sigma0, 1.05)
  expect_gt(mean(res$reject[is_sig]), mean(oracle$reject[is_sig]) - 0.03)
  expect_lt(sum(res$reject & !is_sig) / max(1, sum(res$reject)), 0.08)

  set.seed(111)
  res <- spatial_fdr_sabha(z = rnorm(5e4, 0.5, 1.2), group = rep(1:50, each = 1000))
  expect_equal(res$mu0, 0.5, tolerance = 0.1)
  expect_equal(res$sigma0, 1.2, tolerance = 0.05)
})

test_that("pi0 and weights follow sorted group labels, not first appearance", {
  set.seed(1)
  n <- 200
  group <- c(rep(2, n), rep(1, n))
  p <- c(runif(n), rbeta(n, 0.1, 5))
  res <- spatial_fdr_sabha(p = p, group = group, lambda = 0, empirical_null = FALSE)
  expect_identical(as.integer(res$groups), as.integer(factor(group)))
  expect_equal(res$pi0_raw, pi0_reference(p, group, 0.5), tolerance = 1e-12)
  expect_gt(mean(res$weights[group == 1]), mean(res$weights[group == 2]))
  a <- 1 / res$pi0_smooth[as.integer(factor(group))]
  expect_equal(res$weights, a * length(p) / sum(a), tolerance = 1e-12)
})

test_that("results are invariant to relabelling groups and permuting features", {
  set.seed(99)
  G <- 12; per <- 60
  group <- rep(seq_len(G), each = per)
  z <- rnorm(G * per)
  z[group %in% 3:5] <- z[group %in% 3:5] + 3
  neighbors <- lapply(seq_len(G), function(g) setdiff(c(g - 1, g + 1), c(0, G + 1)))
  base <- spatial_fdr_sabha(z = z, group = group, neighbors = neighbors, lambda = 1)
  perm <- sample(G)
  nb2 <- vector("list", G)
  for (g in seq_len(G)) nb2[[perm[g]]] <- perm[neighbors[[g]]]
  rel <- spatial_fdr_sabha(z = z, group = perm[group], neighbors = nb2, lambda = 1)
  expect_identical(rel$reject, base$reject)
  expect_equal(rel$q, base$q, tolerance = 1e-12)
  expect_equal(rel$pi0_smooth[perm], base$pi0_smooth, tolerance = 1e-12)
  o <- sample(length(z))
  pf <- spatial_fdr_sabha(z = z[o], group = group[o], neighbors = neighbors, lambda = 1)
  expect_equal(pf$q[order(o)], base$q, tolerance = 1e-12)
})

test_that("neighbour smoothing uses the groups the neighbours list names", {
  p <- c(rep(0.9, 10), rep(0.001, 10), c(rep(0.9, 5), rep(0.001, 5)))
  group <- c(rep(3, 10), rep(1, 10), rep(2, 10))
  res <- spatial_fdr_sabha(p = p, group = group, neighbors = list(2L, c(1L, 3L), 2L),
                           lambda = 1, empirical_null = FALSE)
  raw <- pi0_reference(p, group, 0.5)
  expect_equal(raw, c(0, 1, 1))
  expected <- c((raw[1] + raw[2]) / 2, (raw[2] + mean(raw[c(1, 3)])) / 2, (raw[3] + raw[2]) / 2)
  expect_equal(res$pi0_smooth, pmax(expected, 0.05), tolerance = 1e-12)
})

test_that("weighted BH kernels match a direct R implementation (with NAs)", {
  set.seed(2026)
  m <- 2000
  p <- c(rbeta(300, 0.2, 6), runif(m - 300))
  w <- rexp(m) + 0.1
  p[sample(m, 40)] <- NA
  for (alpha in c(0.01, 0.05, 0.2)) {
    got <- sabha_weighted_bh_cpp(p, w, alpha = alpha)
    ref <- wbh_reference(p, w, alpha)
    expect_identical(as.logical(got$reject), ref$reject)
    expect_equal(got$k, ref$k)
    expect_equal(got$threshold, ref$threshold, tolerance = 1e-12)
    expect_equal(got$w_norm, ref$w_norm, tolerance = 1e-12)
    q <- sabha_bh_qvalues_scaled_cpp(p / got$w_norm)
    expect_equal(q, ref$q, tolerance = 1e-12)
    expect_identical(ref$reject, !is.na(q) & q <= alpha)
  }
})

test_that("end-to-end equals weighted BH with 1/pi0 weights; NA handling; validation", {
  set.seed(77)
  group <- sample(rep(1:8, each = 100))
  p <- runif(length(group))
  sig <- group %in% c(2, 7)
  p[sig] <- rbeta(sum(sig), 0.15, 5)
  p[c(5, 17, 230)] <- NA
  group[c(40, 41)] <- NA
  res <- spatial_fdr_sabha(p = p, group = group, lambda = 0, empirical_null = FALSE)
  pi0 <- pmin(pmax(pi0_reference(p[!is.na(group)], group[!is.na(group)], 0.5), 0.05), 1)
  expect_equal(res$pi0_smooth, pi0, tolerance = 1e-12)
  ref <- wbh_reference(p, ifelse(is.na(group), 1, 1 / pi0[group]), 0.05)
  expect_identical(res$reject, ref$reject)
  expect_equal(res$q, ref$q, tolerance = 1e-12)
  expect_true(all(is.na(res$q[is.na(p)])))
  expect_true(all(is.na(res$groups[c(40, 41)])))

  expect_error(spatial_fdr_sabha(p = p, group = group[-1]), "group length must match")
  expect_error(spatial_fdr_sabha(p = p, group = group, alpha = 1), "alpha must be in")
  expect_error(spatial_fdr_sabha(p = p, group = group, tau = 0), "tau must be in")
  expect_error(spatial_fdr_sabha(p = p, group = group, neighbors = list(1L)),
               "neighbors must be a list of length")
  expect_error(spatial_fdr_sabha(z = 1, p = 1, group = 1), "not both")
  # theoretical null with fewer than 10 finite z
  z <- c(-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 2.5)
  r <- spatial_fdr_sabha(z = z, group = rep(1:2, each = 4))
  expect_equal(c(r$mu0, r$sigma0), c(0, 1))
  expect_equal(r$p, 2 * pnorm(abs(z), lower.tail = FALSE), tolerance = 1e-12)
})

test_that("matrix input processes columns independently", {
  I <- sabha_case_inputs()
  Z <- cbind(a = I$zA, b = rev(I$zA))
  res <- spatial_fdr_sabha(z = Z, group = I$gA, neighbors = I$nbA)
  for (j in 1:2) {
    one <- spatial_fdr_sabha(z = Z[, j], group = I$gA, neighbors = I$nbA)
    expect_identical(res$q[, j], one$q)
    expect_identical(res$reject[, j], one$reject)
    expect_identical(res$pi0_smooth[, j], one$pi0_smooth)
    expect_identical(unname(res$sigma0[j]), one$sigma0)
  }
  expect_identical(colnames(res$q), c("a", "b"))
})

# ---------------------------------------------------------------------------
# Post-hoc integration
# ---------------------------------------------------------------------------

test_that("posthoc fdr:sabha uses z with the empirical null and equals the plain function", {
  I <- sabha_case_inputs()
  z <- array(c(I$zA, -I$zA * 1.1), c(2000, 1, 2),
             dimnames = list(NULL, "meta", c("c1", "c2")))
  g <- new_gds(list(z = z), space_sample_labels(sprintf("f%04d", 1:2000)),
               subjects = "meta", contrasts = c("c1", "c2"))
  out <- compute(posthoc(g, "fdr:sabha", options = list(group = I$gA, neighbors = I$nbA)))
  for (k in 1:2) {
    ref <- spatial_fdr_sabha(z = z[, 1, k], group = I$gA, neighbors = I$nbA)
    expect_identical(assay(out, "q")[, 1, k], ref$q)
  }
  # z input: no theoretical-null p is persisted next to the empirical-null q
  expect_null(out$assays[["p"]])
  att <- metadata(out)$attachments[["posthoc/fdr:sabha/q"]]
  expect_equal(att$statistic, "z")
  expect_equal(att$sigma0[1, 1], 1.0882101848346413, tolerance = 1e-10)
  expect_equal(dim(att$pi0_smooth), c(20L, 1L, 2L))

  # p mode adjusts the (derived) p-values without an empirical null
  outp <- compute(posthoc(g, "fdr:sabha", options = list(group = I$gA, statistic = "p")))
  refp <- spatial_fdr_sabha(p = 2 * pnorm(-abs(z[, 1, 1])), group = I$gA)
  expect_equal(assay(outp, "q")[, 1, 1], refp$q, tolerance = 1e-15)
})

test_that("posthoc fdr:sabha handles per-term families and row_data grouping", {
  set.seed(5)
  N <- 10; B <- 120
  subj <- paste0("s", 1:N)
  x <- rnorm(N)
  beta <- array(rnorm(B * N), c(B, N, 1), dimnames = list(NULL, subj, "c1"))
  beta[1:30, , 1] <- beta[1:30, , 1] + outer(rep(2, 30), x)
  g <- with_col_data(as_gds(list(beta = beta, var = array(1, dim(beta)))),
                     data.frame(x = x, row.names = subj))
  parcels <- rep(1:6, each = 20)
  fit <- compute(reduce(g, method = "ols:voxelwise", formula = ~ x))
  fit <- with_row_data(fit, data.frame(parcel = parcels))
  out <- compute(posthoc(fit, "fdr:sabha"))
  expect_true(all(c("q_coef:(Intercept)", "q_coef:x") %in% names(assays(out))))
  ref <- spatial_fdr_sabha(p = assay(fit, "p_coef:x")[, 1, 1], group = parcels)
  expect_equal(assay(out, "q_coef:x")[, 1, 1], ref$q, tolerance = 1e-15)
  # fdr:spatial (parcel-level Simes) is a different method and is unchanged
  sp <- compute(posthoc(fit, "fdr:spatial"))
  expect_false(isTRUE(all.equal(assay(sp, "q_coef:x"), assay(out, "q_coef:x"))))
  expect_error(compute(posthoc(compute(reduce(g, method = "ols:voxelwise", formula = ~ x)),
                               "fdr:sabha")), "needs a spatial grouping")
})

test_that("spatial_fdr_sabha validates inputs and falls back to the theoretical null", {
  expect_error(spatial_fdr_sabha(group = 1:3), "either z or p")
  p <- c(0.001, 0.2, 0.5, 0.04, 0.9, 0.03)
  g <- c(1, 1, 1, 2, 2, 2)
  expect_error(spatial_fdr_sabha(p = p, group = g, min_pi0 = 2), "min_pi0")
  expect_error(spatial_fdr_sabha(p = p, group = rep(NA, 6)), "No valid groups")
  expect_warning(neg <- spatial_fdr_sabha(p = p, group = g, lambda = -1), "lambda < 0")
  expect_identical(neg$q, spatial_fdr_sabha(p = p, group = g, lambda = 0)$q)
  # fewer than 10 finite z: the empirical null is N(0, 1), so p = 2 * pnorm(-|z|)
  z <- c(3.1, -0.4, 0.2, 2.5, -1.1, 0.7)
  tz <- spatial_fdr_sabha(z = z, group = g)
  expect_equal(c(tz$mu0, tz$sigma0), c(0, 1))
  expect_identical(tz$q, spatial_fdr_sabha(p = 2 * pnorm(-abs(z)), group = g)$q)
})

test_that("posthoc fdr:sabha resolves sources and per-term z families", {
  I <- sabha_case_inputs()
  z <- array(I$zA, c(2000, 1, 1), dimnames = list(NULL, "meta", "c1"))
  g <- new_gds(list(z = z, z_extra = z * 0.5), space_sample_labels(sprintf("f%04d", 1:2000)),
               subjects = "meta", contrasts = "c1")
  expect_error(compute(posthoc(g, "fdr:sabha", options = list(group = I$gA, statistic = "t"))),
               "statistic must be")
  expect_error(compute(posthoc(g, "fdr:sabha", options = list(group = I$gA, source = "nope"))),
               "source assay not found")
  out <- compute(posthoc(g, "fdr:sabha", options = list(group = I$gA, source = "z_extra")))
  expect_identical(assay(out, "q")[, 1, 1],
                   spatial_fdr_sabha(z = z[, 1, 1] * 0.5, group = I$gA)$q)
  # per-term z families map to q_coef:<term>
  gz <- new_gds(list(`z_coef:a` = z, `z_coef:b` = -z), space_sample_labels(sprintf("f%04d", 1:2000)),
                subjects = "meta", contrasts = "c1")
  oz <- compute(posthoc(gz, "fdr:sabha", options = list(group = I$gA)))
  expect_identical(assay(oz, "q_coef:b")[, 1, 1], spatial_fdr_sabha(z = -z[, 1, 1], group = I$gA)$q)
  gp <- new_gds(list(p = array(runif(2000), c(2000, 1, 1))),
                space_sample_labels(sprintf("f%04d", 1:2000)), subjects = "meta", contrasts = "c1")
  expect_error(compute(posthoc(gp, "fdr:sabha", options = list(group = I$gA, statistic = "z"))),
               "no z or z_coef")
})

test_that("Storey pi0 counts only tests with a finite p-value", {
  # group 1: 10 finite tests, 4 above tau, plus 10 missing; group 2: complete
  p <- c(c(0.9, 0.8, 0.7, 0.6, rep(0.01, 6)), rep(NA, 10),
         c(0.9, 0.95, 0.7, 0.6, 0.55, rep(0.02, 5)))
  g <- rep(c(1, 2), c(20, 10))
  res <- spatial_fdr_sabha(p = p, group = g, tau = 0.5, lambda = 0)
  expect_equal(res$pi0_raw, c(4 / (0.5 * 10), 5 / (0.5 * 10)))
  # the missing tests do not change group 1's estimate
  res_drop <- spatial_fdr_sabha(p = p[!is.na(p)], group = g[!is.na(p)], tau = 0.5, lambda = 0)
  expect_equal(res$pi0_raw, res_drop$pi0_raw)
  expect_equal(res$q[!is.na(p)], res_drop$q)
})
