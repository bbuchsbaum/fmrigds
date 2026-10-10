# Internal: keep derived assays coherent without overwriting existing fields

.sync_derived <- function(arrays) {
  # se <-> var
  if (is.null(arrays$se) && !is.null(arrays$var)) {
    arrays$se <- tryCatch(derive_se(arrays), error = function(e) arrays$se)
  }
  if (is.null(arrays$var) && !is.null(arrays$se)) {
    arrays$var <- tryCatch(derive_var(arrays), error = function(e) arrays$var)
  }
  # t from beta/var
  if (is.null(arrays$t) && all(c("beta", "var") %in% names(arrays))) {
    arrays$t <- tryCatch(derive_t(arrays), error = function(e) arrays$t)
  }
  # z/p from available stats
  if (is.null(arrays$z)) {
    arrays$z <- tryCatch(derive_z(arrays), error = function(e) arrays$z)
  }
  if (is.null(arrays$p)) {
    arrays$p <- tryCatch(derive_p(arrays), error = function(e) arrays$p)
  }
  arrays
}

# Internal: rule for assays after a space transformation (map_to()/align()).
#
# A linear map or alignment changes the sample axis, so every assay must be
# either explicitly propagated to the target space or dropped. Only the assays
# named in `propagated` (computed for the target space by the caller, e.g.
# beta/var and, when correctly aggregated, df) survive. Everything else (se, t,
# z, p, q, n_eff, custom assays, ...) is discarded and the standard derived
# statistics are re-derived from the propagated ones. When no df survives but
# the source carried test statistics (z/p), a Wald z = beta / se is used so
# that p stays consistent with the new effects.
.reset_after_transform <- function(arrays, propagated, had_stats = FALSE) {
  arrays <- arrays[intersect(propagated, names(arrays))]
  arrays <- Filter(Negate(is.null), arrays)
  if (isTRUE(had_stats) && is.null(arrays$df) && all(c("beta", "var") %in% names(arrays))) {
    arrays$z <- arrays$beta / sqrt(arrays$var)
  }
  .sync_derived(arrays)
}

# Internal: carry df across a transformation without Satterthwaite
# aggregation. This is only valid when df is constant across the source
# samples for each subject x contrast (e.g. one first-level model df per
# subject); then it is broadcast to the `n_target` target samples. Otherwise
# NULL is returned and df is dropped.
.broadcast_constant_df <- function(df, n_target) {
  if (is.null(df)) return(NULL)
  d <- dim(df)
  out <- array(NA_real_, dim = c(n_target, d[2], d[3]))
  for (j in seq_len(d[2])) {
    for (k in seq_len(d[3])) {
      vals <- unique(df[, j, k])
      if (length(vals) != 1L) return(NULL)
      out[, j, k] <- vals
    }
  }
  out
}

