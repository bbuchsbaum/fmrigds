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

