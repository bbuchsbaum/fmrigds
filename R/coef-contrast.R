#' Linear contrasts of group-level coefficients
#'
#' Computes per-feature linear contrasts of the coefficients estimated by a
#' regression-type reducer (for example `"ols:voxelwise"`, `"meta:fe_reg"`,
#' `"meta:re_reg"`, or the `lmm:*` reducers) on a realised GDS. The estimates
#' are read from the `coef:<term>` assays and their sampling covariance from
#' the `cov:<a>:<b>` assays (or, failing that, from the packed attachment
#' returned by [coef_cov_tri()]); single-term contrasts need only
#' `se_coef:<term>`.
#'
#' For each contrast \eqn{c} and feature the estimate is \eqn{c^T \hat\beta}
#' with standard error \eqn{\sqrt{c^T V c}}. With `joint = TRUE` the rows of
#' `C` are tested together: \eqn{F = (C\hat\beta)^T (C V C^T)^{-1}
#' (C\hat\beta) / q} on \eqn{(q, \nu)} degrees of freedom for `stat = "t"`, or
#' the Wald statistic \eqn{q F} referred to \eqn{\chi^2_q} for `stat = "z"`.
#'
#' @param x A realised [`gds`] produced by a regression-type reducer.
#' @param C Contrast weights. A numeric vector named by term (one contrast), or
#'   a numeric matrix with one row per contrast and column names giving terms.
#'   Terms not named in `C` receive weight zero; names that are not
#'   coefficient terms of `x` are an error. Row names (or `name`) label the
#'   output assays.
#' @param stat Reference distribution: `"auto"` (default) uses `"t"` when the
#'   reducer reported `t_coef:<term>` assays and `"z"` otherwise (e.g. the
#'   normal-approximation meta-regression reducers).
#' @param df Optional degrees of freedom for `stat = "t"`: a scalar or an array
#'   conformable with the coefficient assays. By default the per-term
#'   `df_coef:<term>` assays are used when present (taking the smallest over
#'   the terms a contrast involves, which is conservative), else `df_res`,
#'   else `df`.
#' @param alternative Alternative hypothesis for single contrasts:
#'   `"two.sided"` (default), `"greater"`, or `"less"`. Joint tests are always
#'   upper-tailed.
#' @param joint Logical; if `TRUE`, test all rows of a matrix `C` jointly with
#'   one F (or chi-square) statistic instead of one t (or z) test per row.
#' @param name Optional label(s) for the output assays: one per row of `C`, or
#'   a single label when `joint = TRUE` (default `"joint"`).
#'
#' @return `x` with additional assays, one family per contrast `<name>`:
#'   - single contrasts: `contrast:<name>` (estimate), `se_contrast:<name>`,
#'     `t_contrast:<name>` and `df_contrast:<name>` (or `z_contrast:<name>`),
#'     and `p_contrast:<name>`;
#'   - joint tests: `F_contrast:<name>` with `df1_contrast:<name>` and
#'     `df2_contrast:<name>` (or `chi2_contrast:<name>` with
#'     `df1_contrast:<name>`), and `p_contrast:<name>`.
#'
#'   Features whose contrast variance is not positive and finite are `NA`.
#'   Adjust the p-values with, e.g.,
#'   `posthoc(y, "fdr:bh", options = list(source = "p_contrast:<name>"))`.
#' @seealso [coef_cov_tri()], [coef_array()], [reducer-ols-voxelwise]
#' @export
#' @examples
#' set.seed(1)
#' subj <- paste0("s", 1:12)
#' beta <- array(rnorm(4 * 12), c(4, 12, 1),
#'               dimnames = list(paste0("v", 1:4), subj, "c1"))
#' g <- as_gds(list(beta = beta, var = array(1, dim(beta), dimnames(beta))))
#' g <- with_col_data(g, data.frame(age = rnorm(12), dose = rnorm(12),
#'                                  row.names = subj))
#' fit <- compute(reduce(g, method = "ols:voxelwise", formula = ~ age + dose,
#'                       options = list(return_cov = "tri")))
#' res <- coef_contrast(fit, c(age = 1, dose = -1), name = "age_vs_dose")
#' assay(res, "t_contrast:age_vs_dose")[, 1, 1]
coef_contrast <- function(x,
                          C,
                          stat = c("auto", "t", "z"),
                          df = NULL,
                          alternative = c("two.sided", "greater", "less"),
                          joint = FALSE,
                          name = NULL) {
  if (!inherits(x, "gds")) {
    stop("`x` must be a realised gds (call compute() on a plan first).", call. = FALSE)
  }
  stat <- match.arg(stat)
  alternative <- match.arg(alternative)
  arrays <- x$assays
  terms <- sub("^coef:", "", grep("^coef:", names(arrays), value = TRUE))
  if (!length(terms)) {
    stop("`x` has no coef:<term> assays; run a regression-type reducer first.", call. = FALSE)
  }
  Cm <- .coef_contrast_matrix(C, terms)
  q <- nrow(Cm)
  if (isTRUE(joint) && q < 1L) stop("`C` must have at least one row.", call. = FALSE)
  labels <- .coef_contrast_labels(Cm, name, joint)

  if (identical(stat, "auto")) {
    stat <- if (any(grepl("^t_coef:", names(arrays))) || !is.null(df)) "t" else "z"
  }

  coef_of <- function(term) arrays[[paste0("coef:", term)]]
  base <- coef_of(terms[1L])
  dims <- dim(base)
  cov_lookup <- .coef_contrast_cov_lookup(x, dims)

  df_for <- function(used) {
    if (!identical(stat, "t")) return(NULL)
    if (!is.null(df)) {
      if (!is.numeric(df) || !(length(df) == 1L || length(df) == length(base))) {
        stop("`df` must be a scalar or conformable with the coefficient assays.", call. = FALSE)
      }
      return(array(as.numeric(df), dim = dims))
    }
    per_term <- lapply(used, function(t) arrays[[paste0("df_coef:", t)]])
    if (length(used) && all(!vapply(per_term, is.null, logical(1)))) {
      return(array(do.call(pmin, lapply(per_term, as.numeric)), dim = dims))
    }
    src <- arrays[["df_res"]] %||% arrays[["df"]]
    if (is.null(src)) {
      stop("stat = \"t\" needs degrees of freedom: supply `df`, or use a reducer ",
           "that reports df_coef:<term>, df_res or df.", call. = FALSE)
    }
    array(as.numeric(src), dim = dims)
  }

  out <- list()
  if (!isTRUE(joint)) {
    for (r in seq_len(q)) {
      w <- Cm[r, ]
      used <- terms[w != 0]
      est <- array(0, dim = dims)
      for (t in used) est <- est + w[[t]] * coef_of(t)
      v <- .coef_contrast_quadform(w, w, used, used, cov_lookup)
      se <- sqrt(v)
      se[!is.finite(v) | v <= 0] <- NA_real_
      tz <- est / se
      lab <- labels[r]
      out[[paste0("contrast:", lab)]] <- est
      out[[paste0("se_contrast:", lab)]] <- se
      if (identical(stat, "t")) {
        dfa <- df_for(used)
        p <- .p_from_t(tz, dfa, alternative)
        out[[paste0("t_contrast:", lab)]] <- tz
        out[[paste0("df_contrast:", lab)]] <- dfa
      } else {
        p <- .p_from_z(tz, alternative)
        out[[paste0("z_contrast:", lab)]] <- tz
      }
      dim(p) <- dims
      out[[paste0("p_contrast:", lab)]] <- p
    }
  } else {
    used <- terms[colSums(Cm != 0) > 0]
    est <- lapply(seq_len(q), function(r) {
      e <- array(0, dim = dims)
      for (t in used) if (Cm[r, t] != 0) e <- e + Cm[r, t] * coef_of(t)
      e
    })
    S <- matrix(list(), q, q)
    for (r in seq_len(q)) for (s in r:q) {
      S[[r, s]] <- .coef_contrast_quadform(Cm[r, ], Cm[s, ], used, used, cov_lookup)
      S[[s, r]] <- S[[r, s]]
    }
    W <- rep(NA_real_, length(base))
    for (e in seq_along(W)) {
      b <- vapply(est, `[`, numeric(1), e)
      M <- matrix(vapply(S, `[`, numeric(1), e), q, q)
      if (!all(is.finite(b)) || !all(is.finite(M))) next
      sol <- tryCatch(solve(M, b), error = function(err) NULL)
      if (is.null(sol)) next
      W[e] <- sum(b * sol)
    }
    W[!is.finite(W) | W < 0] <- NA_real_
    lab <- labels[1L]
    if (identical(stat, "t")) {
      dfa <- df_for(used)
      Fv <- W / q
      out[[paste0("F_contrast:", lab)]] <- array(Fv, dim = dims)
      out[[paste0("df1_contrast:", lab)]] <- array(q, dim = dims)
      out[[paste0("df2_contrast:", lab)]] <- dfa
      p <- stats::pf(Fv, q, as.numeric(dfa), lower.tail = FALSE)
    } else {
      out[[paste0("chi2_contrast:", lab)]] <- array(W, dim = dims)
      out[[paste0("df1_contrast:", lab)]] <- array(q, dim = dims)
      p <- stats::pchisq(W, q, lower.tail = FALSE)
    }
    out[[paste0("p_contrast:", lab)]] <- array(p, dim = dims)
  }
  out <- lapply(out, function(a) {
    a <- array(as.numeric(a), dim = dims)
    dimnames(a) <- dimnames(base)
    a
  })

  new_gds(
    assays = c(arrays[setdiff(names(arrays), names(out))], out),
    space = x$space,
    subjects = x$subjects,
    contrasts = x$contrasts,
    col_data = x$col_data,
    row_data = x$row_data,
    metadata = x$metadata
  )
}

.p_from_t <- function(t, df, alternative = "two.sided") {
  t <- as.numeric(t)
  df <- as.numeric(df)
  if (identical(alternative, "greater")) {
    stats::pt(t, df, lower.tail = FALSE)
  } else if (identical(alternative, "less")) {
    stats::pt(t, df, lower.tail = TRUE)
  } else {
    2 * stats::pt(-abs(t), df)
  }
}

.coef_contrast_matrix <- function(C, terms) {
  if (is.data.frame(C)) C <- as.matrix(C)
  if (!is.numeric(C)) stop("`C` must be numeric.", call. = FALSE)
  if (is.matrix(C)) {
    cn <- colnames(C)
    rn <- rownames(C)
  } else {
    cn <- names(C)
    rn <- NULL
    C <- matrix(C, nrow = 1L)
  }
  if (is.null(cn) || anyNA(cn) || any(!nzchar(cn))) {
    stop("`C` must be named by coefficient term (available terms: ",
         paste(terms, collapse = ", "), ").", call. = FALSE)
  }
  if (anyDuplicated(cn)) stop("`C` names terms more than once.", call. = FALSE)
  unknown <- setdiff(cn, terms)
  if (length(unknown)) {
    stop("Unknown term(s) in `C`: ", paste(unknown, collapse = ", "),
         ". Available terms: ", paste(terms, collapse = ", "), ".", call. = FALSE)
  }
  if (any(!is.finite(C))) stop("`C` must contain finite weights.", call. = FALSE)
  out <- matrix(0, nrow = nrow(C), ncol = length(terms), dimnames = list(rn, terms))
  out[, cn] <- C
  if (any(rowSums(out != 0) == 0)) stop("Every contrast in `C` needs a non-zero weight.", call. = FALSE)
  out
}

.coef_contrast_labels <- function(Cm, name, joint) {
  q <- nrow(Cm)
  if (isTRUE(joint)) {
    lab <- name %||% "joint"
    if (length(lab) != 1L) stop("`name` must be a single label when joint = TRUE.", call. = FALSE)
    return(as.character(lab))
  }
  lab <- name %||% rownames(Cm)
  if (is.null(lab)) {
    lab <- vapply(seq_len(q), function(r) {
      w <- Cm[r, Cm[r, ] != 0]
      paste(sprintf("%s%s", ifelse(w == 1, "", ifelse(w == -1, "-", paste0(format(w), "*"))),
                    names(w)), collapse = "+")
    }, character(1))
  }
  if (length(lab) != q || anyNA(lab) || any(!nzchar(lab)) || anyDuplicated(lab)) {
    stop("Contrast labels (`name` or row names of `C`) must be ", q,
         " unique non-empty strings.", call. = FALSE)
  }
  as.character(lab)
}

# Returns function(a, b) -> covariance array of coef:a and coef:b, read from
# `cov:a:b`/`cov:b:a` assays, `se_coef:a`^2 on the diagonal, or the packed
# reducer attachments.
.coef_contrast_cov_lookup <- function(x, dims) {
  arrays <- x$assays
  att <- .coef_contrast_attachment_cov(x, dims)
  function(a, b) {
    v <- arrays[[paste0("cov:", a, ":", b)]] %||% arrays[[paste0("cov:", b, ":", a)]]
    if (is.null(v)) v <- att(a, b)
    if (is.null(v) && identical(a, b)) {
      se <- arrays[[paste0("se_coef:", a)]]
      if (!is.null(se)) v <- se^2
    }
    if (is.null(v)) {
      stop(sprintf(
        paste0("No covariance between coefficients '%s' and '%s'. Contrasts that combine ",
               "terms need the coefficient covariance: rerun the reducer with ",
               "options = list(return_cov = \"tri\") where supported."),
        a, b), call. = FALSE)
    }
    v
  }
}

.coef_contrast_attachment_cov <- function(x, dims) {
  att <- x$metadata$attachments %||% list()
  keys <- names(att)
  keep <- vapply(att, function(a) is.list(a) && identical(a$type, "cov_tri"), logical(1))
  att <- att[keep]
  keys <- keys[keep]
  contrasts <- as.character(x$contrasts)
  function(a, b) {
    if (!length(att) || dims[2L] != 1L) return(NULL)
    out <- array(NA_real_, dim = dims)
    found <- FALSE
    for (k in seq_len(dims[3L])) {
      idx <- which(endsWith(keys, paste0("/contrast=", k)) |
                     endsWith(keys, paste0("/contrast=", contrasts[k])))
      if (!length(idx)) next
      e <- att[[idx[1L]]]
      i <- match(a, e$terms); j <- match(b, e$terms)
      if (is.na(i) || is.na(j)) next
      lo <- min(i, j); hi <- max(i, j); p <- length(e$terms)
      row <- (lo - 1L) * p - (lo - 1L) * (lo - 2L) / 2 + (hi - lo + 1L)
      out[, 1L, k] <- as.numeric(e$cov_tri[row, ])
      found <- TRUE
    }
    if (found) out else NULL
  }
}

# sum_i sum_j u_i v_j Cov(coef_i, coef_j) over the named terms.
.coef_contrast_quadform <- function(u, v, terms_u, terms_v, cov_lookup) {
  acc <- NULL
  for (a in terms_u) {
    if (u[[a]] == 0) next
    for (b in terms_v) {
      if (v[[b]] == 0) next
      term <- u[[a]] * v[[b]] * cov_lookup(a, b)
      acc <- if (is.null(acc)) term else acc + term
    }
  }
  acc
}
