# Post-hoc registry -----------------------------------------------------------

.gds_posthoc <- new.env(parent = emptyenv())

#' Register a post-hoc method
#'
#' Register or replace a post-hoc handler identified by `name`.
#'
#' The handler `fun` is called as `fun(arrays, opts)` where `arrays` is a
#' named list of assays shaped `sample x subject x contrast`. It should return
#' either (a) a named list of new/updated assays to merge (e.g., `list(q = ...)`),
#' or (b) the full `arrays` list with modifications. Implementations must
#' preserve the input array dimensions.
#'
#' @param name String identifier, e.g., "fdr:bh"
#' @param fun Function(arrays, opts) -> list of arrays (delta or full arrays)
#' @param requires Character vector of required input assays, e.g., c("p")
#' @param provides Character vector of output assay names, e.g., c("q")
#' @param overwrite Logical; if FALSE and a handler with `name` exists, error
#' @param case_deletion Optional declaration of whether and how the post-hoc
#'   conclusion can be recomputed for selected case deletions.
#' @return Invisibly, the registered post-hoc method `name`.
#' @seealso [posthoc-methods] for the built-in methods and their output assays.
#' @export
register_posthoc <- function(name,
                             fun,
                             requires,
                             provides,
                             overwrite = TRUE,
                             case_deletion = NULL) {
  stopifnot(is.character(name), length(name) == 1L, is.function(fun))
  if (!isTRUE(overwrite) && !is.null(.gds_posthoc[[name]])) {
    stop(sprintf("Post-hoc method '%s' already exists; set overwrite = TRUE to replace.", name), call. = FALSE)
  }
  .gds_posthoc[[name]] <- list(
    name = name,
    fun = fun,
    requires = as.character(requires),
    provides = as.character(provides),
    case_deletion = .normalize_posthoc_case_deletion(case_deletion)
  )
  invisible(name)
}

.normalize_posthoc_case_deletion <- function(contract = NULL) {
  defaults <- list(
    supported = FALSE,
    mode = "unavailable",
    requires_full_space = FALSE,
    deterministic = TRUE,
    seed_contract = NULL
  )
  if (is.null(contract)) return(defaults)
  if (!is.list(contract)) {
    stop("case_deletion must be NULL or a list.", call. = FALSE)
  }
  unknown <- setdiff(names(contract), names(defaults))
  if (length(unknown)) {
    stop("Unknown case_deletion fields: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- utils::modifyList(defaults, contract, keep.null = TRUE)
  for (field in c("supported", "requires_full_space", "deterministic")) {
    if (!is.logical(out[[field]]) || length(out[[field]]) != 1L || is.na(out[[field]])) {
      stop("case_deletion ", field, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.character(out$mode) || length(out$mode) != 1L || is.na(out$mode) ||
      !out$mode %in% c("unavailable", "recompute", "selected_recompute")) {
    stop(
      "case_deletion mode must be unavailable, recompute, or selected_recompute.",
      call. = FALSE
    )
  }
  if (isTRUE(out$supported) && identical(out$mode, "unavailable")) {
    stop("A supported case_deletion contract must declare a recomputation mode.", call. = FALSE)
  }
  if (!is.null(out$seed_contract) &&
      (!is.character(out$seed_contract) || length(out$seed_contract) != 1L ||
       is.na(out$seed_contract) || !nzchar(out$seed_contract))) {
    stop("case_deletion seed_contract must be NULL or one non-empty string.", call. = FALSE)
  }
  out
}

.posthoc_case_deletion_preflight <- function(nodes) {
  if (!length(nodes)) {
    return(data.frame(
      method = character(), supported = logical(), mode = character(),
      requires_full_space = logical(), deterministic = logical(),
      seed_contract = character(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(nodes, function(node) {
    method <- get_posthoc(node$method)
    contract <- method$case_deletion %||% .normalize_posthoc_case_deletion()
    data.frame(
      method = as.character(node$method),
      supported = isTRUE(contract$supported),
      mode = as.character(contract$mode),
      requires_full_space = isTRUE(contract$requires_full_space),
      deterministic = isTRUE(contract$deterministic),
      seed_contract = contract$seed_contract %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Get post-hoc method by name
#'
#' @param name Method identifier (e.g., "fdr:bh")
#' @return A list with fields `name`, `fun`, `requires`, `provides`, or `NULL` if not found
#' @export
get_posthoc <- function(name) {
  .gds_posthoc[[name]]
}

#' List registered post-hoc methods
#' @return A sorted character vector of registered post-hoc method names.
#' @seealso [posthoc-methods] for what each method outputs.
#' @export
list_posthoc <- function() {
  sort(ls(.gds_posthoc))
}

#' Built-in post-hoc correction methods and their output assays
#'
#' Methods are selected by name in [posthoc()] and add their adjusted-p assays
#' to the realised GDS. The FDR methods adjust the object's p-value assay(s):
#' a genuine bare `p` (from meta/combine reducers or a tabular `p` column) yields
#' a single `q`; a per-term `p_coef:<term>` family (from regression/LMM reducers
#' such as `ols:voxelwise`) yields one `q_coef:<term>` per term. Pass
#' `options = list(source = "p_coef:<term>")` to adjust one specific term.
#'
#' @section Methods and outputs:
#' - `"fdr:bh"`, `"fdr:by"`: add `q` (or per-term `q_coef:<term>`),
#'   Benjamini-Hochberg / Benjamini-Yekutieli adjusted p-values, computed
#'   column-wise per subject x contrast.
#' - `"fdr:spatial"`: adds `q` (or per-term `q_coef:<term>`), group-wise spatial
#'   FDR. Requires a discrete spatial grouping (one parcel/cluster id per
#'   sample) via `options$group` or a recognised `row_data` column
#'   (`feature_group`, `spatial_group`, `group`, `parcel`, `parcel_id`).
#' - `"nt:tfce_fwer"` (needs the optional `neurothresh` package): adds `q`
#'   (TFCE FWER-corrected p), `sig_mask` (1/0 at `alpha`), and `tfce` (the
#'   TFCE-enhanced statistic map).
#' - `"nt:cluster_fdr_perm"` (needs `neurothresh`): adds `q` (cluster-FDR q
#'   broadcast to voxels) and `sig_mask`.
#'
#' Required input assays are auto-derived (and persisted) when absent: e.g.
#' `fdr:bh` on a z-only GDS adds a derived `p` alongside `q`.
#'
#' @return This is a documentation page (no return value).
#' @seealso [posthoc()], [list_posthoc()], [list_assays()]
#' @name posthoc-methods
NULL

#' Unregister a post-hoc method
#'
#' @param name Method identifier to remove
#' @return Invisibly, TRUE if removed, FALSE if not found
#' @export
unregister_posthoc <- function(name) {
  if (is.null(.gds_posthoc[[name]])) return(invisible(FALSE))
  rm(list = name, envir = .gds_posthoc)
  invisible(TRUE)
}

# Built-ins: FDR BH/BY --------------------------------------------------------

# Resolve the p-value assay(s) an FDR post-hoc should adjust, returning a named
# list mapping each *output* assay name to its source p-array. This avoids the
# fragile `arrays$p` partial match (which silently resolved to a lone
# `p_coef:<term>` for single-term models and errored for multi-term models):
#
#   * `options$source = "<assay>"` -> adjust exactly that assay; output is
#     `q_coef:<term>` when the source is `p_coef:<term>`, else `q`.
#   * a genuine bare `p` assay (meta/combine reducers, tabular `p` columns,
#     or one derived from z/t+df) -> a single `q` (back-compatible).
#   * one or more `p_coef:<term>` assays (regression/LMM reducers) -> one
#     `q_coef:<term>` per term, mirroring the `coef:`/`se_coef:`/`t_coef:`/
#     `p_coef:` family.
.fdr_p_sources <- function(arrays, opts = list()) {
  nms <- names(arrays)
  if (!is.null(opts$source)) {
    src <- opts$source
    if (!src %in% nms) {
      stop("post-hoc FDR: options$source assay not found: ", src, call. = FALSE)
    }
    out <- if (grepl("^p_coef:", src)) sub("^p_coef:", "q_coef:", src) else "q"
    return(stats::setNames(list(arrays[[src]]), out))
  }
  if (!is.null(arrays[["p"]])) {
    return(list(q = arrays[["p"]]))
  }
  pcoef <- grep("^p_coef:", nms, value = TRUE)
  if (length(pcoef)) {
    return(stats::setNames(
      lapply(pcoef, function(n) arrays[[n]]),
      sub("^p_coef:", "q_coef:", pcoef)
    ))
  }
  # Last resort: derive a bare p from z or t+df (e.g. when the handler is called
  # directly without .ensure_required_arrays having added `p`).
  if (!is.null(arrays[["z"]])) {
    return(list(q = derive_p(list(z = arrays[["z"]]))))
  }
  if (!is.null(arrays[["t"]]) && !is.null(arrays[["df"]])) {
    return(list(q = derive_p(list(t = arrays[["t"]], df = arrays[["df"]]))))
  }
  stop("post-hoc FDR requires p values (or z/t with df)", call. = FALSE)
}

.posthoc_fdr_template <- function(method) {
  force(method)
  function(arrays, opts) {
    sources <- .fdr_p_sources(arrays, opts)
    out_list <- vector("list", length(sources))
    names(out_list) <- names(sources)
    for (out_name in names(sources)) {
      p <- sources[[out_name]]
      dims <- dim(p)
      out <- array(NA_real_, dim = dims)
      for (j in seq_len(dims[2])) {
        for (k in seq_len(dims[3])) {
          out[, j, k] <- stats::p.adjust(p[, j, k], method = method)
        }
      }
      out_list[[out_name]] <- out
    }
    out_list
  }
}

.register_builtin_posthoc <- function() {
  fdr_deletion <- list(
    supported = TRUE,
    mode = "recompute",
    requires_full_space = FALSE,
    deterministic = TRUE
  )
  register_posthoc(
    "fdr:bh", .posthoc_fdr_template("BH"),
    requires = c("p"), provides = c("q", "q_coef"),
    case_deletion = fdr_deletion
  )
  register_posthoc(
    "fdr:by", .posthoc_fdr_template("BY"),
    requires = c("p"), provides = c("q", "q_coef"),
    case_deletion = fdr_deletion
  )
  # Register spatial FDR if available
  if (exists(".posthoc_spatial_fdr", mode = "function")) {
    register_posthoc(
      "fdr:spatial", .posthoc_spatial_fdr(),
      requires = c("p"), provides = c("q", "q_coef"),
      case_deletion = list(
        supported = TRUE,
        mode = "recompute",
        requires_full_space = TRUE,
        deterministic = TRUE
      )
    )
  }
  # Register neurothresh-backed TFCE if optional dependency is available
  if (exists(".posthoc_neurothresh_tfce_fwer", mode = "function") &&
      requireNamespace("neurothresh", quietly = TRUE)) {
    register_posthoc(
      "nt:tfce_fwer",
      .posthoc_neurothresh_tfce_fwer(),
      requires = c("z"),
      provides = c("q", "sig_mask", "tfce"),
      case_deletion = list(
        supported = TRUE,
        mode = "selected_recompute",
        requires_full_space = TRUE,
        deterministic = FALSE,
        seed_contract = "options$seed"
      )
    )
    register_posthoc(
      "nt:cluster_fdr_perm",
      .posthoc_neurothresh_cluster_fdr_perm(),
      requires = c("z"),
      provides = c("q", "sig_mask"),
      case_deletion = list(
        supported = TRUE,
        mode = "selected_recompute",
        requires_full_space = TRUE,
        deterministic = FALSE,
        seed_contract = "options$seed"
      )
    )
  }
}

# initialize defaults when package loads
NULL
