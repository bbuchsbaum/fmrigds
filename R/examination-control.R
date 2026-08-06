# Group-examination control -------------------------------------------------

#' Configure group examination
#'
#' `examination_control()` validates the bounded execution, residual-geometry,
#' review, staging, and numerical settings used by [examine_group()]. Keeping
#' these settings in a typed object makes defaults visible and reproducible.
#'
#' @param block_size Positive number of samples read per scan block.
#' @param geometry Named list controlling the deterministic low-rank residual
#'   representation. Supported entries are `rank`, `oversample`, `cap`,
#'   `balance_contrasts`, and `stability_replicates`.
#' @param review Named list controlling review status. Supported entries are
#'   `surprise`, `influence`, `quality`, and `min_stability`.
#' @param exact_refit_n Number of selected random-effects cases eligible for an
#'   exact heterogeneity refit.
#' @param retain_n Number of highest-priority subjects retained for localization
#'   when no explicit `retain` set is supplied.
#' @param staging Named list with optional `tempdir`. Staging stores are always
#'   removed on success or failure.
#' @param tolerance Named list of numerical tolerances: `rank`, `leverage`, and
#'   `degeneracy`.
#'
#' @return An object of class `gds_examination_control`.
#' @export
examination_control <- function(
    block_size = 1024L,
    geometry = list(
      rank = 32L,
      oversample = 8L,
      cap = 8,
      balance_contrasts = TRUE,
      stability_replicates = 2L
    ),
    review = list(
      surprise = list(
        energy_threshold = 2.5,
        tail_threshold = 0.01,
        residual_threshold = 3
      ),
      influence = list(
        energy_threshold = 1,
        max_abs_threshold = 2
      ),
      quality = list(
        coverage_fraction = list(direction = "low", threshold = 0.8)
      ),
      min_stability = 0.7
    ),
    exact_refit_n = 5L,
    retain_n = 5L,
    staging = list(tempdir = NULL),
    tolerance = list(
      rank = sqrt(.Machine$double.eps),
      leverage = 1e-8,
      degeneracy = 1e-12
    )) {
  geometry_defaults <- list(
    rank = 32L,
    oversample = 8L,
    cap = 8,
    balance_contrasts = TRUE,
    stability_replicates = 2L
  )
  review_defaults <- list(
    surprise = list(
      energy_threshold = 2.5,
      tail_threshold = 0.01,
      residual_threshold = 3
    ),
    influence = list(
      energy_threshold = 1,
      max_abs_threshold = 2
    ),
    quality = list(
      coverage_fraction = list(direction = "low", threshold = 0.8)
    ),
    min_stability = 0.7
  )
  staging_defaults <- list(tempdir = NULL)
  tolerance_defaults <- list(
    rank = sqrt(.Machine$double.eps),
    leverage = 1e-8,
    degeneracy = 1e-12
  )

  geometry <- .examination_merge_settings(geometry, geometry_defaults, "geometry")
  review <- .examination_merge_settings(review, review_defaults, "review")
  staging <- .examination_merge_settings(staging, staging_defaults, "staging")
  tolerance <- .examination_merge_settings(tolerance, tolerance_defaults, "tolerance")

  block_size <- .examination_positive_integer(block_size, "block_size")
  geometry$rank <- .examination_positive_integer(geometry$rank, "geometry rank")
  geometry$oversample <- .examination_nonnegative_integer(
    geometry$oversample,
    "geometry oversample"
  )
  geometry$stability_replicates <- .examination_nonnegative_integer(
    geometry$stability_replicates,
    "geometry stability_replicates"
  )
  if (!is.numeric(geometry$cap) || length(geometry$cap) != 1L ||
      !is.finite(geometry$cap) || geometry$cap <= 0) {
    stop("geometry cap must be a positive finite number.", call. = FALSE)
  }
  geometry$balance_contrasts <- .examination_scalar_logical(
    geometry$balance_contrasts,
    "geometry balance_contrasts"
  )

  if (!is.numeric(review$min_stability) || length(review$min_stability) != 1L ||
      is.na(review$min_stability) || review$min_stability < 0 ||
      review$min_stability > 1) {
    stop("review min_stability must be between 0 and 1.", call. = FALSE)
  }
  if (!is.list(review$quality)) {
    stop("review quality must be a list.", call. = FALSE)
  }
  review$surprise <- .validate_surprise_review(review$surprise)
  review$influence <- .validate_influence_review(review$influence)
  review$quality <- .validate_quality_review(review$quality)

  exact_refit_n <- .examination_nonnegative_integer(exact_refit_n, "exact_refit_n")
  retain_n <- .examination_nonnegative_integer(retain_n, "retain_n")
  if (!is.null(staging$tempdir) &&
      (!is.character(staging$tempdir) || length(staging$tempdir) != 1L ||
       is.na(staging$tempdir) || !nzchar(staging$tempdir))) {
    stop("staging tempdir must be NULL or one non-empty path.", call. = FALSE)
  }
  for (field in names(tolerance)) {
    value <- tolerance[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("tolerance ", field, " must be a positive finite number.", call. = FALSE)
    }
  }

  structure(
    list(
      block_size = block_size,
      geometry = geometry,
      review = review,
      exact_refit_n = exact_refit_n,
      retain_n = retain_n,
      staging = staging,
      tolerance = tolerance
    ),
    class = "gds_examination_control"
  )
}

.examination_merge_settings <- function(value, defaults, name) {
  if (!is.list(value)) {
    stop(name, " settings must be a list.", call. = FALSE)
  }
  unknown <- setdiff(names(value), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown ", name, " settings: ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  utils::modifyList(defaults, value, keep.null = TRUE)
}

.examination_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != as.integer(value) || value < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

.examination_nonnegative_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != as.integer(value) || value < 0L) {
    stop(name, " must be a nonnegative integer.", call. = FALSE)
  }
  as.integer(value)
}

.examination_scalar_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.validate_surprise_review <- function(spec) {
  defaults <- list(
    energy_threshold = 2.5,
    tail_threshold = 0.01,
    residual_threshold = 3
  )
  spec <- .examination_merge_settings(spec, defaults, "review surprise")
  for (field in c("energy_threshold", "residual_threshold")) {
    if (!is.numeric(spec[[field]]) || length(spec[[field]]) != 1L ||
        !is.finite(spec[[field]]) || spec[[field]] <= 0) {
      stop("review surprise ", field, " must be positive.", call. = FALSE)
    }
  }
  if (!is.numeric(spec$tail_threshold) || length(spec$tail_threshold) != 1L ||
      !is.finite(spec$tail_threshold) || spec$tail_threshold < 0 ||
      spec$tail_threshold > 1) {
    stop("review surprise tail_threshold must be between 0 and 1.", call. = FALSE)
  }
  spec
}

.validate_influence_review <- function(spec) {
  defaults <- list(energy_threshold = 1, max_abs_threshold = 2)
  spec <- .examination_merge_settings(spec, defaults, "review influence")
  for (field in names(defaults)) {
    if (!is.numeric(spec[[field]]) || length(spec[[field]]) != 1L ||
        !is.finite(spec[[field]]) || spec[[field]] <= 0) {
      stop("review influence ", field, " must be positive.", call. = FALSE)
    }
  }
  spec
}

.validate_quality_review <- function(spec) {
  if (!is.list(spec)) stop("review quality must be a list.", call. = FALSE)
  for (name in names(spec)) {
    item <- spec[[name]]
    if (!is.list(item) ||
        !identical(sort(names(item)), sort(c("direction", "threshold")))) {
      stop(
        "review quality '", name,
        "' must declare exactly direction and threshold.",
        call. = FALSE
      )
    }
    if (!is.character(item$direction) || length(item$direction) != 1L ||
        !item$direction %in% c("high", "low")) {
      stop("review quality direction must be high or low.", call. = FALSE)
    }
    if (!is.numeric(item$threshold) || length(item$threshold) != 1L ||
        !is.finite(item$threshold)) {
      stop("review quality threshold must be finite numeric.", call. = FALSE)
    }
  }
  spec
}
