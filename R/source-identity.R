# Source identity receipts -------------------------------------------------

.source_identity_schema_version <- "1.0.0"

.normalize_source_identity_policy <- function(policy = "sha256") {
  if (is.null(policy)) policy <- "sha256"
  if (!is.character(policy) || length(policy) != 1L || is.na(policy)) {
    stop("`source_identity` must be \"sha256\" or \"none\".", call. = FALSE)
  }
  match.arg(policy, c("sha256", "none"))
}

.source_stat <- function(path) {
  info <- file.info(path)
  if (!nrow(info) || is.na(info$size[[1L]])) return(NULL)
  list(
    byteSize = unname(as.numeric(info$size[[1L]])),
    modifiedTime = format(
      as.POSIXct(info$mtime[[1L]], tz = "UTC"),
      "%Y-%m-%dT%H:%M:%OS6Z",
      tz = "UTC"
    )
  )
}

.source_stat_equal <- function(x, y) {
  if (is.null(x) || is.null(y) ||
      !identical(as.numeric(x$byteSize), as.numeric(y$byteSize))) return(FALSE)
  if (!is.null(x$modifiedTime) && !is.null(y$modifiedTime)) {
    return(identical(as.character(x$modifiedTime), as.character(y$modifiedTime)))
  }
  # Compatibility for an in-memory receipt made by an earlier development
  # build. JSON round-tripping can lose sub-microsecond precision.
  if (!is.null(x$modifiedEpoch) && !is.null(y$modifiedEpoch)) {
    return(abs(as.numeric(x$modifiedEpoch) - as.numeric(y$modifiedEpoch)) <= 1e-6)
  }
  FALSE
}

.source_identity_failure <- function(status, path, reason, stat = NULL) {
  list(
    schemaVersion = .source_identity_schema_version,
    byteSize = stat$byteSize %||% NA_real_,
    digest = NULL,
    identityStatus = status,
    verifiedAt = .iso_utc(Sys.time()),
    reason = reason,
    private = list(
      locator = normalizePath(path, winslash = "/", mustWork = FALSE),
      stat = stat,
      verifyOnRead = FALSE
    )
  )
}

.hash_source_file <- function(locator) {
  digest::digest(locator, algo = "sha256", serialize = FALSE, file = TRUE)
}

# Hash one file exactly once while recording the stat values immediately before
# and after hashing. `on_error = "record"` is used by diagnostics/tests that
# need to distinguish unavailable inputs; adapters use the fail-closed default.
.capture_file_identity <- function(path,
                                   policy = "sha256",
                                   on_error = c("error", "record")) {
  policy <- .normalize_source_identity_policy(policy)
  on_error <- match.arg(on_error)
  path <- as.character(path)
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Source path must be one non-empty string.", call. = FALSE)
  }

  fail <- function(status, reason, stat = NULL) {
    receipt <- .source_identity_failure(status, path, reason, stat)
    if (identical(on_error, "error")) {
      stop("Cannot establish source identity for '", path, "': ", reason, call. = FALSE)
    }
    receipt
  }

  if (!file.exists(path) || dir.exists(path)) {
    return(fail("unavailable", "file is missing or is not a regular file"))
  }
  before <- .source_stat(path)
  if (is.null(before)) {
    return(fail("unavailable", "file metadata could not be read"))
  }

  locator <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (identical(policy, "none")) {
    return(list(
      schemaVersion = .source_identity_schema_version,
      byteSize = before$byteSize,
      digest = NULL,
      identityStatus = "not-requested",
      verifiedAt = .iso_utc(Sys.time()),
      reason = "source_identity policy is none",
      private = list(
        locator = locator,
        stat = before,
        verifyOnRead = TRUE,
        policy = policy
      )
    ))
  }

  value <- tryCatch(
    .hash_source_file(locator),
    error = function(e) e
  )
  if (inherits(value, "error")) {
    return(fail("unavailable", conditionMessage(value), before))
  }
  after <- .source_stat(locator)
  if (!.source_stat_equal(before, after)) {
    return(fail(
      "changed-during-read",
      "file size or modification time changed while its digest was being computed",
      after %||% before
    ))
  }

  list(
    schemaVersion = .source_identity_schema_version,
    byteSize = after$byteSize,
    digest = list(algorithm = "sha256", value = tolower(value)),
    identityStatus = "verified",
    verifiedAt = .iso_utc(Sys.time()),
    reason = NULL,
    private = list(
      locator = locator,
      stat = after,
      verifyOnRead = TRUE,
      policy = policy
    )
  )
}

.source_role_record <- function(role, ordinal, pair = NULL) {
  out <- list(role = as.character(role), ordinal = as.integer(ordinal))
  if (!is.null(pair)) out$pair <- as.integer(pair)
  out
}

.source_entity_id <- function(receipt, roles, ordinal, occurrence) {
  identity <- if (!is.null(receipt$digest$value)) {
    receipt$digest$value
  } else {
    .portable_sha256(list(
      status = receipt$identityStatus,
      roles = roles,
      ordinal = as.integer(ordinal)
    ))$value
  }
  suffix <- .portable_sha256(list(
    identity = identity,
    occurrence = occurrence
  ))$value
  sprintf("source-%04d-%s", as.integer(ordinal), substr(suffix, 1L, 16L))
}

# `records` is an ordered list of list(path, role, ordinal, pair, kind). The
# same physical path is hashed once and retains every logical role.
.source_entities_from_files <- function(records,
                                        policy = "sha256",
                                        on_error = c("error", "record")) {
  policy <- .normalize_source_identity_policy(policy)
  on_error <- match.arg(on_error)
  if (!length(records)) return(list())
  if (!is.list(records)) stop("Source records must be a list.", call. = FALSE)

  normalized <- lapply(seq_along(records), function(i) {
    record <- records[[i]]
    if (!is.list(record) || is.null(record$path) || is.null(record$role)) {
      stop("Each source record must provide `path` and `role`.", call. = FALSE)
    }
    path <- as.character(record$path)
    list(
      path = path,
      key = normalizePath(path, winslash = "/", mustWork = FALSE),
      kind = as.character(record$kind %||% "file"),
      role = as.character(record$role),
      ordinal = as.integer(record$ordinal %||% i),
      pair = if (is.null(record$pair)) NULL else as.integer(record$pair),
      expectedStat = record$expected_stat %||% NULL
    )
  })

  keys <- unique(vapply(normalized, `[[`, character(1L), "key"))
  entities <- lapply(keys, function(key) {
    idx <- which(vapply(normalized, function(x) identical(x$key, key), logical(1L)))
    group <- normalized[idx]
    receipt <- .capture_file_identity(group[[1L]]$path, policy, on_error)
    expected_stats <- Filter(
      Negate(is.null),
      lapply(group, `[[`, "expectedStat")
    )
    changed_before_identity <- length(expected_stats) && any(!vapply(
      expected_stats,
      .source_stat_equal,
      logical(1L),
      y = receipt$private$stat %||% NULL
    ))
    if (changed_before_identity) {
      receipt$identityStatus <- "changed-during-read"
      receipt$reason <- "file size or modification time changed between adapter read and identity capture"
      receipt$digest <- NULL
      receipt$private$preReadStat <- expected_stats[[1L]]
      if (identical(on_error, "error")) {
        stop(
          "Cannot establish source identity for '", group[[1L]]$path,
          "': ", receipt$reason,
          call. = FALSE
        )
      }
    }
    roles <- lapply(group, function(x) .source_role_record(x$role, x$ordinal, x$pair))
    ordinal <- min(vapply(group, `[[`, integer(1L), "ordinal"))
    list(
      id = .source_entity_id(
        receipt,
        roles,
        ordinal,
        .next_portable_occurrence("source-entity")
      ),
      schemaVersion = .source_identity_schema_version,
      kind = group[[1L]]$kind,
      roles = roles,
      byteSize = receipt$byteSize,
      digest = receipt$digest,
      identityStatus = receipt$identityStatus,
      verifiedAt = receipt$verifiedAt,
      reason = receipt$reason,
      private = receipt$private
    )
  })
  .assert_unique_source_entity_ids(entities)
  entities
}

.source_nonfile_entity <- function(kind,
                                   role = "object",
                                   ordinal = 1L,
                                   status = "unavailable",
                                   reason = "source has no stable byte representation") {
  roles <- list(.source_role_record(role, ordinal))
  seed <- .portable_sha256(list(
    kind = kind,
    roles = roles,
    status = status,
    occurrence = .next_portable_occurrence("source-entity")
  ))$value
  list(
    id = sprintf("source-%04d-%s", as.integer(ordinal), substr(seed, 1L, 16L)),
    schemaVersion = .source_identity_schema_version,
    kind = as.character(kind),
    roles = roles,
    byteSize = NA_real_,
    digest = NULL,
    identityStatus = status,
    verifiedAt = .iso_utc(Sys.time()),
    reason = reason,
    private = list(verifyOnRead = FALSE)
  )
}

.assert_unique_source_entity_ids <- function(entities) {
  if (!length(entities)) return(invisible(TRUE))
  ids <- vapply(entities, `[[`, character(1L), "id")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Source entity IDs must be present and unique.", call. = FALSE)
  }
  invisible(TRUE)
}

.portable_source_entity <- function(entity) {
  keep <- c(
    "id", "schemaVersion", "kind", "roles", "byteSize", "digest",
    "identityStatus", "reason"
  )
  out <- entity[intersect(keep, names(entity))]
  # JSON has no portable NA numeric. Unknown sizes are represented by absence.
  if (!is.null(out$byteSize) && (length(out$byteSize) != 1L || !is.finite(out$byteSize))) {
    out$byteSize <- NULL
  }
  out
}

.portable_source_entities <- function(entities) {
  lapply(entities %||% list(), .portable_source_entity)
}

.metadata_with_source_entities <- function(metadata = list(), entities = list()) {
  metadata$provenance <- .normalize_provenance(metadata$provenance %||% NULL)
  .append_source_entities(metadata, entities)
}

.deactivate_source_verification <- function(metadata) {
  entities <- metadata$provenance$entities %||% list()
  if (!length(entities)) return(metadata)
  for (i in seq_along(entities)) {
    if (is.null(entities[[i]]$private)) entities[[i]]$private <- list()
    entities[[i]]$private$verifyOnRead <- FALSE
  }
  metadata$provenance$entities <- entities
  metadata
}

.verify_source_entities <- function(entities,
                                    phase = c("before-read", "after-read"),
                                    fail = TRUE) {
  phase <- match.arg(phase)
  changed <- list()
  for (i in seq_along(entities %||% list())) {
    entity <- entities[[i]]
    private <- entity$private %||% list()
    if (!isTRUE(private$verifyOnRead) || is.null(private$locator)) next
    current <- .source_stat(private$locator)
    if (!.source_stat_equal(private$stat, current)) {
      entity$identityStatus <- "changed-during-read"
      entity$reason <- paste("file size or modification time changed", phase)
      entity$private$currentStat <- current
      changed[[length(changed) + 1L]] <- entity
    }
  }
  if (length(changed) && isTRUE(fail)) {
    locators <- vapply(changed, function(x) x$private$locator, character(1L))
    stop(
      "Source identity changed ", phase, ": ", paste(locators, collapse = ", "),
      ". Re-open the source and rebuild the plan.",
      call. = FALSE
    )
  }
  changed
}
