# =============================================================================
# Image Catalog - File discovery and organization for neuroimaging stat maps
# =============================================================================

# -----------------------------------------------------------------------------
# Constructor and Class Definition
# -----------------------------------------------------------------------------

#' Create a new image catalog object
#'
#' Low-level constructor for image_catalog objects. Users should typically
#' use [image_catalog()] for file discovery instead.
#'
#' @param files Character vector of absolute file paths.
#' @param metadata Data frame with one row per file containing extracted
#'   path components and user-assigned metadata.
#' @param root_dir Root directory used for discovery (for relative path display).
#' @param pattern Original glob pattern used for discovery.
#' @param assay_map Named list mapping assay names to file selection criteria.
#'
#' @return An object of class `c("image_catalog", "list")`.
#' @export
new_image_catalog <- function(files,
                               metadata = NULL,
                               root_dir = NULL,
                               pattern = NULL,
                               assay_map = list()) {

  stopifnot(is.character(files))

  if (is.null(metadata)) {
    metadata <- data.frame(
      file = files,
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      files = files,
      metadata = metadata,
      root_dir = root_dir,
      pattern = pattern,
      assay_map = assay_map,
      created = Sys.time()
    ),
    class = c("image_catalog", "list")
  )
}

#' Print method for image_catalog
#'
#' @param x An image_catalog object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the catalog.
#' @method print image_catalog
#' @rdname print.image_catalog
#' @export
print.image_catalog <- function(x, ...) {
  cat("Image Catalog\n")
  cat("  Files: ", length(x$files), "\n", sep = "")
  if (!is.null(x$root_dir)) cat("  Root: ", x$root_dir, "\n", sep = "")
  if (!is.null(x$pattern)) cat("  Pattern: ", x$pattern, "\n", sep = "")

  # Show metadata columns (excluding internal ones)
  meta_cols <- setdiff(names(x$metadata), c("file", "basename", "relpath"))
  if (length(meta_cols)) {
    cat("  Metadata columns: ", paste(meta_cols, collapse = ", "), "\n", sep = "")
  }

  # Show unique subjects if available

if ("subject" %in% names(x$metadata)) {
    subjs <- unique(x$metadata$subject)
    subjs <- subjs[!is.na(subjs)]
    n_subjs <- length(subjs)
    if (n_subjs > 0) {
      subj_display <- if (n_subjs <= 5) {
        paste0(" (", paste(subjs, collapse = ", "), ")")
      } else {
        paste0(" (", paste(utils::head(subjs, 3), collapse = ", "), ", ...)")
      }
      cat("  Subjects: ", n_subjs, subj_display, "\n", sep = "")
    }
  }

  # Show assay mappings
  if (length(x$assay_map)) {
    cat("  Assay mappings: ", paste(names(x$assay_map), collapse = ", "), "\n", sep = "")
  }

  invisible(x)
}

#' Summary method for image_catalog
#'
#' @param object An image_catalog object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the catalog.
#' @method summary image_catalog
#' @rdname summary.image_catalog
#' @export
summary.image_catalog <- function(object, ...) {
  meta <- object$metadata

  cat("Image Catalog Summary\n")
  cat("=====================\n")
  cat("Total files: ", length(object$files), "\n", sep = "")
  if (!is.null(object$root_dir)) cat("Root directory: ", object$root_dir, "\n", sep = "")
  if (!is.null(object$pattern)) cat("Discovery pattern: ", object$pattern, "\n", sep = "")
  cat("\nMetadata columns:\n")

  # Summarize each metadata column
  for (col in names(meta)) {
    if (col %in% c("file", "relpath")) next
    vals <- meta[[col]]
    n_unique <- length(unique(vals[!is.na(vals)]))
    n_na <- sum(is.na(vals))
    cat("  ", col, ": ", n_unique, " unique",
        if (n_na > 0) paste0(" (", n_na, " NA)") else "", "\n", sep = "")
  }

  # Assay map summary
  if (length(object$assay_map)) {
    cat("\nAssay mappings:\n")
    for (nm in names(object$assay_map)) {
      spec <- object$assay_map[[nm]]
      if (is.character(spec)) {
        cat("  ", nm, ": pattern '", spec, "'\n", sep = "")
      } else if (is.list(spec)) {
        cat("  ", nm, ": ", names(spec)[1], " == '", spec[[1]], "'\n", sep = "")
      }
    }
  }

  invisible(object)
}

# -----------------------------------------------------------------------------
# Discovery Functions
# -----------------------------------------------------------------------------

#' Discover image files and create a catalog
#'
#' Sweeps a directory tree for files matching a glob pattern and extracts
#' metadata from file paths using regex capture groups.
#'
#' @param root Root directory to search.
#' @param pattern Glob pattern for file discovery (e.g., `"sub-*/analysis/*.nii.gz"`).
#'   Supports `*` for single directory level and `**` for recursive matching.
#' @param path_regex Optional regex with named capture groups for metadata extraction.
#'   Named groups like `(?<subject>[^/]+)` become metadata columns.
#' @param recursive Whether to search subdirectories. Default `TRUE`. Only used
#'   if pattern does not contain `**`.
#'
#' @return An image_catalog object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Discover files with subject extraction
#' catalog <- image_catalog(
#'   root = "study",
#'   pattern = "sub-*/analysis/*.nii.gz",
#'   path_regex = "sub-(?<subject>[^/]+)"
#' )
#' }
image_catalog <- function(root,
                          pattern = "**/*.nii*",
                          path_regex = NULL,
                          recursive = TRUE) {
  if (!dir.exists(root)) {
    stop("Root directory does not exist: ", root, call. = FALSE)
  }

  root <- normalizePath(root, mustWork = TRUE)

  # Discover files using glob pattern
  files <- .catalog_discover_files(root, pattern, recursive)

  if (!length(files)) {
    warning("No files found matching pattern '", pattern, "'", call. = FALSE)
    return(new_image_catalog(character(), root_dir = root, pattern = pattern))
  }

  # Extract metadata from paths
  metadata <- .catalog_extract_metadata(files, path_regex, root)

  new_image_catalog(
    files = files,
    metadata = metadata,
    root_dir = root,
    pattern = pattern
  )
}

#' Discover files using glob pattern
#'
#' @param root Root directory
#' @param pattern Glob pattern
#' @param recursive Whether to search recursively
#' @return Character vector of file paths
#' @keywords internal
.catalog_discover_files <- function(root, pattern, recursive) {
  # Build full glob pattern
  full_pattern <- file.path(root, pattern)

  # Use Sys.glob for pattern expansion
  files <- Sys.glob(full_pattern)

  # If no results and pattern doesn't have **, try recursive list.files
  if (!length(files) && !grepl("\\*\\*", pattern) && recursive) {
    # Convert glob to regex for list.files
    regex_pattern <- .glob_to_regex(basename(pattern))
    search_dir <- file.path(root, dirname(pattern))
    if (dir.exists(search_dir)) {
      files <- list.files(
        search_dir,
        pattern = regex_pattern,
        recursive = recursive,
        full.names = TRUE
      )
    }
  }

  # Ensure absolute paths and sort
  files <- normalizePath(files, mustWork = FALSE)
  files <- files[file.exists(files)]
  sort(files)
}

#' Convert simple glob pattern to regex
#' @keywords internal
.glob_to_regex <- function(pattern) {
  # First, mark ** with a placeholder to avoid double-processing
  pattern <- gsub("\\*\\*", "\001DOUBLESTAR\001", pattern)
  # Also mark single * with placeholder
  pattern <- gsub("\\*", "\002SINGLESTAR\002", pattern)
  # Mark ? with placeholder
  pattern <- gsub("\\?", "\003QUESTION\003", pattern)

  # Escape regex special chars
  pattern <- gsub("\\.", "\\\\.", pattern)
  pattern <- gsub("\\+", "\\\\+", pattern)
  pattern <- gsub("\\^", "\\\\^", pattern)
  pattern <- gsub("\\$", "\\\\$", pattern)
  pattern <- gsub("\\{", "\\\\{", pattern)
  pattern <- gsub("\\}", "\\\\}", pattern)
  pattern <- gsub("\\|", "\\\\|", pattern)
  pattern <- gsub("\\(", "\\\\(", pattern)
  pattern <- gsub("\\)", "\\\\)", pattern)
  pattern <- gsub("\\[", "\\\\[", pattern)
  pattern <- gsub("\\]", "\\\\]", pattern)

  # Convert glob wildcards from placeholders
  pattern <- gsub("\001DOUBLESTAR\001", ".*", pattern, fixed = TRUE)       # ** -> .*
  pattern <- gsub("\002SINGLESTAR\002", "[^/]*", pattern, fixed = TRUE)    # * -> [^/]*
  pattern <- gsub("\003QUESTION\003", ".", pattern, fixed = TRUE)          # ? -> .

  paste0("^", pattern, "$")
}

#' Extract metadata from file paths using regex
#'
#' @param files Character vector of file paths
#' @param path_regex Regex pattern with named capture groups
#' @param root Root directory for relative path computation
#' @return Data frame with extracted components
#' @keywords internal
.catalog_extract_metadata <- function(files, path_regex, root = NULL) {
  # Base metadata always includes file path and basename
  meta <- data.frame(
    file = files,
    basename = basename(files),
    stringsAsFactors = FALSE
  )

  if (!is.null(root)) {
    # Compute relative paths for cleaner display
    root_escaped <- gsub("([.+^${}|()\\[\\]\\\\])", "\\\\\\1", root)
    meta$relpath <- sub(paste0("^", root_escaped, "/?"), "", files)
  }

  # Extract using named capture groups if regex provided
  if (!is.null(path_regex)) {
    group_names <- .extract_capture_group_names(path_regex)

    if (length(group_names)) {
      # Use regmatches with perl-style named groups
      matches <- regmatches(files, regexec(path_regex, files, perl = TRUE))

      for (nm in group_names) {
        meta[[nm]] <- vapply(matches, function(m) {
          if (length(m) == 0 || is.na(m[1])) return(NA_character_)
          # Find the named group in match names
          match_names <- names(m)
          if (!is.null(match_names) && nm %in% match_names) {
            val <- m[nm]
            if (is.na(val) || val == "") NA_character_ else val
          } else {
            NA_character_
          }
        }, character(1), USE.NAMES = FALSE)
      }
    }
  }

  meta
}

#' Extract named capture group names from a regex pattern
#'
#' @param pattern Regex pattern with named capture groups
#' @return Character vector of group names
#' @keywords internal
.extract_capture_group_names <- function(pattern) {
  # Match (?<name>...) patterns
  m <- gregexpr("\\(\\?<([^>]+)>", pattern, perl = TRUE)
  if (m[[1]][1] == -1L) return(character())

  # Extract the matched strings
  matched <- regmatches(pattern, m)[[1]]

  # Extract just the names from (?<name>
  gsub("\\(\\?<([^>]+)>", "\\1", matched)
}

# -----------------------------------------------------------------------------
# Metadata Assignment Functions
# -----------------------------------------------------------------------------

#' Assign metadata via regex pattern matching
#'
#' Apply regex patterns to filenames or paths to extract values for new
#' metadata columns.
#'
#' @param catalog An image_catalog object.
#' @param column Name of new column to create.
#' @param pattern Regex pattern to match (should contain a capture group).
#' @param source Which existing column to match against. Default `"basename"`.
#' @param replacement Replacement string using backreferences. Default `"\\1"`.
#'
#' @return Updated image_catalog.
#' @export
#'
#' @examples
#' \dontrun{
#' catalog <- assign_meta(catalog, "run", "run-([0-9]+)", source = "basename")
#' catalog <- assign_meta(catalog, "contrast", "([^_]+)\\.nii", replacement = "\\1")
#' }
assign_meta <- function(catalog,
                        column,
                        pattern,
                        source = "basename",
                        replacement = "\\1") {
  stopifnot(inherits(catalog, "image_catalog"))

  src_col <- catalog$metadata[[source]]
  if (is.null(src_col)) {
    stop("Source column '", source, "' not found in catalog metadata", call. = FALSE)
  }

  # Check which rows match the pattern
  matches <- grepl(pattern, src_col, perl = TRUE)

  # Extract using sub with backreference
  # Pattern should contain capture group(s), replacement uses backrefs like \1
  values <- rep(NA_character_, length(src_col))
  values[matches] <- sub(
    paste0(".*?", pattern, ".*"),
    replacement,
    src_col[matches],
    perl = TRUE
  )

  catalog$metadata[[column]] <- values
  catalog
}

#' Join external metadata to catalog
#'
#' Merge a data frame with the catalog metadata by a key column.
#'
#' @param catalog An image_catalog object.
#' @param data External data.frame to join.
#' @param by Column name(s) for joining (must exist in both catalog$metadata and data).
#' @param suffix Suffix for duplicate column names from data. Default `".external"`.
#'
#' @return Updated image_catalog.
#' @export
#'
#' @examples
#' \dontrun{
#' demographics <- data.frame(subject = c("01", "02"), age = c(25, 30))
#' catalog <- join_meta(catalog, demographics, by = "subject")
#' }
join_meta <- function(catalog,
                      data,
                      by,
                      suffix = ".external") {
  stopifnot(inherits(catalog, "image_catalog"))
  stopifnot(is.data.frame(data))

  # Check join keys exist
  missing_cat <- setdiff(by, names(catalog$metadata))
  missing_dat <- setdiff(by, names(data))

  if (length(missing_cat)) {
    stop("Join key(s) not found in catalog: ", paste(missing_cat, collapse = ", "), call. = FALSE)
  }
  if (length(missing_dat)) {
    stop("Join key(s) not found in data: ", paste(missing_dat, collapse = ", "), call. = FALSE)
  }

  # Store original file order
  orig_files <- catalog$files

  # Merge using base R merge
  new_meta <- merge(
    catalog$metadata,
    data,
    by = by,
    all.x = TRUE,
    suffixes = c("", suffix)
  )

  # Restore original row order
  new_meta <- new_meta[match(orig_files, new_meta$file), ]
  rownames(new_meta) <- NULL

  catalog$metadata <- new_meta
  catalog
}

# -----------------------------------------------------------------------------
# Assay Mapping
# -----------------------------------------------------------------------------

#' Map files to assay types
#'
#' Specify which files correspond to which statistical assays (beta, se, var, t).
#' This mapping is used by [as_gds.image_catalog()] to construct the GDS.
#'
#' @param catalog An image_catalog object.
#' @param ... Named arguments where names are assay types and values are either:
#'   - A regex pattern to match against filenames, OR
#'   - A named list like `list(column = value)` for exact column matching
#'
#' @return Updated image_catalog with assay_map.
#' @export
#'
#' @examples
#' \dontrun{
#' # By filename pattern
#' catalog <- map_assays(catalog,
#'   beta = "cope[0-9]+\\.nii",
#'   se = "varcope[0-9]+\\.nii"
#' )
#'
#' # By metadata column value
#' catalog <- map_assays(catalog,
#'   beta = list(stat_type = "cope"),
#'   se = list(stat_type = "varcope")
#' )
#' }
map_assays <- function(catalog, ...) {
  stopifnot(inherits(catalog, "image_catalog"))

  maps <- list(...)

  if (!length(maps)) {
    stop("No assay mappings provided", call. = FALSE)
  }

  valid_assays <- c("beta", "var", "se", "t", "z", "p", "df")

  unknown <- setdiff(names(maps), valid_assays)
  if (length(unknown)) {
    warning("Non-standard assay types (will still be mapped): ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }

  catalog$assay_map <- maps
  catalog
}

#' Resolve assay map to file indices
#'
#' @param catalog An image_catalog object
#' @param assay_name Name of the assay to resolve
#' @return Integer vector of file indices
#' @keywords internal
.resolve_assay_map <- function(catalog, assay_name) {
  spec <- catalog$assay_map[[assay_name]]
  if (is.null(spec)) return(integer(0))

  if (is.character(spec)) {
    # Regex match on basename
    which(grepl(spec, catalog$metadata$basename, perl = TRUE))
  } else if (is.list(spec)) {
    # Column == value match
    col_nm <- names(spec)[1]
    col_val <- spec[[1]]
    if (!col_nm %in% names(catalog$metadata)) {
      stop("Column '", col_nm, "' not found in catalog metadata", call. = FALSE)
    }
    which(catalog$metadata[[col_nm]] == col_val)
  } else {
    stop("Invalid assay map specification for '", assay_name, "'", call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

#' Validate catalog consistency
#'
#' Check that the catalog has consistent file coverage across subjects
#' and report any missing or extra files.
#'
#' @param x An image_catalog object.
#' @param by Column name to group by for consistency check. Default `"subject"`.
#' @param expect Optional character vector of expected file basenames per subject.
#' @param ... Additional arguments (ignored).
#'
#' @return A `catalog_validation_report` object.
#' @export
#'
#' @examples
#' \dontrun{
#' report <- validate(catalog)
#' print(report)
#'
#' report <- validate(catalog, expect = c("cope1.nii.gz", "varcope1.nii.gz"))
#' }
validate.image_catalog <- function(x,
                                   by = "subject",
                                   expect = NULL,
                                   ...) {
  catalog <- x
  stopifnot(inherits(catalog, "image_catalog"))

  meta <- catalog$metadata
  issues <- list()

  # Check if grouping column exists
  if (!by %in% names(meta)) {
    issues$missing_group_column <- by
    return(.catalog_validation_report(catalog, issues))
  }

  # Remove rows with NA in grouping column
  valid_rows <- !is.na(meta[[by]])
  if (!all(valid_rows)) {
    issues$na_in_group_column <- sum(!valid_rows)
  }
  meta_valid <- meta[valid_rows, , drop = FALSE]

  if (nrow(meta_valid) == 0) {
    issues$no_valid_groups <- TRUE
    return(.catalog_validation_report(catalog, issues))
  }

  groups <- split(meta_valid, meta_valid[[by]])

  # Count files per group
  counts <- vapply(groups, nrow, integer(1))

  # Check for uneven counts
  if (length(unique(counts)) > 1L) {
    issues$uneven_counts <- data.frame(
      group = names(counts),
      n_files = counts,
      stringsAsFactors = FALSE
    )
  }

  # Check for consistent file types across groups
  if ("basename" %in% names(meta_valid) && length(groups) > 1L) {
    # Get unique basenames per group
    basenames_per_group <- lapply(groups, function(g) {
      sort(unique(g$basename))
    })

    # Use first group as reference
    ref_pattern <- basenames_per_group[[1]]

    inconsistent <- list()
    for (i in seq_along(basenames_per_group)[-1]) {
      grp_name <- names(basenames_per_group)[i]
      grp_basenames <- basenames_per_group[[i]]

      if (!setequal(grp_basenames, ref_pattern)) {
        inconsistent[[grp_name]] <- list(
          missing = setdiff(ref_pattern, grp_basenames),
          extra = setdiff(grp_basenames, ref_pattern)
        )
      }
    }

    if (length(inconsistent)) {
      issues$inconsistent_files <- inconsistent
    }
  }

  # Check assay map coverage
  if (length(catalog$assay_map)) {
    empty_assays <- character()
    for (nm in names(catalog$assay_map)) {
      idx <- .resolve_assay_map(catalog, nm)
      if (!length(idx)) {
        empty_assays <- c(empty_assays, nm)
      }
    }
    if (length(empty_assays)) {
      issues$empty_assay_map <- empty_assays
    }
  }

  # Check for missing expected files
  if (!is.null(expect)) {
    missing_expected <- list()
    for (grp_name in names(groups)) {
      grp <- groups[[grp_name]]
      present <- grp$basename
      missing <- setdiff(expect, present)
      if (length(missing)) {
        missing_expected[[grp_name]] <- missing
      }
    }
    if (length(missing_expected)) {
      issues$missing_expected <- missing_expected
    }
  }

  .catalog_validation_report(catalog, issues)
}

#' Create validation report object
#' @keywords internal
.catalog_validation_report <- function(catalog, issues) {
  n_subjects <- NA_integer_
  if ("subject" %in% names(catalog$metadata)) {
    subjs <- catalog$metadata$subject
    n_subjects <- length(unique(subjs[!is.na(subjs)]))
  }

  structure(
    list(
      valid = length(issues) == 0L,
      n_files = length(catalog$files),
      n_subjects = n_subjects,
      issues = issues
    ),
    class = "catalog_validation_report"
  )
}

#' Print validation report
#'
#' @param x A catalog_validation_report object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the report.
#' @method print catalog_validation_report
#' @rdname print.catalog_validation_report
#' @export
print.catalog_validation_report <- function(x, ...) {
  cat("Catalog Validation Report\n")
  cat("  Files: ", x$n_files, "\n", sep = "")
  if (!is.na(x$n_subjects)) {
    cat("  Subjects: ", x$n_subjects, "\n", sep = "")
  }

  if (x$valid) {
    cat("  Status: VALID\n")
  } else {
    cat("  Status: ISSUES FOUND\n")
    for (nm in names(x$issues)) {
      issue <- x$issues[[nm]]
      if (nm == "uneven_counts") {
        cat("    - Uneven file counts across groups:\n")
        for (i in seq_len(nrow(issue))) {
          cat("        ", issue$group[i], ": ", issue$n_files[i], " files\n", sep = "")
        }
      } else if (nm == "inconsistent_files") {
        cat("    - Inconsistent files across groups:\n")
        for (grp in names(issue)) {
          info <- issue[[grp]]
          if (length(info$missing)) {
            cat("        ", grp, " missing: ", paste(info$missing, collapse = ", "), "\n", sep = "")
          }
          if (length(info$extra)) {
            cat("        ", grp, " extra: ", paste(info$extra, collapse = ", "), "\n", sep = "")
          }
        }
      } else if (nm == "empty_assay_map") {
        cat("    - Assay mappings with no matching files: ", paste(issue, collapse = ", "), "\n", sep = "")
      } else if (nm == "missing_expected") {
        cat("    - Missing expected files:\n")
        for (grp in names(issue)) {
          cat("        ", grp, ": ", paste(issue[[grp]], collapse = ", "), "\n", sep = "")
        }
      } else if (nm == "missing_group_column") {
        cat("    - Grouping column '", issue, "' not found in metadata\n", sep = "")
      } else if (nm == "na_in_group_column") {
        cat("    - ", issue, " files have NA in grouping column\n", sep = "")
      } else {
        cat("    - ", nm, "\n", sep = "")
      }
    }
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# Serialization
# -----------------------------------------------------------------------------

#' Write catalog to JSON file
#'
#' Save an image catalog to a JSON file for later retrieval with [read_catalog()].
#'
#' @param catalog An image_catalog object.
#' @param file Output file path (should end in `.json`).
#'
#' @return Invisibly returns the file path.
#' @export
#'
#' @examples
#' \dontrun{
#' write_catalog(catalog, "my_study_manifest.json")
#' }
write_catalog <- function(catalog, file) {
  stopifnot(inherits(catalog, "image_catalog"))

  # Convert to serializable structure
  data <- list(
    version = "1.0.0",
    root_dir = catalog$root_dir,
    pattern = catalog$pattern,
    created = format(catalog$created, "%Y-%m-%dT%H:%M:%S"),
    assay_map = catalog$assay_map,
    metadata = as.list(catalog$metadata),
    files = catalog$files
  )

  jsonlite::write_json(data, file, auto_unbox = TRUE, pretty = TRUE)
  invisible(file)
}

#' Read catalog from JSON file
#'
#' Load an image catalog previously saved with [write_catalog()].
#'
#' @param file Path to catalog JSON file.
#'
#' @return An image_catalog object.
#' @export
#'
#' @examples
#' \dontrun{
#' catalog <- read_catalog("my_study_manifest.json")
#' }
read_catalog <- function(file) {
  if (!file.exists(file)) {
    stop("Catalog file does not exist: ", file, call. = FALSE)
  }

  data <- jsonlite::read_json(file, simplifyVector = TRUE)

  # Reconstruct metadata data.frame
  metadata <- as.data.frame(data$metadata, stringsAsFactors = FALSE)

  # Handle NULL assay_map
  assay_map <- data$assay_map
  if (is.null(assay_map)) assay_map <- list()

  new_image_catalog(
    files = data$files,
    metadata = metadata,
    root_dir = data$root_dir,
    pattern = data$pattern,
    assay_map = assay_map
  )
}

# -----------------------------------------------------------------------------
# Subsetting and Query
# -----------------------------------------------------------------------------

#' Subset an image catalog
#'
#' Filter catalog files based on metadata conditions.
#'
#' @param x An image_catalog object.
#' @param subset Logical expression evaluated in the context of metadata columns.
#' @param ... Additional arguments (ignored).
#'
#' @return A filtered image_catalog.
#' @export
#'
#' @examples
#' \dontrun{
#' # Keep only files from subject "01"
#' sub_catalog <- subset(catalog, subject == "01")
#'
#' # Multiple conditions
#' sub_catalog <- subset(catalog, subject %in% c("01", "02") & run == "1")
#' }
subset.image_catalog <- function(x, subset, ...) {
  meta <- x$metadata

  if (missing(subset)) {
    return(x)
  }

  # Evaluate subset expression in metadata context
  e <- substitute(subset)
  r <- eval(e, meta, parent.frame())

  if (!is.logical(r)) {
    stop("'subset' must evaluate to a logical vector", call. = FALSE)
  }

  # Handle NAs - treat as FALSE
  r[is.na(r)] <- FALSE

  # Filter
  new_files <- x$files[r]
  new_meta <- meta[r, , drop = FALSE]
  rownames(new_meta) <- NULL

  new_image_catalog(
    files = new_files,
    metadata = new_meta,
    root_dir = x$root_dir,
    pattern = x$pattern,
    assay_map = x$assay_map
  )
}

#' Get unique values of a metadata column
#'
#' @param x An image_catalog object.
#' @param incomparables Passed to base unique (ignored for catalog method).
#' @param column Name of metadata column to get unique values from.
#' @param ... Additional arguments (ignored).
#'
#' @return Unique values from the specified column.
#' @export
unique.image_catalog <- function(x, incomparables = FALSE, column = "subject", ...) {
  if (!column %in% names(x$metadata)) {
    stop("Column '", column, "' not found in catalog metadata", call. = FALSE)
  }
  unique(x$metadata[[column]])
}

# -----------------------------------------------------------------------------
# GDS Integration
# -----------------------------------------------------------------------------
#' Convert image_catalog to GDS plan
#'
#' Create a GDS plan from an image catalog. The catalog's assay mappings
#' determine which files are used for beta and se assays, and subject-level
#' metadata from the catalog is attached as col_data.
#'
#' @param x An image_catalog object with assay mappings set via [map_assays()].
#' @param mask Optional mask file path or NeuroVol object.
#' @param ... Additional arguments passed to [gds()].
#'
#' @return A gds_plan object.
#' @export
#'
#' @examples
#' \dontrun{
#' catalog <- image_catalog("study/", "sub-*/stats/*.nii.gz",
#'                          path_regex = "sub-(?<subject>[^/]+)") |>
#'   map_assays(beta = "cope", se = "varcope")
#'
#' plan <- as_gds(catalog, mask = "mask.nii.gz")
#' result <- compute(plan)
#' }
as_gds.image_catalog <- function(x, mask = NULL, ...) {
  stopifnot(inherits(x, "image_catalog"))

  # Validate assay map exists
  if (!length(x$assay_map)) {
    stop("Catalog must have assay mappings. Use map_assays() first.", call. = FALSE)
  }

  # Resolve files for each assay
  assay_files <- list()
  for (nm in names(x$assay_map)) {
    idx <- .resolve_assay_map(x, nm)
    if (!length(idx)) {
      stop("No files match assay '", nm, "' mapping", call. = FALSE)
    }
    assay_files[[nm]] <- x$files[idx]
  }

  # Require at least beta
  if (!"beta" %in% names(assay_files)) {
    stop("Catalog must map at least 'beta' assay", call. = FALSE)
  }

  # Build nifti_source specification
  src_args <- list(beta = assay_files$beta)
  if ("se" %in% names(assay_files)) {
    src_args$se <- assay_files$se
  }

  src <- do.call(nifti_source, src_args)

  # Build plan
  plan <- gds(src, format = "nifti", mask = mask, ...)

  # Attach col_data from catalog metadata if subject column exists
  if ("subject" %in% names(x$metadata)) {
    # Aggregate to one row per subject (use beta files' subjects)
    beta_idx <- .resolve_assay_map(x, "beta")
    beta_meta <- x$metadata[beta_idx, , drop = FALSE]

    # Keep unique subjects
    subject_meta <- beta_meta[!duplicated(beta_meta$subject), , drop = FALSE]
    rownames(subject_meta) <- subject_meta$subject

    # Remove file-specific columns for col_data
    drop_cols <- c("file", "basename", "relpath")
    keep_cols <- setdiff(names(subject_meta), drop_cols)
    subject_meta <- subject_meta[, keep_cols, drop = FALSE]

    # Only attach if there are useful columns beyond subject
    if (ncol(subject_meta) > 1) {
      plan <- with_col_data(plan, subject_meta)
    }
  }

  plan
}
