# Neuropublish adapter -----------------------------------------------------

.np_local_id_pattern <- "^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$"

.np_slug <- function(x, fallback = "value") {
  x <- enc2utf8(as.character(x %||% fallback))
  x <- tolower(gsub("[^A-Za-z0-9._-]+", "-", x))
  x <- gsub("(^[-._]+|[-._]+$)", "", x)
  if (!nzchar(x)) fallback else x
}

.np_stable_id <- function(prefix, label, key = label) {
  slug <- substr(.np_slug(label), 1L, 54L)
  suffix <- substr(.portable_sha256(list(prefix = prefix, key = key))$value, 1L, 12L)
  substr(paste(prefix, slug, suffix, sep = "-"), 1L, 128L)
}

.np_safe_existing_id <- function(id, prefix = "provenance") {
  id <- as.character(id)
  if (length(id) == 1L && !is.na(id) && grepl(.np_local_id_pattern, id)) return(id)
  .np_stable_id(prefix, id, id)
}

.np_assert_unique_ids <- function(records, what) {
  ids <- if (length(records)) vapply(records, `[[`, character(1L), "id") else character()
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Generated duplicate or empty ", what, " id(s).", call. = FALSE)
  }
  invisible(TRUE)
}

.publication_reduction_receipt <- function(x) {
  receipts <- metadata(x)$reduction_receipts %||% list()
  if (length(receipts) != 1L) {
    stop(
      "A publishable group GDS must contain exactly one portable reduction receipt; found ",
      length(receipts), ". Recompute the realised group result with one reduce() node.",
      call. = FALSE
    )
  }
  receipt <- receipts[[1L]]
  if (!identical(receipt$status, "portable") ||
      is.null(receipt$planNodeId) || is.null(receipt$nInputSubjects)) {
    stop(
      "The reduction receipt is incomplete or non-portable; recompute it with the current fmrigds version.",
      call. = FALSE
    )
  }
  receipt
}

.admit_neuropublish_gds <- function(x) {
  if (!inherits(x, "gds")) stop("`x` must inherit from gds.", call. = FALSE)
  if (!inherits(space(x), "space_voxel")) {
    stop("Neuropublish currently accepts only realised voxel-space GDS results.", call. = FALSE)
  }
  if (!identical(as.character(subjects(x)), "meta")) {
    stop(
      "Neuropublish currently accepts one realised group output, not subject-level arrays.",
      call. = FALSE
    )
  }
  dimensions <- lapply(assays(x), dim)
  if (!length(dimensions) || any(vapply(dimensions, function(d) {
    length(d) != 3L || d[[2L]] != 1L
  }, logical(1L)))) {
    stop("Every published assay must be a realised 3-D GDS array with group subject dimension 1.", call. = FALSE)
  }
  .publication_reduction_receipt(x)
}

.known_ordinary_p_reducers <- c(
  "meta:fe", "meta:re", "meta:fe_reg", "meta:re_reg",
  "ols:voxelwise", "perm:onesample", "perm:twosample",
  "combine:stouffer", "combine:fisher", "combine:lancaster"
)

.fmrigds_measure_descriptor <- function(assay_name,
                                        reducer_id,
                                        measure_overrides = NULL,
                                        label_overrides = NULL) {
  override <- if (!is.null(measure_overrides)) measure_overrides[[assay_name]] %||% NULL else NULL
  label_override <- if (!is.null(label_overrides)) label_overrides[[assay_name]] %||% NULL else NULL
  if (!is.null(override)) {
    return(list(measure = as.character(override), label = label_override %||% assay_name, trusted = FALSE))
  }
  effect <- identical(assay_name, "beta_g") || identical(assay_name, "beta") || startsWith(assay_name, "coef:")
  standard_error <- identical(assay_name, "se_g") || startsWith(assay_name, "se_coef:")
  t_stat <- identical(assay_name, "t_g") || startsWith(assay_name, "t_coef:")
  z_stat <- identical(assay_name, "z_g") || startsWith(assay_name, "z_coef:")
  ordinary_p <- (identical(assay_name, "p_g") || startsWith(assay_name, "p_coef:")) &&
    reducer_id %in% .known_ordinary_p_reducers

  if (effect) return(list(
    measure = "org.neuropublish.measure/effect",
    label = label_override %||% "Effect",
    trusted = TRUE
  ))
  if (standard_error) return(list(
    measure = "org.neuropublish.measure/standard-error",
    label = label_override %||% "Standard error",
    trusted = TRUE
  ))
  if (t_stat) return(list(
    measure = "org.neuropublish.measure/t-statistic",
    label = label_override %||% "t statistic",
    trusted = TRUE
  ))
  if (z_stat) return(list(
    measure = "org.neuropublish.measure/z-statistic",
    label = label_override %||% "z statistic",
    trusted = TRUE
  ))
  if (ordinary_p) return(list(
    measure = "org.neuropublish.measure/p-value",
    label = label_override %||% "p value",
    trusted = TRUE
  ))

  known <- list(
    var_g = c("org.bbuchsbaum.fmrigds.measure/sampling-variance", "Sampling variance"),
    var = c("org.bbuchsbaum.fmrigds.measure/sampling-variance", "Sampling variance"),
    tau2 = c(
      "org.bbuchsbaum.fmrigds.measure/between-study-heterogeneity-variance",
      "Between-study heterogeneity (\u03c4\u00b2)"
    ),
    Q = c("org.bbuchsbaum.fmrigds.measure/cochran-q", "Cochran Q"),
    I2 = c("org.bbuchsbaum.fmrigds.measure/i-squared", "I\u00b2"),
    n_eff = c("org.bbuchsbaum.fmrigds.measure/effective-sample-size", "Effective sample size"),
    p_perm = c("org.bbuchsbaum.fmrigds.measure/permutation-p-value", "Permutation p value"),
    p_fwer = c("org.bbuchsbaum.fmrigds.measure/fwer-adjusted-p-value", "FWER-adjusted p value")
  )
  if (!is.null(known[[assay_name]])) {
    return(list(
      measure = known[[assay_name]][[1L]],
      label = label_override %||% known[[assay_name]][[2L]],
      trusted = FALSE
    ))
  }
  # Unknown names never gain trusted semantics through resemblance. They stay
  # renderable as an explicitly untrusted producer-namespaced scalar map.
  list(
    measure = "org.bbuchsbaum.fmrigds.measure/custom-scalar-map",
    label = label_override %||% assay_name,
    trusted = FALSE
  )
}

.default_neuropublish_assays <- function(x) {
  available <- names(assays(x))
  coefficient <- grep("^coef:", available, value = TRUE)
  selected <- character()
  if (length(coefficient)) {
    selected <- coefficient
    terms <- sub("^coef:", "", coefficient)
    se <- paste0("se_coef:", terms)
    t <- paste0("t_coef:", terms)
    z <- paste0("z_coef:", terms)
    selected <- c(selected, se[se %in% available])
    primary <- ifelse(t %in% available, t, ifelse(z %in% available, z, NA_character_))
    selected <- c(selected, primary[!is.na(primary)])
  } else {
    effect <- if ("beta_g" %in% available) "beta_g" else if ("beta" %in% available) "beta" else character()
    se <- if ("se_g" %in% available) "se_g" else character()
    statistic <- if ("t_g" %in% available) "t_g" else if ("z_g" %in% available) "z_g" else character()
    selected <- c(effect, se, statistic)
  }
  if ("tau2" %in% available) selected <- c(selected, "tau2")
  unique(selected)
}

.neuropublish_categorical_assays <- function(x) {
  info <- metadata(x)$examination %||% list()
  unique(c(
    metadata(x)$categorical_assays %||% character(),
    metadata(x)$non_interpolable_assays %||% character(),
    info$categorical_assays %||% character(),
    info$non_interpolable_assays %||% character()
  ))
}

.neuropublish_provenance_spec <- function(x, field_ids) {
  raw_provenance <- metadata(x)$provenance %||% NULL
  provenance <- .normalize_provenance(raw_provenance)
  graph_problems <- validate_provenance_graph(raw_provenance, error = FALSE)
  if (!isTRUE(graph_problems) && length(graph_problems) &&
      !identical(provenance$status, "legacy-incomplete")) {
    stop(
      "Cannot publish an invalid provenance graph: ",
      paste(graph_problems, collapse = "; "),
      call. = FALSE
    )
  }
  entities <- .portable_source_entities(provenance$entities %||% list())
  activities <- provenance$graph %||% list()
  internal_ids <- c(
    if (length(entities)) vapply(entities, `[[`, character(1L), "id") else character(),
    if (length(activities) && all(vapply(activities, .is_portable_activity, logical(1L)))) {
      vapply(activities, `[[`, character(1L), "id")
    } else character()
  )
  mapped <- vapply(internal_ids, .np_safe_existing_id, character(1L))
  names(mapped) <- internal_ids
  if (anyDuplicated(mapped)) {
    for (i in seq_along(mapped)) {
      mapped[[i]] <- .np_stable_id("provenance", internal_ids[[i]], list(internal_ids[[i]], i))
    }
  }

  entities_out <- lapply(entities, function(entity) {
    entity$id <- unname(mapped[[entity$id]] %||% entity$id)
    entity
  })
  activities_out <- list()
  edges <- list()
  if (length(activities) && all(vapply(activities, .is_portable_activity, logical(1L)))) {
    activities_out <- lapply(activities, function(activity) {
      list(
        id = unname(mapped[[activity$id]] %||% activity$id),
        label = activity$op,
        semanticId = activity$semanticId,
        params = activity$params,
        inputs = unname(vapply(
          as.character(activity$inputs),
          function(id) mapped[[id]] %||% id,
          character(1L)
        )),
        timestamp = activity$timestamp,
        software = activity$software,
        digest = activity$digest
      )
    })
    edges <- unlist(lapply(activities_out, function(activity) {
      lapply(activity$inputs, function(input) list(from = input, to = activity$id, role = "used"))
    }), recursive = FALSE)
    head_ids <- vapply(as.character(provenance$heads %||% character()), function(id) mapped[[id]] %||% id, character(1L))
    if (length(head_ids)) {
      edges <- c(edges, unlist(lapply(head_ids, function(head) {
        lapply(field_ids, function(field) list(from = head, to = field, role = "generated"))
      }), recursive = FALSE))
    }
  }
  list(
    status = provenance$status,
    graphReceipt = provenance$graphReceipt,
    entities = unname(entities_out),
    activities = unname(activities_out),
    edges = unname(edges)
  )
}

.neuropublish_warning_spec <- function(x, provenance) {
  out <- list()
  add <- function(id, message, concerns = NULL) {
    out[[length(out) + 1L]] <<- list(id = id, message = message, concerns = concerns)
  }
  if (isTRUE(metadata(x)$synthetic_var)) {
    add("synthetic-variance", "Input variance was a synthetic placeholder rather than measured uncertainty.")
  }
  if (isTRUE(metadata(x)$sample_labels_synthetic)) {
    add("synthetic-sample-labels", "Sample labels were generated and do not carry anatomical identity.")
  }
  statuses <- if (length(provenance$entities)) {
    vapply(provenance$entities, function(entity) entity$identityStatus %||% "unavailable", character(1L))
  } else character()
  if (!length(statuses) || any(statuses != "verified")) {
    add(
      "source-identity-incomplete",
      "One or more input sources lack a verified SHA-256 byte identity; local paths were not published."
    )
  }
  if (!identical(provenance$status, "complete")) {
    add(
      "provenance-incomplete",
      paste0("The fmrigds provenance state is '", provenance$status, "'; no historical edges were inferred from list order.")
    )
  }
  .np_assert_unique_ids(out, "warning")
  out
}

# Pure, producer-owned mapping. This function deliberately constructs no
# Neuropublish objects, writes no assets, and performs no staging/network work.
.neuropublish_spec_gds <- function(x,
                                   assays = NULL,
                                   diagnostics = character(),
                                   measure_overrides = NULL,
                                   label_overrides = NULL,
                                   analysis_label = NULL) {
  receipt <- .admit_neuropublish_gds(x)
  available <- names(fmrigds::assays(x))
  requested <- assays %||% .default_neuropublish_assays(x)
  requested <- unique(c(as.character(requested), as.character(diagnostics)))
  if (!length(requested)) stop("No publishable assays were selected.", call. = FALSE)
  missing <- setdiff(requested, available)
  if (length(missing)) stop("Unknown requested assay(s): ", paste(missing, collapse = ", "), call. = FALSE)
  categorical <- intersect(requested, .neuropublish_categorical_assays(x))
  if (length(categorical)) {
    stop(
      "Categorical/non-interpolable assay(s) cannot be published as scalar volumes: ",
      paste(categorical, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(measure_overrides) && (is.null(names(measure_overrides)) || any(!nzchar(names(measure_overrides))))) {
    stop("`measure_overrides` must be a named character vector/list.", call. = FALSE)
  }
  if (!is.null(label_overrides) && (is.null(names(label_overrides)) || any(!nzchar(names(label_overrides))))) {
    stop("`label_overrides` must be a named character vector/list.", call. = FALSE)
  }

  analysis_id <- .np_stable_id("analysis", receipt$reducerId, list(
    receipt$reducerId, receipt$planNodeId, receipt$planDigest$value %||% receipt$planDigest$status
  ))
  domain_id <- .np_stable_id("domain", "volume", list(
    dim = as.integer(space(x)$dim),
    affine = as.numeric(t(unname(space(x)$affine))),
    template = space(x)$template_id %||% "unknown"
  ))

  fields <- list()
  estimand_keys <- list()
  field_order <- list()
  for (assay_name in requested) {
    array <- fmrigds::assay(x, assay_name)
    if (!is.numeric(array)) stop("Assay '", assay_name, "' is not a numeric scalar field.", call. = FALSE)
    term <- if (grepl("^(coef|se_coef|t_coef|z_coef|p_coef):", assay_name)) {
      sub("^[^:]+:", "", assay_name)
    } else NULL
    descriptor <- .fmrigds_measure_descriptor(
      assay_name,
      receipt$reducerId,
      measure_overrides,
      label_overrides
    )
    for (contrast_index in seq_along(contrasts(x))) {
      contrast <- contrasts(x)[[contrast_index]]
      estimand_key <- paste(contrast, term %||% "", sep = "\r")
      estimand_id <- .np_stable_id(
        "estimand",
        paste(contrast, term %||% "effect", sep = "-"),
        list(contrast = contrast, term = term)
      )
      estimand_keys[[estimand_key]] <- estimand_id
      field_id <- .np_stable_id(
        "field",
        paste(assay_name, contrast, sep = "-"),
        list(assay = assay_name, contrast = contrast, term = term)
      )
      asset_id <- .np_stable_id("asset", field_id, field_id)
      field_order[[estimand_id]] <- (field_order[[estimand_id]] %||% 0L) + 1L
      label <- if (identical(assay_name, "tau2")) {
        "Between-study heterogeneity (\u03c4\u00b2)"
      } else {
        paste(c(contrast, term, descriptor$label), collapse = " \u2014 ")
      }
      fields[[length(fields) + 1L]] <- list(
        id = field_id,
        assetId = asset_id,
        assay = assay_name,
        contrast = contrast,
        contrastIndex = contrast_index,
        term = term,
        estimand = estimand_id,
        measure = descriptor$measure,
        trustedMeasure = descriptor$trusted,
        label = label,
        domain = domain_id,
        order = field_order[[estimand_id]]
      )
    }
  }
  .np_assert_unique_ids(fields, "result field")
  asset_records <- lapply(fields, function(field) list(id = field$assetId))
  .np_assert_unique_ids(asset_records, "asset")

  estimands <- lapply(seq_along(estimand_keys), function(i) {
    key <- names(estimand_keys)[[i]]
    split_at <- regexpr("\r", key, fixed = TRUE)[[1L]]
    contrast <- substr(key, 1L, split_at - 1L)
    term <- substr(key, split_at + 1L, nchar(key))
    list(
      id = unname(estimand_keys[[i]]),
      label = if (nzchar(term)) paste(contrast, term, sep = " \u2014 ") else contrast,
      order = i
    )
  })
  .np_assert_unique_ids(estimands, "estimand")

  provenance <- .neuropublish_provenance_spec(x, vapply(fields, `[[`, character(1L), "id"))
  warnings <- .neuropublish_warning_spec(x, provenance)
  spec <- list(
    schemaVersion = "1.0.0",
    sensitivity = "group-level",
    domain = list(id = domain_id),
    analysis = list(
      id = analysis_id,
      label = analysis_label %||% paste("fmrigds", receipt$reducerId),
      sampleSize = receipt$nInputSubjects,
      method = receipt,
      estimands = estimands
    ),
    fields = fields,
    warnings = warnings,
    provenance = provenance
  )

  # This is both a privacy assertion and a regression tripwire for future
  # mapper edits. The pure spec is exactly the information boundary used below.
  encoded <- .canonical_portable_json(spec)
  private_locators <- vapply(
    metadata(x)$provenance$entities %||% list(),
    function(entity) entity$private$locator %||% "",
    character(1L)
  )
  leaked <- private_locators[nzchar(private_locators) & vapply(private_locators, function(path) {
    grepl(path, encoded, fixed = TRUE)
  }, logical(1L))]
  if (length(leaked)) stop("Portable publication spec contains a private source locator.", call. = FALSE)
  spec
}

.neuropublish_display <- function(display, field) {
  candidate <- display[[field$id]] %||% display[[field$assay]] %||% NULL
  if (is.null(candidate)) return(NULL)
  if (!is.list(candidate)) stop("Display recommendations must be named lists.", call. = FALSE)
  threshold <- candidate$threshold
  window <- candidate$window
  if (!is.list(threshold) || !is.list(window) || is.null(candidate$colormap)) {
    stop("Display needs threshold, window, and colormap members.", call. = FALSE)
  }
  threshold <- do.call(neuropublish::np_threshold, threshold)
  window <- do.call(neuropublish::np_window, window)
  neuropublish::np_display(
    threshold = threshold,
    window = window,
    colormap = candidate$colormap,
    opacity = candidate$opacity %||% NULL
  )
}

#' Convert a realised GDS group result to a Neuropublish result
#'
#' This producer-owned adapter publishes a deliberately small set of scalar
#' volume fields from one realised voxel-space group result. A portable
#' reduction receipt is required; subject-level data, local source paths,
#' design matrices, and custom weight values are never placed in the manifest.
#'
#' The `neuropublish` package is optional and its S3 method is registered only
#' after that namespace loads, in either package load order.
#'
#' @param x A realised group-level [`gds`] in voxel space.
#' @param title,synopsis Manifest title and synopsis.
#' @param assays Optional character vector of assays to publish. The default is
#'   effect, standard error, one primary t/z statistic, and random-effects
#'   heterogeneity when available.
#' @param diagnostics Additional assay names to include as diagnostic fields.
#' @param measure_overrides Optional named measure IDs for custom fields. These
#'   remain producer-authored, untrusted extension measures.
#' @param labels Optional named human labels keyed by assay name.
#' @param display Optional named Neuropublish display recommendations keyed by
#'   field ID or assay name. Recommendations are validated by Neuropublish and
#'   no threshold is invented when they are absent.
#' @param staging Directory used by Neuropublish to stage NIfTI assets.
#' @param space_name Optional human name for the voxel coordinate space.
#' @param analysis_label Optional label for the emitted analysis.
#' @param ... Reserved for forward-compatible adapter options.
#'
#' @return A `neuropublish_result` from the optional `neuropublish` package.
#' @examples
#' \dontrun{
#' library(neuropublish)
#' result <- as_neuropublish(compute(reduce(plan, method = "meta:re")))
#' }
#' @name as_neuropublish.gds
#'
# Delayed S3 method: fmrigds owns the mapping while neuropublish remains an
# optional package and is never loaded merely because fmrigds starts.
as_neuropublish.gds <- function(x,
                                title = "fmrigds group result",
                                synopsis = "A realised group-level fmrigds result.",
                                assays = NULL,
                                diagnostics = character(),
                                measure_overrides = NULL,
                                labels = NULL,
                                display = list(),
                                staging = NULL,
                                space_name = NULL,
                                analysis_label = NULL,
                                ...) {
  if (!requireNamespace("neuropublish", quietly = TRUE)) {
    stop(
      "The optional 'neuropublish' package is required for as_neuropublish.gds().",
      call. = FALSE
    )
  }
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("The optional 'neuroim2' package is required to stage volume assets.", call. = FALSE)
  }
  staging <- staging %||% neuropublish::np_staging_dir()
  spec <- .neuropublish_spec_gds(
    x,
    assays = assays,
    diagnostics = diagnostics,
    measure_overrides = measure_overrides,
    label_overrides = labels,
    analysis_label = analysis_label
  )

  volume_cache <- new.env(parent = emptyenv())
  volume_for <- function(field) {
    if (is.null(volume_cache[[field$assay]])) {
      volumes <- as_neurovol_list(x, assay = field$assay, by = "contrast")
      if (length(volumes) != length(contrasts(x))) {
        stop("Could not align exported volumes to the contrast axis.", call. = FALSE)
      }
      volume_cache[[field$assay]] <- volumes
    }
    volume_cache[[field$assay]][[field$contrastIndex]]
  }
  first_volume <- volume_for(spec$fields[[1L]])
  manifest <- neuropublish::np_manifest(title, synopsis, sensitivity = "group-level")
  manifest <- neuropublish::np_add(
    manifest,
    "domains",
    neuropublish::np_domain_volume(
      neuroim2::space(first_volume),
      id = spec$domain$id,
      space_name = space_name %||% space(x)$template_id %||% "unknown"
    )
  )

  for (field in spec$fields) {
    asset <- neuropublish::np_asset_volume(field$assetId, volume_for(field), staging = staging)
    manifest <- neuropublish::np_add(manifest, "assets", asset)
  }
  method <- neuropublish::np_record(
    "org.bbuchsbaum.fmrigds/reduction",
    "1.0",
    spec$analysis$method
  )
  estimands <- lapply(spec$analysis$estimands, function(estimand) {
    neuropublish::np_estimand(estimand$id, estimand$label, estimand$order)
  })
  manifest <- neuropublish::np_add(
    manifest,
    "analyses",
    neuropublish::np_analysis(
      spec$analysis$id,
      spec$analysis$label,
      method = method,
      sample_size = spec$analysis$sampleSize,
      estimands = estimands
    )
  )
  for (field in spec$fields) {
    manifest <- neuropublish::np_add(
      manifest,
      "resultFields",
      neuropublish::np_field(
        field$id,
        field$estimand,
        field$measure,
        field$domain,
        selection = list(level = "group"),
        representations = list(neuropublish::np_volume_rep(field$assetId)),
        order = field$order,
        published_display = .neuropublish_display(display, field),
        label = field$label
      )
    )
  }
  for (warning in spec$warnings) {
    manifest <- neuropublish::np_add(
      manifest,
      "warnings",
      neuropublish::np_warning(warning$id, warning$message, warning$concerns)
    )
  }

  entities <- lapply(spec$provenance$entities, function(entity) {
    record <- neuropublish::np_entity(
      entity$id,
      label = paste("fmrigds", entity$kind %||% "source"),
      hosted = FALSE
    )
    record$schema <- list(id = "org.bbuchsbaum.fmrigds/source-entity", version = "1.0")
    payload <- entity
    payload$id <- NULL
    record$payload <- payload
    record
  })
  activities <- lapply(spec$provenance$activities, function(activity) {
    payload <- activity
    payload$id <- NULL
    payload$label <- NULL
    payload$inputs <- NULL
    neuropublish::np_activity(
      activity$id,
      neuropublish::np_record("org.bbuchsbaum.fmrigds/activity", "1.0", payload),
      label = activity$label
    )
  })
  edges <- lapply(spec$provenance$edges, function(edge) {
    neuropublish::np_edge(edge$from, edge$to, edge$role)
  })
  manifest$provenance <- neuropublish::np_provenance(entities, activities, edges)
  neuropublish::np_result(manifest)
}

.register_neuropublish_method <- function() {
  register <- function() {
    if (!"neuropublish" %in% loadedNamespaces()) return(invisible(FALSE))
    registerS3method(
      "as_neuropublish",
      "gds",
      as_neuropublish.gds,
      envir = asNamespace("neuropublish")
    )
    invisible(TRUE)
  }
  if ("neuropublish" %in% loadedNamespaces()) register()
  if (!isTRUE(.gds_state$neuropublish_hook_registered)) {
    setHook(packageEvent("neuropublish", "onLoad"), function(...) register(), action = "append")
    .gds_state$neuropublish_hook_registered <- TRUE
  }
  invisible(NULL)
}
