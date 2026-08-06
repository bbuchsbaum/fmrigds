# Reusable direct scanner for terminal analyses ---------------------------

.scan_compiled_plan <- function(compiled,
                                assays,
                                block_size = 1024L,
                                initialize,
                                update,
                                finalize,
                                on_error = NULL) {
  if (!inherits(compiled, "gds_examination_plan")) {
    stop("compiled must be a gds_examination_plan.", call. = FALSE)
  }
  if (!identical(compiled$scan_strategy, "direct")) {
    stop(
      "This examination prefix requires ", compiled$scan_strategy,
      " execution: ", paste(compiled$scan_reasons, collapse = "; "),
      call. = FALSE
    )
  }
  callbacks <- list(initialize = initialize, update = update, finalize = finalize)
  if (!all(vapply(callbacks, is.function, logical(1)))) {
    stop("initialize, update, and finalize must be functions.", call. = FALSE)
  }
  if (!is.null(on_error) && !is.function(on_error)) {
    stop("on_error must be NULL or a function.", call. = FALSE)
  }
  block_size <- as.integer(block_size)
  if (length(block_size) != 1L || is.na(block_size) || block_size < 1L) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }

  plan <- compiled$plan
  adapter <- get_adapter(plan$source$adapter)
  selected <- compiled$axis_selection
  source_samples <- selected$sample
  if (!length(source_samples)) stop("Examination selection contains no samples.", call. = FALSE)

  context <- list(
    compiled = compiled,
    subjects = plan$source$probe$subjects[selected$subject],
    contrasts = plan$source$probe$contrasts[selected$contrast],
    source_samples = source_samples,
    sample_labels = .plan_sample_labels(plan)[source_samples],
    col_data = .scan_subset_frame(
      col_data(plan),
      plan$source$probe$subjects[selected$subject]
    ),
    contrast_data = .scan_subset_frame(
      contrast_data(plan),
      plan$source$probe$contrasts[selected$contrast]
    )
  )
  receipt <- list(
    scan_strategy = "direct",
    adapter = plan$source$adapter,
    adapter_reads = 0L,
    bytes_read = 0,
    blocks = 0L,
    samples_requested = length(source_samples),
    samples_emitted = 0L,
    peak_rss_bytes = .process_rss_bytes(),
    started_at = Sys.time(),
    elapsed_seconds = NA_real_
  )

  handle <- adapter$open(plan$source$source)
  on.exit({
    if (!is.null(adapter$close)) {
      try(adapter$close(handle), silent = TRUE)
    }
  }, add = TRUE)

  state <- initialize(context)
  run <- tryCatch(
    {
      starts <- seq.int(1L, length(source_samples), by = block_size)
      for (block_number in seq_along(starts)) {
        start <- starts[[block_number]]
        stop_at <- min(start + block_size - 1L, length(source_samples))
        ordinal <- seq.int(start, stop_at)
        source_index <- source_samples[ordinal]
        block <- list(
          number = as.integer(block_number),
          ordinal = ordinal,
          source_sample = source_index,
          source_subject = selected$subject,
          source_contrast = selected$contrast,
          sample_label = context$sample_labels[ordinal]
        )
        arrays <- .read_examination_block(
          adapter,
          handle,
          plan,
          assays,
          list(
            sample = source_index,
            subject = selected$subject,
            contrast = selected$contrast
          )
        )
        receipt$adapter_reads <- receipt$adapter_reads + 1L
        receipt$blocks <- receipt$blocks + 1L
        receipt$bytes_read <- receipt$bytes_read + .arrays_bytes(arrays)
        if (block_number %% 32L == 0L) {
          receipt$peak_rss_bytes <- max(
            receipt$peak_rss_bytes,
            .process_rss_bytes(),
            na.rm = TRUE
          )
        }

        prefix_result <- .apply_examination_block_prefix(
          arrays,
          compiled,
          source_index,
          context
        )
        arrays <- prefix_result$arrays
        block$source_sample <- prefix_result$source_sample
        block$sample_label <- prefix_result$sample_label
        receipt$samples_emitted <- receipt$samples_emitted + dim(arrays[[1L]])[1L]
        state <- update(state, arrays, block, context)
      }
      receipt$elapsed_seconds <- as.numeric(
        difftime(Sys.time(), receipt$started_at, units = "secs")
      )
      receipt$peak_rss_bytes <- max(
        receipt$peak_rss_bytes,
        .process_rss_bytes(),
        na.rm = TRUE
      )
      value <- finalize(state, receipt, context)
      list(value = value, receipt = receipt)
    },
    error = function(e) {
      receipt$elapsed_seconds <- as.numeric(
        difftime(Sys.time(), receipt$started_at, units = "secs")
      )
      if (!is.null(on_error)) {
        try(on_error(state, e, receipt, context), silent = TRUE)
      }
      stop(e)
    }
  )

  structure(run, class = "gds_scan_result")
}

.process_rss_bytes <- function() {
  if (requireNamespace("ps", quietly = TRUE)) {
    value <- tryCatch(
      as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]),
      error = function(e) NA_real_
    )
    return(value)
  }
  value <- tryCatch(
    suppressWarnings(as.numeric(system2(
      "ps",
      c("-o", "rss=", "-p", as.character(Sys.getpid())),
      stdout = TRUE,
      stderr = FALSE
    ))),
    error = function(e) NA_real_
  )
  if (!length(value) || !is.finite(value[1L])) NA_real_ else value[1L] * 1024
}

.read_examination_block <- function(adapter, handle, plan, assays, block) {
  columns <- plan$meta$adapter_columns %||% list(
    effect_cols = NULL,
    subject_col = NULL,
    sample_col = NULL,
    contrast_col = NULL
  )
  adapter$read(
    handle,
    assays = assays,
    block = block,
    effect_cols = columns$effect_cols,
    subject_col = columns$subject_col,
    sample_col = columns$sample_col,
    contrast_col = columns$contrast_col,
    mask_idx = plan$source$probe$mask_idx,
    spatial_dim = plan$source$probe$spatial_dim,
    n_contrasts = plan$source$probe$n_contrasts,
    temporal_policy = plan$meta$temporal_policy %||% NULL,
    contrast_matrix = plan$meta$contrast_matrix %||% NULL,
    contrast_names = plan$meta$contrast_names %||% NULL
  )
}

.apply_examination_block_prefix <- function(arrays,
                                            compiled,
                                            source_sample,
                                            context) {
  if (!length(compiled$execution_prefix)) {
    return(list(
      arrays = arrays,
      source_sample = source_sample,
      sample_label = .plan_sample_labels(compiled$plan)[source_sample]
    ))
  }
  plan <- compiled$plan
  block_plan <- plan
  block_plan$nodes <- compiled$execution_prefix
  block_plan$source$probe$subjects <- context$subjects
  block_plan$source$probe$contrasts <- context$contrasts
  block_plan$source$probe$dims <- gds_dims(
    sample = length(source_sample),
    subject = length(context$subjects),
    contrast = length(context$contrasts)
  )
  block_space <- .subset_space(plan$source$probe$space, source_sample)
  block_plan$source$probe$space <- block_space
  block_rows <- row_data(plan)
  if (!is.null(block_rows)) block_rows <- block_rows[source_sample, , drop = FALSE]

  result <- .apply_plan_nodes(
    arrays,
    block_plan,
    block_space,
    context$subjects,
    col_data = context$col_data,
    row_data = block_rows,
    contrast_data = context$contrast_data
  )
  local <- result$subset$samples %||% seq_along(source_sample)
  list(
    arrays = result$arrays,
    source_sample = source_sample[local],
    sample_label = .plan_sample_labels(plan)[source_sample[local]]
  )
}

.scan_subset_frame <- function(data, ids) {
  if (is.null(data)) return(NULL)
  if (!is.null(rownames(data)) && all(ids %in% rownames(data))) {
    return(data[ids, , drop = FALSE])
  }
  data
}

.arrays_bytes <- function(arrays) {
  sum(vapply(arrays, function(x) as.numeric(utils::object.size(x)), numeric(1)))
}
