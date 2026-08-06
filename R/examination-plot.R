# Static group-examination views ------------------------------------------

utils::globalVariables(c(
  "display", "display_alpha", "display_group", "feature",
  "influence_energy", "metric_label", "review_label", "review_status",
  "subject", "surprise_energy", "validity_concern", "value", "x", "y"
))

.validate_gds_examination <- function(x) {
  if (!inherits(x, "gds_examination")) {
    stop("x must be a gds_examination.", call. = FALSE)
  }
  invisible(x)
}

.require_examination_plotting <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Static group-examination plots require the suggested ggplot2 package.",
      call. = FALSE
    )
  }
}

#' Plot a group examination
#'
#' The default action plane places model surprise on the horizontal axis and
#' influence on the requested group estimand on the vertical axis. Secondary
#' views show subject-by-diagnostic summaries or the residual geometry. When
#' `subject` is supplied, the retained observed, expected, predictive-residual,
#' and deletion maps are shown for that subject.
#'
#' @param x A `gds_examination`.
#' @param subject Optional retained subject ID. Supplying it selects the linked
#'   subject view regardless of `type`.
#' @param type One of `"action"`, `"heatmap"`, or `"embedding"`.
#' @param contrast Optional contrast subset.
#' @param estimand Optional estimand subset.
#' @param group_by Optional column in `subject_data` mapped to point shape.
#' @param assays Optional selected-subject assays to display.
#' @param slice Optional one-based axial slice for voxel-space subject maps.
#' @param ... Reserved for future plot controls.
#'
#' @return A `ggplot` object.
#' @export
plot.gds_examination <- function(x,
                                 subject = NULL,
                                 type = c("action", "heatmap", "embedding"),
                                 contrast = NULL,
                                 estimand = NULL,
                                 group_by = NULL,
                                 assays = NULL,
                                 slice = NULL,
                                 ...) {
  .validate_gds_examination(x)
  .require_examination_plotting()
  if (!is.null(subject)) {
    return(.plot_examination_subject(
      x,
      subject = subject,
      contrast = contrast,
      estimand = estimand,
      assays = assays,
      slice = slice
    ))
  }
  type <- match.arg(type)
  switch(
    type,
    action = .plot_examination_action(
      x, contrast = contrast, estimand = estimand, group_by = group_by
    ),
    heatmap = .plot_examination_heatmap(
      x, contrast = contrast, estimand = estimand
    ),
    embedding = .plot_examination_embedding(x, group_by = group_by)
  )
}

.examination_select_values <- function(requested, available, name) {
  available <- unique(as.character(available))
  if (is.null(requested)) return(available)
  requested <- as.character(requested)
  unknown <- setdiff(requested, available)
  if (length(unknown)) {
    stop("Unknown ", name, ": ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  requested
}

.examination_screening_estimands <- function(x) {
  data <- x$estimand_data
  if ("ranking_stage" %in% names(data)) {
    data <- data[data$ranking_stage == "screening", , drop = FALSE]
  }
  data
}

.examination_action_data <- function(x,
                                     contrast = NULL,
                                     estimand = NULL,
                                     group_by = NULL) {
  contrasts <- .examination_select_values(
    contrast, x$contrast_data$contrast, "contrast"
  )
  influence <- .examination_screening_estimands(x)
  estimands <- .examination_select_values(
    estimand, influence$estimand, "estimand"
  )
  surprise <- x$contrast_data[
    x$contrast_data$contrast %in% contrasts,
    c(
      "subject", "contrast", "surprise_energy", "tail_extent",
      "surprise_eligible_n", "surprise_status", "surprise_stability",
      "zero_intercept_gain", "gain_status"
    ),
    drop = FALSE
  ]
  influence <- influence[
    influence$contrast %in% contrasts & influence$estimand %in% estimands,
    c(
      "subject", "contrast", "estimand", "mode", "influence_energy",
      "max_abs_delta_stat", "eligible_n", "status", "stability", "stable"
    ),
    drop = FALSE
  ]
  out <- merge(
    surprise,
    influence,
    by = c("subject", "contrast"),
    all = TRUE,
    sort = FALSE
  )
  out <- merge(out, x$subject_data, by = "subject", all.x = TRUE, sort = FALSE)
  out$available <- is.finite(out$surprise_energy) &
    is.finite(out$influence_energy) &
    out$surprise_status == "available" & out$status == "available"
  out$validity_concern <- factor(
    ifelse(
      !is.na(out$review_source) & out$review_source == "quality",
      "validity concern",
      "no validity concern"
    ),
    levels = c("no validity concern", "validity concern")
  )
  out$review_label <- ifelse(out$review_status == "review", out$subject, NA_character_)
  out$display_stability <- pmin(out$surprise_stability, out$stability, na.rm = TRUE)
  out$display_stability[!is.finite(out$display_stability)] <- NA_real_
  out$display_alpha <- ifelse(
    is.finite(out$display_stability),
    0.25 + 0.75 * pmax(0, pmin(1, out$display_stability)),
    0.35
  )
  if (!is.null(group_by)) {
    if (!is.character(group_by) || length(group_by) != 1L ||
        !group_by %in% names(out)) {
      stop("group_by must name one subject_data column.", call. = FALSE)
    }
    group_value <- as.character(out[[group_by]])
    group_value[is.na(group_value)] <- "(missing)"
    out$display_group <- factor(group_value)
  } else {
    out$display_group <- factor("cohort")
  }
  attr(out, "group_by") <- group_by
  out
}

.examination_empty_plot <- function(title, message, caption = NULL) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
    ggplot2::annotate("text", x = 0, y = 0, label = message, colour = "#52606d") +
    ggplot2::labs(title = title, caption = caption) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#172b4d"),
      plot.caption = ggplot2::element_text(colour = "#52606d", hjust = 0)
    )
}

.examination_caption <- function(data, embedding = NULL) {
  available <- data$available %||% rep(TRUE, nrow(data))
  modes <- unique(stats::na.omit(data$mode %||% character()))
  eligible <- c(
    data$surprise_eligible_n %||% numeric(),
    data$eligible_n %||% numeric()
  )
  eligible <- eligible[is.finite(eligible)]
  pieces <- c(
    if (length(modes)) paste0("Influence mode: ", paste(modes, collapse = ", ")),
    if (length(eligible)) paste0("eligible features: ", min(eligible), "-", max(eligible)),
    paste0("unavailable points: ", sum(!available))
  )
  if (!is.null(embedding)) {
    pieces <- c(
      pieces,
      paste0("residual geometry captured energy: ",
             formatC(100 * embedding$captured_energy, digits = 1, format = "f"), "%")
    )
  }
  paste(pieces, collapse = " | ")
}

.plot_examination_action <- function(x,
                                     contrast = NULL,
                                     estimand = NULL,
                                     group_by = NULL) {
  data <- .examination_action_data(x, contrast, estimand, group_by)
  shown <- data[data$available, , drop = FALSE]
  if (!nrow(shown)) {
    return(.examination_empty_plot(
      "Group examination: action plane",
      "Surprise or influence is unavailable for this selection.",
      .examination_caption(data)
    ))
  }
  p <- ggplot2::ggplot(
    shown,
    ggplot2::aes(
      x = surprise_energy,
      y = influence_energy,
      fill = review_status,
      colour = validity_concern,
      shape = display_group,
      alpha = display_alpha
    )
  ) +
    ggplot2::geom_vline(
      xintercept = x$config$control$review$surprise$energy_threshold,
      colour = "#9fb3c8", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = x$config$control$review$influence$energy_threshold,
      colour = "#9fb3c8", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_point(size = 3.2, stroke = 1.05, na.rm = TRUE) +
    ggplot2::geom_text(
      data = shown[!is.na(shown$review_label), , drop = FALSE],
      ggplot2::aes(label = review_label),
      inherit.aes = TRUE,
      nudge_y = 0.025 * max(1, diff(range(shown$influence_energy, finite = TRUE))),
      show.legend = FALSE,
      check_overlap = TRUE,
      alpha = 1,
      colour = "#172b4d",
      size = 3.1,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = c(none = "white", review = "#ffb000", insufficient = "#d9e2ec"),
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c("no validity concern" = "#334e68", "validity concern" = "#d1495b"),
      drop = TRUE
    ) +
    ggplot2::scale_shape_manual(values = rep(21:25, length.out = nlevels(shown$display_group))) +
    ggplot2::scale_alpha_identity() +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21, colour = "#334e68", alpha = 1)),
      colour = ggplot2::guide_legend(override.aes = list(shape = 21, fill = "white", alpha = 1)),
      shape = if (is.null(attr(data, "group_by"))) {
        "none"
      } else {
        ggplot2::guide_legend(
          override.aes = list(fill = "white", colour = "#334e68", alpha = 1)
        )
      }
    ) +
    ggplot2::labs(
      title = "Group examination: action plane",
      subtitle = "Unexpectedness and influence are separate review signals",
      x = "Model surprise (capped residual energy)",
      y = "Group influence (capped deletion-statistic energy)",
      fill = "Review status",
      colour = "Data validity",
      shape = attr(data, "group_by") %||% NULL,
      caption = .examination_caption(data)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", colour = "#172b4d"),
      plot.subtitle = ggplot2::element_text(colour = "#52606d"),
      plot.caption = ggplot2::element_text(colour = "#52606d", hjust = 0),
      legend.position = "bottom"
    )
  if (length(unique(shown$contrast)) > 1L || length(unique(shown$estimand)) > 1L) {
    p <- p + ggplot2::facet_grid(
      rows = ggplot2::vars(estimand),
      cols = ggplot2::vars(contrast),
      scales = "free"
    )
  }
  p
}

.metric_percentile <- function(value) {
  value <- as.numeric(value)
  out <- rep(NA_real_, length(value))
  ok <- is.finite(value)
  if (!any(ok)) return(out)
  if (sum(ok) == 1L || diff(range(value[ok])) == 0) {
    out[ok] <- 0
  } else {
    out[ok] <- (rank(value[ok], ties.method = "average") - 1) / (sum(ok) - 1)
  }
  out
}

.examination_heatmap_data <- function(x, contrast = NULL, estimand = NULL) {
  contrasts <- .examination_select_values(contrast, x$contrast_data$contrast, "contrast")
  influence <- .examination_screening_estimands(x)
  estimands <- .examination_select_values(estimand, influence$estimand, "estimand")
  rows <- list()
  index <- 1L
  for (contrast_name in contrasts) {
    cd <- x$contrast_data[x$contrast_data$contrast == contrast_name, , drop = FALSE]
    rows[[index]] <- data.frame(
      subject = cd$subject,
      contrast = contrast_name,
      estimand = NA_character_,
      metric = "model surprise",
      metric_label = paste(contrast_name, "model surprise", sep = " / "),
      value = cd$surprise_energy,
      display = .metric_percentile(cd$surprise_energy),
      status = cd$surprise_status,
      mode = "predictive",
      stable = cd$surprise_stability,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
    gain_display <- tanh(cd$zero_intercept_gain - 1)
    rows[[index]] <- data.frame(
      subject = cd$subject,
      contrast = contrast_name,
      estimand = NA_character_,
      metric = "map gain",
      metric_label = paste(contrast_name, "map gain minus one", sep = " / "),
      value = cd$zero_intercept_gain,
      display = gain_display,
      status = cd$gain_status,
      mode = "predictive",
      stable = cd$surprise_stability,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
    for (estimand_name in estimands) {
      ed <- influence[
        influence$contrast == contrast_name & influence$estimand == estimand_name,
        , drop = FALSE
      ]
      rows[[index]] <- data.frame(
        subject = ed$subject,
        contrast = contrast_name,
        estimand = estimand_name,
        metric = "group influence",
        metric_label = paste(
          contrast_name, estimand_name, "group influence", sep = " / "
        ),
        value = ed$influence_energy,
        display = .metric_percentile(ed$influence_energy),
        status = ed$status,
        mode = ed$mode,
        stable = ed$stability,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  out <- do.call(rbind, rows)
  subject_order <- x$subject_data$subject[
    order(-x$subject_data$review_priority, x$subject_data$subject, na.last = TRUE)
  ]
  out$subject <- factor(out$subject, levels = subject_order)
  out$metric_label <- factor(out$metric_label, levels = rev(unique(out$metric_label)))
  out
}

.plot_examination_heatmap <- function(x, contrast = NULL, estimand = NULL) {
  data <- .examination_heatmap_data(x, contrast, estimand)
  available <- is.finite(data$display) & data$status == "available"
  if (!any(available)) {
    return(.examination_empty_plot(
      "Subject by diagnostic heatmap",
      "The selected summaries are unavailable."
    ))
  }
  data$display[!available] <- NA_real_
  ggplot2::ggplot(data, ggplot2::aes(x = subject, y = metric_label, fill = display)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, na.value = "#d9e2ec", limits = c(-1, 1)
    ) +
    ggplot2::labs(
      title = "Subject by contrast and estimand",
      subtitle = "Ranks for surprise/influence; signed deviation from unit gain",
      x = NULL,
      y = NULL,
      fill = "Display score",
      caption = paste0(
        "Unavailable cells: ", sum(!available),
        " | influence modes: ", paste(unique(stats::na.omit(data$mode[data$metric == "group influence"])), collapse = ", ")
      )
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", colour = "#172b4d"),
      plot.caption = ggplot2::element_text(colour = "#52606d", hjust = 0),
      legend.position = "bottom"
    )
}

.plot_examination_embedding <- function(x, group_by = NULL) {
  embedding <- x$embedding
  coordinates <- embedding$coordinates %||% data.frame()
  dimensions <- grep("^dimension", names(coordinates), value = TRUE)
  if (!identical(embedding$status, "available") || !length(dimensions)) {
    return(.examination_empty_plot(
      "Model-adjusted residual geometry",
      paste0("Residual geometry is ", embedding$status %||% "unavailable", ".")
    ))
  }
  data <- merge(coordinates, x$subject_data, by = "subject", all.x = TRUE, sort = FALSE)
  data$x <- data[[dimensions[1L]]]
  data$y <- if (length(dimensions) >= 2L) data[[dimensions[2L]]] else 0
  data$review_label <- ifelse(data$review_status == "review", data$subject, NA_character_)
  data$display_alpha <- ifelse(
    is.finite(data$stability),
    0.25 + 0.75 * pmax(0, pmin(1, data$stability)),
    0.35
  )
  if (!is.null(group_by)) {
    if (!is.character(group_by) || length(group_by) != 1L || !group_by %in% names(data)) {
      stop("group_by must name one subject_data column.", call. = FALSE)
    }
    group_value <- as.character(data[[group_by]])
    group_value[is.na(group_value)] <- "(missing)"
    data$display_group <- factor(group_value)
  } else {
    data$display_group <- factor("cohort")
  }
  explained <- embedding$explained_energy
  labels <- paste0(
    "Residual dimension ", seq_len(min(2L, length(explained))),
    " (", formatC(100 * explained[seq_len(min(2L, length(explained)))], digits = 1, format = "f"), "%)"
  )
  if (length(labels) == 1L) labels <- c(labels, "No second retained dimension")
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = x, y = y, shape = display_group, alpha = display_alpha)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#d9e2ec", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#d9e2ec", linewidth = 0.3) +
    ggplot2::geom_point(
      ggplot2::aes(fill = review_status), colour = "#334e68",
      size = 3.2, stroke = 0.9
    ) +
    ggplot2::geom_text(
      data = data[!is.na(data$review_label), , drop = FALSE],
      ggplot2::aes(label = review_label),
      nudge_y = 0.02 * max(1, diff(range(data$y, finite = TRUE))),
      show.legend = FALSE,
      check_overlap = TRUE,
      alpha = 1,
      colour = "#172b4d",
      size = 3.1
    ) +
    ggplot2::scale_fill_manual(
      values = c(none = "white", review = "#ffb000", insufficient = "#d9e2ec"),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(values = rep(21:25, length.out = nlevels(data$display_group))) +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(
      title = "Model-adjusted residual geometry",
      subtitle = "A navigational low-rank view, not a subject classification",
      x = labels[1L],
      y = labels[2L],
      fill = "Review status",
      shape = group_by %||% NULL,
      caption = paste0(
        "Captured residual energy: ",
        formatC(100 * embedding$captured_energy, digits = 1, format = "f"),
        "% | stability: ", embedding$stability_method,
        " | sketch rank: ", embedding$sketch_rank
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", colour = "#172b4d"),
      plot.caption = ggplot2::element_text(colour = "#52606d", hjust = 0),
      legend.position = "bottom"
    )
}

.examination_subject_assays <- function(x,
                                        estimand = NULL,
                                        requested_assays = NULL) {
  available <- names(x$subject_maps$assays)
  if (!is.null(requested_assays)) {
    unknown <- setdiff(requested_assays, available)
    if (length(unknown)) {
      stop("Unknown selected-subject assays: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    return(requested_assays)
  }
  screening <- .examination_screening_estimands(x)
  estimand <- .examination_select_values(estimand, screening$estimand, "estimand")[1L]
  first_available <- function(candidates) candidates[candidates %in% available][1L]
  selected <- c(
    first_available("observed"),
    first_available(c("expected_exact", "expected")),
    first_available(c("predictive_residual_exact", "predictive_residual")),
    first_available(c(
      paste0("delta_stat_exact:", estimand),
      paste0("delta_stat:", estimand)
    ))
  )
  unique(stats::na.omit(selected))
}

.examination_subject_map_data <- function(x,
                                          subject,
                                          contrast = NULL,
                                          estimand = NULL,
                                          requested_assays = NULL,
                                          slice = NULL) {
  if (is.null(x$subject_maps)) {
    stop("This examination retained no selected-subject maps.", call. = FALSE)
  }
  subject <- as.character(subject)
  if (length(subject) != 1L || !subject %in% subjects(x$subject_maps)) {
    stop("subject must identify one retained subject.", call. = FALSE)
  }
  map_contrasts <- x$subject_maps$contrasts
  selected_contrasts <- .examination_select_values(
    contrast, map_contrasts, "contrast"
  )
  assay_names <- .examination_subject_assays(x, estimand, requested_assays)
  subject_index <- match(subject, subjects(x$subject_maps))
  contrast_index <- match(selected_contrasts, map_contrasts)
  values <- x$subject_maps$assays[assay_names]
  sp <- space(x$subject_maps)
  if (inherits(sp, "space_voxel")) {
    if (is.null(slice)) slice <- as.integer(ceiling(sp$dim[3L] / 2))
    slice <- as.integer(slice)
    if (length(slice) != 1L || is.na(slice) || slice < 1L || slice > sp$dim[3L]) {
      stop("slice must be a valid one-based axial slice.", call. = FALSE)
    }
    rows <- list()
    index <- 1L
    sample_index <- sp$mask_idx %||% seq_len(prod(sp$dim))
    for (assay_name in assay_names) {
      for (k in seq_along(selected_contrasts)) {
        full <- rep(NA_real_, prod(sp$dim))
        full[sample_index] <- values[[assay_name]][, subject_index, contrast_index[k]]
        volume <- array(full, dim = sp$dim)
        grid <- expand.grid(x = seq_len(sp$dim[1L]), y = seq_len(sp$dim[2L]))
        grid$value <- as.vector(volume[, , slice])
        grid$assay <- assay_name
        grid$contrast <- selected_contrasts[k]
        grid$slice <- slice
        rows[[index]] <- grid
        index <- index + 1L
      }
    }
    out <- do.call(rbind, rows)
    out$assay <- factor(out$assay, levels = assay_names)
    out$contrast <- factor(out$contrast, levels = selected_contrasts)
    attr(out, "display") <- "voxel"
    return(out)
  }

  labels <- if (inherits(sp, c("space_parcels", "space_sample_labels"))) {
    sp$labels
  } else {
    as.character(seq_len(dim(values[[1L]])[1L]))
  }
  rows <- list()
  index <- 1L
  for (assay_name in assay_names) {
    for (k in seq_along(selected_contrasts)) {
      rows[[index]] <- data.frame(
        feature = seq_along(labels),
        feature_label = labels,
        value = values[[assay_name]][, subject_index, contrast_index[k]],
        assay = assay_name,
        contrast = selected_contrasts[k],
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  out <- do.call(rbind, rows)
  out$assay <- factor(out$assay, levels = assay_names)
  out$contrast <- factor(out$contrast, levels = selected_contrasts)
  attr(out, "display") <- "feature"
  out
}

.plot_examination_subject <- function(x,
                                      subject,
                                      contrast = NULL,
                                      estimand = NULL,
                                      assays = NULL,
                                      slice = NULL) {
  data <- .examination_subject_map_data(
    x, subject, contrast, estimand, assays, slice
  )
  reason <- .examination_review_reason(x, subject)
  modes <- unique(metadata(x$subject_maps)$examination$assay_modes[data$assay])
  caption <- paste(
    strwrap(
      paste0(reason, " | map modes: ", paste(unique(modes), collapse = ", ")),
      width = 135L
    ),
    collapse = "\n"
  )
  if (identical(attr(data, "display"), "voxel")) {
    limit <- max(abs(data$value), na.rm = TRUE)
    if (!is.finite(limit) || limit == 0) limit <- 1
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, fill = value)) +
        ggplot2::geom_raster(na.rm = TRUE) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(assay),
          cols = ggplot2::vars(contrast),
          switch = "y"
        ) +
        ggplot2::scale_fill_gradient2(
          low = "#2166ac", mid = "white", high = "#b2182b",
          midpoint = 0, limits = c(-limit, limit), na.value = "#eef2f6"
        ) +
        ggplot2::labs(
          title = paste("Selected subject", subject),
          subtitle = paste("Axial slice", unique(data$slice)),
          x = NULL, y = NULL, fill = "value", caption = caption
        ) +
        ggplot2::theme_void(base_size = 10) +
        ggplot2::theme(
          strip.text = ggplot2::element_text(colour = "#172b4d"),
          strip.placement = "outside",
          plot.title = ggplot2::element_text(face = "bold", colour = "#172b4d"),
          plot.caption = ggplot2::element_text(colour = "#52606d", hjust = 0),
          legend.position = "bottom"
        )
    )
  }
  ggplot2::ggplot(data, ggplot2::aes(x = feature, y = value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#d9e2ec", linewidth = 0.35) +
    ggplot2::geom_line(colour = "#2166ac", linewidth = 0.55, na.rm = TRUE) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(assay),
      cols = ggplot2::vars(contrast),
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::labs(
      title = paste("Selected subject", subject),
      subtitle = "Observed, leave-one-out expected, residual, and group-statistic deletion profiles",
      x = "Canonical sample index", y = NULL, caption = caption
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.placement = "outside",
      plot.title = ggplot2::element_text(face = "bold", colour = "#172b4d"),
      plot.caption = ggplot2::element_text(colour = "#52606d", hjust = 0)
    )
}
