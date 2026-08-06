# Self-contained group-examination report --------------------------------

.examination_format_number <- function(value, digits = 3L) {
  ifelse(
    is.finite(value),
    formatC(value, digits = digits, format = "fg", flag = "#"),
    "unavailable"
  )
}

.examination_review_reason <- function(x, subject) {
  .validate_gds_examination(x)
  subject <- as.character(subject)
  if (length(subject) != 1L || !subject %in% x$subject_data$subject) {
    stop("subject must identify one examined subject.", call. = FALSE)
  }
  row <- x$subject_data[x$subject_data$subject == subject, , drop = FALSE]
  if (identical(row$review_status, "insufficient")) {
    return("Insufficient eligible data; no review interpretation is assigned.")
  }
  if (identical(row$review_status, "none")) {
    source <- c(
      data_validity = row$quality_percentile,
      model_surprise = row$surprise_percentile,
      group_influence = row$influence_percentile
    )
    source[!is.finite(source)] <- -Inf
    component <- if (all(source == -Inf)) "available diagnostic" else
      names(source)[which.max(source)]
    return(paste0(
      "No absolute review criterion met. ",
      gsub("_", " ", component, fixed = TRUE),
      " provided this subject's highest relative ordering, but the absolute discrepancy was modest."
    ))
  }
  if (identical(row$review_source, "quality")) {
    return(paste0(
      "Priority review from the data-validity component. ",
      row$review_reason,
      " Model surprise and group influence remain separate diagnostics."
    ))
  }
  if (identical(row$review_source, "surprise")) {
    candidates <- x$contrast_data[x$contrast_data$subject == subject, , drop = FALSE]
    candidates <- candidates[
      order(-candidates$surprise_energy, candidates$contrast, na.last = TRUE),
      , drop = FALSE
    ]
    chosen <- candidates[1L, , drop = FALSE]
    screening <- .examination_screening_estimands(x)
    screening <- screening[
      screening$subject == subject & screening$contrast == chosen$contrast,
      , drop = FALSE
    ]
    mode <- if (nrow(screening)) screening$mode[1L] else "unavailable"
    gain <- if (identical(chosen$gain_status, "available")) {
      paste0("; signed map gain ", .examination_format_number(chosen$zero_intercept_gain))
    } else {
      paste0("; map gain ", chosen$gain_status)
    }
    return(paste0(
      "Priority review from model surprise in contrast ", chosen$contrast,
      " (capped residual energy ", .examination_format_number(chosen$surprise_energy),
      ", tail extent ", .examination_format_number(chosen$tail_extent), gain,
      "; status ", chosen$surprise_status, "; mode ", mode,
      "). Group influence was not the component that set review status."
    ))
  }
  candidates <- .examination_screening_estimands(x)
  candidates <- candidates[candidates$subject == subject, , drop = FALSE]
  candidates <- candidates[
    order(
      -candidates$influence_energy,
      candidates$contrast,
      candidates$estimand,
      na.last = TRUE
    ),
    , drop = FALSE
  ]
  chosen <- candidates[1L, , drop = FALSE]
  paste0(
    "Priority review from group influence for estimand ", chosen$estimand,
    " in contrast ", chosen$contrast,
    " (capped deletion-statistic energy ",
    .examination_format_number(chosen$influence_energy),
    ", maximum absolute change ",
    .examination_format_number(chosen$max_abs_delta_stat),
    "; status ", chosen$status, "; mode ", chosen$mode,
    "). Model surprise was not the component that set review status."
  )
}

.examination_ranked_regions <- function(x, subject, contrast = NULL, n = 10L) {
  if (is.null(x$subject_maps)) return(data.frame())
  subject_index <- match(subject, x$subject_maps$subjects)
  if (is.na(subject_index)) return(data.frame())
  map_contrasts <- x$subject_maps$contrasts
  contrast <- contrast %||% map_contrasts[1L]
  contrast_index <- match(contrast, map_contrasts)
  if (is.na(contrast_index)) return(data.frame())
  available <- names(x$subject_maps$assays)
  residual_name <- c("predictive_residual_exact", "predictive_residual")
  residual_name <- residual_name[residual_name %in% available][1L]
  estimand <- unique(.examination_screening_estimands(x)$estimand)[1L]
  deletion_name <- c(
    paste0("delta_stat_exact:", estimand),
    paste0("delta_stat:", estimand)
  )
  deletion_name <- deletion_name[deletion_name %in% available][1L]
  if (is.na(residual_name) && is.na(deletion_name)) return(data.frame())

  n_sample <- dim(x$subject_maps$assays[[1L]])[1L]
  rows <- x$subject_maps$row_data
  candidates <- c("region", "parcel", "parcel_id", "atlas", "label", "name")
  label_column <- if (is.null(rows)) character() else intersect(candidates, names(rows))[1L]
  kind <- "region"
  if (length(label_column) && !is.na(label_column)) {
    labels <- as.character(rows[[label_column]])
  } else if (inherits(space(x$subject_maps), "space_parcels")) {
    labels <- space(x$subject_maps)$labels
  } else if (inherits(space(x$subject_maps), "space_sample_labels")) {
    labels <- space(x$subject_maps)$labels
    kind <- "feature"
  } else {
    return(data.frame())
  }
  if (length(labels) != n_sample) return(data.frame())
  residual <- if (!is.na(residual_name)) {
    x$subject_maps$assays[[residual_name]][, subject_index, contrast_index]
  } else {
    rep(NA_real_, n_sample)
  }
  deletion <- if (!is.na(deletion_name)) {
    x$subject_maps$assays[[deletion_name]][, subject_index, contrast_index]
  } else {
    rep(NA_real_, n_sample)
  }
  groups <- unique(labels[!is.na(labels) & nzchar(labels)])
  output <- lapply(groups, function(label) {
    use <- labels == label
    data.frame(
      label = label,
      mean_abs_predictive_residual = if (any(is.finite(residual[use]))) {
        mean(abs(residual[use]), na.rm = TRUE)
      } else {
        NA_real_
      },
      max_abs_delta_stat = if (any(is.finite(deletion[use]))) {
        max(abs(deletion[use]), na.rm = TRUE)
      } else {
        NA_real_
      },
      eligible_n = sum(is.finite(residual[use]) | is.finite(deletion[use])),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, output)
  score <- pmax(
    ifelse(is.finite(output$mean_abs_predictive_residual), output$mean_abs_predictive_residual, -Inf),
    ifelse(is.finite(output$max_abs_delta_stat), output$max_abs_delta_stat, -Inf)
  )
  output <- output[order(-score, output$label), , drop = FALSE]
  output <- head(output, as.integer(n))
  names(output)[1L] <- kind
  rownames(output) <- NULL
  output
}

.html_escape <- function(value) {
  value <- as.character(value)
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub("\"", "&quot;", value, fixed = TRUE)
  gsub("'", "&#39;", value, fixed = TRUE)
}

.html_id <- function(value) {
  stem <- gsub("[^A-Za-z0-9_-]+", "-", as.character(value))
  paste0(stem, "-", substr(digest::digest(value, algo = "xxhash64"), 1L, 8L))
}

.html_table <- function(data, class = "data-table", empty = "None") {
  if (is.null(data) || !nrow(data)) {
    return(paste0("<p class=\"muted\">", .html_escape(empty), "</p>"))
  }
  format_cell <- function(value) {
    if (length(value) != 1L || is.na(value)) return("unavailable")
    if (is.numeric(value)) return(.examination_format_number(value))
    as.character(value)
  }
  head <- paste0(
    "<thead><tr>",
    paste0("<th>", .html_escape(names(data)), "</th>", collapse = ""),
    "</tr></thead>"
  )
  body <- vapply(seq_len(nrow(data)), function(i) {
    cells <- vapply(data[i, , drop = FALSE], function(value) {
      paste0("<td>", .html_escape(format_cell(value)), "</td>")
    }, character(1))
    paste0("<tr>", paste0(cells, collapse = ""), "</tr>")
  }, character(1))
  paste0(
    "<div class=\"table-wrap\"><table class=\"", class, "\">",
    head, "<tbody>", paste0(body, collapse = ""), "</tbody></table></div>"
  )
}

.svg_scale <- function(value, domain, range) {
  if (!length(value)) return(numeric())
  if (!all(is.finite(domain)) || diff(domain) <= 0) return(rep(mean(range), length(value)))
  range[1L] + (value - domain[1L]) / diff(domain) * diff(range)
}

.html_action_plane <- function(x, width = 760L, height = 430L) {
  examination <- x
  data <- examination$subject_data
  ok <- is.finite(data$surprise_score) & is.finite(data$influence_score)
  if (!any(ok)) {
    return("<div class=\"empty\">Action-plane metrics are unavailable.</div>")
  }
  data <- data[ok, , drop = FALSE]
  margin <- c(left = 70, right = 24, top = 30, bottom = 58)
  x_range <- range(c(0, data$surprise_score), finite = TRUE)
  y_range <- range(c(0, data$influence_score), finite = TRUE)
  point_x <- .svg_scale(
    data$surprise_score, x_range, c(margin["left"], width - margin["right"])
  )
  point_y <- .svg_scale(
    data$influence_score, y_range, c(height - margin["bottom"], margin["top"])
  )
  x_threshold <- examination$config$control$review$surprise$energy_threshold
  y_threshold <- examination$config$control$review$influence$energy_threshold
  threshold_x <- .svg_scale(x_threshold, x_range, c(margin["left"], width - margin["right"]))
  threshold_y <- .svg_scale(y_threshold, y_range, c(height - margin["bottom"], margin["top"]))
  threshold_x <- max(margin["left"], min(width - margin["right"], threshold_x))
  threshold_y <- max(margin["top"], min(height - margin["bottom"], threshold_y))
  retained <- examination$config$retained_subjects %||% character()
  points <- vapply(seq_len(nrow(data)), function(i) {
    fill <- if (data$review_status[i] == "review") "#ffb000" else
      if (data$review_status[i] == "insufficient") "#d9e2ec" else "#ffffff"
    stroke <- if (identical(data$review_source[i], "quality")) "#d1495b" else "#334e68"
    anchor <- .html_id(data$subject[i])
    title <- paste0(
      data$subject[i], ": surprise ", .examination_format_number(data$surprise_score[i]),
      ", influence ", .examination_format_number(data$influence_score[i]),
      ", ", data$review_status[i]
    )
    label <- if (data$review_status[i] == "review") {
      label_x <- if (point_x[i] > width - 100) point_x[i] - 7 else point_x[i] + 7
      label_anchor <- if (point_x[i] > width - 100) "end" else "start"
      label_y <- if (point_y[i] < margin["top"] + 14) point_y[i] + 15 else point_y[i] - 7
      paste0(
        "<text x=\"", label_x, "\" y=\"", label_y,
        "\" text-anchor=\"", label_anchor,
        "\" class=\"point-label\">", .html_escape(data$subject[i]), "</text>"
      )
    } else {
      ""
    }
    circle <- paste0(
      "<circle cx=\"", point_x[i], "\" cy=\"", point_y[i],
      "\" r=\"5.5\" fill=\"", fill, "\" stroke=\"", stroke,
      "\" stroke-width=\"2\"><title>", .html_escape(title),
      "</title></circle>", label
    )
    if (data$subject[i] %in% retained) {
      paste0("<a href=\"#subject-", anchor, "\">", circle, "</a>")
    } else {
      circle
    }
  }, character(1))
  paste0(
    "<svg class=\"board-svg\" viewBox=\"0 0 ", width, " ", height,
    "\" role=\"img\" aria-label=\"Action plane\">",
    "<rect x=\"", margin["left"], "\" y=\"", margin["top"], "\" width=\"",
    width - margin["left"] - margin["right"], "\" height=\"",
    height - margin["top"] - margin["bottom"], "\" fill=\"#fbfcfe\"/>",
    "<line class=\"threshold\" x1=\"", threshold_x, "\" x2=\"", threshold_x,
    "\" y1=\"", margin["top"], "\" y2=\"", height - margin["bottom"], "\"/>",
    "<line class=\"threshold\" y1=\"", threshold_y, "\" y2=\"", threshold_y,
    "\" x1=\"", margin["left"], "\" x2=\"", width - margin["right"], "\"/>",
    "<line class=\"axis\" x1=\"", margin["left"], "\" x2=\"", width - margin["right"],
    "\" y1=\"", height - margin["bottom"], "\" y2=\"", height - margin["bottom"], "\"/>",
    "<line class=\"axis\" x1=\"", margin["left"], "\" x2=\"", margin["left"],
    "\" y1=\"", margin["top"], "\" y2=\"", height - margin["bottom"], "\"/>",
    paste0(points, collapse = ""),
    "<text class=\"axis-label\" x=\"", width / 2, "\" y=\"", height - 16,
    "\">model surprise</text>",
    "<text class=\"axis-label\" transform=\"translate(18 ", height / 2,
    ") rotate(-90)\">group influence</text>",
    "</svg>"
  )
}

.html_embedding <- function(x, width = 540L, height = 330L) {
  embedding <- x$embedding
  coordinates <- embedding$coordinates %||% data.frame()
  dimensions <- grep("^dimension", names(coordinates), value = TRUE)
  if (!identical(embedding$status, "available") || !length(dimensions)) {
    return(paste0(
      "<div class=\"empty\">Residual geometry is ",
      .html_escape(embedding$status %||% "unavailable"), ".</div>"
    ))
  }
  data <- merge(coordinates, x$subject_data, by = "subject", all.x = TRUE, sort = FALSE)
  x_value <- data[[dimensions[1L]]]
  y_value <- if (length(dimensions) > 1L) data[[dimensions[2L]]] else rep(0, nrow(data))
  px <- .svg_scale(x_value, range(x_value, finite = TRUE), c(38, width - 18))
  py <- .svg_scale(y_value, range(y_value, finite = TRUE), c(height - 38, 18))
  retained <- x$config$retained_subjects %||% character()
  points <- vapply(seq_len(nrow(data)), function(i) {
    anchor <- .html_id(data$subject[i])
    fill <- if (data$review_status[i] == "review") "#ffb000" else "#ffffff"
    opacity <- if (is.finite(data$stability[i])) 0.25 + 0.75 * data$stability[i] else 0.35
    circle <- paste0(
      "<circle cx=\"", px[i], "\" cy=\"", py[i],
      "\" r=\"5\" fill=\"", fill, "\" stroke=\"#334e68\" opacity=\"", opacity,
      "\"><title>", .html_escape(data$subject[i]), "</title></circle>"
    )
    if (data$subject[i] %in% retained) {
      paste0("<a href=\"#subject-", anchor, "\">", circle, "</a>")
    } else {
      circle
    }
  }, character(1))
  paste0(
    "<svg class=\"board-svg\" viewBox=\"0 0 ", width, " ", height,
    "\" role=\"img\" aria-label=\"Residual geometry\">",
    "<line class=\"axis\" x1=\"38\" x2=\"", width - 18,
    "\" y1=\"", height - 38, "\" y2=\"", height - 38, "\"/>",
    "<line class=\"axis\" x1=\"38\" x2=\"38\" y1=\"18\" y2=\"",
    height - 38, "\"/>", paste0(points, collapse = ""),
    "<text class=\"small-label\" x=\"", width / 2, "\" y=\"", height - 8,
    "\">dimension 1</text></svg>",
    "<p class=\"caption\">Captured residual energy: ",
    formatC(100 * embedding$captured_energy, digits = 1, format = "f"),
    "%; stability: ", .html_escape(embedding$stability_method), ".</p>"
  )
}

.html_colour <- function(value) {
  if (!is.finite(value)) return("#d9e2ec")
  value <- max(-1, min(1, value))
  if (value >= 0) {
    red <- c(247, 247, 247) + value * (c(178, 24, 43) - c(247, 247, 247))
  } else {
    red <- c(247, 247, 247) + (-value) * (c(33, 102, 172) - c(247, 247, 247))
  }
  sprintf("#%02x%02x%02x", as.integer(red[1]), as.integer(red[2]), as.integer(red[3]))
}

.html_heatmap <- function(x) {
  data <- .examination_heatmap_data(x)
  subjects_order <- levels(data$subject)
  metrics <- levels(data$metric_label)
  header <- paste0(
    "<tr><th>diagnostic</th>",
    paste0("<th class=\"vertical\">", .html_escape(subjects_order), "</th>", collapse = ""),
    "</tr>"
  )
  body <- vapply(metrics, function(metric) {
    cells <- vapply(subjects_order, function(subject) {
      row <- data[as.character(data$subject) == subject & as.character(data$metric_label) == metric, , drop = FALSE]
      if (!nrow(row)) return("<td class=\"heat unavailable\"></td>")
      title <- paste0(
        subject, ": ", row$metric, " = ", .examination_format_number(row$value),
        "; status ", row$status, "; mode ", row$mode
      )
      paste0(
        "<td class=\"heat\" style=\"background:", .html_colour(row$display),
        "\"><title>", .html_escape(title), "</title></td>"
      )
    }, character(1))
    paste0("<tr><th>", .html_escape(metric), "</th>", paste0(cells, collapse = ""), "</tr>")
  }, character(1))
  paste0(
    "<div class=\"table-wrap\"><table class=\"heatmap\"><thead>", header,
    "</thead><tbody>", paste0(body, collapse = ""), "</tbody></table></div>"
  )
}

.html_polyline_panel <- function(data, width = 330L, height = 120L) {
  finite <- is.finite(data$value)
  if (!any(finite)) return("<div class=\"empty compact\">Unavailable</div>")
  x <- .svg_scale(data$feature, range(data$feature, finite = TRUE), c(8, width - 8))
  y <- .svg_scale(data$value, range(data$value[finite], finite = TRUE), c(height - 12, 12))
  points <- paste0(x[finite], ",", y[finite], collapse = " ")
  zero <- if (range(data$value[finite])[1L] <= 0 && range(data$value[finite])[2L] >= 0) {
    .svg_scale(0, range(data$value[finite], finite = TRUE), c(height - 12, 12))
  } else {
    NA_real_
  }
  paste0(
    "<svg class=\"map-svg\" viewBox=\"0 0 ", width, " ", height, "\">",
    if (is.finite(zero)) paste0(
      "<line x1=\"8\" x2=\"", width - 8, "\" y1=\"", zero,
      "\" y2=\"", zero, "\" stroke=\"#d9e2ec\"/>"
    ) else "",
    "<polyline points=\"", points,
    "\" fill=\"none\" stroke=\"#2166ac\" stroke-width=\"1.4\"/></svg>"
  )
}

.html_voxel_panel <- function(data, width = 230L, height = 200L) {
  x_values <- sort(unique(data$x))
  y_values <- sort(unique(data$y))
  cell_w <- width / length(x_values)
  cell_h <- height / length(y_values)
  limit <- max(abs(data$value), na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) limit <- 1
  cells <- vapply(seq_len(nrow(data)), function(i) {
    value <- data$value[i] / limit
    paste0(
      "<rect x=\"", (match(data$x[i], x_values) - 1) * cell_w,
      "\" y=\"", (length(y_values) - match(data$y[i], y_values)) * cell_h,
      "\" width=\"", cell_w + 0.2, "\" height=\"", cell_h + 0.2,
      "\" fill=\"", .html_colour(value), "\"><title>",
      .html_escape(.examination_format_number(data$value[i])), "</title></rect>"
    )
  }, character(1))
  paste0(
    "<svg class=\"map-svg\" viewBox=\"0 0 ", width, " ", height, "\">",
    paste0(cells, collapse = ""), "</svg>"
  )
}

.html_subject_maps <- function(x, subject) {
  data <- .examination_subject_map_data(x, subject)
  keys <- unique(data[, c("assay", "contrast"), drop = FALSE])
  cards <- vapply(seq_len(nrow(keys)), function(i) {
    selected <- data[data$assay == keys$assay[i] & data$contrast == keys$contrast[i], , drop = FALSE]
    plot <- if (identical(attr(data, "display"), "voxel")) {
      .html_voxel_panel(selected)
    } else {
      .html_polyline_panel(selected)
    }
    paste0(
      "<div class=\"map-card\"><h4>", .html_escape(keys$assay[i]),
      "</h4><p class=\"eyebrow\">", .html_escape(keys$contrast[i]),
      "</p>", plot, "</div>"
    )
  }, character(1))
  paste0("<div class=\"map-grid\">", paste0(cards, collapse = ""), "</div>")
}

.html_subject_section <- function(x, subject) {
  id <- .html_id(subject)
  subject_row <- x$subject_data[x$subject_data$subject == subject, , drop = FALSE]
  core <- c(
    "subject", "coverage_fraction", "surprise_score", "influence_score",
    "review_priority", "review_status", "review_source", "review_reason", "retained",
    "quality_score", "quality_percentile", "surprise_percentile", "influence_percentile"
  )
  qc_columns <- setdiff(names(subject_row), core)
  qc <- if (length(qc_columns)) {
    data.frame(metric = qc_columns, value = unlist(subject_row[qc_columns], use.names = FALSE), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  contrast_profile <- x$contrast_data[
    x$contrast_data$subject == subject,
    c(
      "contrast", "coverage_fraction", "surprise_energy", "tail_extent",
      "surprise_status", "surprise_stability", "weighted_correlation",
      "correlation_status", "zero_intercept_gain", "gain_status"
    ),
    drop = FALSE
  ]
  influence_profile <- x$estimand_data[
    x$estimand_data$subject == subject,
    c(
      "contrast", "estimand", "mode", "ranking_stage", "influence_energy",
      "max_abs_delta_stat", "eligible_n", "status", "stability"
    ),
    drop = FALSE
  ]
  regions <- .examination_ranked_regions(x, subject)
  conclusion <- x$conclusion$results
  if (!is.null(conclusion)) {
    conclusion <- conclusion[conclusion$subject == subject, , drop = FALSE]
  }
  paste0(
    "<section class=\"subject-section\" id=\"subject-", id, "\">",
    "<div class=\"section-heading\"><div><p class=\"eyebrow\">selected subject</p><h2>",
    .html_escape(subject), "</h2></div><a class=\"back\" href=\"#top\">back to board</a></div>",
    "<p class=\"reason\">", .html_escape(.examination_review_reason(x, subject)), "</p>",
    .html_subject_maps(x, subject),
    "<div class=\"two-col compact-grid\"><div><h3>Contrast profile</h3>",
    .html_table(contrast_profile), "</div><div><h3>Estimand sensitivity</h3>",
    .html_table(influence_profile), "</div></div>",
    "<div class=\"two-col compact-grid\"><div><h3>First-level and model covariates</h3>",
    .html_table(qc, empty = "No additional covariates were supplied."),
    "</div><div><h3>Ranked localization</h3>",
    .html_table(regions, empty = "No parcel or sample labels are available for localization."),
    "</div></div><div class=\"conclusion-box\"><h3>Threshold-dependent conclusion sensitivity</h3>",
    "<p class=\"muted\">This selected recomputation is separate from the continuous influence diagnostics above.</p>",
    .html_table(
      conclusion,
      empty = "No inherited post-hoc conclusion was requested or recomputation is unavailable."
    ),
    "</div></section>"
  )
}

.examination_report_css <- function() {
  paste0(
    ":root{--ink:#172b4d;--muted:#52606d;--line:#d9e2ec;--paper:#f4f7fb;",
    "--panel:#fff;--accent:#2166ac;--warn:#d1495b;--review:#ffb000}",
    "*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);",
    "font-family:ui-sans-serif,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;line-height:1.45}",
    "main{max-width:1320px;margin:auto;padding:28px}h1,h2,h3,h4{margin:.2rem 0 .7rem;line-height:1.15}",
    "h1{font-size:2rem}h2{font-size:1.4rem}h3{font-size:1rem}h4{font-size:.88rem}",
    ".eyebrow{text-transform:uppercase;letter-spacing:.09em;color:var(--muted);font-size:.72rem;font-weight:700;margin:0}",
    ".lede,.muted,.caption{color:var(--muted)}.caption{font-size:.78rem}.board,.subject-section{",
    "background:var(--panel);border:1px solid var(--line);border-radius:14px;box-shadow:0 8px 28px rgba(23,43,77,.06)}",
    ".board{padding:20px}.board-grid{display:grid;grid-template-columns:minmax(0,1.55fr) minmax(280px,.8fr);gap:18px}",
    ".panel{border:1px solid var(--line);border-radius:10px;padding:14px;overflow:hidden}.panel.wide{grid-column:1/-1}",
    ".summary-strip{display:flex;flex-wrap:wrap;gap:8px;margin:16px 0}.pill{background:#eaf0f7;border-radius:999px;padding:6px 10px;font-size:.8rem}",
    ".board-svg,.map-svg{display:block;width:100%;height:auto}.axis{stroke:#829ab1;stroke-width:1}.threshold{stroke:#9fb3c8;stroke-dasharray:5 4}",
    ".axis-label,.small-label,.point-label{fill:var(--muted);font-size:12px}.point-label{font-size:10px;font-weight:700}",
    ".table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.78rem}th,td{padding:6px 8px;border-bottom:1px solid #edf1f5;text-align:left;white-space:nowrap}",
    "th{color:var(--muted);font-weight:650;background:#fbfcfe;position:sticky;top:0}.heatmap th.vertical{writing-mode:vertical-rl;transform:rotate(180deg);height:70px;padding:4px}",
    ".heat{min-width:18px;height:18px;padding:0}.unavailable{background:#d9e2ec}.review-list{list-style:none;padding:0;margin:0}.review-list li{border-bottom:1px solid var(--line);padding:10px 0}",
    ".review-list a{color:var(--ink);font-weight:700;text-decoration:none}.reason{padding:11px 13px;background:#fff8df;border-left:4px solid var(--review);border-radius:5px}",
    ".subject-section{margin-top:24px;padding:22px}.section-heading{display:flex;align-items:start;justify-content:space-between}.back{color:var(--accent);text-decoration:none;font-size:.8rem}",
    ".map-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:12px}.map-card{border:1px solid var(--line);border-radius:8px;padding:10px;min-width:0}",
    ".two-col{display:grid;grid-template-columns:1fr 1fr;gap:16px}.compact-grid{margin-top:16px}.empty{display:grid;place-items:center;min-height:160px;color:var(--muted);background:#fbfcfe;border-radius:8px}.empty.compact{min-height:90px}",
    ".conclusion-box{margin-top:16px;padding:14px;border:1px solid #b8cbe0;border-radius:8px;background:#f8fbff}",
    "details{margin-top:18px}summary{cursor:pointer;font-weight:700}.footer{color:var(--muted);font-size:.75rem;margin:24px 0}",
    "@media(max-width:850px){main{padding:12px}.board-grid,.two-col{grid-template-columns:1fr}.panel.wide{grid-column:auto}}",
    "@media print{body{background:white}.board,.subject-section{box-shadow:none;break-inside:avoid}a{color:inherit}}"
  )
}

#' Write a self-contained group-examination report
#'
#' The report uses inline CSS and SVG only. Its action plane and residual
#' geometry link to retained subject drill-down sections; it has no external
#' JavaScript, font, image, or stylesheet dependency.
#'
#' @param x A `gds_examination`.
#' @param file Output HTML path.
#' @param title Report title.
#' @param subjects Optional retained subject IDs to include. By default, all
#'   subjects with retained maps are included.
#' @param overwrite Whether an existing file may be replaced.
#' @param ... Reserved for future report controls.
#'
#' @return The normalized output path, invisibly.
#' @export
write_report <- function(x,
                         file,
                         title = "Group examination",
                         subjects = NULL,
                         overwrite = FALSE,
                         ...) {
  .validate_gds_examination(x)
  if (!is.character(title) || length(title) != 1L || is.na(title) || !nzchar(title)) {
    stop("title must be one non-empty string.", call. = FALSE)
  }
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be one non-empty path.", call. = FALSE)
  }
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Report already exists; set overwrite = TRUE to replace it.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  }
  retained <- if (is.null(x$subject_maps)) character() else x$subject_maps$subjects
  if (is.null(subjects)) {
    selected <- retained
  } else {
    selected <- unique(as.character(subjects))
    if (anyNA(selected) || any(!nzchar(selected))) {
      stop("subjects must contain non-empty retained subject IDs.", call. = FALSE)
    }
    unknown <- setdiff(selected, retained)
    if (length(unknown)) {
      stop(
        "Reports can include only subjects with retained maps: ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  queue <- x$subject_data[x$subject_data$review_status == "review", , drop = FALSE]
  queue <- queue[order(-queue$review_priority, queue$subject), , drop = FALSE]
  queue_html <- if (!nrow(queue)) {
    "<p class=\"muted\">No subject met the absolute review criteria. The ordering remains available for inspection.</p>"
  } else {
    items <- vapply(seq_len(nrow(queue)), function(i) {
      anchor <- .html_id(queue$subject[i])
      linked <- if (queue$subject[i] %in% selected) {
        paste0("<a href=\"#subject-", anchor, "\">", .html_escape(queue$subject[i]), "</a>")
      } else {
        paste0("<strong>", .html_escape(queue$subject[i]), "</strong>")
      }
      paste0(
        "<li>", linked, "<br><span class=\"muted\">",
        .html_escape(.examination_review_reason(x, queue$subject[i])),
        "</span></li>"
      )
    }, character(1))
    paste0("<ul class=\"review-list\">", paste0(items, collapse = ""), "</ul>")
  }
  subject_sections <- if (length(selected)) {
    paste0(vapply(selected, function(subject) {
      .html_subject_section(x, subject)
    }, character(1)), collapse = "")
  } else {
    ""
  }
  availability <- x$availability
  unavailable <- availability[availability$status != "available", , drop = FALSE]
  conclusion_summary <- x$conclusion$results
  html <- paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>", .html_escape(title), "</title><style>", .examination_report_css(),
    "</style></head><body><main id=\"top\"><header><p class=\"eyebrow\">fmrigds</p><h1>",
    .html_escape(title), "</h1><p class=\"lede\">How expected is each subject under the intended group model, what explains any discrepancy, and how much does that subject change the group conclusion?</p>",
    "<div class=\"summary-strip\"><span class=\"pill\">N = ", x$cohort$n_included,
    "</span><span class=\"pill\">model: ", .html_escape(x$cohort$model),
    "</span><span class=\"pill\">formula: ", .html_escape(x$cohort$formula %||% "none"),
    "</span><span class=\"pill\">review cases: ", x$cohort$review_n,
    "</span><span class=\"pill\">variance: ", .html_escape(x$cohort$variance_mode),
    "</span></div></header><section class=\"board\"><div class=\"board-grid\">",
    "<div class=\"panel\"><h2>Action plane</h2><p class=\"muted\">Horizontal position is model surprise; vertical position is influence on the group statistic. Threshold lines are review gates, not exclusion rules.</p>",
    .html_action_plane(x), "</div><div class=\"panel\"><h2>Review queue</h2>",
    queue_html, "</div><div class=\"panel\"><h2>Residual geometry</h2>",
    .html_embedding(x), "</div><div class=\"panel\"><h2>Availability</h2>",
    .html_table(unavailable, empty = "All requested diagnostics are available."),
    "</div><div class=\"panel wide\"><h2>Subject x contrast x estimand</h2>",
    .html_heatmap(x), "</div><div class=\"panel wide\"><h2>Threshold-dependent conclusion sensitivity</h2>",
    "<p class=\"muted\">Full and selected-deletion results are stored separately from continuous influence. Only methods with a declared case-deletion contract are recomputed.</p>",
    .html_table(
      conclusion_summary,
      empty = "No inherited post-hoc conclusion was requested."
    ),
    "</div></div><details><summary>Execution and provenance</summary>",
    .html_table(data.frame(
      field = c("source plan digest", "examination digest", "adapter reads", "bytes read", "elapsed seconds", "peak RSS bytes", "retained map assays"),
      value = c(
        x$provenance$source_plan_digest,
        x$provenance$examination_digest,
        x$provenance$scan_receipt$adapter_reads,
        x$provenance$scan_receipt$bytes_read,
        x$provenance$scan_receipt$elapsed_seconds,
        x$provenance$scan_receipt$peak_rss_bytes,
        x$provenance$scan_receipt$retained_map_count
      ),
      stringsAsFactors = FALSE
    )), "</details></section>", subject_sections,
    "<p class=\"footer\">This examination supports inspection and sensitivity analysis. It does not classify subjects.</p>",
    "</main></body></html>"
  )
  directory <- dirname(file)
  if (!dir.exists(directory)) {
    stop("Report directory does not exist: ", directory, call. = FALSE)
  }
  writeLines(html, con = file, useBytes = TRUE)
  invisible(normalizePath(file, mustWork = TRUE))
}
