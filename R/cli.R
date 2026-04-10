# nocov start
.is_flag <- function(x) {
  is.character(x) && length(x) == 1L && nzchar(x) && grepl("^-", x)
}

.split_csv <- function(x) {
  if (is.null(x)) return(character())
  x <- as.character(x)
  x <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

.parse_kv <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "=", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) {
    stop("Expected KEY=VALUE but got: ", x, call. = FALSE)
  }
  key <- trimws(parts[1L])
  val <- trimws(paste(parts[-1L], collapse = "="))
  if (!nzchar(key)) stop("Empty key in KEY=VALUE: ", x, call. = FALSE)
  list(key = key, value = val)
}

.cli_parse_args <- function(args) {
  opts <- list()
  pos <- character()
  i <- 1L
  while (i <= length(args)) {
    tok <- args[[i]]
    if (identical(tok, "--")) {
      if (i < length(args)) pos <- c(pos, args[(i + 1L):length(args)])
      break
    }
    if (startsWith(tok, "--")) {
      kv <- substring(tok, 3L)
      if (grepl("=", kv, fixed = TRUE)) {
        parts <- strsplit(kv, "=", fixed = TRUE)[[1L]]
        key <- parts[[1L]]
        val <- paste(parts[-1L], collapse = "=")
        opts[[key]] <- c(opts[[key]], val)
        i <- i + 1L
        next
      }
      key <- kv
      if (i == length(args) || .is_flag(args[[i + 1L]])) {
        opts[[key]] <- c(opts[[key]], TRUE)
        i <- i + 1L
        next
      }
      opts[[key]] <- c(opts[[key]], args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (identical(tok, "-h")) {
      opts$help <- c(opts$help, TRUE)
      i <- i + 1L
      next
    }
    if (identical(tok, "-v")) {
      opts$version <- c(opts$version, TRUE)
      i <- i + 1L
      next
    }
    pos <- c(pos, tok)
    i <- i + 1L
  }
  list(opts = opts, pos = pos)
}

.cli_print <- function(...) {
  cat(..., sep = "")
}

.cli_warn <- function(...) {
  cat(..., "\n", sep = "", file = stderr())
}

.cli_quit <- function(status = 0L) {
  quit(save = "no", status = as.integer(status), runLast = FALSE)
}

.cli_version <- function() {
  .cli_print("fmrigds ", as.character(utils::packageVersion("fmrigds")), "\n")
}

.cli_get_opt <- function(opts, name, default = NULL) {
  val <- opts[[name]]
  if (is.null(val) || !length(val)) default else val[[1L]]
}

.cli_has_flag <- function(opts, name) {
  isTRUE(.cli_get_opt(opts, name, FALSE))
}

.cli_cast_scalar <- function(x) {
  x <- trimws(as.character(x)[1L])
  if (!nzchar(x)) return("")
  if (identical(x, "TRUE")) return(TRUE)
  if (identical(x, "FALSE")) return(FALSE)
  if (identical(x, "NA")) return(NA)
  if (grepl("^[+-]?[0-9]+$", x)) return(as.integer(x))
  if (grepl("^[+-]?([0-9]*\\.[0-9]+|[0-9]+\\.?)([eE][+-]?[0-9]+)?$", x)) {
    return(as.numeric(x))
  }
  x
}

.cli_cast_value <- function(x) {
  parts <- .split_csv(x)
  if (!length(parts)) return(character())
  parsed <- lapply(parts, .cli_cast_scalar)
  if (length(parsed) == 1L) return(parsed[[1L]])

  if (all(vapply(parsed, is.logical, logical(1)))) {
    return(as.logical(unlist(parsed, use.names = FALSE)))
  }
  if (all(vapply(parsed, function(val) is.integer(val) || is.numeric(val), logical(1)))) {
    return(as.numeric(unlist(parsed, use.names = FALSE)))
  }
  if (all(vapply(parsed, is.character, logical(1)))) {
    return(as.character(unlist(parsed, use.names = FALSE)))
  }

  as.character(unlist(lapply(parsed, as.character), use.names = FALSE))
}

.cli_parse_options_kv <- function(x) {
  if (is.null(x) || !length(x)) return(list())
  out <- list()
  for (tok in x) {
    kv <- .parse_kv(tok)
    out[[kv$key]] <- .cli_cast_value(kv$value)
  }
  out
}

.cli_parse_subset_spec <- function(x) {
  kv <- .parse_kv(x)
  axis <- kv$key
  if (!axis %in% c("sample", "subject", "contrast")) {
    stop("Subset axis must be one of sample, subject, contrast", call. = FALSE)
  }
  values <- .cli_cast_value(kv$value)
  if (is.numeric(values) && all(!is.na(values)) && all(abs(values - round(values)) < .Machine$double.eps^0.5)) {
    values <- as.integer(values)
  }
  list(axis = axis, values = values)
}

.cli_read_table <- function(path, kind) {
  if (!file.exists(path)) stop(kind, " file not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
    csv = data.table::fread(path, data.table = FALSE),
    tsv = data.table::fread(path, sep = "\t", data.table = FALSE),
    parquet = {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("arrow package required to read parquet metadata files", call. = FALSE)
      }
      as.data.frame(arrow::read_parquet(path))
    },
    data.table::fread(path, data.table = FALSE)
  )
  if (!is.data.frame(df) || !nrow(df)) stop(kind, " file has no rows: ", path, call. = FALSE)
  df
}

.cli_read_keyed_data <- function(path, id_col = NULL, kind = "metadata") {
  df <- .cli_read_table(path, kind = kind)
  id_col <- if (is.null(id_col) || !nzchar(id_col)) names(df)[1L] else as.character(id_col)[1L]
  if (!id_col %in% names(df)) {
    id_col <- names(df)[1L]
  }
  ids <- as.character(df[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids))) stop(kind, " id column contains missing/empty identifiers", call. = FALSE)
  if (anyDuplicated(ids)) stop(kind, " id column must be unique", call. = FALSE)
  rownames(df) <- ids
  df[[id_col]] <- NULL
  df
}

.cli_read_matrix <- function(path) {
  mat <- as.matrix(.cli_read_table(path, kind = "contrast-matrix"))
  storage.mode(mat) <- "double"
  mat
}

.cli_infer_format_from_path <- function(path) {
  if (is.null(path) || isTRUE(is.na(path))) return("auto")
  p <- tolower(as.character(path))
  if (grepl("\\.(csv|tsv|parquet)$", p)) return("tabular")
  if (grepl("\\.nii(\\.gz)?$", p)) return("nifti")
  if (grepl("\\.(h5|hdf5)$", p)) return("h5")
  "auto"
}

.cli_infer_out_format <- function(path) {
  if (is.null(path) || isTRUE(is.na(path))) return(NULL)
  p <- tolower(as.character(path))
  if (grepl("\\.(h5|hdf5)$", p)) return("h5")
  if (grepl("\\.csv$", p)) return("csv")
  if (grepl("\\.parquet$", p)) return("parquet")
  if (grepl("\\.nii(\\.gz)?$", p)) return("nifti")
  NULL
}

.cli_source_inputs <- function(opts, pos) {
  if (!is.null(opts[["input"]])) {
    vals <- unlist(lapply(opts[["input"]], .split_csv), use.names = FALSE)
    vals <- vals[nzchar(vals)]
    if (length(vals)) return(vals)
  }
  pos[nzchar(pos)]
}

.cli_source_flag_names <- function() {
  c(
    "input", "format", "prefer", "mask",
    "catalog-root", "catalog-pattern", "catalog-regex", "assay",
    "subject-col", "sample-col", "contrast-col", "beta-col", "var-col", "se-col",
    "contrast-data-cols", "temporal-policy", "contrast-matrix-file", "contrast-names",
    "source-option"
  )
}

.cli_assert_plan_source_conflicts <- function(opts, pos) {
  offenders <- .cli_source_flag_names()
  offenders <- offenders[vapply(offenders, function(name) !is.null(opts[[name]]), logical(1))]
  if (length(pos)) offenders <- c(offenders, "<positional-input>")
  if (length(offenders)) {
    stop(
      "--load-plan cannot be combined with source construction arguments: ",
      paste(offenders, collapse = ", "),
      call. = FALSE
    )
  }
}

.cli_source_options <- function(opts) {
  out <- .cli_parse_options_kv(opts[["source-option"]] %||% NULL)

  subject_col <- .cli_get_opt(opts, "subject-col", NULL)
  sample_col <- .cli_get_opt(opts, "sample-col", NULL)
  contrast_col <- .cli_get_opt(opts, "contrast-col", NULL)
  beta_col <- .cli_get_opt(opts, "beta-col", NULL)
  var_col <- .cli_get_opt(opts, "var-col", NULL)
  se_col <- .cli_get_opt(opts, "se-col", NULL)

  if (!is.null(subject_col)) out$subject_col <- subject_col
  if (!is.null(sample_col)) out$sample_col <- sample_col
  if (!is.null(contrast_col)) out$contrast_col <- contrast_col
  if (!is.null(beta_col) || !is.null(var_col) || !is.null(se_col)) {
    effect_cols <- out$effect_cols %||% list()
    if (!is.list(effect_cols)) effect_cols <- as.list(effect_cols)
    if (is.null(effect_cols$beta)) effect_cols$beta <- beta_col %||% "beta"
    if (is.null(effect_cols$var)) effect_cols$var <- var_col %||% "var"
    if (!is.null(se_col)) effect_cols$se <- se_col
    out$effect_cols <- effect_cols
  }

  if (!is.null(opts[["contrast-data-cols"]])) {
    out$contrast_data_cols <- unlist(lapply(opts[["contrast-data-cols"]], .split_csv), use.names = FALSE)
  }
  if (!is.null(opts[["temporal-policy"]])) {
    out$temporal_policy <- .cli_get_opt(opts, "temporal-policy")
  }
  if (!is.null(opts[["contrast-matrix-file"]])) {
    out$contrast_matrix <- .cli_read_matrix(.cli_get_opt(opts, "contrast-matrix-file"))
  }
  if (!is.null(opts[["contrast-names"]])) {
    out$contrast_names <- unlist(lapply(opts[["contrast-names"]], .split_csv), use.names = FALSE)
  }

  out
}

.cli_base_plan_from_opts <- function(opts, pos) {
  plan_file <- .cli_get_opt(opts, "load-plan", .cli_get_opt(opts, "plan", NULL))
  if (!is.null(plan_file)) {
    .cli_assert_plan_source_conflicts(opts, pos)
    return(load_plan(plan_file))
  }

  mask_path <- .cli_get_opt(opts, "mask", NULL)
  if (!is.null(opts[["catalog-root"]])) {
    root <- .cli_get_opt(opts, "catalog-root")
    pattern <- .cli_get_opt(opts, "catalog-pattern", "**/*.nii*")
    pregex <- .cli_get_opt(opts, "catalog-regex", NULL)
    catalog <- image_catalog(root = root, pattern = pattern, path_regex = pregex)
    if (!is.null(opts[["assay"]])) {
      assay_specs <- lapply(opts[["assay"]], .parse_kv)
      assay_args <- list()
      for (spec in assay_specs) assay_args[[spec$key]] <- spec$value
      catalog <- do.call(map_assays, c(list(catalog), assay_args))
    }
    return(as_gds(catalog, mask = mask_path))
  }

  inputs <- .cli_source_inputs(opts, pos)
  if (!length(inputs)) stop("Missing --input <path> (or use --load-plan)", call. = FALSE)
  source <- if (length(inputs) == 1L) inputs[[1L]] else inputs
  format <- .cli_get_opt(opts, "format", .cli_infer_format_from_path(inputs[[1L]]))
  prefer <- .cli_get_opt(opts, "prefer", NULL)
  gds_args <- c(
    list(source = source, format = format),
    if (!is.null(prefer)) list(prefer = prefer) else list(),
    if (!is.null(mask_path)) list(mask = mask_path) else list(),
    .cli_source_options(opts)
  )
  do.call(gds, gds_args)
}

.cli_attach_metadata <- function(plan, opts) {
  if (!is.null(opts[["col-data"]])) {
    id_col <- .cli_get_opt(opts, "col-data-id", "subject")
    plan <- with_col_data(plan, .cli_read_keyed_data(.cli_get_opt(opts, "col-data"), id_col, kind = "col-data"))
  }
  if (!is.null(opts[["row-data"]])) {
    id_col <- .cli_get_opt(opts, "row-data-id", "sample")
    plan <- with_row_data(plan, .cli_read_keyed_data(.cli_get_opt(opts, "row-data"), id_col, kind = "row-data"))
  }
  if (!is.null(opts[["contrast-data"]])) {
    id_col <- .cli_get_opt(opts, "contrast-data-id", "contrast")
    plan <- with_contrast_data(plan, .cli_read_keyed_data(.cli_get_opt(opts, "contrast-data"), id_col, kind = "contrast-data"))
  }
  plan
}

.cli_apply_pipeline_opts <- function(plan, opts) {
  subset_specs <- opts[["subset"]] %||% character()
  if (!is.null(opts[["subset-sample"]])) {
    subset_specs <- c(subset_specs, paste0("sample=", unlist(opts[["subset-sample"]], use.names = FALSE)))
  }
  if (!is.null(opts[["subset-subject"]])) {
    subset_specs <- c(subset_specs, paste0("subject=", unlist(opts[["subset-subject"]], use.names = FALSE)))
  }
  if (!is.null(opts[["subset-contrast"]])) {
    subset_specs <- c(subset_specs, paste0("contrast=", unlist(opts[["subset-contrast"]], use.names = FALSE)))
  }
  for (spec in subset_specs) {
    parsed <- .cli_parse_subset_spec(spec)
    plan <- do.call(subset, c(list(plan), stats::setNames(list(parsed$values), parsed$axis)))
  }

  if (!is.null(opts[["mask-policy"]]) || !is.null(opts[["mask-scope"]]) || !is.null(opts[["mask-threshold"]])) {
    rule <- .cli_get_opt(opts, "mask-policy", "intersection")
    if (identical(rule, "custom")) {
      stop("CLI does not support custom mask callbacks; use plan JSON or R for that workflow", call. = FALSE)
    }
    scope <- .cli_get_opt(opts, "mask-scope", "group")
    threshold <- as.numeric(.cli_get_opt(opts, "mask-threshold", 0.95))
    plan <- mask(plan, MaskPolicy(scope = scope, rule = rule, threshold = threshold))
  }

  if (!is.null(opts[["derive"]])) {
    what <- unlist(lapply(opts[["derive"]], .split_csv), use.names = FALSE)
    if (length(what)) plan <- derive(plan, what)
  }

  if (!is.null(opts[["align"]])) {
    for (family in opts[["align"]]) {
      plan <- align(plan, family = as.character(family))
    }
  }

  if (!is.null(opts[["reduce"]])) {
    reduce_opts <- .cli_parse_options_kv(opts[["reduce-option"]] %||% NULL)
    plan <- reduce(
      plan,
      method = as.character(.cli_get_opt(opts, "reduce")),
      weights = as.character(.cli_get_opt(opts, "weights", "1/var")),
      by = as.character(.cli_get_opt(opts, "reduce-by", "contrast")),
      formula = .cli_get_opt(opts, "formula", NULL),
      options = reduce_opts
    )
  }

  if (!is.null(opts[["posthoc"]])) {
    posthoc_opts <- .cli_parse_options_kv(opts[["posthoc-option"]] %||% NULL)
    for (method in opts[["posthoc"]]) {
      plan <- posthoc(plan, method = as.character(method), options = posthoc_opts)
    }
  }

  plan
}

.cli_node_summary <- function(node) {
  switch(node$op,
    reduce = paste0("reduce:", node$method),
    map = paste0("map:", .space_brief(node$target_space)),
    mask_policy = paste0("mask:", node$policy$rule),
    subset_axis = "subset",
    align_to_group = paste0("align:", node$family_name %||% node$family$name %||% "(anonymous)"),
    write = paste0("write:", node$format),
    posthoc = paste0("posthoc:", node$method),
    derive = paste0("derive:", paste(node$what %||% character(), collapse = ",")),
    node$op
  )
}

.cli_plan_summary <- function(plan) {
  probe <- plan$source$probe %||% list()
  dims <- probe$dims %||% NULL
  list(
    adapter = plan$source$adapter,
    dims = if (is.null(dims)) NULL else as.list(as.integer(dims)),
    assays = unname(as.character(probe$assays %||% character())),
    subjects = unname(as.character(probe$subjects %||% character())),
    contrasts = unname(as.character(probe$contrasts %||% character())),
    space = .space_brief(probe$space %||% NULL),
    has_col_data = !is.null(col_data(plan)),
    has_row_data = !is.null(row_data(plan)),
    has_contrast_data = !is.null(contrast_data(plan)),
    map_families = names(.map_registry(plan)),
    nodes = lapply(plan$nodes, function(node) {
      list(
        op = node$op,
        summary = .cli_node_summary(node)
      )
    })
  )
}

.cli_gds_summary <- function(g) {
  dims <- dim(assays(g)[[1L]])
  ai <- lapply(names(assays(g)), assay_info)
  roles <- stats::setNames(vapply(ai, function(x) x$role %||% "?", character(1)), names(ai))
  list(
    dims = as.list(as.integer(dims)),
    space = .space_brief(g$space),
    subjects = as.character(g$subjects),
    contrasts = as.character(g$contrasts),
    assays = as.list(roles),
    has_col_data = !is.null(g$col_data),
    has_row_data = !is.null(g$row_data),
    map_families = names(g$metadata$map_families %||% list())
  )
}

.cli_array_summary <- function(arrays) {
  lapply(arrays, function(arr) {
    list(
      dim = as.list(as.integer(dim(arr))),
      values = as.vector(arr)
    )
  })
}

.cli_emit_json <- function(x) {
  .cli_print(jsonlite::toJSON(x, pretty = TRUE, auto_unbox = TRUE, null = "null", dataframe = "rows"), "\n")
}

.cli_save_plan_if_requested <- function(plan, opts) {
  path <- .cli_get_opt(opts, "save-plan", NULL)
  if (!is.null(path)) save_plan(plan, path)
  invisible(plan)
}

.cli_plan_has_write <- function(plan) {
  any(vapply(plan$nodes, function(node) identical(node$op, "write"), logical(1)))
}

.cli_render_table <- function(df, max_rows = 20L) {
  if (!is.data.frame(df) || !nrow(df)) {
    .cli_print("(no rows)\n")
    return(invisible(df))
  }
  shown <- utils::head(df, max_rows)
  print(shown, row.names = FALSE)
  if (nrow(df) > max_rows) {
    .cli_print("... ", nrow(df) - max_rows, " more row(s)\n")
  }
  invisible(df)
}

.cli_help_main <- function() {
  .cli_print(
    "fmrigds -- ergonomic command-line access to fmrigds plans and analyses\n\n",
    "Usage:\n",
    "  fmrigds <command> [options]\n",
    "  fmrigds help [command]\n\n",
    "Core commands:\n",
    "  run       Build or load a plan, execute it, and write outputs\n",
    "  plan      Build or load a plan, inspect it, validate it, or save it\n",
    "  probe     Inspect an input source or saved plan\n",
    "  preview   Execute a small preview block through a plan\n",
    "  list      List reducers, posthoc methods, or adapters\n\n",
    "Registry aliases:\n",
    "  reducers  Same as: fmrigds list reducers\n",
    "  posthoc   Same as: fmrigds list posthoc\n",
    "  adapters  Same as: fmrigds list adapters\n\n",
    "Global options:\n",
    "  -h, --help        Show help\n",
    "  -v, --version     Print version\n\n",
    "Workflow:\n",
    "  1. fmrigds probe --input group.csv\n",
    "  2. fmrigds plan  --input group.csv --derive t,p --reduce fixed --save-plan plan.json\n",
    "  3. fmrigds run   --load-plan plan.json --out results.h5\n\n",
    "Try:\n",
    "  fmrigds help run\n",
    "  fmrigds help plan\n"
  )
}

.cli_help_shared_source <- function() {
  .cli_print(
    "Source options:\n",
    "  --input <path[,path...]>       Input file(s); positional paths also work\n",
    "  --load-plan <file>             Load a saved plan JSON instead of probing a source\n",
    "  --format <auto|tabular|nifti|h5|fmristore|nftab>\n",
    "  --prefer <adapter>             Preferred adapter when auto-detect is ambiguous\n",
    "  --source-option <k=v>          Repeatable adapter-specific passthrough\n\n",
    "Metadata options:\n",
    "  --col-data <path>              Subject metadata CSV/TSV/Parquet\n",
    "  --col-data-id <column>         Key column; default: subject\n",
    "  --row-data <path>              Sample metadata CSV/TSV/Parquet\n",
    "  --row-data-id <column>         Key column; default: sample\n",
    "  --contrast-data <path>         Contrast metadata CSV/TSV/Parquet\n",
    "  --contrast-data-id <column>    Key column; default: contrast\n\n",
    "Catalog mode:\n",
    "  --catalog-root <dir>\n",
    "  --catalog-pattern <glob>       Default: **/*.nii*\n",
    "  --catalog-regex <regex>        Optional named capture groups\n",
    "  --assay <name=pattern>         Repeatable, e.g. beta=cope, se=varcope\n",
    "  --mask <path>                  Optional mask for NIfTI catalog ingestion\n\n",
    "Tabular/fmristore helpers:\n",
    "  --subject-col <name>           Default: subject\n",
    "  --sample-col <name>            Default: sample\n",
    "  --contrast-col <name>          Default: contrast\n",
    "  --beta-col <name>              Default: beta\n",
    "  --var-col <name>               Default: var\n",
    "  --se-col <name>                Optional standard error column\n",
    "  --contrast-data-cols <cols>    Comma-separated columns to lift into contrast_data\n",
    "  --temporal-policy <as_is|mean|design>\n",
    "  --contrast-matrix-file <path>  Numeric CSV/TSV matrix for fmristore design mode\n",
    "  --contrast-names <names>       Comma-separated contrast labels\n"
  )
}

.cli_help_shared_pipeline <- function() {
  .cli_print(
    "Pipeline options:\n",
    "  --subset <axis=values>         Repeatable; axis is sample|subject|contrast\n",
    "  --subset-sample <values>       Shortcut for --subset sample=...\n",
    "  --subset-subject <values>      Shortcut for --subset subject=...\n",
    "  --subset-contrast <values>     Shortcut for --subset contrast=...\n",
    "  --mask-policy <rule>           intersection, union, or threshold\n",
    "  --mask-scope <group|subject>   Default: group\n",
    "  --mask-threshold <num>         Default: 0.95\n",
    "  --derive <stat[,stat...]>      e.g. t,z,p\n",
    "  --align <family>               Repeatable registered alignment family\n",
    "  --reduce <method>              e.g. fixed, random, meta:fe_reg, lmm:ri\n",
    "  --weights <1/var|n_eff|equal|custom>  Default: 1/var\n",
    "  --reduce-by <axis>             Default: contrast\n",
    "  --formula <R formula>          Quote it, e.g. \"~ 1 + age + group\"\n",
    "  --reduce-option <k=v>          Repeatable reducer options\n",
    "  --posthoc <method>             Repeatable, e.g. fdr:bh\n",
    "  --posthoc-option <k=v>         Applied to each posthoc step\n",
    "  --save-plan <file>             Write the resulting plan JSON\n"
  )
}

.cli_help_run <- function() {
  .cli_print(
    "fmrigds run -- execute a group analysis plan\n\n",
    "Usage:\n",
    "  fmrigds run --input group.csv --reduce fixed --out results.h5\n",
    "  fmrigds run --load-plan plan.json\n",
    "  fmrigds run --catalog-root study_dir --assay beta=cope --assay se=varcope \\\n",
    "    --reduce fixed --derive t,p --posthoc fdr:bh --out group_results.h5\n\n"
  )
  .cli_help_shared_source()
  .cli_print("\n")
  .cli_help_shared_pipeline()
  .cli_print(
    "\nOutput options:\n",
    "  --out <path>                   Output path; optional if loaded plan already writes\n",
    "  --out-format <h5|csv|parquet|nifti>\n",
    "  --write-option <k=v>           Repeatable write_out options\n",
    "  --dry-run                      Print the plan instead of computing it\n",
    "  --no-validate                  Skip validate(plan) before compute\n",
    "  --verbose                      Print plan details before compute\n",
    "  --json                         Emit machine-readable summary on dry-run or completion\n"
  )
}

.cli_help_plan <- function() {
  .cli_print(
    "fmrigds plan -- build, inspect, validate, and save plans\n\n",
    "Usage:\n",
    "  fmrigds plan --input group.csv --derive t,p --reduce fixed\n",
    "  fmrigds plan --load-plan plan.json --posthoc fdr:bh --save-plan new-plan.json\n\n"
  )
  .cli_help_shared_source()
  .cli_print("\n")
  .cli_help_shared_pipeline()
  .cli_print(
    "\nPlan options:\n",
    "  --no-validate                  Skip validate(plan)\n",
    "  --json                         Emit plan summary as JSON\n"
  )
}

.cli_help_probe <- function() {
  .cli_print(
    "fmrigds probe -- inspect a source or saved plan\n\n",
    "Usage:\n",
    "  fmrigds probe --input group.csv\n",
    "  fmrigds probe --catalog-root study_dir --assay beta=cope --assay se=varcope\n",
    "  fmrigds probe --load-plan plan.json\n\n"
  )
  .cli_help_shared_source()
  .cli_print(
    "\nProbe options:\n",
    "  --json                         Emit probe summary as JSON\n"
  )
}

.cli_help_preview <- function() {
  .cli_print(
    "fmrigds preview -- execute a small preview block through a plan\n\n",
    "Usage:\n",
    "  fmrigds preview --input group.csv --reduce fixed --n 3\n",
    "  fmrigds preview --load-plan plan.json --show-assay beta --include-col-data\n\n"
  )
  .cli_help_shared_source()
  .cli_print("\n")
  .cli_help_shared_pipeline()
  .cli_print(
    "\nPreview options:\n",
    "  --n <int>                      Number of samples to preview; default: 5\n",
    "  --show-assay <name[,name...]>  Restrict printed assay columns\n",
    "  --include-col-data             Join col_data columns into the preview table\n",
    "  --raw                          Preview adapter-level arrays instead of final GDS rows\n",
    "  --json                         Emit preview data as JSON\n"
  )
}

.cli_help_list <- function() {
  .cli_print(
    "fmrigds list -- inspect CLI registries\n\n",
    "Usage:\n",
    "  fmrigds list reducers\n",
    "  fmrigds list posthoc\n",
    "  fmrigds list adapters\n\n",
    "Options:\n",
    "  --json                         Emit registry contents as JSON\n"
  )
}

.cli_registry_rows <- function(kind) {
  switch(kind,
    reducers = {
      ids <- list_reducers()
      lapply(ids, function(id) {
        reducer <- get_reducer(id)
        list(
          name = id,
          requires = as.character(reducer$requires %||% character()),
          provides = as.character(reducer$provides %||% character()),
          input_shape = reducer$input_shape %||% "contrastwise"
        )
      })
    },
    posthoc = {
      ids <- list_posthoc()
      lapply(ids, function(id) {
        handler <- get_posthoc(id)
        list(
          name = id,
          requires = as.character(handler$requires %||% character()),
          provides = as.character(handler$provides %||% character())
        )
      })
    },
    adapters = {
      ids <- sort(ls(.adapter_registry))
      lapply(ids, function(id) list(name = id))
    },
    stop("Unknown registry: ", kind, call. = FALSE)
  )
}

.cli_print_registry <- function(kind, rows) {
  title <- switch(kind,
    reducers = "Reducers",
    posthoc = "Post-hoc methods",
    adapters = "Adapters",
    kind
  )
  .cli_print(title, ":\n")
  if (!length(rows)) {
    .cli_print("  (none)\n")
    return(invisible(rows))
  }
  for (row in rows) {
    .cli_print("  ", row$name, "\n")
    if (!is.null(row$requires) && length(row$requires)) {
      .cli_print("    requires: ", paste(row$requires, collapse = ", "), "\n")
    }
    if (!is.null(row$provides) && length(row$provides)) {
      .cli_print("    provides: ", paste(row$provides, collapse = ", "), "\n")
    }
    if (!is.null(row$input_shape) && nzchar(row$input_shape)) {
      .cli_print("    input-shape: ", row$input_shape, "\n")
    }
  }
  invisible(rows)
}

.cli_cmd_run <- function(args) {
  parsed <- .cli_parse_args(args)
  opts <- parsed$opts
  pos <- parsed$pos

  if (.cli_has_flag(opts, "help")) {
    .cli_help_run()
    return(invisible(NULL))
  }

  plan <- .cli_base_plan_from_opts(opts, pos)
  plan <- .cli_attach_metadata(plan, opts)
  plan <- .cli_apply_pipeline_opts(plan, opts)

  out_path <- .cli_get_opt(opts, "out", NULL)
  if (!is.null(out_path)) {
    out_format <- .cli_get_opt(opts, "out-format", NULL)
    if (is.null(out_format) || !nzchar(out_format)) out_format <- .cli_infer_out_format(out_path)
    if (is.null(out_format)) stop("Could not infer --out-format from: ", out_path, call. = FALSE)
    write_opts <- .cli_parse_options_kv(opts[["write-option"]] %||% NULL)
    plan <- write_out(plan, path = out_path, format = out_format, options = write_opts)
  } else if (!.cli_plan_has_write(plan) && !.cli_has_flag(opts, "dry-run")) {
    stop("--out <path> is required unless the loaded/built plan already contains a write step", call. = FALSE)
  }

  .cli_save_plan_if_requested(plan, opts)

  if (.cli_has_flag(opts, "dry-run")) {
    if (.cli_has_flag(opts, "json")) {
      .cli_emit_json(.cli_plan_summary(plan))
    } else {
      explain(plan)
      explain_plan(plan)
    }
    return(invisible(plan))
  }

  if (!.cli_has_flag(opts, "no-validate")) validate(plan)
  if (.cli_has_flag(opts, "verbose") && !.cli_has_flag(opts, "json")) {
    explain(plan)
    explain_plan(plan)
  }

  g <- compute(plan)
  if (.cli_has_flag(opts, "json")) {
    .cli_emit_json(.cli_gds_summary(g))
  } else {
    explain(g)
  }
  invisible(g)
}

.cli_cmd_plan <- function(args) {
  parsed <- .cli_parse_args(args)
  opts <- parsed$opts
  pos <- parsed$pos

  if (.cli_has_flag(opts, "help")) {
    .cli_help_plan()
    return(invisible(NULL))
  }

  plan <- .cli_base_plan_from_opts(opts, pos)
  plan <- .cli_attach_metadata(plan, opts)
  plan <- .cli_apply_pipeline_opts(plan, opts)
  if (!.cli_has_flag(opts, "no-validate")) validate(plan)
  .cli_save_plan_if_requested(plan, opts)

  if (.cli_has_flag(opts, "json")) {
    .cli_emit_json(.cli_plan_summary(plan))
  } else {
    explain(plan)
    explain_plan(plan)
  }
  invisible(plan)
}

.cli_cmd_probe <- function(args) {
  parsed <- .cli_parse_args(args)
  opts <- parsed$opts
  pos <- parsed$pos

  if (.cli_has_flag(opts, "help")) {
    .cli_help_probe()
    return(invisible(NULL))
  }

  plan <- .cli_base_plan_from_opts(opts, pos)
  plan <- .cli_attach_metadata(plan, opts)
  if (.cli_has_flag(opts, "json")) {
    .cli_emit_json(.cli_plan_summary(plan))
  } else {
    explain(plan)
  }
  invisible(plan)
}

.cli_cmd_preview <- function(args) {
  parsed <- .cli_parse_args(args)
  opts <- parsed$opts
  pos <- parsed$pos

  if (.cli_has_flag(opts, "help")) {
    .cli_help_preview()
    return(invisible(NULL))
  }

  plan <- .cli_base_plan_from_opts(opts, pos)
  plan <- .cli_attach_metadata(plan, opts)
  plan <- .cli_apply_pipeline_opts(plan, opts)
  .cli_save_plan_if_requested(plan, opts)

  n <- as.integer(.cli_get_opt(opts, "n", 5L))
  show_assays <- if (is.null(opts[["show-assay"]])) character() else unlist(lapply(opts[["show-assay"]], .split_csv), use.names = FALSE)

  if (.cli_has_flag(opts, "raw")) {
    assays_req <- if (length(show_assays)) show_assays else plan$source$probe$assays
    arrays <- preview(plan, n = n, assays = assays_req)
    if (.cli_has_flag(opts, "json")) {
      .cli_emit_json(.cli_array_summary(arrays))
    } else {
      for (name in names(arrays)) {
        .cli_print(name, " dim=", paste(dim(arrays[[name]]), collapse = "x"), "\n")
      }
    }
    return(invisible(arrays))
  }

  g <- preview(plan, n = n)
  include_col_data <- .cli_has_flag(opts, "include-col-data")
  df <- gds_to_tibble(
    g,
    assays = if (length(show_assays)) show_assays else NULL,
    include_col_data = include_col_data
  )

  if (.cli_has_flag(opts, "json")) {
    .cli_emit_json(list(summary = .cli_gds_summary(g), data = df))
  } else {
    explain(g)
    .cli_render_table(df)
  }
  invisible(g)
}

.cli_cmd_list <- function(args, alias_kind = NULL) {
  parsed <- .cli_parse_args(args)
  opts <- parsed$opts
  pos <- parsed$pos

  if (.cli_has_flag(opts, "help")) {
    .cli_help_list()
    return(invisible(NULL))
  }

  kind <- if (!is.null(alias_kind)) alias_kind else if (length(pos)) pos[[1L]] else ""
  if (!kind %in% c("reducers", "posthoc", "adapters")) {
    stop("list expects one of: reducers, posthoc, adapters", call. = FALSE)
  }

  rows <- .cli_registry_rows(kind)
  if (.cli_has_flag(opts, "json")) {
    .cli_emit_json(rows)
  } else {
    .cli_print_registry(kind, rows)
  }
  invisible(rows)
}

.cli_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- .cli_parse_args(args)
  if (.cli_has_flag(parsed$opts, "version")) {
    .cli_version()
    return(invisible(NULL))
  }

  if (!length(args) || .cli_has_flag(parsed$opts, "help")) {
    .cli_help_main()
    return(invisible(NULL))
  }

  cmd <- args[[1L]]
  rest <- if (length(args) > 1L) args[-1L] else character()

  if (identical(cmd, "help")) {
    sub <- if (length(rest)) rest[[1L]] else ""
    switch(sub,
      run = .cli_help_run(),
      plan = .cli_help_plan(),
      probe = .cli_help_probe(),
      preview = .cli_help_preview(),
      list = .cli_help_list(),
      reducers = .cli_help_list(),
      posthoc = .cli_help_list(),
      adapters = .cli_help_list(),
      .cli_help_main()
    )
    return(invisible(NULL))
  }

  switch(cmd,
    run = .cli_cmd_run(rest),
    plan = .cli_cmd_plan(rest),
    probe = .cli_cmd_probe(rest),
    preview = .cli_cmd_preview(rest),
    list = .cli_cmd_list(rest),
    reducers = .cli_cmd_list(rest, alias_kind = "reducers"),
    posthoc = .cli_cmd_list(rest, alias_kind = "posthoc"),
    adapters = .cli_cmd_list(rest, alias_kind = "adapters"),
    stop("Unknown command: ", cmd, call. = FALSE)
  )
}

#' CLI entrypoint for the `fmrigds` command
#'
#' Intended to be called from a small wrapper script in `inst/bin/fmrigds`.
#' The CLI mirrors the package's lazy plan grammar:
#'
#' - `probe` inspects a source or saved plan
#' - `plan` builds, validates, and saves plans without executing them
#' - `run` executes plans and writes outputs
#' - `preview` materialises a small sample block for quick inspection
#' - `list` exposes reducer, post-hoc, and adapter registries
#'
#' Power-user control is available through repeatable `--source-option`,
#' `--reduce-option`, `--posthoc-option`, and `--write-option` flags, plus
#' `--load-plan` / `--save-plan` for plan-first workflows.
#'
#' @keywords internal
fmrigds_cli_exec <- function(args = commandArgs(trailingOnly = TRUE)) {
  tryCatch(
    .cli_main(args),
    error = function(e) {
      .cli_warn("fmrigds: error: ", conditionMessage(e))
      .cli_warn("Run `fmrigds --help` for usage.")
      .cli_quit(1L)
    }
  )
}
# nocov end
