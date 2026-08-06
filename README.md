# fmrigds

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/fmrigds/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/fmrigds/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/fmrigds/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/bbuchsbaum/fmrigds/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

**Format-agnostic group-level analysis for fMRI**

`fmrigds` provides a unified framework for working with first-level fMRI statistical maps across multiple data formats and spatial representations. The package implements a lazy evaluation model that enables efficient, reproducible group and meta-analyses without requiring data format conversions or space transformations upfront.

## Key Features

### Universal Data Access
- **Multiple input formats**: CSV/TSV, Parquet, NIfTI, HDF5, and fmristore layouts
- **Automatic format detection**: Point to a file path and the appropriate adapter loads automatically
- **Multiple spatial representations**: Voxels, parcels/ROIs, basis components (PCA/ICA), or surface vertices

### Lazy Evaluation
Build analysis pipelines declaratively without loading data into memory:
- `gds()` creates a lazy plan from any supported source
- Chain operations (`subset()`, `derive()`, `map_to()`, `reduce()`, `mask()`, `write_out()`)
- `compute()` executes the entire pipeline, with automatic optimization and streaming for large datasets

### Statistical Methods
- **Derivation engine**: Automatically compute derived statistics (t, z, p-values, effect sizes)
- **Variance propagation**: Correct uncertainty handling through transformations and mappings
- **Group-level reducers**: Fixed-effects, random-effects (DerSimonian-Laird), meta-regression
- **Repeated-measures LMMs**: Fast Gaussian mixed models for shared-design multiresponse data with random intercepts or one random slope
- **Evidence combiners**: Stouffer's method, Fisher's method, Lancaster's method for combining statistics

### Data Export
- **HDF5**: Native `/gds` format with full provenance and metadata
- **Tabular**: CSV/Parquet export for ROI or parcel-level results
- **NIfTI**: Voxel-space maps compatible with standard neuroimaging tools
- All dependencies for specific formats are optional (install only what you need)

### Reproducibility
- **Full provenance tracking**: Every operation is recorded with timestamps and parameters
- **Persistent metadata**: Provenance, spatial alignments, and analysis parameters saved with outputs
- **Computational digests**: Unique identifiers for analysis pipelines enable exact reproducibility

## Quick Start

### Basic workflow

```r
library(fmrigds)

# Load data from any format (auto-detected)
plan <- gds("roi_stats.csv")

# Build a lazy analysis pipeline
plan <- plan |>
  subset(contrast = "Faces>Places") |>
  derive("t") |>
  reduce(method = "fixed")

# Execute and get results
result <- compute(plan)

# Access results
assays(result)         # Statistical maps
subjects(result)       # "meta" (group-level)
space(result)          # Spatial representation
```

## Command-line interface

This package bundles an `fmrigds` command wrapper under `inst/bin/`. After
installing `fmrigds`, resolve its installed path via `system.file()` and then
use the CLI to probe sources, assemble plans, preview small blocks, and run
full analyses:

```bash
# Show help
"$(Rscript -e 'cat(system.file(\"bin\", \"fmrigds\", package = \"fmrigds\"))')" --help

# Inspect a source
"$(Rscript -e 'cat(system.file(\"bin\", \"fmrigds\", package = \"fmrigds\"))')" probe \
  --input group.csv

# Build and save a plan without executing it
"$(Rscript -e 'cat(system.file(\"bin\", \"fmrigds\", package = \"fmrigds\"))')" plan \
  --input group.csv \
  --derive t,p \
  --reduce fixed \
  --posthoc fdr:bh \
  --save-plan fixed-plan.json

# Run a fixed-effects analysis from a CSV and write an HDF5 GDS
"$(Rscript -e 'cat(system.file(\"bin\", \"fmrigds\", package = \"fmrigds\"))')" run \
  --input group.csv --reduce fixed --out results.h5

# Reuse a saved plan later
"$(Rscript -e 'cat(system.file(\"bin\", \"fmrigds\", package = \"fmrigds\"))')" run \
  --load-plan fixed-plan.json --out results.h5
```

Main commands:

- `probe`: summarize a source or saved plan
- `plan`: build, validate, inspect, and save plans
- `preview`: execute a small sample block and print a tidy preview table
- `run`: execute a plan and write outputs
- `list`: inspect reducers, post-hoc methods, and adapters

The CLI is designed around the same lazy grammar as the R API. Common paths are
first-class (`--derive`, `--reduce`, `--posthoc`, `--out`), while more advanced
workflows can still be steered through repeatable passthrough flags such as
`--source-option`, `--reduce-option`, `--posthoc-option`, and `--write-option`.
For plan-first workflows, `--save-plan` and `--load-plan` let you separate plan
construction from execution.

## Interop Surface (for external packages)

These small, public helpers make it easy to “just use GDS” from other R
packages without re-implementing adapters or boilerplate.

- Coercion front door
  - `as_gds(x)`: Turn in-memory structures into a validated GDS
  - Methods for `list` (named 3D arrays), `array` (3D), and `data.frame`
  - `as_gds.data.frame()` mirrors CSV/Parquet ingestion via the tabular adapter

- Covariates and models
  - `with_col_data(x, df)`: Attach/align subject-level covariates to a plan or GDS
  - `model_matrix(x, ~ formula)`: Build a model matrix from attached col_data

- Tidy export
  - `gds_to_tibble(g, assays, drop_na)`: Long-form table for downstream analysis

- Introspection and validation
  - `explain(x)`: Human-readable summary for a plan or realised GDS
  - `validate(x)`: Structural checks with actionable messages
  - `explain_plan(plan)`: Tidy table of pending plan nodes
  - `preview(plan, n, assays=NULL)`: Execute a tiny block for a quick peek

- Post-hoc registry (FDR, etc.)
  - `register_posthoc()`, `list_posthoc()`, `get_posthoc()`
  - Built-ins: `"fdr:bh"`, `"fdr:by"`
  - `posthoc(plan, method)`: Lazy verb to add a post-hoc step

- Compatibility helpers
  - `assert_compatible_spaces(g1, g2)`; `common_mask(g1, g2, rule)`
  - `harmonise_contrasts(g, map)`; `relabel_subjects(g, mapping)`

- Alignment helpers (map families)
  - `make_linear_family()`, `make_warp_family()` for common alignment inputs
  - Sugar: `register_alignment()`, `list_alignments()`, `get_alignment()`
  - Families persist through HDF5 and are used by `align()`

- Weight hooks for reducers
  - `attach_weight(g, name, array)`; `use_weight(g, name)` → `reduce(weights="custom")`

- Space utilities
  - `space_from_nifti(path, mask=NULL)`; `space_subset(space, idx)`

### Meta-analysis with covariates

```r
# Subject-level covariates for meta-regression
col_data <- data.frame(
  age = c(25, 32, 28, 35),
  group = c("control", "patient", "control", "patient"),
  row.names = c("sub-01", "sub-02", "sub-03", "sub-04")
)

# Fixed-effects meta-regression
plan <- gds("roi_stats.csv", col_data = col_data) |>
  reduce(method = "meta:fe_reg", formula = ~ age + group)

# Or attach covariates later
plan <- with_col_data(plan, col_data)

result <- compute(plan)

# Results include coefficients and SEs for each parameter
names(assays(result))
# "coef:(Intercept)", "coef:age", "coef:grouppatient",
# "se_coef:(Intercept)", "se_coef:age", "se_coef:grouppatient"
```

### Examine the cohort before interpreting the group fit

`examine_group()` branches from the subject-level prefix of an analysis. It
uses the reducer's frozen formula, estimands, variance rules, and subject order
to report three different questions: whether the data have a validity concern,
whether each subject is unexpected under the fitted model, and how much deleting
that subject changes the requested group statistic. It does not change the
analysis or classify subjects for removal.

```r
analysis <- gds(paths, col_data = participants) |>
  subset(contrast = c("faces>places", "tools>faces")) |>
  reduce(
    method = "meta:re_reg",
    formula = ~ group + age + sex + site
  ) |>
  posthoc("fdr:bh")

exam <- examine_group(
  analysis,
  estimands = "grouppatient",
  quality = c("mean_fd", "tsnr")
)

plot(exam)                         # surprise versus influence
plot(exam, type = "embedding")    # model-adjusted residual geometry
plot(exam, subject = "sub-017")   # retained subject maps
write_report(exam, "group-examination.html")

fit <- compute(analysis)
```

The review queue uses absolute, stability-aware gates. The percentile-based
`review_priority` only orders inspection. For random-effects models, the first
pass holds the full-data heterogeneity estimate fixed; retained subjects then
receive exact refits, with both modes recorded in the result. See
`vignette("group-examination")` for the result contract, availability states,
and interpretation of the plots.

### Repeated-measures mixed models

For common neuroimaging repeated-measures workflows, `reduce()` also supports a
restricted but efficient Gaussian mixed-model family. The fast path assumes the
same observation layout and the same fixed/random design for every sample:

```r
result <- gds("roi_long.csv", contrast_data_cols = "time") |>
  reduce(
    method = "lmm:ri_slope1",
    formula = ~ time,
    options = list(
      slope = "time",
      covariance = "diag",
      fit = "REML",
      theta_mode = "voxelwise"
    )
  ) |>
  compute()

names(assays(result))
# "coef:(Intercept)", "coef:time", "vc_intercept", "vc_slope",
# "lambda_intercept", "lambda_slope", ...
```

Available repeated-measures reducers:

- `method = "lmm:ri"` for a random-intercept model
- `method = "lmm:ri_slope1"` for a random intercept plus one within-subject slope
- `options$theta_mode = "pooled"` to share variance parameters across samples
- `options$theta_mode = "voxelwise"` to fit variance parameters separately per sample

This is intentionally narrower than `lmer()`: one grouping factor only,
Gaussian outcomes only, and no general random-effects formula parser.

### Working with different formats

```r
# Tabular data (CSV, TSV, Parquet)
plan1 <- gds("results.csv")

# NIfTI volumes
plan2 <- gds("cope1.nii.gz")

# HDF5 with /gds group
plan3 <- gds("group_data.h5")

# fmristore layouts (auto-detected)
plan4 <- gds("subject01_fmristore.h5")
```

### Efficient streaming for large datasets

```r
# Process in blocks without loading everything into memory
result <- compute(plan, block = list(sample = 10000))

# Or write directly to disk
plan_with_write <- write_out(plan, "output.h5", format = "h5")
compute(plan_with_write)
```

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("bbuchsbaum/fmrigds")

# Optional packages for specific workflows
install.packages(c("arrow", "neuroim2", "neurotabs", "neurothresh"))
```

## Development

```bash
# Run test coverage (requires the covr package)
make coverage
```

## Documentation

Browse the complete [fmrigds documentation website](https://bbuchsbaum.github.io/fmrigds/),
read the source vignettes directly on GitHub, or run
`browseVignettes("fmrigds")` after installation to open the built articles locally.

- [Getting started](vignettes/fmrigds.Rmd) --- core pipeline tutorial
- [Command-line workflows](vignettes/cli-workflows.Rmd) --- probe, plan, preview, and run analyses from the shell
- [Repeated-measures mixed models](vignettes/repeated-measures-lmm.Rmd) --- restricted Gaussian LMM workflow
- [Spatial operations](vignettes/spatial-operations.Rmd) --- masking, alignment, and space mapping
- [Post-hoc corrections and spatial FDR](vignettes/as-plan-and-spatial-fdr.Rmd) --- standard and spatial FDR
- [Group examination](vignettes/group-examination.Rmd) --- model surprise, influence, residual geometry, and review reports
- [fmristore HDF5 ingestion](vignettes/fmristore-ingestion.Rmd) --- reading fmristore HDF5 files
- [Technical details](notes/TECHNICAL_SPECIFICATION.md) --- full design specification
- **Function reference:** `?gds`, `?compute`, `?reduce`, `?align`, `?mask`, and `?map_to`

## Release Focus

The current development branch is focused on a clean path to `1.0.0`:

- polished end-user workflows for tabular ROI data, voxelwise NIfTI inputs, and HDF5/fmristore sources
- a stable core grammar built around `gds()`, the lazy verbs, and `compute()`
- cross-format consistency, provenance-aware export, and predictable object summaries

## Contributing

Bug reports and feature requests are welcome on the [issue tracker](https://github.com/bbuchsbaum/fmrigds/issues).

## Citation

If you use fmrigds in your research, please cite:

```
Buchsbaum, B. R. (2025). fmrigds: Format-agnostic group-level analysis for fMRI.
R package. https://github.com/bbuchsbaum/fmrigds
```

## License

MIT

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The defaults are configured via `params$family` and `params$preset` (family = 'red', preset = 'interaction'). The pkgdown site uses `template: { package: albersdown }` together with generated `pkgdown/extra.css` and `pkgdown/extra.js` so the theme is linked and activated on site pages.
<!-- albersdown:theme-note:end -->
