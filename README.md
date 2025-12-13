# fmrigds

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
plan <- plan %>%
  subset(contrast = "Faces>Places") %>%
  derive("t") %>%
  reduce(method = "fixed")

# Execute and get results
result <- compute(plan)

# Access results
assays(result)         # Statistical maps
subjects(result)       # "meta" (group-level)
space(result)          # Spatial representation
```

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
plan <- gds("roi_stats.csv", col_data = col_data) %>%
  reduce(method = "meta:fe_reg", formula = ~ age + group)

# Or attach covariates later
plan <- with_col_data(plan, col_data)

result <- compute(plan)

# Results include coefficients and SEs for each parameter
names(assays(result))
# "coef:(Intercept)", "coef:age", "coef:grouppatient",
# "se_coef:(Intercept)", "se_coef:age", "se_coef:grouppatient"
```

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

# Optional dependencies for specific formats
install.packages(c("arrow", "hdf5r", "RNifti"))
```

## Development

```bash
# Run test coverage (requires the covr package)
make coverage
```

## Documentation

- **Getting Started**: See `vignette("intro-gdsfmri")` for a complete tutorial
- **Working with realized GDS and spatial FDR**: See `vignette("as-plan-and-spatial-fdr")`
- **fmristore Integration**: See `vignette("fmristore-ingestion")` for reading fmristore files
- **Technical Details**: See `TECHNICAL_SPECIFICATION.md` for the full design specification
- **Package Documentation**: `?gds`, `?compute`, `?reduce` for function references

## Development Status

The package is functional and tested with comprehensive test coverage (224+ tests).
Core APIs are stabilizing for the 0.1.0 release.

### Planned Enhancements
- Enhanced visualization tools for provenance graphs
- Additional storage backends (TileDB, Arrow IPC)
- Optimized spatial mapping kernels
- pkgdown documentation site

## Contributing

Bug reports and feature requests are welcome on the [issue tracker](https://github.com/bbuchsbaum/fmrigds/issues).

## Citation

If you use fmrigds in your research, please cite:

```
Buchsbaum, B. R. (2025). fmrigds: Format-agnostic group-level analysis for fMRI.
R package version 0.1.0. https://github.com/bbuchsbaum/fmrigds
```

## License

GPL-3
