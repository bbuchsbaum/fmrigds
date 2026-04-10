# fmrigds 0.1.0.9000

## CLI overhaul
- Expanded the bundled `fmrigds` CLI into a plan-oriented interface with `probe`, `plan`, `preview`, `run`, and `list` commands.
- Added saved-plan workflows via `--save-plan` / `--load-plan`, plus repeatable passthrough flags for adapter, reducer, post-hoc, and writer options.
- Added CLI support for attaching subject/sample/contrast metadata tables, mask policies, preview tables, and registry inspection.
- Added a new `cli-workflows` vignette covering the main command patterns and advanced control flags.
- Added focused CLI tests covering plan construction, JSON output, raw previews, and reducer registry listing.

## Sprint 10 / Repeated-measures LMM polish
- Added restricted repeated-measures Gaussian LMM reducers: `lmm:ri` and `lmm:ri_slope1`.
- Added `theta_mode = "pooled"` and `theta_mode = "voxelwise"` for shared- vs sample-specific variance parameters.
- Added long-table repeated-measures ingestion via `contrast_data_cols`, so contrast metadata can be extracted directly from tabular sources.
- Hardened the low-rank LMM kernels against numerical asymmetry in intermediate SPD matrices.
- Expanded mixed-model tests to cover `lme4` parity, pooled/voxelwise equivalence for identical responses, and sample-specific theta recovery.
- Added a dedicated repeated-measures vignette and reducer-specific documentation for the supported LMM contract.
- Added CI coverage for `lme4`-backed restricted LMM parity tests plus a scheduled benchmark workflow.
- Warm-started voxelwise theta optimization from pooled and previous-sample fits to reduce repeated optimizer cold starts.

## Release-track cleanup
- Aligned package metadata, docs, and provenance fields around the `fmrigds` name.
- Replaced scaffold-era package documentation with a user-facing package overview.
- Added console print methods for realised GDS objects and plans.
- Tightened README examples and installation guidance to reflect the current package surface.

## Sprint 9 / GDS-first integration
- Added in-memory tabular support: `gds(<data.frame>)` now works via the tabular adapter (no temporary files needed).
- Added `roi_col` alias in tabular ingestion to treat ROI/parcel labels as the sample axis.
- NIfTI source normalization: adapters now accept `list(beta=..., se=...)`, directories, or vectors; classification by filename patterns for beta vs se.
- Introduced an in-memory adapter and `as_plan.gds()` so verbs like `reduce()`/`posthoc()` accept realized GDS directly.
- Added ergonomic alias `plan(x)` for `as_plan(x)`.
- Relaxed `as_gds()` to accept 2-D inputs and standardized all assays to 3-D `[sample × subject × contrast]` shapes.
- Auto-derivation of `var` from `se` (and vice versa when required) in reducer/posthoc preflight.
- Implemented spatial FDR (`fdr:spatial`) posthoc with group-wise Simes + weighted BH; documented usage and options.
- Promoted OLS covariance triangles to assays (`cov:<term_i>:<term_j>`) while retaining attachments for metadata.
- Added S3 accessors for plans: `subjects.gds_plan()` and `contrasts.gds_plan()`.
- New helper `nifti_source(beta=NULL, se=NULL)` for explicit NIfTI list sources.

Guidance for integrators (fmrireg):
- Prefer `gds(df, format="tabular", roi_col="roi")` or `as_gds(df, mapping=list(roi="roi", ...))` over temporary CSV writes.
- Use `plan(gds_obj)` or `as_plan(gds_obj)` to chain `reduce()`/`posthoc()` on realized results.
- Expect consistent 3-D assay shapes; 2-D inputs are upcast automatically.

## Sprint 6
- Added an HDF5 storage adapter (`gds-h5`) with probe and block-read support so plans can lazily stream assays from disk.
- Documented the HDF5 layout under `inst/schemas/hdf5.md` and wired helpers for writing simple gds stores.
- Expanded tests to cover the new adapter and serialization pipeline.

## Sprint 7
- `compute()` now supports `sink = "h5"` and `write_out()` queues, emitting provenance-aware `gds-h5` stores (with tests).
- Map families registered on plans/GDS are serialised under `/gds/alignments`, rehydrated at load time, and available via `align(plan, "<name>")`.
