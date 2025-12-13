# fmrigds (development)

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
