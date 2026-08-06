# fmrigds 0.1.0.9000

## Documentation

- Standardized every vignette and pkgdown article on the red `interaction`
  Albers theme, restored theme activation in the fmristore article, and added
  a post-build check that prevents article theme drift.

## Model-conditioned group examination

- Added `examine_group()` as a terminal branch from subject-level plans. It
  reports data validity, cross-fitted model surprise, deletion influence on
  named estimands, and model-adjusted residual geometry without changing the
  source analysis.
- Added exact fixed-effect, meta-regression, and voxelwise OLS diagnostics,
  plus fixed-`tau2` random-effects screening and exact retained-subject refits.
- Added deterministic block scanning and staging, action-plane and drill-down
  plots, selected post-hoc sensitivity analyses, and self-contained HTML reports.
- Corrected the DerSimonian-Laird meta-regression denominator to use
  `tr(W - W X (X' W X)^-1 X' W)` and aligned R/C++ variance-validity rules.
- Hardened CI and pkgdown provisioning for the GitHub-only `neuroim2`,
  `neurotabs`, and `neurothresh` development dependencies.

## Bug fixes (reported issues #1, #5, #6, #7)
- `write_nifti_assays()` and `write_out(format = "nifti")` now preserve the
  input spatial affine (spacing, origin, orientation) and write `scl_slope = 1`
  for 3D outputs, instead of emitting an identity NeuroSpace. Exported group
  maps are now geometrically correct without a manual re-stamp (#6).
- `ols:voxelwise` (and therefore `one_sample()`/`group_ols()`) now applies
  per-sample listwise deletion for non-finite subjects: a single `NaN` subject
  no longer poisons a voxel and an all-`NaN` subject no longer poisons the whole
  map. Samples with too few finite observations are returned as `NA`, a new
  `n_obs` assay reports the effective sample size, and a summary warning is
  emitted when any sample is reduced (#7).
- `one_sample()` / `group_ols(~ 1)` now build the intercept-only design
  automatically and no longer require `col_data` (#1).
- Variance-weighted reducers (`fixed`/`random`/`meta:*`, or any `weights =
  "1/var"` reduction) now error with an actionable message when applied to a
  beta-only GDS whose `var` assay is the synthetic unit-variance placeholder,
  rather than silently producing meaningless group standard errors (#5).
- Documented reducer output assay names (`reduce()` and the `ols:voxelwise`
  topic), the single-contrast scope of `gds_from_nifti_maps()`, and the
  filename pairing contract for `nifti_source(beta=, se=)` (#2, #3, #4).

## NIfTI raw-map ingestion
- Beta-only NIfTI sources now materialise as realised GDS objects by adding a synthetic `var = 1` assay with a warning, matching `gds_from_neurovols()` behavior for raw maps without uncertainty images.
- Documented beta-only `gds_from_nifti_maps()` workflows and the recommendation to use `ols:voxelwise` rather than fixed/random-effects meta-analysis when the variance assay is synthetic.
- Tightened `gds_from_nifti_maps()` contrast relabeling for one-contrast raw-map layouts.
- Added scalar-map workflow helpers: `gds_from_scalar_maps()` / `as_scalar_map_gds()`, `group_ols()`, `one_sample()`, `two_sample()`, and `write_nifti_assays()` with a per-file manifest.

## CLI overhaul
- Expanded the bundled `fmrigds` CLI into a plan-oriented interface with `probe`, `plan`, `preview`, `run`, and `list` commands.
- Added saved-plan workflows via `--save-plan` / `--load-plan`, plus repeatable passthrough flags for adapter, reducer, post-hoc, and writer options.
- Added CLI support for attaching subject/sample/contrast metadata tables, mask policies, preview tables, and registry inspection.
- Added a new `cli-workflows` vignette covering the main command patterns and advanced control flags.
- Added focused CLI tests covering plan construction, JSON output, raw previews, and reducer registry listing.

## Sprint 10 / Repeated-measures LMM polish
- Added restricted repeated-measures Gaussian LMM reducers: `lmm:ri` and `lmm:ri_slope1`.
- Added variance-aware `lmm:ri_knownvar` and `lmm:ri_slope1_knownvar`
  reducers. They fit `diag(var) + vc_resid * I + Z G Z'`, reject synthetic or
  nonpositive variances, and are numerically checked against
  `metafor::rma.mv()`.
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
