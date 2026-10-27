# Issues 36-39: initial local repair evidence

These pre-rebase checks are historical. See integration.md for the current
upstream-integrated validation and remaining limitations.

Base commit: `d58e3581609aec2b0fb2c482a0ddb35006e50f5d`.
Final source and test hashes: `source-sha256.txt`. Changes were uncommitted when these initial checks ran.
R 4.5.1, macOS arm64, system Accelerate BLAS; full sessionInfo in reproduction logs.

All four original reports reproduced against a fresh native build in a temporary
source copy. The checkout's tracked native binaries were not modified.

- #36: Nonzero map support governs missingness for effects and independent
  variances, including sparse inputs. Missing contributors remain missing.
  Dense fully observed inputs retain matrix multiplication. All-zero independent
  map rows retain zero effect/variance. Satterthwaite df propagates missingness
  instead of throwing for missing contributing variance/df.
- #37: Default voxel compatibility checks the affine exactly, alongside existing
  dimensions and optional identity fields. Explicit fields can request weaker
  checks. This is not evidence of anatomical registration or complete ordered
  support equivalence.
- #38: Validate all resulting contrast IDs before assignment; reject collisions
  and malformed maps. Align and rename contrast metadata and assay dimnames,
  including simultaneous swaps. Distinct estimates are never aggregated.
- #39: Default eps is zero in R and C++; positive finite variances are not
  floored by default. Explicit positive eps retains regularization and warns
  when applied. Invalid variances remain excluded. Remove dimensionful floors
  from DL denominators, and correct the regression denominator to
  tr(W) - tr(W X (X' W X)^-1 X' W). Intercept-only regression agrees with
  ordinary meta-analysis. Scaling factors 1e-8, 1, and 1e8 preserve inference;
  tests include heterogeneous weights and analytic equal-variance and DL cases.
  These checks do not establish stability at all floating-point extremes.

## Checks

- Fresh native baseline reproduction: confirms all four reports.
- Focused tests: 368 passed, 0 failures/warnings/skips.
- Full testthat suite with RGL_USE_NULL=true and options(rgl.useNULL=TRUE):
  1707 passed, 0 failures/warnings/skips.
- A final diagnostic-formatting edit in compat.R was followed by the affected
  compatibility and issue regression tests: 258 passed, 0 failures/warnings/skips.
- git diff --check passed.
- Roxygen Rd files regenerated. Existing unresolved gds_plan/plan topic links
  were reported by load_source; only four intended Rd files changed.

The first full-suite run aborted with SIGABRT in rgl/macOS OpenGL context
creation. Loading fmristore in a separate process without fmrigds reproduced
that abort. The headless suite ran every test without skips. Arrow emitted
sandbox CPU-info diagnostics; testthat emitted a package-build-version startup
warning. Neither produced a test warning or failure.

Full R CMD check was not run: 17 GiB free was below the inherited workspace
20 GiB precheck threshold. No hosted CI, publication, or deployment was run.
Existing unrelated untracked UPSTREAM-ISSUES-from-fmrireg.md was preserved.

Reproduction: run scripts/repro-issues-36-39.R with the path to a temporary source
copy whose src/*.o and src/*.so have been omitted; pkgload builds fresh native code.
Run the suite there with RGL_USE_NULL=true and options(rgl.useNULL=TRUE).
Raw logs and exit-status/timing sidecars are retained beside this note.
