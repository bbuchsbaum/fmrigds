# Sprint 8 Plan — Performance, Exporters & Release Prep (of 8)

**Sprint:** 8 of 8  
**Duration:** 2–3 weeks  
**Prerequisites:** Sprints 1–7 complete (core plan verbs, adapters, HDF5 I/O, map persistence).

---

## Objectives

1. Hardening & performance: profile plan execution, improve block streaming, and tighten memory use.
2. Export surface polish: finish remaining exporters (`write_out()` formats, fmristore shim if deferred) and add summary utilities.
3. Package readiness: documentation pass (pkgdown, README, vignettes), finalize tests, and prep for initial release tag.

---

## Scope & Tasks

### 1. Performance & Robustness

**Files:** `R/compute.R`, `R/plan-optimizer.R`, `R/h5-write.R`, `tests/testthat/test-optimizer.R`

- Implement additional plan optimizations (subset pushdown, derive fusion, map/reduce coalescing) in `optimize_plan()`.
- Benchmark block sizes for HDF5 reads/writes; expose tuning knobs via `compute(block = ...)` and document suggested defaults.
- Add guardrails for large plans: sanity-check memory, warn on missing map families, improve error messaging.
- Tests: ensure optimizer rewrites expected node orders; add regression tests for block-size handling.

### 2. Exporters & fmristore Integration

**Files:** `R/verb-write.R`, `R/write-exporters.R` (new), `R/adapter-fmristore.R` (if still pending), `inst/schemas/hdf5.md`

- Add CSV/Parquet exporters for parcel/basis spaces (`write_out(..., format = "csv"|"parquet")`).
- Provide a minimal NIfTI writer for voxel spaces (using RNifti or existing utilities).
- Finalise fmristore adapter: detection, space mapping, and fallback to `/gds` layout when present.
- Update schema docs to reference exporter expectations and fmristore interop.
- Tests: round-trip exporters on toy data; fmristore adapter smoke tests.

### 3. Documentation & Release Prep

**Files:** `README.md` (new), `NEWS.md`, `pkgdown/`, `vignettes/`

- Write README and pkgdown home page summarising GDS concepts and quick-start examples.
- Add vignette covering lazy plans + HDF5 write/read workflow.
- Ensure roxygen docs complete and run `devtools::document()`.
- Prepare `DESCRIPTION`/`NEWS.md` for 0.1.0 release; confirm license headers.
- CI: ensure `devtools::check()` clean (R CMD check, lintr optional).

---

## Deliverables

- ✅ Optimised executor with documented performance guidelines.
- ✅ `write_out()` supports CSV/Parquet/NIfTI; fmristore adapter available (or explicitly scoped out).
- ✅ Documentation suite ready (README, vignette, pkgdown config, NEWS).
- ✅ Package passes `R CMD check` and is tag-ready.

---

## Out-of-Scope

- TileDB/Arrow backends (post-0.1 roadmap).
- Advanced provenance viewers or GUI tools.
- Distributed execution engines.

---

## Acceptance Criteria

1. Benchmarked plan runs show measurable improvement (documented in sprint notes).
2. Exporters tested and produce correct output for representative spaces.
3. fmristore adapter (if included) successfully opens example stores.
4. Documentation builds without warnings; `devtools::check()` passes.

---

## After Sprint 8

- Plan post-release roadmap (TileDB backend, richer provenance UI, OT kernels).
- Gather user feedback and prioritise feature backlog for 0.2 cycle.

