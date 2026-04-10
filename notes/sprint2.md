# Sprint 2 Plan — Storage Adapters & Execution Scaffold (of 8)

**Sprint:** 2 of 8
**Duration:** 3 weeks
**Prerequisites:** Sprint 1 complete (core foundations)
**Target:** ≥95% test coverage for new modules (adapters, compute scaffold)

## Objectives
- Introduce the adapter layer so `gds()` can open real data sources (tabular + NIfTI).
- Provide an execution scaffold (`compute()`/`execute_stream()`) that walks plans and streams blocks with functional block iteration.
- Ensure adapter detection/registration is predictable and tested.
- Deliver end-to-end pipeline: source → plan → compute → realized GDS

## Scope & Tasks

### 1. Adapter Registry & Front Door
- Implement the adapter registry (`register_adapter`, `get_adapter`, `detect_adapter`) as defined in TECHNICAL_SPECIFICATION.md §7.1–7.2.
- Implement `gds()` front door: detect adapter, open/probe source, return `gds_plan` with populated `gds_source`.
- Provide `close` handling (on plan finalization or via `compute()` context).

### 2. Tabular Adapter (CSV/Parquet)
- Implement adapter functions per TECHNICAL_SPECIFICATION.md §7.3:
  - `detect`: file extension (`.csv`, `.tsv`, `.parquet`) + lightweight file existence validation.
  - `open`: read data frame (CSV via `data.table::fread`, Parquet via `arrow::read_parquet`).
  - `probe`: infer axes (samples/subjects/contrasts), build `space_parcels`, produce metadata (`schema_version = "0.1.0"`).
    - Auto-detect columns or accept explicit mapping via `subject_col`, `sample_col`, `contrast_col`, `effect_cols` arguments.
  - `read`: reshape from long format into `[sample × subject × contrast]` arrays.
    - Support `block` argument (sample window) for chunked reading.
    - For Sprint 2, block size can be generous (full dataset or 10K samples).
- Add tests using temporary CSV with small synthetic dataset (multiple subjects, ROIs, contrasts).
- Document assumptions: expects long format (one row per sample×subject×contrast) per spec.

### 3. NIfTI Adapter (read-only MVP)
- Implement adapter skeleton using `RNifti` per TECHNICAL_SPECIFICATION.md §7.4:
  - `detect`: handles single `.nii`/`.nii.gz` files, vectors of files, or directories.
  - `open`: store file list/layout metadata (BIDS-like or flat structure).
  - `probe`:
    - Determine voxel dims and affine from first file header.
    - Extract subjects/contrasts from filenames (BIDS pattern or sequential numbering).
    - Compute mask (union over first 3-5 files, or load if provided).
    - Return `space_voxel` with `storage = "packed"` and `mask_idx`.
  - `read`:
    - Apply mask and return arrays in `[sample × subject × contrast]` format.
    - Support `block` argument (sample window) to enable chunked reading.
    - For Sprint 2, implement basic blocking (read full volumes, then slice by sample indices).
- Tests with synthetic NIfTI (create via `RNifti::writeNifti` to temp files):
  - Verify probe extracts correct dims and mask length.
  - Verify read applies mask correctly.
  - Test block reading with sample indices.

### 4. Execution Scaffold
- Implement `compute()` per TECHNICAL_SPECIFICATION.md §3.9 and §5:
  - Call `optimize_plan()` (still identity for Sprint 2; real implementation in Sprint 4).
  - Iterate blocks via `execute_stream()` with functional block iteration on sample axis.
  - Return realized `gds` using existing constructors.
- Implement `execute_stream()` structure (§5.4):
  - Open source handle via adapter.
  - Iterate sample blocks (configurable block size, default 100K samples).
  - For each block:
    - Read assays from adapter.
    - Apply operations sequentially (see operation status below).
    - Write block to sink (memory sink functional; HDF5 sink stubbed).
  - Finalize sink and construct realized GDS.
- **Operation execution status for Sprint 2:**
  - ✅ `subset()`: **Functional** (filter samples/subjects/contrasts)
  - ✅ `derive()`: **Basic implementation** (var ↔ se only; t/z deferred to Sprint 3)
  - ❌ `align()`: Pass-through with TODO marker (Sprint 4)
  - ❌ `mask()`: Pass-through with TODO marker (Sprint 4)
  - ❌ `map_to()`: Pass-through with TODO marker (Sprint 3)
  - ❌ `reduce()`: Pass-through with TODO marker (Sprint 5)
  - ❌ `write_out()`: Pass-through with TODO marker (Sprint 7)
- Add provenance recording per §4.4:
  - `plan_digest` using `digest_plan()`
  - `computed_at` timestamp
  - Append to `metadata$provenance$log`

### 5. Testing & Tooling
- Expand `testthat` suite with ≥95% coverage target:
  - **Adapter registry tests** (`test-adapter-registry.R`):
    - Registration, retrieval, detection scoring.
    - Probe result validation.
  - **Tabular adapter tests** (`test-adapter-tabular.R`):
    - CSV/Parquet detection.
    - Probe metadata extraction with auto-detect and explicit columns.
    - Read reshaping to 3D arrays.
    - Block reading correctness.
  - **NIfTI adapter tests** (`test-adapter-nifti.R`):
    - File/directory detection.
    - Probe with mask creation.
    - Read with masking applied.
    - Block reading with sample indices.
  - **Compute tests** (`test-compute.R`):
    - Basic compute() on tabular adapter.
    - Basic compute() on NIfTI adapter.
    - Derive execution (var ↔ se).
    - Subset execution (dimension reduction).
    - Provenance recording.
  - **Integration tests** (`test-integration.R`):
    - End-to-end: CSV → plan → compute → realized GDS.
    - End-to-end: NIfTI → plan → compute → realized GDS.
    - Multi-step pipeline: subset → derive → compute.
- Update GitHub Actions workflow to install optional packages (`arrow`, `data.table`, `RNifti`).
- Run `devtools::document()` and ensure NAMESPACE reflects new exports.
- Verify `R CMD check` passes with zero errors, warnings, or notes.

## Deliverables

**Core Functionality:**
- ✅ Functional `gds()` front door returning plans bound to adapters
- ✅ Working tabular adapter (CSV/Parquet) with auto-detection and explicit column mapping
- ✅ Working NIfTI adapter with mask creation and packed storage
- ✅ `compute()` pipeline with block streaming and memory sink
- ✅ Basic operation execution: `subset()` and `derive()` (var ↔ se)
- ✅ Provenance recording (digest, timestamp, log)

**Testing & Quality:**
- ✅ ≥95% test coverage for adapters and compute modules
- ✅ Integration tests demonstrating end-to-end workflows
- ✅ CI pipeline passing with new dependencies (`arrow`, `data.table`, `RNifti`)
- ✅ `R CMD check` clean (zero errors/warnings/notes)

**Documentation:**
- ✅ Roxygen2 documentation for all exported functions
- ✅ Updated NAMESPACE
- ✅ Examples in function documentation
- ✅ NEWS.md entry for Sprint 2

## Technical Debt & Deferred Work

**Deferred to Sprint 3 (Statistical Operations):**
- Full `derive()` execution (t from beta/var, z from t+df, F-statistics)
- `map_to()` implementation with variance propagation
- Statistical correctness enforcement (refuse to map t/z without effect scale)

**Deferred to Sprint 4 (Spatial Operations):**
- `align()` execution with MapFamily
- `mask()` execution with MaskPolicy
- Complete optimizer implementation (pushdown, coalescing, fusion)
- `explain()` function for plan visualization

**Deferred to Sprint 5 (Meta-Analysis):**
- `reduce()` execution (fixed-effects, Stouffer, Fisher)
- df aggregation (Satterthwaite)

**Deferred to Sprint 6 (Persistence):**
- HDF5 sink (currently stubbed; writes to memory)
- Plan serialization (save_plan/load_plan)

**Deferred to Sprint 7 (Integration):**
- `write_out()` exporters (NIfTI, CSV, Parquet)
- fmristore adapter integration

## Dependencies & Notes

**Technical Specification References:**
- Adapter interface: §7.1–7.2
- Tabular adapter: §7.3
- NIfTI adapter: §7.4
- Compute execution: §3.9, §5.4
- Provenance: §4.4

**Package Dependencies:**
- **Required imports:** Matrix, digest, jsonlite
- **Suggested packages:** data.table (CSV), arrow (Parquet), RNifti (NIfTI)

**Implementation Notes:**
- Keep adapters **read-only**; write support deferred to Sprint 7
- Optimizer is **pass-through** (identity); real implementation in Sprint 4
- HDF5 sink is **stubbed**; falls back to memory sink
- Focus on **correctness and test coverage** over performance optimization

## Acceptance Criteria

**Sprint 2 is complete when:**
1. ✅ User can load CSV data via `gds()` and get a plan
2. ✅ User can load NIfTI data via `gds()` and get a plan
3. ✅ `compute(plan)` materializes data to a realized GDS
4. ✅ `subset()` correctly filters dimensions
5. ✅ `derive(c("var", "se"))` correctly converts between var and se
6. ✅ Block streaming handles datasets with >10K samples without memory issues
7. ✅ Provenance is recorded (digest, timestamp, log)
8. ✅ All unit tests pass (≥95% coverage for new modules)
9. ✅ Integration tests demonstrate realistic workflows
10. ✅ CI pipeline passes cleanly
11. ✅ Documentation is complete and accurate

## Example Usage (End of Sprint 2)

```r
library(gdsfmri)

# Example 1: Tabular (CSV) workflow
plan <- gds("roi_stats.csv") %>%
  subset(subject = c("sub-01", "sub-02", "sub-03")) %>%
  derive(c("var", "se"))

result <- compute(plan, sink = "memory")
print(result)
# GDS object
# ==========
# Dimensions: 5 samples × 3 subjects × 2 contrasts
# Space: parcels
# Assays: beta, var, se

# Example 2: NIfTI workflow
plan <- gds(c("sub-01/beta.nii.gz", "sub-02/beta.nii.gz")) %>%
  subset(contrast = "faces_vs_places")

result <- compute(plan, sink = "memory", block = list(sample = 10000))
print(result)
# GDS object
# ==========
# Dimensions: 45231 samples × 2 subjects × 1 contrast
# Space: voxel (packed)
# Assays: beta
```

## Upcoming Sprints (3–8)

**Sprint 3:** Statistical Operations & Derivations (3 weeks)
- Complete `derive()` execution (all derivation rules)
- Variance propagation (independent + covariance modes)
- `map_to()` implementation with uncertainty propagation

**Sprint 4:** Spatial Operations (3 weeks)
- `align()` execution with MapFamily
- `mask()` execution with MaskPolicy
- Complete optimizer (pushdown, coalescing, fusion)
- `explain()` function

**Sprint 5:** Meta-Analysis & Reduction (2 weeks)
- `reduce()` execution (fixed-effects, Stouffer, Fisher)
- df aggregation (Satterthwaite)

**Sprint 6:** Persistence & Serialization (3 weeks)
- HDF5 adapter with GDS layout
- Plan serialization (save_plan/load_plan)
- Provenance persistence

**Sprint 7:** Integration & Export (3 weeks)
- fmristore adapter integration
- `write_out()` exporters (NIfTI, CSV, Parquet)
- Factorial contrast utilities

**Sprint 8:** Production Ready & Release (3 weeks)
- Comprehensive documentation and vignettes
- Performance optimization
- v0.1.0 release preparation

---

**Sprint 2 Version:** 2.0 (Enhanced)
**Last Updated:** 2025-01-XX
**Status:** Ready for implementation

