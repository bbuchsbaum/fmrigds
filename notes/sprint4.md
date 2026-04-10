# Sprint 4 Plan — Alignment, Masking, and Optimizer (of 8)

**Sprint:** 4 of 8  
**Duration:** 3 weeks  
**Prerequisites:** Sprints 1–3 complete (core GDS, adapters, derivations, variance propagation, map_to).

---

## Objectives

1. Execute subject→group alignment via `MapFamily` (orthogonal/linear, deformable, OT).
2. Implement mask policies (`mask()` verb) with storage repacking.
3. Introduce optimizer rewrites (subset pushdown, derive coalescing, early align/mask) and integrate with `compute()`.
4. Expand plan serialization/provenance to capture alignment/mask operations.

---

## Scope & Tasks

### 1. Alignment Execution (`align()`)

**Files:** `R/verb-align.R` (new), `R/align-exec.R` (new), `R/compute.R`

- Implement `align()` verb returning plan nodes.
- Runtime implementation:
  * Linear/orthogonal/OT map families via per-subject operators (reuse variance propagation from Sprint 3).
  * Support both dense matrices and functional descriptors (via callbacks).
  * Record new space (group space) in plan metadata.
- Integrate with `compute()` pipeline (apply align before mask/map).
- Tests: synthetic per-subject rotations, transport plans verifying β/var results and df aggregation.

### 2. Mask Policies (`mask()`)

**Files:** `R/verb-mask.R` (new), `R/mask-exec.R` (new), `R/space.R`

- Implement `MaskPolicy` object and `mask()` verb (group/subject-level, intersection/union/threshold/custom).
- Execution: derive analysis mask, repack packed/dense storage, update space metadata, ensure map/mask interplay.
- Tests: per-subject voxel masks, group threshold; ensure repacking yields expected sample counts.

### 3. Optimizer Enhancements

**File:** `R/plan-optimizer.R` (rename from placeholder), `R/compute.R`

- Implement rewrite rules:
  1. Push `subset` before `align`, `mask`, `derive`, `map`.
  2. Coalesce chained `derive` nodes.
  3. Reorder `align` before `mask` and `map`.
  4. Fuse consecutive `mask` operations.
- Update `compute()` to call new optimizer and rely on output order.
- Tests: plan-level assertions verifying rewrite invariants.

### 4. Provenance & Serialization

- Extend provenance nodes to include `align` and `mask` parameters.
- Update plan digest/canonicalization to include optimizer results.
- Serialization stub (JSON) for plan capture (introduce `save_plan`, `load_plan` returning plan objects).

### 5. Documentation & Tests

- Add roxygen for new verbs and helpers.
- Expand DESCRIPTION Imports if needed (e.g., for matrix ops).
- Tests:
  * Alignment unit/integration tests.
  * Mask policy unit/integration tests.
  * Optimizer rewrite tests.
  * Plan serialization round-trip test.

---

## Deliverables

- ✅ `align()` & execution covering orthogonal/linear/OT families.
- ✅ `mask()` verb with policies and storage repacking.
- ✅ Optimizer rewrite engine with core transformations.
- ✅ Provenance/digest updates for align/mask.
- ✅ `save_plan`/`load_plan` stubs for serialization.
- ✅ Comprehensive tests covering new features.

---

## Out-of-Scope / Deferred

- Kernel-based covariance propagation (Sprint 5).
- Random-effects meta-analysis (`reduce_random`) (Sprint 6).
- Advanced regridding (surface ↔ voxel) (Sprint 6/7).

---

## Acceptance Criteria

1. `align()` supports orthogonal, linear, and OT families with correct β/var/df outputs.
2. `mask()` policies compute analysis masks and update storage without breaking adapters.
3. Optimizer rewrites produce canonical plan order; compute relies on optimized plan.
4. Tests demonstrate alignment+mask pipeline (subset → align → mask → map → derive → compute) end-to-end.
5. Plan digest/serialization captures new nodes.
6. CI passes with tests and documentation updated.

---

## After Sprint 4

- Sprint 5: Reducers (`reduce()` fixed/random), Stouffer/Fisher integration, df management.
- Sprint 6: HDF5 backend (gds-h5), plan serialization to disk, provenance append.
- Sprint 7: fmristore adapter integration, exporters (`write_out()`), factorial utilities.
- Sprint 8: Performance tuning, pkgdown, release prep.

