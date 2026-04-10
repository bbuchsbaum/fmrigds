# Sprint 6 Plan — HDF5 Backend & Plan Serialization (of 8)

**Sprint:** 6 of 8  
**Duration:** 3 weeks  
**Prerequisites:** Sprints 1–5 complete (alignment, masking, reducers, adapters, derivations).

---

## Objectives

1. Integrate an HDF5 backend (gds-h5) enabling persistence and lazy loading.
2. Finalize plan serialization (`save_plan`, `load_plan`) with full node coverage and provenance digests.
3. Extend provenance logging to include reducer/mask/align operations with file-backed metadata.

---

## Scope & Tasks

### 1. HDF5 Adapter (`gds-h5`)

**Files:** `R/adapter-h5.R` (new), `inst/schemas/hdf5.md`

- Implement adapter functions:
  * `detect`: check for `/gds` group (`gds-h5` mini spec).
  * `open`: return handle with `hdf5r::H5File`.
  * `probe`: read axes (`/gds/axes/subjects`, `/gds/axes/contrasts`), space info, assays.
  * `read`: support block reads along sample axis using chunked datasets.
  * `close`: close file handle.
- Add writer helpers (`write_gds_h5()`, `append_provenance_h5()`).
- Update `register_builtin_adapters()` to include h5.
- Tests: generate temporary HDF5 file, round-trip via `gds()`/`compute()`, ensure subset/block reads.

### 2. Plan Serialization Enhancements

**Files:** `R/plan-serialization.R`

- Extend serialization to handle `align_to_group`, `mask_policy`, `map` with `map_linear` descriptors, `reduce` nodes.
- Include provenance digest in serialized JSON.
- Implement `load_plan()` to reconstruct nodes (falling back to placeholder for unsupported custom functions).
- Tests: plan ↔ JSON round-trip with all node types.

### 3. Provenance & Metadata

**Files:** `R/gds-class.R`, `R/compute.R`

- Extend provenance logs to include `align`, `mask`, `map`, `reduce`, `write_out` nodes with parameters.
- When `compute()` writes to HDF5 (future Sprint 7), ensure provenance append is ready.
- Tests: validate provenance graph entries after compute with mixed nodes.

### 4. Documentation & Tooling

- Update `DESCRIPTION` (Imports: hdf5r).
- Add `inst/schemas/hdf5.md` describing layout.
- NEWS entry summarizing HDF5 support.
- roxygen docs for new adapter and serialization helpers.
- Ensure new files covered by tests.

---

## Deliverables

- ✅ HDF5 adapter with read support and tests.
- ✅ Serialization covering all plan nodes and provenance digest in JSON.
- ✅ Provenance logging expanded to new operations.
- ✅ Updated docs/NEWS.

---

## Out-of-scope (later sprints)

- HDF5 write support integrated with `write_out()` (Sprint 7).
- Full provenance append to HDF5 (Sprint 7).
- TileDB/Arrow backends (future).

---

## Acceptance Criteria

1. `gds()` can open HDF5 files and `compute()` realizes them with subsetting.
2. `save_plan`/`load_plan` round-trip plans containing subset/align/mask/map/reduce nodes.
3. Provenance graph includes all verb parameters.
4. Tests cover HDF5 adapter and serialization features.
5. CI passes.

---

## Next

- Sprint 7: fmristore adapter integration, exporters (`write_out()`), factorial utilities.
- Sprint 8: Performance tuning, pkgdown, release prep.

