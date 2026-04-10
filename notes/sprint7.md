# Sprint 7 Plan — Output, fmristore, and Alignment Persistence (of 8)

**Sprint:** 7 of 8  
**Duration:** 3 weeks  
**Prerequisites:** Sprints 1–6 complete (core plan verbs, adapters, HDF5 backend, serialization, provenance).

---

## Objectives

1. Implement write-side support so `compute(..., sink = "h5")` and `write_out()` can materialise results in the `gds-h5` layout (including provenance append).
2. Introduce an fmristore adapter that reads existing fmristore files and prefers `/gds` groups when present.
3. Persist alignment/map families to disk and expose loader helpers so downstream plans can reuse registered transforms.

---

## Scope & Tasks

### 1. HDF5 Write-Out & Provenance

**Files:** `R/compute.R`, `R/h5-write.R`, `R/adapter-h5.R`, `inst/schemas/hdf5.md`

- Extend `compute()` with `sink = c("memory","h5")` (accept path + options).
- Implement block writer (`.h5_write_dataset`) mirroring reader semantics.
- Append provenance log entries to `/gds/provenance` after compute.
- Update schema doc with write expectations and provenance layout.
- Tests: round-trip `compute(..., sink="h5")`, verify provenance append, ensure partial-block write paths.

### 2. fmristore Adapter

**Files:** `R/adapter-fmristore.R` (new), `NAMESPACE`

- Implement detect/probe/read via `fmristore::detect_h5_type()`; fall back to `gds-h5` when `/gds` present.
- Map voxel/parcellated/latent fmristore structures to `space_*` descriptors and LinearMaps.
- Register adapter in `.onLoad`; add tests covering voxel + parcel ingestion, ensuring compatibility with plan verbs.

### 3. Alignment Persistence Helpers

**Files:** `R/verb-align.R`, `R/align-exec.R`, `R/h5-write.R`, `tests/testthat/test-align.R`

- Extend `register_map()`/`align()` to serialize map families under `/gds/alignments/<name>`.
- Provide `list_maps()`/`load_map_family()` helpers (used by adapter probe).
- Tests: register + persist orthogonal matrices, reload via plan constructed from file, ensure compute uses stored transforms.

### 4. Documentation & Release Notes

- Update `NEWS.md` (Sprint 7 section) and add schema/persistence notes.
- Roxygen docs for new arguments (`compute`, `write_out`, fmristore adapter).
- (Optional) Brief README paragraph describing gds-h5 writer & fmristore integration.

---

## Deliverables

- ✅ `compute()` and `write_out()` can create `gds-h5` stores with provenance logs.
- ✅ fmristore adapter supports voxel/parcel/latent inputs (prefers `/gds` layout).
- ✅ Alignment families persist and reload from disk.
- ✅ Updated schema docs, NEWS, and reference documentation.

---

## Out-of-Scope (Sprint 8)

- TileDB/Arrow backends.
- Streaming writes to alternative sinks (e.g., Parquet).
- Performance profiling & pkgdown build (targeted for Sprint 8).

---

## Acceptance Criteria

1. `compute(plan, sink="h5", path=tmp)` produces a readable `gds-h5` file with provenance appended.
2. `gds(fmristore_file)` returns a plan with correct space/assay metadata for voxel/parcel/latent variants.
3. A plan that registers and persists map families can be saved, reloaded, and computed on a new session.
4. All tests (including new adapter + persistence suites) pass.

---

## Next

- Sprint 8: performance tuning, fmristore exporter polish, pkgdown site, release packaging.

