# Sprint 1 Plan — gdsfmri Implementation (of 8)

## Objectives
- Establish concrete foundations from the technical specification for the core nouns (GDS, Space, Map, UncertaintyRule) and the lazy plan infrastructure.
- Deliver a runnable package with meaningful functionality beyond the current stubs, including validation and basic derivations.
- Set up automated testing and CI hooks to support subsequent sprints.

## Scope & Tasks

### 1. Core Data Structures
- Implement `new_gds()` with full validation (dims, assay requirements, metadata defaults) — see §2.1 & §2.6.
- Add accessor functions (`assay`, `assays`, `space`, `subjects`, `contrasts`, `col_data`, `row_data`, `metadata`).
- Implement `gds_metadata()` (spec §4.2) and provenance helpers (`provenance_node`, `add_provenance_node`).

### 2. Space Constructors
- Implement `space_voxel()`/`space_voxels()` with mask handling and storage modes (§2.2.1).
- Implement `space_parcels()`, `space_surface()`, `space_basis()` including projector conventions (§§2.2.2–2.2.4).
- Add basic validation utilities (`validate_space_*`) to confirm sample counts vs. assays.

### 3. Map & Uncertainty Primitives
- Implement `map_linear()` with validation of operator dimensions and traits (§2.3.1).
- Implement `MapFamily()` plus helper constructors (`OrthogonalFamily`, `OTFamily`, `WarpFamily`) (§2.3.2).
- Implement `UncertaintyRule()` with argument checks (§2.4).

### 4. Plan Infrastructure
- Implement `gds_plan()`, `gds_source()`, operation node constructors (`op_subset_axis`, `op_derive`, `op_align_to_group`, `op_mask_policy`, `op_map`, `op_reduce`, `op_write`) (§§2.5, 8.4).
- Implement `add_op()`, `as_plan()` helper, and `subset()` verb (lazy only) with basic validation.
- Stub optimizer entry point `optimize_plan()` returning the plan unchanged (detailed rewrites arrive later).

### 5. Derivation & Registry Foundations
- Implement assay registry (`register_assay`, `assay_info`, `can_map_linear`) with default registrations from §4.1.
- Implement `derive()` verb that records the operation (execution deferred).

### 6. Package & Tooling Setup
- Replace placeholder stubs with real implementations; ensure exported functions use roxygen2 tags.
- Configure `devtools::document()` workflow (DESCRIPTION, NAMESPACE) and set up `lintr` defaults.
- Establish `testthat` structure with unit tests for items above (constructor validations, space constructors, map/uncertainty creation, plan node formation).
- Configure GitHub Actions (or equivalent CI) to run `R CMD check` and unit tests.

## Deliverables
- `gdsfmri` package able to construct validated GDS objects, space definitions, and basic plan scaffolding.
- Unit tests covering new functionality (≥70% coverage target for modules touched in this sprint).
- CI workflow validating check/test on push.
- Updated documentation (`NEWS`, README stub, pkgdown scaffold optional) capturing Sprint 1 functionality.

## Dependencies & Notes
- All implementations should adhere to TECHNICAL_SPECIFICATION.md rev 0.1.0 and fmrigds_blueprint.md.
- Leave adapters, executor logic, variance propagation, and meta-analysis calculations for future sprints.
- Maintain feature flags or TODO markers for functionality planned in later sprints to avoid regressions.

## Looking Ahead (Sprints 2–8)
- Sprint 2: Storage adapters (tabular & NIfTI) + execute_stream scaffolding.
- Sprint 3: Derivation executor, variance propagation, and map_to computation.
- Sprint 4: mask(), align(), MapFamily execution, and optimizer pushdown rules.
- Sprint 5: meta-analysis reducers, Stouffer/Fisher combiners, df handling.
- Sprint 6: HDF5 adapter (Path B), provenance persistence, plan serialization.
- Sprint 7: fmristore adapter integration, write_out exporters, factorial contrast utilities.
- Sprint 8: Polish, documentation, v0.1.0 release readiness (pkgdown, examples, migration guide).

