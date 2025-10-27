# gdsfmri Implementation Sprint Plan (Sprints 1-8)

**Package:** `gdsfmri`
**Version Target:** 0.1.0
**Timeline:** 8 sprints (2-3 weeks each, ~4-6 months total)
**Testing:** testthat, minimum 90% coverage by Sprint 8
**Documentation:** roxygen2, pkgdown site by Sprint 8

---

## Sprint Overview

| Sprint | Focus Area | Key Deliverables | Dependencies |
|--------|-----------|------------------|--------------|
| 1 | Core Foundations | GDS, Space, Map, Plan infrastructure | None |
| 2 | Storage Adapters | Tabular & NIfTI adapters, compute scaffolding | Sprint 1 |
| 3 | Statistical Operations | Derivation execution, variance propagation | Sprints 1-2 |
| 4 | Spatial Operations | Alignment, masks, optimizer | Sprints 1-3 |
| 5 | Meta-Analysis | Fixed-effects, Stouffer, Fisher | Sprints 1-4 |
| 6 | Persistence | HDF5 adapter, plan serialization | Sprints 1-5 |
| 7 | Integration | fmristore adapter, exporters | Sprints 1-6 |
| 8 | Production Ready | Polish, docs, examples, release | Sprints 1-7 |

---

## Sprint 1: Core Foundations ✓

**Status:** COMPLETED (per sprint1.md)

**Objectives:**
- Establish core data structures (GDS, Space, Map, UncertaintyRule, Plan)
- Build validation framework
- Set up package infrastructure and testing

**Key Deliverables:**
- ✓ `new_gds()` with complete validation
- ✓ All Space constructors (voxel, parcels, surface, basis)
- ✓ `map_linear()`, `MapFamily()` with helper constructors
- ✓ `UncertaintyRule()` implementation
- ✓ Plan infrastructure (`gds_plan`, operation nodes, `add_op`, `as_plan`)
- ✓ Assay registry with default registrations
- ✓ Basic `subset()` and `derive()` verbs (lazy only)
- ✓ Package structure with testthat, CI workflow
- ✓ Unit tests (≥70% coverage for touched modules)

**Technical Debt:** None (foundation sprint)

---

## Sprint 2: Storage Adapters & Compute Scaffolding

**Objectives:**
- Implement first working adapters (tabular, NIfTI)
- Build compute execution framework
- Enable end-to-end lazy→realized pipeline

**Key Deliverables:**

### 2.1 Adapter Infrastructure
- Adapter registry (`register_adapter`, `get_adapter`, `detect_adapter`)
- Adapter interface specification enforcement
- Probe result validation utilities

### 2.2 Tabular Adapter (CSV/Parquet)
- CSV/TSV reader with auto-detection
- Column mapping for beta, var, subject, sample, contrast
- Parcel space auto-construction
- Block reading support (for large CSVs)
- Probe implementation returning complete metadata

### 2.3 NIfTI Adapter
- Single file and multi-file (BIDS-like) layout detection
- Mask creation/loading utilities
- Voxel space construction from NIfTI headers
- Subject/contrast extraction from file paths
- On-disk lazy access via RNifti or neuroim2

### 2.4 Compute Framework
- `compute()` function with sink selection (memory, HDF5 stub)
- `init_executor()` with adapter handle management
- `execute_stream()` block processor (basic version)
- Block iteration with sample-axis chunking
- Simple operation application (subset, derive stubs)

### 2.5 Testing & Integration
- Unit tests for adapters (tabular, NIfTI)
- Integration test: CSV → compute() → realized GDS
- Integration test: NIfTI → compute() → realized GDS
- Edge cases: missing files, malformed CSVs, dimension mismatches

**Technical Debt:**
- HDF5 sink is stubbed (writes to memory for now)
- Optimizer is pass-through (no actual optimization yet)
- Limited derive execution (only var↔se, basic t/z)

**Dependencies:** Sprint 1

**Estimated Effort:** 3 weeks

---

## Sprint 3: Statistical Operations & Derivations

**Objectives:**
- Implement complete derivation engine
- Build variance propagation machinery
- Enable `map_to()` with uncertainty propagation

**Key Deliverables:**

### 3.1 Derivation Engine
- `execute_derive()` with all derivation rules (§3.3, §6.3)
- var ↔ se, t from beta/var, z from t/df, p from F
- Overwrite control and validation
- df handling for t→z conversions

### 3.2 Variance Propagation
- `propagate_variance_independent()` implementation (§6.1.1)
- `propagate_variance_covariance()` implementation (§6.1.2)
- Satterthwaite df aggregation (§6.2.1)
- Test suite for variance propagation correctness

### 3.3 map_to() Execution
- Assay-aware mapping logic (location vs. test statistics)
- Integration with UncertaintyRule modes
- Recomputation of t/z after mapping beta/var
- Validation: refuse to map t/z without effect scale
- Stouffer/Fisher combiners for z-only mapping (basic implementation)

### 3.4 Compute Integration
- Wire `op_derive` into executor
- Wire `op_map` into executor with variance propagation
- Block-wise application of mappings
- Memory efficiency checks

### 3.5 Testing
- Unit tests for all derivation paths
- Variance propagation correctness tests (independence, covariance)
- map_to() integration tests with synthetic data
- Edge cases: NA handling, dimension mismatches, zero variance

**Technical Debt:**
- Full covariance mode may need optimization (large blocks)
- Kernel-based covariance deferred

**Dependencies:** Sprints 1-2

**Estimated Effort:** 3 weeks

---

## Sprint 4: Spatial Operations (Alignment & Masks)

**Objectives:**
- Implement alignment transforms (subject→group)
- Add mask operations and packing
- Build plan optimizer with pushdown rules

**Key Deliverables:**

### 4.1 Alignment Execution
- `execute_align()` for subject-specific transforms
- Support for OrthogonalFamily, OTFamily, WarpFamily
- Per-subject variance propagation
- Space updating (from→to transition)

### 4.2 Mask Operations
- `compute_mask_union()`, `compute_mask_intersection()`, `compute_mask_threshold()`
- `pack_assays()` and `unpack_assays()` utilities
- `execute_mask_policy()` with space updates
- Storage mode transitions (dense→packed)

### 4.3 Late-Binding Map Registration
- `register_map()` for attaching transforms to GDS/Plan
- Validation of subject lists and dimensions
- Space graph metadata structure
- Provenance tracking for registered maps

### 4.4 Plan Optimizer (Real Implementation)
- `pushdown_subset()` - push to adapter
- `coalesce_derives()` - merge consecutive derives
- `push_mask_early()` - after align, before heavy ops
- `fuse_map_reduce()` - kernel fusion markers
- Optimizer integration tests

### 4.5 Compute Integration
- Wire `op_align_to_group` into executor
- Wire `op_mask_policy` into executor
- Optimizer execution in `compute()`
- `explain()` function for readable plan output

### 4.6 Testing
- Alignment correctness (orthogonal, OT, synthetic warps)
- Mask policy tests (intersection, threshold)
- Pack/unpack round-trip tests
- Optimizer correctness (subset pushdown, fusion)

**Technical Debt:**
- Deformable warp families need external registration tool support
- OT plan storage/serialization deferred

**Dependencies:** Sprints 1-3

**Estimated Effort:** 3 weeks

---

## Sprint 5: Meta-Analysis & Reduction

**Objectives:**
- Implement complete `reduce()` functionality
- Add fixed-effects meta-analysis
- Implement Stouffer and Fisher combiners

**Key Deliverables:**

### 5.1 Fixed-Effects Meta-Analysis
- `combine_fixed()` implementation (§6.4.1)
- Inverse-variance weighting
- Pooled variance computation
- n_eff weighting option
- Dimensional reduction (subjects→1)

### 5.2 Evidence Combiners
- `combine_stouffer()` with equal and weighted modes (§6.4.2)
- `combine_fisher()` implementation (§6.4.3)
- p-value to z-score conversions
- Sign preservation logic

### 5.3 reduce() Execution
- `execute_reduce()` dispatcher
- Method selection (fixed, stouffer, fisher)
- By-contrast vs. by-sample aggregation
- df handling for pooled statistics

### 5.4 Random-Effects Scaffold
- DerSimonian-Laird placeholder
- Between-study variance estimation stub
- Documentation of future implementation

### 5.5 Compute Integration
- Wire `op_reduce` into executor
- Block-wise reduction support
- Memory management for large reductions

### 5.6 Testing
- Fixed-effects correctness (synthetic data with known answers)
- Stouffer combination validation
- Fisher combination validation
- Edge cases: single subject, all NA, heterogeneous df

**Technical Debt:**
- Random-effects is stubbed (Sprint 6 or post-v0.1.0)

**Dependencies:** Sprints 1-4

**Estimated Effort:** 2 weeks

---

## Sprint 6: Persistence & Serialization

**Objectives:**
- Implement HDF5 adapter with GDS layout
- Add plan serialization (JSON/YAML)
- Enable provenance persistence

**Key Deliverables:**

### 6.1 HDF5 GDS Adapter
- Writer: `write_gds_h5()` with spec layout (§7.5)
- Reader: `read_gds_h5()` with probe
- Chunking strategy (sample-axis primary)
- Compression (gzip or zstd)
- Metadata persistence (JSON in attributes)

### 6.2 HDF5 Sink for compute()
- `sink = "HDF5"` implementation
- Block-wise writing to HDF5 datasets
- Space serialization to /gds/space/*
- Provenance appending to /gds/provenance

### 6.3 Plan Serialization
- `save_plan()` with JSON/YAML support (§5.2)
- `load_plan()` with deserialization
- Operator serialization (matrices, descriptors)
- Round-trip tests

### 6.4 Provenance System
- `explain()` implementation (§5.3)
- `digest_plan()` stable hashing (§5.4)
- Provenance graph persistence
- Human-readable log generation

### 6.5 Testing
- HDF5 write/read round-trip
- Plan save/load round-trip
- Provenance correctness (node graph, hashes)
- Large file handling (>1GB datasets)

**Technical Debt:**
- Zarr/TileDB backends deferred
- Plan migration (schema versioning) deferred

**Dependencies:** Sprints 1-5

**Estimated Effort:** 3 weeks

---

## Sprint 7: Integration & Export

**Objectives:**
- Integrate with fmristore
- Implement write_out() exporters
- Add factorial contrast utilities

**Key Deliverables:**

### 7.1 fmristore Adapter
- Detect fmristore HDF5 files
- Path A: read existing layouts (labeled_volume, parcellated, latent)
- Path B: read/write /gds shim group
- Space extraction from fmristore metadata
- Map construction (voxel→parcel, latent→voxel)

### 7.2 write_out() Exporters
- NIfTI exporter (requires voxel space)
- CSV exporter (good for parcels/ROI)
- Parquet exporter (Arrow integration)
- Format-specific options (compression, precision)

### 7.3 Contrast Info & Factorial Designs
- `contrast_info()` constructor (§4.3)
- Factorial design specification
- Dependency tracking (main vs. interaction)
- Integration with metadata

### 7.4 Compute Integration
- Wire `op_write` into executor
- Multi-format writing
- Lazy write execution (deferred until compute)

### 7.5 Migration Utilities
- `as_gds()` for legacy group_data
- `as_group_data()` for backward compat
- Shim layer in fmrireg (if applicable)

### 7.6 Testing
- fmristore integration tests (read/write)
- NIfTI export correctness (header preservation)
- CSV/Parquet export validation
- Factorial contrast metadata tests

**Technical Debt:**
- fmristore write utilities need fmristore package updates
- Surface exporters (GIFTI) deferred

**Dependencies:** Sprints 1-6

**Estimated Effort:** 3 weeks

---

## Sprint 8: Production Ready & Release

**Objectives:**
- Comprehensive documentation
- Real-world examples and vignettes
- Performance optimization
- v0.1.0 release

**Key Deliverables:**

### 8.1 Documentation
- Complete roxygen2 documentation for all exported functions
- pkgdown website with rendered docs
- 5+ vignettes:
  - Introduction to GDS
  - Working with Spaces
  - Lazy Execution & Plans
  - Alignment & Spatial Operations
  - fmristore Integration
- Migration guide from group_data

### 8.2 Examples & Real Data
- Example datasets in inst/extdata/
- ROI analysis example (tabular)
- Voxel analysis example (NIfTI)
- Multi-subject alignment example
- Meta-analysis example

### 8.3 Performance Optimization
- Profiling with real fMRI datasets
- Memory optimization in executors
- Sparse matrix optimizations
- Block size tuning

### 8.4 Testing Completeness
- Achieve ≥90% code coverage
- Integration tests with real fMRI data
- Stress tests (large N subjects, whole-brain voxels)
- Edge case coverage

### 8.5 Release Preparation
- CRAN submission checks
- NEWS.md with complete changelog
- README with badges and quick start
- DESCRIPTION with all dependencies
- LICENSE file
- Code of Conduct
- Contributing guidelines

### 8.6 Refinement
- User feedback incorporation (if beta testers available)
- Bug fixes from testing
- API polish (consistency, naming)
- Error message improvement

**Technical Debt:**
- Mark future features (random-effects, advanced surfaces) in TODO
- Document extension points for community

**Dependencies:** Sprints 1-7

**Estimated Effort:** 3 weeks

---

## Post-v0.1.0 Roadmap

**Future Sprints (0.2.0+):**
- Random-effects meta-analysis (DerSimonian-Laird, REML)
- Surface operations (GIFTI I/O, surface mapping)
- Advanced covariance structures (AR1, compound symmetry)
- Multi-level models (crossed random effects)
- Bayesian meta-analysis hooks
- GPU acceleration for variance propagation
- Distributed computing (future/clustermq integration)

---

## Risk Management

**Technical Risks:**
1. **Memory scaling:** Large voxel × subject arrays may exceed RAM
   - Mitigation: Aggressive blocking, HDF5 sink by Sprint 6

2. **Adapter complexity:** NIfTI layouts are heterogeneous
   - Mitigation: Start with simple layouts, iterate

3. **Statistical correctness:** Variance propagation is subtle
   - Mitigation: Extensive unit tests with known answers, Sprint 3 focus

4. **Performance:** R loops may be slow for large datasets
   - Mitigation: Vectorization, Rcpp if needed (Sprint 8)

**Process Risks:**
1. **Scope creep:** Blueprint is comprehensive
   - Mitigation: Strict sprint boundaries, defer non-critical features

2. **Testing burden:** 90% coverage is aggressive
   - Mitigation: Test-driven development from Sprint 1

3. **Documentation lag:** Complex API needs good docs
   - Mitigation: Roxygen from day 1, Sprint 8 dedicated to docs

---

## Success Criteria

**Sprint-Level:**
- Each sprint delivers working, tested code
- No regressions in prior sprint functionality
- CI passes on all commits

**Release-Level (v0.1.0):**
- ≥90% test coverage
- All core workflows (ROI, voxel, alignment, meta-analysis) functional
- Complete documentation (man pages, vignettes, pkgdown site)
- CRAN-ready package structure
- Zero critical bugs in issue tracker

---

## Dependencies & Tools

**R Packages (Imports):**
- Matrix (sparse matrices)
- digest (hashing)
- jsonlite (JSON)
- hdf5r (HDF5 I/O)
- RNifti or neuroim2 (NIfTI I/O)

**R Packages (Suggests):**
- testthat (testing)
- knitr, rmarkdown (vignettes)
- pkgdown (website)
- arrow (Parquet)
- data.table (fast CSV)
- fmristore (integration)

**Development Tools:**
- devtools, usethis
- roxygen2
- lintr, styler
- covr (coverage)
- GitHub Actions (CI)

---

## Sprint Timeline (Estimated)

| Sprint | Start | Duration | End | Cumulative |
|--------|-------|----------|-----|------------|
| 1 | Week 0 | 2 weeks | Week 2 | 2 weeks |
| 2 | Week 2 | 3 weeks | Week 5 | 5 weeks |
| 3 | Week 5 | 3 weeks | Week 8 | 8 weeks |
| 4 | Week 8 | 3 weeks | Week 11 | 11 weeks |
| 5 | Week 11 | 2 weeks | Week 13 | 13 weeks |
| 6 | Week 13 | 3 weeks | Week 16 | 16 weeks |
| 7 | Week 16 | 3 weeks | Week 19 | 19 weeks |
| 8 | Week 19 | 3 weeks | Week 22 | 22 weeks |

**Total:** ~22 weeks (5.5 months)

---

## Notes for Implementation

1. **Incremental Development:** Each sprint should build incrementally; avoid large refactors
2. **Test-Driven:** Write tests before/during implementation, not after
3. **Documentation:** Roxygen comments as you code, not at the end
4. **Code Review:** If multiple developers, review before merging
5. **Continuous Integration:** CI must pass before sprint completion
6. **User Feedback:** Gather feedback early (Sprint 4+) from potential users

---

**Document Version:** 1.0
**Last Updated:** 2025-01-XX
**Next Review:** After Sprint 2 completion
