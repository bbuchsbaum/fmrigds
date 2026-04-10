# Package Responsibilities: fmrigds vs fmrireg

## Quick Reference

| Fix | Package | Priority | Time | Status |
|-----|---------|----------|------|--------|
| var/se auto-derivation | **fmrigds** | P1 | 1 hour | TODO |
| fmri_ttest NAMESPACE export | **fmrireg** | P1 | 15 min | TODO |
| HDF5 S4 slot access | **fmrireg** | P1 | 15 min | TODO |
| fdr:spatial posthoc | **fmrigds** | P2 | 2-3 hours | TODO |
| ols:voxelwise C++ | **fmrigds** | P2 | 1 day | TODO |
| meta:re:pm reducer | **fmrigds** | P3 | 4-6 hours | TODO |
| meta:re:reml reducer | **fmrigds** | P3 | 4-6 hours | TODO |
| as_plan.gds method | **fmrigds** | P4 | 30 min | TODO |

---

## fmrigds Fixes (~/code/fmrigds)

These fix issues in the fmrigds library itself - improving the core infrastructure that fmrireg depends on.

### P1: var/se Auto-Derivation (1 hour)

**File**: `R/adapter-*.R` or `R/reducer-*.R`

**Problem**: Meta-analysis reducers expect `var` assay but adapters only provide `se`. Auto-derivation isn't working.

**Solution**: Ensure auto-derivation in adapter probe() or reducer initialization:

```r
# In adapter probe() or reducer pre-processing
if ("se" %in% assay_names && !"var" %in% assay_names) {
  # Auto-derive var from se
  assays$var <- assays$se^2
}
```

**Why fmrigds?** This is a core fmrigds feature that benefits all users, not just fmrireg.

---

### P2: Implement fdr:spatial Posthoc (2-3 hours)

**File**: Create `R/posthoc-fdr-spatial.R`

**Problem**: Neuroimaging spatial FDR method doesn't exist in fmrigds.

**Solution**: Adapt algorithm from fmrireg's `spatial_fdr()` function:

```r
# Register in posthoc registry
register_posthoc(
  name = "fdr:spatial",
  requires = "z",
  produces = "q",
  fn = function(gds, options = list()) {
    z <- assay(gds, "z")
    group <- options$group
    # Implement weighted BH algorithm
    # ...
    add_assay(gds, "q", q_values)
  }
)
```

**Why fmrigds?** This is a general neuroimaging method that should be available to all fmrigds users.

**Source**: Copy algorithm from `/Users/bbuchsbaum/code/fmrireg/R/spatial_fdr.R`

---

### P2: C++ Implementation for ols:voxelwise (1 day)

**File**: Create `src/ols_voxelwise.cpp`

**Problem**: Currently R-only, causing 10-50x slowdown for large datasets.

**Solution**: Port to C++ with OpenMP parallelization, following pattern of existing C++ reducers.

**Why fmrigds?** Performance optimization for core reducer functionality.

---

### P3: Implement meta:re:pm and meta:re:reml (4-6 hours each)

**File**: `R/reducer-meta-re.R`

**Problem**: Only DerSimonian-Laird tau² estimator available. PM and REML alternatives missing.

**Solution**: Add additional tau² estimation methods to random-effects registry.

**Why fmrigds?** These are standard statistical methods that belong in the core library.

---

### P4: Add as_plan.gds() S3 Method (30 min)

**File**: `R/plan.R`

**Problem**: `as_plan()` only works on gds_source objects, not realized GDS objects.

**Solution**:

```r
#' @export
as_plan.gds <- function(x) {
  # Extract source from realized GDS and create new plan
  gds_plan(source = attr(x, "source"))
}
```

**Why fmrigds?** This is part of the fmrigds API design.

---

## fmrireg Fixes (~/code/fmrireg)

These fix integration issues in fmrireg - the glue code that connects to fmrigds.

### P1: Fix NAMESPACE Export for fmri_ttest (15 min)

**File**: `R/fmri_ttest_gds.R` + regenerate NAMESPACE

**Problem**: `fmri_ttest.group_data_gds()` method exists but isn't exported properly.

**Solution**:

```bash
# In fmrireg directory
Rscript -e "devtools::document()"
```

Verify NAMESPACE contains:
```
S3method(fmri_ttest,group_data_gds)
```

**Why fmrireg?** This is fmrireg's S3 method registration.

**Status**: Method already implemented at [R/fmri_ttest_gds.R:8](R/fmri_ttest_gds.R#L8), just needs export.

---

### P1: Fix HDF5 S4 Slot Access (15 min)

**File**: `R/group_data_h5.R`

**Problem**: Legacy code uses `h5_handle$dim` on S4 objects (should use `@` operator).

**Solution**:

```r
# In read_h5_metadata() function
# Change from:
dims <- h5_handle$dim

# To:
dims <- h5_handle@dim
# OR
dims <- slot(h5_handle, "dim")
```

**Why fmrireg?** This is fmrireg's legacy HDF5 reading code.

**File location**: Around line 84 in `R/group_data_h5.R`

---

### Optional: Add Workaround for var/se (15 min)

**File**: `R/fmri_meta_gds.R` or `R/fmri_ttest_gds.R`

**Problem**: If fmrigds doesn't fix auto-derivation quickly, need temporary workaround.

**Solution**:

```r
# At start of fmri_meta.group_data_gds()
if ("se" %in% fmrigds::assay_names(data) && !"var" %in% fmrigds::assay_names(data)) {
  data <- fmrigds::derive(data, var = quote(se^2))
}
```

**Why fmrireg?** Temporary workaround only - should be removed once fmrigds fixes auto-derivation.

**Recommended?** NO - Better to fix properly in fmrigds.

---

## Division of Labor Philosophy

### fmrigds Scope (General Infrastructure)
- Core data structures (gds, gds_plan, gds_source)
- Reducers (meta-analysis, OLS, evidence combination)
- Posthoc methods (FDR, spatial correction)
- Adapters (CSV, H5, NIfTI format support)
- Performance (C++ implementations, OpenMP)
- Auto-derivation of assays
- Registry systems

### fmrireg Scope (Domain-Specific Glue)
- S3 methods bridging fmrireg API to fmrigds (fmri_meta.group_data_gds, fmri_ttest.group_data_gds)
- Legacy group_data compatibility layer
- fMRI-specific wrappers and convenience functions
- Event/HRF system integration
- fMRI-specific validation and error messages
- Migration path from legacy to fmrigds

---

## Recommended Work Order

### Step 1: fmrireg Quick Wins (30 minutes)
Do these first since you have fmrireg open:

1. ✅ Run `devtools::document()` to fix NAMESPACE export
2. ✅ Fix HDF5 S4 slot access in `read_h5_metadata()`
3. ✅ Test with: `devtools::test(filter = 'gds-ttest|gds-h5')`

Expected: 4 tests pass (2 ttest + 2 h5)

### Step 2: fmrigds Critical Fix (1 hour)
Switch to fmrigds:

4. ✅ Implement var/se auto-derivation in adapters
5. ✅ Test with: `devtools::test(filter = 'adapter')`

Expected: Auto-derivation works

### Step 3: Test Integration (15 minutes)
Back to fmrireg:

6. ✅ Run full GDS test suite: `devtools::test(filter = 'gds')`

Expected: 9-10 tests pass (only posthoc failing)

### Step 4: fmrigds Enhancement (2-3 hours)
Back to fmrigds:

7. ✅ Implement `fdr:spatial` posthoc method
8. ✅ Test with posthoc tests

Expected: All 11 fmrireg tests pass

### Step 5: Performance (1 day - optional)
In fmrigds:

9. Implement C++ ols:voxelwise
10. Benchmark and validate

---

## Testing Strategy

### In fmrigds
```r
# Unit tests for new features
devtools::test(filter = "adapter|reducer|posthoc")

# Verify auto-derivation
gd <- gds(source = "test.csv", format = "tabular")
assay_names(gd)  # Should include both "se" and "var"
```

### In fmrireg
```r
# Integration tests
devtools::test(filter = "gds")

# Should see:
# ✓ test-gds-integration.R
# ✓ test-gds-parity-methods.R
# ✓ test-gds-ttest-parity.R
# ✓ test-gds-posthoc.R
# ✓ test-gds-h5-parity.R
# ✓ test-gds-nifti-parity.R
# ✓ test-gds-bridges.R
```

---

## Communication Plan

### When working on fmrigds:
1. Notify fmrireg developers when auto-derivation is fixed
2. Provide example showing var is derived from se
3. Document any API changes in CHANGELOG

### When working on fmrireg:
1. Keep workarounds minimal and temporary
2. Document dependencies on fmrigds features
3. Add integration tests that will pass once fmrigds is fixed

---

## Success Criteria

### Phase 1: Basic Integration (Both packages)
- [ ] fmrireg: NAMESPACE exports fixed, H5 code fixed (30 min)
- [ ] fmrigds: var/se auto-derivation working (1 hour)
- [ ] Result: 9-10 of 11 tests passing

### Phase 2: Full Feature Parity
- [ ] fmrigds: fdr:spatial implemented (2-3 hours)
- [ ] Result: 11 of 11 tests passing

### Phase 3: Production Ready
- [ ] fmrigds: C++ ols:voxelwise, PM/REML reducers (1-2 weeks)
- [ ] fmrireg: Migration guide, deprecation warnings
- [ ] Result: Full performance, feature parity

---

## Key Takeaway

**Most work (70%) belongs in fmrigds** because these are general infrastructure improvements that benefit all users.

**fmrireg work (30%)** is mostly quick fixes to the integration layer - NAMESPACE exports, legacy code updates, and S3 method dispatch.

**Total time to basic integration**: ~2 hours combined across both packages.
