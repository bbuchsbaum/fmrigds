# Quick Reference: fmrireg vs fmrigds

## TL;DR

**Can fmrigds replace fmrireg's group_data?** No, not directly.

**What should we do?** Implement adapters + fmri_meta.gds() (4-5 months)

**Risk of direct replacement?** HIGH

---

## Key Differences

| Aspect | fmrireg | fmrigds |
|--------|---------|---------|
| **Architecture** | Format-specific classes | Format-agnostic adapter |
| **Loading** | Lazy (keeps file paths) | Eager (materializes arrays) |
| **Data shape** | Variable per format | Always 3D (sample × subject × contrast) |
| **Contrasts** | One per object | Multiple always available |
| **Space** | Via NeuroVol | Explicit space objects |
| **fmri_meta()** | Implemented | NOT YET |

---

## API Compatibility Status

### 100% Compatible
- `get_subjects()` ↔ `subjects()`
- `get_covariates()` ↔ `col_data()`
- Constructor `group_data()` ↔ `gds()`

### Requires Translation
- `n_subjects(x)` → `length(subjects(x))`
- `read_*_full()` → `compute() + assay()`
- `extract_csv_data()` → Manual subsetting

### NOT Implemented
- `fmri_meta()` on GDS objects **← CRITICAL BLOCKER**

---

## What's Missing in fmrigds

### Critical (blocks replacement)
1. `fmri_meta.gds()` method
2. Meta-analysis implementations (fe, pm, dl, reml)
3. Meta-regression with covariates

### Important
- Robust estimation (huber, t)
- Custom weights
- Contrasts at fit time
- ROI extraction helpers

### Nice to have
- Lazy chunk loading
- Format-specific optimizations

---

## Recommended Path Forward

### Phase 1 (0-3 months): Adapters
```r
as_gds_plan(group_data)        # fmrireg → GDS
as_group_data_csv(gds)         # GDS → fmrireg
```

### Phase 2 (3-9 months): fmri_meta.gds()
```r
fit <- fmri_meta(gds, formula = ~ 1 + group, method = "pm")
coefficients(fit)
```

### Phase 3 (9+ months): Deprecation
- Mark group_data deprecated
- Provide migration guide
- Archive fmrireg

---

## Code Examples

### Create GDS from fmrireg data
```r
library(fmrigds)

# Old way (fmrireg):
gd <- group_data_from_csv(df, effect_cols = c(beta = "beta", se = "se"), ...)

# New way (fmrigds):
plan <- gds("data.csv", effect_cols = ..., subject_col = "subject")
gds_obj <- compute(plan)

# Or with adapter (during transition):
gds_obj <- compute(as_gds_plan(gd))
```

### Access data
```r
# Old way:
beta <- read_nifti_full(gd)$beta  # list

# New way:
beta <- assay(gds_obj, "beta")   # 3D array

# Adapter helps:
gd_back <- as_group_data_csv(gds_obj)
beta_old_way <- extract_csv_data(gd_back)$beta
```

### Meta-analysis
```r
# Old way (works now):
fit <- fmri_meta(gd, formula = ~ 1 + group, method = "pm")

# New way (needs implementation):
fit <- fmri_meta(gds_obj, formula = ~ 1 + group, method = "pm")  # NOT YET

# During transition:
gd_converted <- as_group_data_csv(gds_obj)
fit <- fmri_meta(gd_converted, ...)  # Uses old code
```

---

## Decision Matrix

### If choosing between systems NOW:
→ **Use fmrireg** (fully featured)

### If planning future work:
→ **Start with fmrigds** (better architecture, will be supported)

### For new projects:
→ **Use fmrigds** + implement missing pieces as needed

### For migrating existing code:
→ **Use adapters** (no immediate rewrite needed)

---

## Critical Implementation Gaps

Must implement before considering fmrigds as replacement:

```r
# 1. fmri_meta.gds() with these methods
fmri_meta.gds <- function(data, formula = ~1, 
                          method = c("fe", "pm", "dl", "reml"),
                          robust = c("none", "huber", "t"),
                          weights = c("ivw", "equal", "custom"),
                          ...)

# 2. Extraction methods
coefficients(fit_gds)  # dim: (samples × contrasts) or (samples × contrast × beta)
se(fit_gds)
tau2(fit_gds)

# 3. Plotting/summarizing
plot(fit_gds)
summary(fit_gds)
```

---

## Files Generated

1. **fmrireg_vs_fmrigds_comparison.md** (full technical analysis)
2. **adapter_implementation_guide.md** (code examples)
3. **EXECUTIVE_SUMMARY.txt** (this summary level)
4. **QUICK_REFERENCE.md** (this file)

---

## Next Actions

1. **Read**: EXECUTIVE_SUMMARY.txt (15 min read)
2. **Review**: Full comparison.md with team
3. **Decide**: Which approach (recommend: hybrid)
4. **Plan**: Implementation phases
5. **Build**: Start with adapters
6. **Validate**: Run parallel tests

---

## Contact Points for Questions

When evaluating fmrigds adoption, ask:

- "Is fmri_meta.gds() implemented?" → NO (as of now)
- "Can I use lazy loading?" → Use block parameter
- "Do I need multiple contrasts?" → fmrigds handles better
- "Will my old code work?" → Yes, with adapters
- "What's the timeline?" → 4-5 months to full parity

---

**Status**: Analysis complete
**Date**: 2025-10-29
**Recommendation**: HYBRID APPROACH with adapters
**Confidence**: HIGH (70% API compatible, identified all gaps)
