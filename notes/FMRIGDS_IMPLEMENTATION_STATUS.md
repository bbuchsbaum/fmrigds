# fmrigds Implementation Status Audit

**Date:** 2025-10-29
**Audited by:** Claude Code
**Purpose:** Comprehensive inventory of reducers, posthoc methods, and adapters

---

## 1. Reducer Registry

### API Overview

The reducer registry is implemented in `/Users/bbuchsbaum/code/fmrigds/R/reducer-registry.R`:

- **Storage:** Environment `.gds_reducers` (package-private)
- **Registration:** `register_reducer(name, fun, requires, provides, options_schema)`
- **Retrieval:** `get_reducer(name)` returns list with `name`, `fun`, `requires`, `provides`, `options_schema`
- **Discovery:** `list_reducers()` returns sorted character vector of registered names
- **Legacy mapping:** `.normalize_reducer_name()` maps old names like "fixed" → "meta:fe"

### Kernel Signature

All reducer kernels follow this signature:
```r
function(beta, var, X, z, p, df, df1, df2, opts) -> named list
```

Where:
- `beta`, `var`, `z`, `p`, `df`, `df1`, `df2`: Input assay matrices (subjects × samples)
- `X`: Design matrix (subjects × predictors) for regression reducers
- `opts`: Named list of options
- **Returns:** Named list of output arrays to be written as assays

---

## 2. Registered Reducers

### Complete Table

| Reducer Name | Status | C++ Impl | R Fallback | Requires | Provides | Options |
|--------------|--------|----------|------------|----------|----------|---------|
| `meta:fe` | ✓ Complete | ✓ `meta_fe_cpp` | ✓ `.colwise_fe` | beta, var | beta_g, var_g, se_g, z_g, p_g, Q, I2 | eps, alternative, min_subjects |
| `meta:re` | ✓ Complete | ✓ `meta_re_dl_cpp` | ✓ `.colwise_tau2_dl` + RE logic | beta, var | beta_g, var_g, se_g, z_g, p_g, tau2, Q, I2 | eps, alternative, min_subjects, tau2="DL" |
| `meta:fe_reg` | ✓ Complete | ✓ `meta_fe_reg_cpp` | ✓ R loop | beta, var, X | coef, se_coef, Q, df_res | eps, min_subjects |
| `meta:re_reg` | ✓ Complete | ✓ `meta_re_reg_dl_cpp` | ✓ R loop | beta, var, X | coef, se_coef, tau2, Q, df_res | eps, min_subjects |
| `combine:stouffer` | ✓ Complete | ✓ `stouffer_combine_cpp` | ✓ R | z | z_g, p_g | weights, min_subjects |
| `combine:fisher` | ✓ Complete | ✓ `fisher_combine_cpp` | ✓ R | p | p_g, chi2, df | min_subjects |
| `combine:lancaster` | ✓ Complete | ✓ `lancaster_combine_cpp` | ✓ R | p | p_g, chi2, df | dfw (auto-derived from df), min_subjects |
| `ols:voxelwise` | ✓ Complete | ✗ None | ✓ R only | beta, X | coef, se_coef, t_coef, p_coef, sigma2, df_res | return_cov ("none"/"tri") |

**Legend:**
- ✓ Complete: Fully implemented and tested
- ✗ None: Not implemented
- C++ Impl: Optimized C++ version available (with OpenMP parallelization)
- R Fallback: Pure R implementation available

### Implementation Details

#### Meta-Analysis Reducers

**meta:fe (Fixed-Effects)**
- **Location:** `R/reducers-core.R:33-41`, `src/reducers_core.cpp:21-90`
- **Method:** Inverse-variance weighted mean
- **Formula:** `μ = Σ(wi·yi) / Σwi` where `wi = 1/var_i`
- **Outputs:**
  - `beta_g`: pooled effect
  - `var_g`: pooled variance = 1/Σwi
  - `se_g`: √var_g
  - `z_g`: μ/se (z-score)
  - `p_g`: p-value (two-sided, less, or greater based on `alternative`)
  - `Q`: Cochran's Q heterogeneity statistic
  - `I2`: I² heterogeneity index
- **Tail options:** "two.sided", "less", "greater"
- **Performance:** Parallel across samples via OpenMP

**meta:re (Random-Effects, DerSimonian-Laird)**
- **Location:** `R/reducers-core.R:43-65`, `src/reducers_core.cpp:95-184`
- **Method:** DerSimonian-Laird tau² estimation
- **Formula:**
  1. Compute FE estimate and Q statistic
  2. `τ² = max(0, (Q - (k-1)) / C)` where `C = Σwi - Σwi²/Σwi`
  3. Re-weight: `wi* = 1/(var_i + τ²)`
  4. `μ = Σ(wi*·yi) / Σwi*`
- **Outputs:** Same as FE plus `tau2`
- **Note:** Currently only DL method implemented; PM, REML not yet available

**meta:fe_reg (Fixed-Effects Meta-Regression)**
- **Location:** `R/reducers-core.R:120-149`, `src/reducers_core.cpp:189-247`
- **Method:** Weighted least squares with `W = diag(1/var_i)`
- **Formula:** `β = (X'WX)^(-1) X'Wy`
- **Outputs:**
  - `coef`: [p × samples] matrix of regression coefficients
  - `se_coef`: standard errors
  - `Q`: residual heterogeneity
  - `df_res`: residual degrees of freedom
- **Expansion:** Per-term assays created during execution (e.g., `coef:(Intercept)`, `coef:age`)

**meta:re_reg (Random-Effects Meta-Regression)**
- **Location:** `R/reducers-core.R:151-190`, `src/reducers_core.cpp:252-337`
- **Method:** Two-stage DL with meta-regression
- **Formula:**
  1. WLS fit to get residuals
  2. `τ² = max(0, (Q - df) / (Σwi - tr(H)))` where `H = X(X'WX)^(-1)X'W`
  3. Re-fit with `W* = diag(1/(var_i + τ²))`
- **Outputs:** Same as `fe_reg` plus `tau2`

#### Evidence Combination Reducers

**combine:stouffer (Stouffer's Z-method)**
- **Location:** `R/reducers-core.R:67-82`, `src/reducers_core.cpp:342-371`
- **Method:** Weighted Z-score combination
- **Formula:** `Z = Σ(wi·zi) / √(Σwi²)`
- **Outputs:** `z_g`, `p_g`
- **Weights:** Optional; defaults to equal weighting

**combine:fisher (Fisher's Method)**
- **Location:** `R/reducers-core.R:84-91`, `src/reducers_core.cpp:383-402`
- **Method:** Combine p-values via chi-squared
- **Formula:** `χ² = -2Σln(pi)`, df = 2k
- **Outputs:** `p_g`, `chi2`, `df`
- **P-value clamping:** (1e-300, 1-1e-16) to avoid log(0)

**combine:lancaster (Lancaster/Weighted Fisher)**
- **Location:** `R/reducers-core.R:93-103`, `src/reducers_core.cpp:407-434`
- **Method:** Weighted Fisher method
- **Formula:** Transform each pi via χ²(2wi) quantile, sum, test against χ²(2Σwi)
- **Outputs:** `p_g`, `chi2`, `df`
- **Weights:** Auto-derived from `df` or `df1` if not provided

#### Regression Reducers

**ols:voxelwise (Voxelwise OLS)**
- **Location:** `R/reducers-core.R:211-229` (registration), `R/reducers-core.R:233-278` (R impl)
- **Method:** Ordinary least squares across subjects, per sample
- **Formula:** `β = (X'X)^(-1)X'y` per sample column
- **Outputs:**
  - `coef`: [p × samples] coefficient matrix
  - `se_coef`: standard errors (√(diag(X'X)^(-1) · σ²))
  - `t_coef`: t-statistics
  - `p_coef`: p-values
  - `sigma2`: residual variance per sample
  - `df_res`: residual df per sample
- **Optional:** `cov_tri` (packed upper-triangular covariance per sample) when `return_cov = "tri"`
- **Performance:** R-only implementation (no C++ version yet)
- **Expansion:** Creates per-term assays like `coef:(Intercept)`, `coef:age`, etc.

---

## 3. Missing Reducers

### Meta-Analysis Variants

| Name | Status | Priority | Complexity | Notes |
|------|--------|----------|------------|-------|
| `meta:re:pm` | Missing | Medium | Low | Paule-Mandel tau² estimator; alternative to DL |
| `meta:re:reml` | Missing | Medium | Medium | REML tau² estimator via iterative optimization |
| `meta:re:ml` | Missing | Low | Medium | Maximum likelihood tau² estimator |
| `meta:re:eb` | Missing | Low | Medium | Empirical Bayes estimator |

**Implementation Path for PM:**
- Similar structure to DL
- Different tau² formula: iterative solution to `Σwi*(yi - μ)² / (1 + τ²wi) = k - 1`
- Complexity: ~50 lines R, ~80 lines C++
- Effort: 2-3 hours

**Implementation Path for REML:**
- Requires iterative optimization (optim or Newton-Raphson)
- Profile likelihood maximization
- Complexity: ~100 lines R, ~200 lines C++ (with optimizer)
- Effort: 1-2 days

### Evidence Combiners

All planned evidence combiners are implemented.

### Regression Reducers

All planned regression reducers are implemented.

---

## 4. Posthoc Registry

### API Overview

The posthoc registry is implemented in `/Users/bbuchsbaum/code/fmrigds/R/posthoc-registry.R`:

- **Storage:** Environment `.gds_posthoc` (package-private)
- **Registration:** `register_posthoc(name, fun, requires, provides)`
- **Retrieval:** `get_posthoc(name)` returns list with `name`, `fun`, `requires`, `provides`
- **Discovery:** `list_posthoc()` returns sorted character vector of registered names

### Kernel Signature

All posthoc kernels follow this signature:
```r
function(arrays, opts) -> named list of arrays
```

Where:
- `arrays`: Named list of 3D arrays [samples × subjects × contrasts]
- `opts`: Named list of options
- **Returns:** Named list of output arrays (typically `q` for FDR-adjusted p-values)

### Registered Posthoc Methods

| Method | Status | Requires | Produces | Method Details |
|--------|--------|----------|----------|----------------|
| `fdr:bh` | ✓ Complete | p | q | Benjamini-Hochberg FDR correction via `p.adjust(..., method="BH")` |
| `fdr:by` | ✓ Complete | p | q | Benjamini-Yekutieli FDR correction via `p.adjust(..., method="BY")` |
| `fdr:spatial` | Missing | p, space | q | Spatial FDR accounting for voxel neighborhoods |

**Implementation Details:**

**fdr:bh / fdr:by**
- **Location:** `R/posthoc-registry.R:37-62`
- **Method:** Standard FDR correction applied independently per (subject, contrast) slice
- **Auto-derivation:** If `p` is missing, attempts to derive from `z` (→ pnorm) or `t + df` (→ pt)
- **Loop structure:** `for j in subjects: for k in contrasts: q[,j,k] = p.adjust(p[,j,k])`
- **Tested:** Yes (`tests/testthat/test-posthoc.R`)

### Missing Posthoc Methods

| Name | Status | Priority | Complexity | Notes |
|------|--------|----------|------------|-------|
| `fdr:spatial` | Missing | High | Medium-High | Requires spatial neighbor graph construction from space info |
| `cluster:tfce` | Missing | Medium | High | Threshold-free cluster enhancement; requires spatial integration |
| `permute:maxT` | Missing | Low | Very High | Permutation-based family-wise error control; requires resampling |

**Implementation Path for fdr:spatial:**
- Requires `space` metadata (voxel/parcels)
- Build adjacency graph from `space_voxel$mask_idx` or `space_parcels$labels`
- Use spatial FDR methods (e.g., Li & Barber 2019, or BH with spatial pooling)
- Complexity: ~150 lines R (graph construction + spatial FDR)
- Effort: 1-2 days
- Considerations: Different strategies for voxel vs parcel spaces

---

## 5. Adapter Registry

### API Overview

The adapter registry is implemented in `/Users/bbuchsbaum/code/fmrigds/R/adapter-registry.R`:

- **Storage:** Environment `.adapter_registry` (package-private)
- **Registration:** `register_adapter(name, detect, open, probe, read, close, ...)`
- **Retrieval:** `get_adapter(name)` returns adapter list
- **Auto-detection:** `detect_adapter(source, prefer=NULL)` scores all adapters and picks best

### Adapter Interface

Each adapter provides:
- `detect(source)`: Returns score in [0,1] or FALSE (higher = better match)
- `open(source, mode, ...)`: Returns handle
- `probe(handle, ...)`: Returns metadata list (assays, dims, subjects, contrasts, space, maps, metadata, columns)
- `read(handle, assays, block, ...)`: Returns named list of 3D arrays
- `close(handle)`: Cleanup

### Registered Adapters

| Adapter | Status | Formats | Assays | Space Types | Notes |
|---------|--------|---------|--------|-------------|-------|
| `tabular` | ✓ Complete | CSV, TSV, Parquet | User-defined via `effect_cols` | sample_labels, voxel (if provided) | Read-only; requires `sample`, `subject`, `contrast` columns |
| `h5` | ✓ Complete | HDF5 with `/gds` group | All in `/gds/assays/*` | voxel, parcels, sample_labels | Full GDS schema v0.1.0 |
| `nifti` | ✓ Complete | .nii, .nii.gz (3D/4D) | beta (+ derived se) | voxel | Uses neuroim2; supports mask; single-subject or multi-file |
| `fmristore` | ✓ Complete | HDF5 (legacy fmristore) | beta, var, se (layout-dependent) | voxel, parcels, basis (latent) | Multi-layout adapter (LabeledVolumeSet, LatentNeuroVec, cluster/ROI) |

### Adapter Details

#### tabular
- **Location:** `R/adapter-tabular.R`
- **Detection score:** 0.8 for .csv/.tsv/.parquet
- **Capabilities:**
  - Flexible column mapping via `effect_cols = list(beta="beta", var="var", ...)`
  - Auto-detects subjects, samples, contrasts from columns
  - Handles missing data (NA in arrays)
  - Block selection on sample/subject/contrast axes
- **Limitations:**
  - Read-only
  - Requires well-formed tabular data
  - No spatial structure (unless `space` explicitly provided)
- **Testing:** `tests/testthat/test-adapter-tabular.R`

#### h5 (Native GDS)
- **Location:** `R/adapter-h5.R`
- **Detection score:** 1.0 if `/gds` group exists
- **Schema:**
  ```
  /gds/assays/<name>         [samples × subjects × contrasts]
  /gds/axes/subjects         [subjects]
  /gds/axes/contrasts        [contrasts]
  /gds/axes/subjects_table/  (optional col_data)
  /gds/space/                (type, voxel/parcels/sample_labels metadata)
  /gds/alignments/           (optional map_family objects)
  ```
- **Capabilities:**
  - Full roundtrip support (read + write)
  - Alignments/maps preserved
  - Subject-level covariates in `/gds/axes/subjects_table`
- **Testing:** `tests/testthat/test-h5-adapter.R`

#### nifti
- **Location:** `R/adapter-nifti.R`
- **Detection score:** 0.9 for .nii/.nii.gz files
- **Capabilities:**
  - 3D (single contrast) or 4D (multi-contrast) NIfTI
  - Optional mask file (separate .nii)
  - Affine transformation extracted from header
  - Multi-subject via vector of file paths
  - Subject IDs from filenames
- **Limitations:**
  - Currently only provides `beta` assay
  - Placeholder `se` filled with NA
  - No variance/statistics from NIfTI alone
- **Testing:** `tests/testthat/test-adapters.R`
- **Dependencies:** Requires `neuroim2` package

#### fmristore
- **Location:** `R/adapter-fmristore.R`
- **Detection score:** 1.0 for native `/gds`, 0.95 for legacy layouts
- **Multi-Layout Support:**
  - **Path A (Legacy fmristore layouts):**
    - LabeledVolumeSet: `/data/<label>`, `/mask`, `/labels`, `/header`
    - LatentNeuroVec: `/basis/basis_matrix`, `/scans/<subject>/embedding`
    - Cluster/ROI: `/cluster_map`, `/scans/<subject>` or `/data`
  - **Path B (Native GDS):** Delegates to `h5` adapter if `/gds` present
- **Capabilities:**
  - Automatic layout detection
  - Temporal policy for time-series data: "as_is", "mean", "design" (with contrast_matrix)
  - Sparse voxel→parcel maps for cluster layouts
  - Multi-subject composition from vector of files
  - Assay detection: beta, var, se (suffix-based for LabeledVolumeSet)
- **Testing:** `tests/testthat/test-adapter-fmristore.R`

### Missing Adapters

| Name | Status | Priority | Complexity | Notes |
|------|--------|----------|------------|-------|
| `gifti` | Missing | Low | Medium | Surface-based neuroimaging (.gii); requires geometry handling |
| `zarr` | Missing | Low | Medium | Cloud-native array storage; requires zarr package |
| `rds` | Missing | Very Low | Low | Native R serialization; trivial to implement |

---

## 6. Assay Registry

The assay registry provides semantic metadata about assays:

- **Location:** `R/assay-registry.R`
- **Storage:** Environment `.assay_registry`
- **API:** `register_assay(name, role, units, variance_of, derive_from)`

### Built-in Assays

| Assay | Role | Units | Variance Of | Derive From | Can Map Linear |
|-------|------|-------|-------------|-------------|----------------|
| beta | location | %BOLD | - | - | ✓ |
| var | variance | %BOLD² | beta | - | ✓ |
| se | stdev | %BOLD | - | var | ✓ |
| t | t | - | - | beta, var, df | ✗ |
| z | z | - | - | t, df | ✗ |
| F | F | - | - | beta, var, df1, df2 | ✗ |
| df | df | - | - | - | ✗ |
| df1 | df | - | - | - | ✗ |
| df2 | df | - | - | - | ✗ |
| n_eff | n_eff | - | - | - | ✗ |
| p | p | - | - | - | ✗ |
| chi2 | chi2 | - | - | - | ✗ |
| logBF | log_evidence | log BF | - | - | ✗ |

**Purpose:** Determines which assays can be linearly mapped through spatial transformations (location, variance, stdev) vs. those requiring re-computation (test statistics).

---

## 7. C++ Performance Layer

### Available C++ Kernels

All located in `src/reducers_core.cpp` with OpenMP parallelization:

1. `meta_fe_cpp(beta, var, min_subj, eps, tail)`
2. `meta_re_dl_cpp(beta, var, min_subj, eps, tail)`
3. `meta_fe_reg_cpp(beta, var, X, min_subj, eps)`
4. `meta_re_reg_dl_cpp(beta, var, X, min_subj, eps)`
5. `stouffer_combine_cpp(z, weights, min_subj)`
6. `fisher_combine_cpp(p, min_subj)`
7. `lancaster_combine_cpp(p, dfw, min_subj)`

**Performance:** Parallel across samples (columns) via `#pragma omp parallel for`

### Missing C++ Implementations

| Function | Priority | Complexity | Speedup Potential | Notes |
|----------|----------|------------|-------------------|-------|
| `ols_voxelwise_cpp` | High | Medium | 10-50x | Pure R currently; batched matrix ops would benefit greatly |
| `fdr_spatial_cpp` | Medium | High | 5-10x | Graph operations; adjacency matrix multiplication |

---

## 8. Testing Coverage

### Test Files

- `test-reduce.R`: Reducer execution (FE, RE, Stouffer, Lancaster, alternative tails)
- `test-posthoc.R`: FDR BH/BY posthoc methods
- `test-ols-voxelwise.R`: OLS reducer with covariance extraction, coef_array
- `test-adapter-tabular.R`: Tabular CSV/TSV/Parquet ingestion
- `test-adapter-h5.R`: Native GDS HDF5 roundtrip
- `test-adapters.R`: NIfTI adapter
- `test-adapter-fmristore.R`: Legacy fmristore layouts
- `test-reduce-regression.R`: Meta-regression reducers
- `test-reducer-kernels.R`: Low-level kernel correctness

**Coverage:** All implemented reducers, posthoc methods, and adapters have tests.

---

## 9. Gap Analysis

### High Priority Gaps

1. **C++ ols_voxelwise**
   - Current: Pure R implementation
   - Impact: ~10-50x speedup for large datasets
   - Effort: 1-2 days (matrix ops in Armadillo)

2. **fdr:spatial posthoc method**
   - Current: Missing
   - Impact: Critical for neuroimaging spatial FDR
   - Effort: 1-2 days (graph construction + spatial FDR logic)
   - Considerations: Voxel vs. parcel space handling

3. **meta:re:pm (Paule-Mandel tau²)**
   - Current: Only DL available
   - Impact: Alternative heterogeneity estimator (more robust in some cases)
   - Effort: 2-3 hours (iterative solver)

### Medium Priority Gaps

4. **meta:re:reml**
   - Current: Missing
   - Impact: Gold-standard heterogeneity estimator
   - Effort: 1-2 days (optimization routine)

5. **Variance propagation for cluster:tfce**
   - Current: No cluster-based methods
   - Impact: Threshold-free cluster enhancement for spatial inference
   - Effort: 3-5 days (spatial integration, permutation support)

### Low Priority Gaps

6. **Additional adapters (GIFTI, Zarr)**
   - Current: Adequate coverage for most use cases
   - Impact: Extended format support
   - Effort: 1-2 days each

---

## 10. Summary

### What Works

✓ **8 reducers** fully implemented:
  - Meta-analysis: FE, RE (DL), FE regression, RE regression
  - Evidence combination: Stouffer, Fisher, Lancaster
  - Regression: OLS voxelwise

✓ **2 posthoc methods**: FDR BH, FDR BY

✓ **4 adapters**: tabular (CSV/TSV/Parquet), h5 (native GDS), nifti, fmristore (multi-layout)

✓ **7 C++ kernels** for performance (all but ols_voxelwise)

✓ **Comprehensive testing** for all implemented features

### What's Missing

**Reducers:**
- meta:re:pm (Paule-Mandel)
- meta:re:reml (REML)
- meta:re:ml, meta:re:eb (lower priority)

**Posthoc:**
- fdr:spatial (high priority)
- cluster:tfce (medium priority)
- permute:maxT (low priority)

**C++ Performance:**
- ols_voxelwise_cpp (high impact)

**Adapters:**
- gifti, zarr, rds (low priority)

### Architecture Strengths

1. **Extensible registries**: Easy to add new reducers/posthoc/adapters without core changes
2. **Hybrid R/C++ strategy**: R fallbacks ensure correctness; C++ provides performance
3. **Flexible adapter system**: Auto-detection and multiple format support
4. **Rich metadata**: Assay registry enables intelligent spatial mapping and derivation
5. **Comprehensive testing**: All implemented features have test coverage

### Recommended Next Steps

1. **Immediate (1-2 days):**
   - Implement `ols_voxelwise_cpp` for performance
   - Implement `fdr:spatial` for neuroimaging workflows

2. **Short-term (1 week):**
   - Add `meta:re:pm` for alternative heterogeneity estimation
   - Add `meta:re:reml` for gold-standard meta-analysis

3. **Medium-term (2-4 weeks):**
   - Implement `cluster:tfce` for spatial inference
   - Add GIFTI adapter for surface-based analysis

---

**End of Audit**
