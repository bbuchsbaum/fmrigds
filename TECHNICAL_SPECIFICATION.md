# fmrigds: Complete Technical Specification for R Implementation

**Version:** 0.1.0
**Package Name:** `gdsfmri`
**Target:** Production-ready R implementation
**Testing Framework:** testthat

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Core Architecture](#2-core-architecture)
3. [API Specification](#3-api-specification)
4. [Data Model & Schema](#4-data-model--schema)
5. [Lazy Execution System](#5-lazy-execution-system)
6. [Statistical Foundations](#6-statistical-foundations)
7. [Storage & I/O](#7-storage--io)
8. [Spatial Operations](#8-spatial-operations)
9. [Implementation Modules](#9-implementation-modules)
10. [Test Specifications](#10-test-specifications)
11. [Migration Strategy](#11-migration-strategy)

---

## 1. Executive Summary

### 1.1 Purpose

GDS (Group Data Set) provides a format-agnostic, space-aware, statistically principled container for first-level fMRI results that flows cleanly into group-level and meta-analytic workflows.

**Core Design Principles:**
- Fixed dimensionality: `[sample × subject × contrast]` across all spaces
- Lazy-first execution with plan-based optimization
- Statistical correctness: proper variance propagation, never average standard errors
- Space-agnostic: works with voxels, parcels, surfaces, latent bases
- Minimal dependencies: no Bioconductor requirement

### 1.2 Five Core Nouns

1. **GDS** - Realized group dataset with assays, axes, space, and metadata
2. **Space** - Defines what the sample axis represents (voxel/parcel/surface/basis)
3. **Map** - Linear mapping between spaces (global or per-subject)
4. **UncertaintyRule** - Variance propagation strategy through mappings
5. **Plan** - Lazy pipeline DAG; serializable, auditable, optimizable

### 1.3 Eight Core Verbs

All verbs return a Plan; only `compute()` touches data:

```r
gds()        # Create Plan from adapter (front door)
subset()     # Select samples/subjects/contrasts
derive()     # Compute var, se, t, z from existing assays
align()      # Apply subject-specific transforms to group space
mask()       # Derive/apply group mask policy
map_to()     # Change spaces with uncertainty propagation
reduce()     # Meta-analysis across subjects
write_out()  # Export to file (lazy; executes at compute)
compute()    # Execute optimized plan; return realized GDS
```

### 1.4 Minimal Example

```r
# ROI analysis with factorial contrasts
plan <- gds("roi_stats.csv") %>%
  subset(contrast = c("A_main", "B_main", "A:B")) %>%
  derive(c("var", "t")) %>%
  reduce(method = "fixed", weights = "1/var") %>%
  write_out("group_roi.csv", format = "csv")

result <- compute(plan, sink = "memory")
```

```r
# Voxel analysis with alignment and spatial mapping
plan <- gds(Sys.glob("sub-*/beta.nii.gz")) %>%
  subset(subject = paste0("sub-", sprintf("%02d", 1:30))) %>%
  derive(c("var", "t")) %>%
  align(map = ortho_family) %>%
  mask(policy = MaskPolicy(scope = "group", rule = "threshold", threshold = 0.95)) %>%
  map_to(target_space = parcels, map = vox_to_parcel_map,
         uncertainty = UncertaintyRule("independent")) %>%
  reduce(method = "fixed") %>%
  write_out("results.h5", format = "h5")

result <- compute(plan, sink = "HDF5", path = "results.h5", block = list(sample = 1e5))
```

---

## 2. Core Architecture

### 2.1 GDS Object (Realized)

A realized GDS is an S3 object with class `c("gds", "list")`.

**Required Slots:**

```r
structure(
  list(
    assays = list(             # Named list of 3D arrays
      beta = array(...),       # [sample × subject × contrast]
      var = array(...),        # same dims
      # ... other assays
    ),
    space = <Space>,           # S3 space object
    subjects = character(),    # Length dim(assays[[1]])[2]
    contrasts = character(),   # Length dim(assays[[1]])[3]
    col_data = data.frame(),   # Rows = subjects (optional)
    row_data = data.frame(),   # Rows = samples (optional)
      metadata = list(           # Schema, units, provenance, etc.
      schema_version = "0.1.0",
      units = list(beta = "% BOLD", var = "% BOLD^2"),
      provenance = list(
        graph = list(),        # DAG nodes
        log = character()      # Human-readable log
      )
    )
  ),
  class = c("gds", "group_data")
)
```

**Constructor:**

```r
new_gds <- function(assays, space, subjects, contrasts,
                    col_data = NULL, row_data = NULL,
                    metadata = list()) {
  # Validation (see §2.6)
  validate_gds_assays(assays)
  validate_gds_dims(assays, subjects, contrasts)
  validate_gds_space(space, nrow(assays[[1]]))

  structure(
    list(
      assays = assays,
      space = space,
      subjects = subjects,
      contrasts = contrasts,
      col_data = col_data %||% data.frame(subject = subjects),
      row_data = row_data %||% data.frame(sample = seq_len(dim(assays[[1]])[1])),
      metadata = c(list(schema_version = "0.1.0"), metadata)
    ),
    class = c("gds", "group_data")
  )
}
```

**Accessors:**

```r
# Extract single assay
assay <- function(x, name = "beta") UseMethod("assay")
assay.gds <- function(x, name = "beta") x$assays[[name]]

# Extract all assays
assays <- function(x) UseMethod("assays")
assays.gds <- function(x) x$assays

# Extract space
space <- function(x) UseMethod("space")
space.gds <- function(x) x$space

# Extract subjects
subjects <- function(x) UseMethod("subjects")
subjects.gds <- function(x) x$subjects

# Extract contrasts
contrasts <- function(x) UseMethod("contrasts")
contrasts.gds <- function(x) x$contrasts

# Extract column data (subject covariates)
col_data <- function(x) UseMethod("col_data")
col_data.gds <- function(x) x$col_data

# Extract row data (sample metadata)
row_data <- function(x) UseMethod("row_data")
row_data.gds <- function(x) x$row_data

# Extract metadata
metadata <- function(x) UseMethod("metadata")
metadata.gds <- function(x) x$metadata
```

### 2.2 Space Objects

All spaces inherit from `gds_space` and define what a "sample" represents.

#### 2.2.1 Voxel Space

```r
space_voxel <- function(dim, affine,
                        mask_bitmap = NULL,
                        mask_idx = NULL,
                        storage = c("dense", "packed"),
                        template_id = NULL) {
  storage <- match.arg(storage)

  # Derive mask_idx from bitmap if needed
  if (!is.null(mask_bitmap) && is.null(mask_idx)) {
    mask_idx <- which(as.vector(mask_bitmap))
  }

  # Validate
  stopifnot(
    length(dim) == 3,
    all(dim > 0),
    is.matrix(affine),
    dim(affine) == c(4, 4)
  )

  if (!is.null(mask_idx)) {
    stopifnot(all(mask_idx > 0), all(mask_idx <= prod(dim)))
  }

  structure(
    list(
      type = "voxel",
      dim = as.integer(dim),
      affine = affine,
      mask_bitmap = mask_bitmap,
      mask_idx = as.integer(mask_idx),
      storage = storage,
      template_id = template_id
    ),
    class = c("space_voxel", "gds_space")
  )
}
```

**Aliases:** `space_voxels()` (same as `space_voxel()`)

**Storage Modes:**
- `dense`: Assays are full `[X*Y*Z × subject × contrast]` (sample axis = all voxels)
- `packed`: Assays store only `length(mask_idx)` samples (sample axis = masked voxels only)

**Dimensional Contract:**
- If `storage = "dense"`: `nrow(assays) == prod(dim)`
- If `storage = "packed"`: `nrow(assays) == length(mask_idx)`

#### 2.2.2 Parcel Space

```r
space_parcels <- function(labels, lookup = NULL, membership = NULL) {
  stopifnot(is.character(labels), length(labels) > 0)

  structure(
    list(
      type = "parcels",
      labels = labels,              # Character vector of ROI/parcel names
      lookup = lookup,              # Optional: named list or data.frame with metadata
      membership = membership       # Optional: list or sparse matrix [samples × voxels]
    ),
    class = c("space_parcels", "gds_space")
  )
}
```

**membership** encoding (for voxel → parcel mapping):
- **List:** `membership[[i]]` = integer vector of voxel indices for parcel i
- **Sparse matrix:** `[parcels × voxels]` with weights (often binary)

#### 2.2.3 Surface Space

```r
space_surface <- function(vertices, faces, hemi, template_id = NULL) {
  stopifnot(
    is.matrix(vertices),
    ncol(vertices) == 3,
    is.matrix(faces),
    hemi %in% c("L", "R", "LR")
  )

  structure(
    list(
      type = "surface",
      vertices = vertices,          # [n_vertices × 3] (x, y, z)
      faces = faces,                # [n_faces × 3] (vertex indices)
      hemi = hemi,                  # Hemisphere
      template_id = template_id     # e.g., "fsaverage", "fslr32k"
    ),
    class = c("space_surface", "gds_space")
  )
}
```

#### 2.2.4 Basis Space

```r
space_basis <- function(k, basis_name = NULL, projector = NULL,
                       voxel_space = NULL) {
  stopifnot(k > 0, is.numeric(k))

  structure(
    list(
      type = "basis",
      k = as.integer(k),            # Number of components
      basis_name = basis_name,      # e.g., "PCA", "ICA", "NMF"
      projector = projector,        # Matrix [k × voxels] or function
      voxel_space = voxel_space     # Optional reference space
    ),
    class = c("space_basis", "gds_space")
  )
}
```

**projector** encoding (normative):
- Operator matrices are `[target × source]`.
- Matrix form: store `B` as `[k × n_voxels]` (basis rows over voxels)
  - forward to basis: `y_basis = B %*% x_voxel`
  - backward to voxel: `x_voxel ≈ t(B) %*% y_basis`
- Function: custom projector mapping between basis and voxel space

### 2.3 Map Objects (Linear Mappings)

#### 2.3.1 map_linear

```r
map_linear <- function(from_space, to_space, operator,
                      by_subject = NULL,
                      traits = list(orthogonal = FALSE, mass_preserving = FALSE),
                      uncertainty = UncertaintyRule("independent"),
                      name = NULL) {
  # Validate
  stopifnot(
    inherits(from_space, "gds_space") || is.character(from_space),
    inherits(to_space, "gds_space") || is.character(to_space)
  )

  structure(
    list(
      name = name %||% paste0(from_space$type, "_to_", to_space$type),
      from = from_space,
      to = to_space,
      operator = operator,          # Matrix or function (global)
      by_subject = by_subject,      # Named list of subject-specific operators
      traits = traits,
      uncertainty = uncertainty
    ),
    class = c("map_linear", "gds_map")
  )
}
```

**Operator Types:**

1. **Global matrix:** Single sparse or dense matrix `[n_target × n_source]`
2. **Subject-specific matrices:** `by_subject` is a named list; each element is a descriptor:
   ```r
   list(
     kind = "matrix",            # "matrix", "affine", "warp", "transport"
     dim = c(n_target, n_source),
     storage = "dense",          # "dense", "csr", "csc"
     data = ...                  # Actual matrix or serialized form
   )
   ```
3. **Function:** Custom mapping function with signature `function(x, ...) { ... }`

**Traits:**
- `orthogonal = TRUE`: Indicates isometry (Euclidean norm preservation)
- `mass_preserving = TRUE`: Row sums = 1 (e.g., OT/Wasserstein)

#### 2.3.2 MapFamily (Subject-Aware Transform Sets)

```r
MapFamily <- function(name, from_space, to_space,
                     type = c("linear", "orthogonal", "affine3d", "deform3d", "ot"),
                     by_subject,
                     traits = list(orthogonal = FALSE, mass_preserving = FALSE),
                     uncertainty = UncertaintyRule("independent")) {
  type <- match.arg(type)

  structure(
    list(
      name = name,
      from = from_space,
      to = to_space,
      type = type,
      by_subject = by_subject,      # Named list: subject_id -> operator descriptor
      traits = traits,
      uncertainty = uncertainty
    ),
    class = c("gds_map_family", "gds_map")
  )
}
```

**Common Families (Helper Constructors):**

```r
# Orthogonal/Procrustes alignment
OrthogonalFamily <- function(name, from_space, to_space, matrices_by_subject,
                            uncertainty = UncertaintyRule("independent")) {
  MapFamily(name, from_space, to_space, type = "orthogonal",
           by_subject = matrices_by_subject,
           traits = list(orthogonal = TRUE, mass_preserving = FALSE),
           uncertainty = uncertainty)
}

# Optimal Transport alignment
OTFamily <- function(name, from_space, to_space, plans_by_subject,
                    uncertainty = UncertaintyRule("independent")) {
  MapFamily(name, from_space, to_space, type = "ot",
           by_subject = plans_by_subject,
           traits = list(orthogonal = FALSE, mass_preserving = TRUE),
           uncertainty = uncertainty)
}

# Affine/deformable warps
WarpFamily <- function(name, from_space, to_space, warps_by_subject,
                      uncertainty = UncertaintyRule("independent")) {
  MapFamily(name, from_space, to_space, type = "deform3d",
           by_subject = warps_by_subject,
           traits = list(orthogonal = FALSE, mass_preserving = FALSE),
           uncertainty = uncertainty)
}
```

### 2.4 UncertaintyRule

```r
UncertaintyRule <- function(mode = c("independent", "cov_provider", "kernel", "none"),
                           df_rule = c("satterthwaite", "none"),
                           cov_provider = NULL,
                           kernel = NULL) {
  mode <- match.arg(mode)
  df_rule <- match.arg(df_rule)

  # Validate
  if (mode == "cov_provider" && is.null(cov_provider)) {
    stop("cov_provider must be specified when mode = 'cov_provider'")
  }
  if (mode == "kernel" && is.null(kernel)) {
    stop("kernel must be specified when mode = 'kernel'")
  }

  structure(
    list(
      mode = mode,
      df_rule = df_rule,
      cov_provider = cov_provider,
      kernel = kernel
    ),
    class = "gds_uncertainty_rule"
  )
}
```

**Modes:**

1. **independent:** Treat variance as diagonal; fast path
   - Formula: `Var(y) = (M ⊙ M) %*% Var(x)`

2. **cov_provider:** User-supplied covariance function
   - Signature: `cov_provider(indices) -> Sigma_block`
   - Applied per target row: `Var(y[i]) = t(w) %*% Sigma %*% w` where `w = M[i,]`

3. **kernel:** Build covariance from spatial correlation model
   - Example: `kernel(dist) = exp(-dist/length_scale)`

4. **none:** Skip variance propagation (discouraged; for debugging only)

**DF Rules:**

1. **satterthwaite:** Approximate combined df using Welch-Satterthwaite:
   ```
   df' ≈ (Σ a_i² s_i²)² / Σ(a_i⁴ s_i⁴ / df_i)
   ```

2. **none:** Do not compute combined df

### 2.5 Plan Object (Lazy Pipeline)

```r
gds_plan <- function(source, nodes = list(), meta = list()) {
  structure(
    list(
      source = source,              # gds_source object
      nodes = nodes,                # List of op_* nodes
      meta = meta                   # Metadata: version, creation time, etc.
    ),
    class = c("gds_plan", "list")
  )
}
```

**Source Object:**

```r
gds_source <- function(adapter_name, source_spec, probe_result = NULL) {
  structure(
    list(
      adapter = adapter_name,       # Adapter identifier
      source = source_spec,         # Path, list of paths, or connection
      probe = probe_result,         # Cached probe result (optional)
      hash = digest::digest(list(adapter_name, source_spec))
    ),
    class = c("gds_source", "list")
  )
}
```

**Operation Nodes:**

All nodes are lists with `op` field specifying the operation type.

```r
# Subset operation
op_subset_axis <- function(sample = NULL, subject = NULL, contrast = NULL) {
  list(op = "subset_axis", sample = sample, subject = subject, contrast = contrast)
}

# Derivation operation
op_derive <- function(what, options = list()) {
  list(op = "derive", what = as.character(what), options = options)
}

# Mapping operation
op_map <- function(target_space, map, uncertainty, combine = NULL) {
  list(op = "map", target_space = target_space, map = map,
       uncertainty = uncertainty, combine = combine)
}

# Alignment operation
op_align_to_group <- function(family) {
  list(op = "align_to_group", family = family)
}

# Mask operation
op_mask_policy <- function(policy) {
  list(op = "mask_policy", policy = policy)
}

# Reduction operation
op_reduce <- function(method, weights, by) {
  list(op = "reduce", method = method, weights = weights, by = by)
}

# Write operation
op_write <- function(path, format, options) {
  list(op = "write", path = path, format = format, options = options)
}
```

**Plan Manipulation:**

```r
# Add operation to plan
add_op <- function(plan, node) {
  stopifnot(inherits(plan, "gds_plan"))
  plan$nodes <- append(plan$nodes, list(node))
  plan
}

# Convert GDS or source to plan
as_plan <- function(x) {
  if (inherits(x, "gds_plan")) return(x)
  if (inherits(x, "gds_source")) return(gds_plan(source = x))
  if (inherits(x, "gds")) {
    # Convert realized GDS to in-memory source
    src <- gds_source("memory", x)
    return(gds_plan(source = src))
  }
  stop("Cannot convert to plan: ", class(x))
}
```

### 2.6 Validation Functions

```r
# Validate assay list
validate_gds_assays <- function(assays) {
  stopifnot(
    is.list(assays),
    length(assays) > 0,
    !is.null(names(assays)),
    all(nzchar(names(assays)))
  )

  # Check dimensions
  dims <- lapply(assays, dim)
  ref_dim <- dims[[1]]
  stopifnot(all(sapply(dims, identical, ref_dim)))

  # Check for required assays
  if ("beta" %in% names(assays)) {
    has_uncertainty <- any(c("var", "se") %in% names(assays))
    if (!has_uncertainty) {
      stop("If 'beta' is present, must have either 'var' or 'se'")
    }
  }

  if ("t" %in% names(assays) && !"df" %in% names(assays)) {
    warning("'t' present without 'df'; some operations may fail")
  }

  invisible(TRUE)
}

# Validate dimensions against axes
validate_gds_dims <- function(assays, subjects, contrasts) {
  dims <- dim(assays[[1]])
  stopifnot(
    length(dims) == 3,
    dims[2] == length(subjects),
    dims[3] == length(contrasts)
  )
  invisible(TRUE)
}

# Validate space against data
validate_gds_space <- function(space, n_samples) {
  if (inherits(space, "space_voxel")) {
    if (space$storage == "dense") {
      expected <- prod(space$dim)
    } else {
      expected <- length(space$mask_idx)
    }
    stopifnot(n_samples == expected)
  }
  invisible(TRUE)
}
```

---

## 3. API Specification

### 3.1 gds() - Front Door Constructor

```r
gds <- function(source, format = c("auto", "nifti", "h5", "csv", "parquet", "fmristore"),
                prefer = NULL, ...) {
  format <- match.arg(format)

  # Detect adapter
  adapter_name <- if (format == "auto") {
    detect_adapter(source, prefer = prefer)
  } else {
    format
  }

  # Get adapter
  adapter <- get_adapter(adapter_name)

  # Open and probe
  handle <- adapter$open(source, ...)
  probe_result <- adapter$probe(handle)

  # Create source
  src <- gds_source(adapter_name, source, probe_result)

  # Return plan
  gds_plan(source = src)
}
```

**Parameters:**
- `source`: File path(s), directory, or connection
- `format`: Force specific adapter or auto-detect
- `prefer`: Preferred adapter if multiple match
- `...`: Adapter-specific options

**Returns:** `gds_plan` object

**Behavior:**
- Does NOT load data into memory
- Probes source to validate structure and extract metadata
- Returns lazy plan that can be composed with verbs

### 3.2 subset() - Select Data Slices

```r
subset.gds_plan <- function(x, sample = NULL, subject = NULL, contrast = NULL, ...) {
  add_op(x, op_subset_axis(sample, subject, contrast))
}
```

**Parameters:**
- `sample`: Integer indices, logical vector, or character names (for parcels/basis)
- `subject`: Character subject IDs or integer indices
- `contrast`: Character contrast names or integer indices

**Returns:** Modified `gds_plan`

**Dimensional Contract:**
- Input: `[S × J × K]`
- Output: `[S' × J' × K']` where `S' ≤ S`, `J' ≤ J`, `K' ≤ K`

**Optimizer Behavior:**
- Pushed down to adapter when possible (avoids loading unused data)
- Multiple subsets are coalesced

### 3.3 derive() - Compute Derived Assays

```r
derive <- function(x, what = c("var", "se", "t", "z"), options = list()) {
  plan <- as_plan(x)
  add_op(plan, op_derive(what, options))
}
```

**Parameters:**
- `what`: Character vector of assays to derive
- `options`: List of derivation options (e.g., df aggregation method)

**Returns:** Modified `gds_plan`

**Derivation Rules:**

| From | To | Formula |
|------|-----|---------|
| `se` | `var` | `var = se^2` |
| `var` | `se` | `se = sqrt(var)` |
| `beta, var` | `t` | `t = beta / se` (requires `df` separately for inference) |
| `t, df` | `z` | `z = qnorm(1 - (2*pt(-abs(t), df))/2) * sign(t)` |
| `z, df` | `t` | `t = sign(z) * qt(pnorm(abs(z)), df)` |
| `F, df1, df2` | `p` | `p = pf(F, df1, df2, lower.tail = FALSE)` |
| `p` | `z` | `z = qnorm(1 - p/2)` for two-sided p (use sign separately if available) |

**Overwrite Behavior:**
- By default, does NOT overwrite existing assays
- Set `options$overwrite = TRUE` to force recomputation

### 3.4 align() - Apply Subject-Specific Transforms

```r
align <- function(x, map, name = NULL) {
  plan <- as_plan(x)

  # Handle MapFamily or registered map name
  family <- if (inherits(map, "gds_map_family")) {
    map
  } else if (is.character(map)) {
    # Look up in plan meta space graph (populated by adapters or register_map())
    plan$meta$space_graph$maps[[map]]
  } else {
    stop("map must be a MapFamily or character name")
  }

  add_op(plan, op_align_to_group(family))
}
```

**Parameters:**
- `map`: `MapFamily` object or registered map name
- `name`: Optional override name

**Returns:** Modified `gds_plan`

**Execution Semantics:**
- For each subject `s` and contrast `k`:
  - `beta'[·, s, k] = M_s %*% beta[·, s, k]`
  - `var'[·, s, k]` = propagated via `family$uncertainty` rule
- After alignment, space is updated to `family$to`
- Recomputes `se`, `t`, `z` if they exist (never map test statistics directly)

**Statistical Contract:**
- MUST have `beta` and `var` (or `se`)
- Refuses to operate on `t/z/F` only (no effect scale)

### 3.5 mask() - Derive and Apply Group Mask

```r
mask <- function(x, policy = MaskPolicy(scope = "group", rule = "intersection", threshold = 0.9)) {
  plan <- as_plan(x)
  add_op(plan, op_mask_policy(policy))
}
```

**MaskPolicy Constructor:**

```r
MaskPolicy <- function(scope = c("subject", "group"),
                      rule = c("intersection", "union", "threshold", "custom"),
                      threshold = 0.95,
                      custom = NULL) {
  scope <- match.arg(scope)
  rule <- match.arg(rule)

  structure(
    list(scope = scope, rule = rule, threshold = threshold, custom = custom),
    class = "gds_mask_policy"
  )
}
```

**Parameters:**
- `scope`: Apply to individual subjects or group level
- `rule`: Mask computation strategy
  - `intersection`: Voxels present in all subjects
  - `union`: Voxels present in any subject
  - `threshold`: Voxels present in ≥ threshold fraction of subjects
  - `custom`: User-supplied function
- `threshold`: Fraction (0-1) for threshold rule
- `custom`: Function with signature `function(mask_list) -> logical_vector`

**Returns:** Modified `gds_plan`

**Execution Behavior:**
- Computes analysis mask from per-subject masks
- Updates `space$mask_idx`
- May repack assays to `storage = "packed"` for efficiency
- Dense storage: samples not in the analysis mask become NA (index preserved)
- Packed storage: samples outside the mask are removed; index mapping retained via `space$mask_idx`

**Optimizer Pushdown:**
- Executed after `align()` but before `map_to()` to minimize data movement

### 3.6 map_to() - Change Spaces with Uncertainty Propagation

```r
map_to <- function(x, target_space, map,
                  uncertainty = UncertaintyRule("independent"),
                  combine = NULL) {
  plan <- as_plan(x)
  add_op(plan, op_map(target_space, map, uncertainty, combine))
}
```

**Parameters:**
- `target_space`: Target `gds_space` object or name
- `map`: `map_linear` object or matrix
- `uncertainty`: `UncertaintyRule` for variance propagation
- `combine`: Optional combiner for z-scores (if no effect scale)

**Returns:** Modified `gds_plan`

**Dimensional Contract:**
- Input: `[S_source × J × K]` in `source_space`
- Map: `M` is `[S_target × S_source]`
- Output: `[S_target × J × K]` in `target_space`

**Assay-Aware Behavior:**

1. **Has `beta` and `var`:**
   - Map location: `beta' = M %*% beta`
   - Propagate variance via `uncertainty` rule
   - Derive `se' = sqrt(var')`
   - Recompute `t' = beta' / se'` if `t` and `df` exist
   - Recompute `z'` from `t'` if `z` exists

2. **Has only `t/z/F` (no effect scale):**
  - Refuse linear mapping unless `combine` is specified (mapping test stats is invalid)
  - If `combine = "stouffer"`: Apply Stouffer's method per target sample (assumes independence)
    ```
    z'[i] = (Σ_j w[i,j] z[j]) / sqrt(Σ_j w[i,j]^2)
    ```
  - If `combine = "fisher"`: Apply Fisher's method per target sample (assumes independence)
    ```
    p_j = 2 * pnorm(-abs(z_j))     # or from t with df: p_j = 2*pt(-abs(t_j), df_j)
    chi2'[i] = -2 Σ_j log(p_j)
    ```
    Note: Fisher’s combination has no natural weighting; prefer Stouffer for weighted evidence.

**Uncertainty Propagation (Independent Mode):**

```r
# Pseudocode for independent propagation
.propagate_independent <- function(M, beta_block, var_block) {
  beta_out <- M %*% beta_block
  var_out <- (M^2) %*% var_block  # Element-wise square of M, then matrix multiply
  list(beta = beta_out, var = var_out)
}
```

**Uncertainty Propagation (Covariance Mode):**

```r
# Pseudocode for covariance provider
.propagate_cov_provider <- function(M, beta_block, var_block, cov_provider) {
  beta_out <- M %*% beta_block
  var_out <- numeric(nrow(M))

  for (i in seq_len(nrow(M))) {
    idx <- which(M[i, ] != 0)
    w <- M[i, idx]
    Sigma_block <- cov_provider(idx)
    var_out[i] <- drop(t(w) %*% Sigma_block %*% w)
  }

  list(beta = beta_out, var = var_out)
}
```

### 3.9 Eager Wrappers (Convenience)

For interactive workflows, eager variants wrap lazy verbs and immediately call `compute()`.

```r
subset_eager <- function(x, ...) compute(subset(x, ...))
derive_eager <- function(x, ...) compute(derive(x, ...))
map_to_eager <- function(x, ...) compute(map_to(x, ...))
```

Use lazy-first in production pipelines; eager wrappers are optional conveniences.

### 3.7 reduce() - Meta-Analysis Across Subjects

```r
reduce <- function(x, method = c("fixed", "random", "stouffer", "fisher"),
                  weights = c("1/var", "equal", "n_eff"),
                  by = "contrast") {
  plan <- as_plan(x)
  method <- match.arg(method)
  weights <- match.arg(weights)

  add_op(plan, op_reduce(method, weights, by))
}
```

**Parameters:**
- `method`: Meta-analysis method
  - `fixed`: Fixed-effects (inverse-variance weighted mean)
  - `random`: Random-effects (DerSimonian-Laird, future)
  - `stouffer`: Stouffer's z-score combination
  - `fisher`: Fisher's combined probability
- `weights`: Weighting scheme
- `by`: Grouping variable (default: within each contrast)

**Returns:** Modified `gds_plan`

**Dimensional Contract:**
- Input: `[S × J × K]`
- Output: `[S × 1 × K]` (or `[S × 1 × 1]` if aggregating contrasts)

**Fixed-Effects Formula:**

```r
# Weights
w[i,j,k] = 1 / var[i,j,k]
W[i,k] = Σ_j w[i,j,k]

# Pooled effect
beta_pooled[i,k] = (Σ_j w[i,j,k] * beta[i,j,k]) / W[i,k]

# Pooled variance
var_pooled[i,k] = 1 / W[i,k]

# Pooled SE and t-statistic
se_pooled[i,k] = sqrt(var_pooled[i,k])
t_pooled[i,k] = beta_pooled[i,k] / se_pooled[i,k]
```

**Stouffer's Method (for z-scores):**

```r
# Equal weights
z_combined[i,k] = Σ_j z[i,j,k] / sqrt(J)

# Weighted (e.g., by n_eff)
z_combined[i,k] = (Σ_j w[i,j,k] * z[i,j,k]) / sqrt(Σ_j w[i,j,k]^2)
```

**Fisher's Method:**

```r
# Convert z to p-values
p[i,j,k] = 2 * pnorm(-abs(z[i,j,k]))

# Combine
chi2[i,k] = -2 * Σ_j log(p[i,j,k])
df_combined[i,k] = 2 * J

# Convert back to p
p_combined[i,k] = pchisq(chi2[i,k], df_combined[i,k], lower.tail = FALSE)
z_combined[i,k] = qnorm(p_combined[i,k] / 2, lower.tail = FALSE) * sign(mean(z[i,,k]))
```

### 3.8 write_out() - Export to File (Lazy)

```r
write_out <- function(x, path, format = c("h5", "nifti", "csv", "parquet"),
                     options = list()) {
  plan <- as_plan(x)
  format <- match.arg(format)

  add_op(plan, op_write(path, format, options))
}
```

**Parameters:**
- `path`: Output file path
- `format`: Output format
  - `h5`: HDF5 with GDS layout (`options$overwrite`, etc.)
  - `nifti`: NIfTI-1 (requires voxel space or mapping; `options$stat` selects assay)
  - `csv`: CSV (tabular; good for parcels/ROI; `options$stats`, `options$drop_na`)
  - `parquet`: Parquet (columnar; good for large tables; requires `arrow`)
- `options`: Format-specific options (examples above)

**Returns:** Modified `gds_plan`

**Execution Behavior:**
- Does NOT write until `compute()` is called
- Can be final sink for `compute()` or intermediate output
- Provenance is written to output metadata and persisted in `/gds/provenance` (HDF5)
- For non-HDF5 formats, exports are derived from the realized `gds` returned by `compute()`

### 3.9 compute() - Execute Plan

```r
compute <- function(plan,
                   sink = c("memory", "HDF5"),
                   path = NULL,
                   block = list(sample = 100000),
                   scheduler = c("sequential", "multicore", "future"),
                   cache = TRUE,
                   verbose = FALSE) {
  sink <- match.arg(sink)
  scheduler <- match.arg(scheduler)

  # Optimize plan
  plan_opt <- optimize_plan(plan)

  if (verbose) {
    message("Optimized plan:")
    print(explain(plan_opt))
  }

  # Initialize executor
  executor <- init_executor(plan_opt, sink, path, block, scheduler, cache)

  # Execute
  result <- execute_plan(executor, plan_opt)

  # Return realized GDS
  result
}
```

**Parameters:**
- `sink`: Where to materialize results
  - `memory`: In-memory R arrays
  - `HDF5`: On-disk HDF5 with delayed access (use `write_out()` for CSV/Parquet/NIfTI)
- `path`: Output path (required for `HDF5` sink)
- `block`: Block size for streaming (list with axis names); forwarded to adapter `read()` for chunked processing
- `scheduler`: Parallelization strategy
- `cache`: Enable plan-level caching
- `verbose`: Print execution details

**Returns:** Realized `gds` object

**Execution Flow:**

1. **Optimize:** Apply rewrite rules to plan
2. **Open:** Open source handle via adapter
3. **Stream:** Process data in blocks
   - Read block from source
   - Apply fused operations (subset → derive → map → reduce)
   - Write block to sink
4. **Finalize:** Close handles, write provenance, return GDS

**Block Streaming:**

```r
# Pseudocode for block execution
execute_stream <- function(source, plan_nodes, sink, block_size) {
  n_samples <- source$probe$dims["sample"]
  block_starts <- seq(1, n_samples, by = block_size)

  for (start in block_starts) {
    end <- min(start + block_size - 1, n_samples)

    # Read block
    block_data <- source$adapter$read_assay(source$handle,
                                            assays = c("beta", "var"),
                                            i = start:end,
                                            s = NULL,
                                            k = NULL)

    # Apply operations
    for (node in plan_nodes) {
      block_data <- apply_op(node, block_data)
    }

    # Write block
    sink$write_assay(sink$handle, block_data, i = start:end)
  }
}
```

---

## 4. Data Model & Schema

### 4.1 Assay Registry

Assays have semantic roles that determine how they're mapped and combined.

```r
# Global assay registry
.assay_registry <- new.env(parent = emptyenv())

# Register assay
register_assay <- function(name, role, units = NULL, variance_of = NULL,
                          derive_from = NULL) {
  role <- match.arg(role, c("location", "variance", "stdev", "z", "t", "F",
                           "df", "n_eff", "p", "chi2"))

  .assay_registry[[name]] <- list(
    name = name,
    role = role,
    units = units,
    variance_of = variance_of,
    derive_from = derive_from
  )

  invisible(NULL)
}

# Get assay info
assay_info <- function(name) {
  .assay_registry[[name]]
}

# Check if assay can be linearly mapped
can_map_linear <- function(name) {
  info <- assay_info(name)
  if (is.null(info)) return(FALSE)
  info$role %in% c("location", "variance", "stdev")
}
```

**Default Registrations:**

```r
# During package load
.onLoad <- function(libname, pkgname) {
  register_assay("beta", role = "location", units = "% BOLD")
  register_assay("var", role = "variance", units = "% BOLD^2", variance_of = "beta")
  register_assay("se", role = "stdev", units = "% BOLD", derive_from = "var")
  register_assay("t", role = "t", derive_from = c("beta", "var", "df"))
  register_assay("z", role = "z", derive_from = c("t", "df"))
  register_assay("F", role = "F", derive_from = c("beta", "var", "df1", "df2"))
  register_assay("df", role = "df")
  register_assay("df1", role = "df")
  register_assay("df2", role = "df")
  register_assay("n_eff", role = "n_eff")
  register_assay("p", role = "p")
  register_assay("chi2", role = "chi2")
}
```

### 4.2 Metadata Structure

```r
# Standard metadata schema
gds_metadata <- function(schema_version = "0.1.0",
                        units = list(),
                        provenance = list(graph = list(), log = character()),
                        software = list(
                          package = "gdsfmri",
                          version = packageVersion("gdsfmri"),
                          R_version = getRversion()
                        ),
                        alignment = NULL,
                        mask_info = NULL,
                        contrast_info = NULL,
                        design_mats = NULL,
                        notes = NULL) {
  list(
    schema_version = schema_version,
    units = units,
    provenance = provenance,
    software = software,
    alignment = alignment,
    mask_info = mask_info,
    contrast_info = contrast_info,
    design_mats = design_mats,
    notes = notes,
    created = Sys.time()
  )
}
```

### 4.3 Contrast Info (Factorial Designs)

```r
contrast_info <- function(contrasts,
                         design_spec = NULL,
                         factors = NULL,
                         C = NULL,
                         type = NULL,
                         depends_on = NULL,
                         scale = NULL,
                         notes = NULL) {
  structure(
    list(
      contrasts = contrasts,
      design_spec = design_spec,      # Formula or list
      factors = factors,              # Named list: factor -> levels
      C = C,                          # Contrast matrix [K × nparams]
      type = type,                    # "main", "interaction", etc.
      depends_on = depends_on,        # List: contrast -> parent contrasts
      scale = scale,                  # Per-contrast units
      notes = notes
    ),
    class = "gds_contrast_info"
  )
}
```

**Example (2×2 Factorial):**

```r
cinfo <- contrast_info(
  contrasts = c("A_main", "B_main", "A:B"),
  design_spec = ~ A * B,
  factors = list(A = c("a1", "a2"), B = c("b1", "b2")),
  type = c("main", "main", "interaction"),
  depends_on = list("A:B" = c("A_main", "B_main"))
)
```

### 4.4 Provenance Graph

```r
# Provenance node
provenance_node <- function(op_name, params, inputs, timestamp = Sys.time(),
                           software = list(), hash = NULL) {
  structure(
    list(
      op = op_name,
      params = params,                # Named list (canonicalized)
      inputs = inputs,                # Input node IDs or hashes
      timestamp = timestamp,
      software = software,            # Package versions
      hash = hash %||% digest::digest(list(op_name, params, inputs))
    ),
    class = "gds_provenance_node"
  )
}

# Add node to provenance graph
add_provenance_node <- function(metadata, op_name, params, inputs = list()) {
  node <- provenance_node(
    op_name,
    params,
    inputs,
    software = list(
      package = "gdsfmri",
      version = packageVersion("gdsfmri")
    )
  )

  metadata$provenance$graph <- append(metadata$provenance$graph, list(node))

  # Human-readable log
  log_entry <- sprintf("[%s] %s(%s)",
                      format(node$timestamp, "%Y-%m-%d %H:%M:%S"),
                      op_name,
                      paste(names(params), params, sep = "=", collapse = ", "))
  metadata$provenance$log <- c(metadata$provenance$log, log_entry)

  metadata
}
```

---

## 5. Lazy Execution System

### 5.1 Plan Optimization

The optimizer rewrites the plan to minimize data movement and computation.

```r
optimize_plan <- function(plan) {
  plan <- pushdown_subset(plan)
  plan <- coalesce_derives(plan)
  plan <- push_mask_early(plan)
  plan <- fuse_map_reduce(plan)
  plan
}
```

#### 5.1.1 Pushdown Subset

Push subset operations to the source adapter when possible.

```r
pushdown_subset <- function(plan) {
  # Find first subset_axis operation
  subset_idx <- which(sapply(plan$nodes, function(n) n$op == "subset_axis"))

  if (length(subset_idx) == 0) return(plan)

  # Collect all subset operations before first non-pushable op
  subsets <- list()
  i <- 1
  while (i <= length(plan$nodes) && plan$nodes[[i]]$op == "subset_axis") {
    subsets <- append(subsets, list(plan$nodes[[i]]))
    i <- i + 1
  }

  if (length(subsets) == 0) return(plan)

  # Merge subsets
  merged <- merge_subsets(subsets)

  # Update source probe with subset info
  plan$source$probe$subset <- merged

  # Remove subset nodes from plan
  plan$nodes <- plan$nodes[-seq_len(length(subsets))]

  plan
}

merge_subsets <- function(subsets) {
  sample <- subject <- contrast <- NULL

  for (sub in subsets) {
    sample <- intersect_index(sample, sub$sample)
    subject <- intersect_index(subject, sub$subject)
    contrast <- intersect_index(contrast, sub$contrast)
  }

  list(sample = sample, subject = subject, contrast = contrast)
}
```

#### 5.1.2 Coalesce Derives

Merge consecutive derive operations into a single pass.

```r
coalesce_derives <- function(plan) {
  i <- 1
  while (i < length(plan$nodes)) {
    if (plan$nodes[[i]]$op == "derive" && plan$nodes[[i+1]]$op == "derive") {
      # Merge
      what_merged <- unique(c(plan$nodes[[i]]$what, plan$nodes[[i+1]]$what))
      options_merged <- c(plan$nodes[[i]]$options, plan$nodes[[i+1]]$options)

      plan$nodes[[i]]$what <- what_merged
      plan$nodes[[i]]$options <- options_merged

      # Remove next node
      plan$nodes <- plan$nodes[-(i+1)]
    } else {
      i <- i + 1
    }
  }

  plan
}
```

#### 5.1.3 Push Mask Early

Execute mask operation after alignment but before heavy operations.

```r
push_mask_early <- function(plan) {
  mask_idx <- which(sapply(plan$nodes, function(n) n$op == "mask_policy"))
  if (length(mask_idx) == 0) return(plan)

  # Find alignment operation
  align_idx <- which(sapply(plan$nodes, function(n) n$op == "align_to_group"))

  if (length(align_idx) > 0) {
    # Reorder: alignment -> mask -> everything else
    target_pos <- max(align_idx) + 1
    if (mask_idx[1] > target_pos) {
      mask_node <- plan$nodes[[mask_idx[1]]]
      plan$nodes <- plan$nodes[-mask_idx[1]]
      plan$nodes <- append(plan$nodes, list(mask_node), after = target_pos - 1)
    }
  }

  plan
}
```

#### 5.1.4 Fuse Map and Reduce

When mapping to a single space followed by reduction, fuse operations.

```r
fuse_map_reduce <- function(plan) {
  i <- 1
  while (i < length(plan$nodes)) {
    if (plan$nodes[[i]]$op == "map" && plan$nodes[[i+1]]$op == "reduce") {
      # Can fuse if both operate on same axis
      # Mark for kernel fusion
      plan$nodes[[i]]$fuse_with_reduce <- TRUE
      plan$nodes[[i+1]]$fused_with_map <- TRUE
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  plan
}
```

### 5.2 Plan Serialization

Plans can be saved and loaded for reproducibility.

```r
save_plan <- function(plan, file, format = c("json", "yaml")) {
  format <- match.arg(format)

  # Serialize to list
  plan_list <- list(
    version = "1.0.0",
    source = serialize_source(plan$source),
    nodes = lapply(plan$nodes, serialize_node),
    meta = plan$meta
  )

  # Write
  if (format == "json") {
    jsonlite::write_json(plan_list, file, pretty = TRUE, auto_unbox = TRUE)
  } else {
    yaml::write_yaml(plan_list, file)
  }

  invisible(file)
}

load_plan <- function(file, format = c("auto", "json", "yaml")) {
  format <- match.arg(format)

  if (format == "auto") {
    format <- if (grepl("\\.json$", file)) "json" else "yaml"
  }

  # Read
  plan_list <- if (format == "json") {
    jsonlite::read_json(file, simplifyVector = FALSE)
  } else {
    yaml::read_yaml(file)
  }

  # Deserialize
  source <- deserialize_source(plan_list$source)
  nodes <- lapply(plan_list$nodes, deserialize_node)

  gds_plan(source, nodes, plan_list$meta)
}
```

### 5.3 Plan Explanation

```r
explain <- function(plan, verbose = FALSE) {
  cat("GDS Plan\n")
  cat("========\n\n")

  cat("Source:\n")
  cat(sprintf("  Adapter: %s\n", plan$source$adapter))
  cat(sprintf("  Files: %s\n", format_source(plan$source$source)))

  if (!is.null(plan$source$probe)) {
    cat(sprintf("  Dimensions: %s\n", paste(plan$source$probe$dims, collapse = " × ")))
    cat(sprintf("  Assays: %s\n", paste(plan$source$probe$assays, collapse = ", ")))
  }

  cat("\nOperations:\n")
  for (i in seq_along(plan$nodes)) {
    node <- plan$nodes[[i]]
    cat(sprintf("  %d. %s\n", i, format_node(node, verbose)))
  }

  cat("\nEstimated output:\n")
  cat(sprintf("  Dimensions: %s\n", estimate_output_dims(plan)))
  cat(sprintf("  Assays: %s\n", paste(estimate_output_assays(plan), collapse = ", ")))

  invisible(plan)
}

format_node <- function(node, verbose = FALSE) {
  switch(node$op,
    subset_axis = {
      parts <- c()
      if (!is.null(node$sample)) parts <- c(parts, sprintf("sample=%s", format_index(node$sample)))
      if (!is.null(node$subject)) parts <- c(parts, sprintf("subject=%s", format_index(node$subject)))
      if (!is.null(node$contrast)) parts <- c(parts, sprintf("contrast=%s", format_index(node$contrast)))
      sprintf("subset(%s)", paste(parts, collapse = ", "))
    },
    derive = sprintf("derive(%s)", paste(node$what, collapse = ", ")),
    align_to_group = sprintf("align(family=%s)", node$family$name),
    mask_policy = sprintf("mask(rule=%s)", node$policy$rule),
    map = sprintf("map_to(space=%s)", format_space(node$target_space)),
    reduce = sprintf("reduce(method=%s)", node$method),
    write = sprintf("write_out(path=%s, format=%s)", node$path, node$format),
    sprintf("<%s>", node$op)
  )
}
```

### 5.4 Plan Digest (Stable Hash)

```r
digest_plan <- function(plan) {
  # Canonical representation
  canonical <- list(
    source_hash = plan$source$hash,
    nodes = lapply(plan$nodes, canonicalize_node)
  )

  digest::digest(canonical, algo = "xxhash64")
}

canonicalize_node <- function(node) {
  # Sort params alphabetically and normalize values
  params <- node[names(node) != "op"]
  params <- params[order(names(params))]

  list(op = node$op, params = params)
}

# Realized objects carry the plan digest for traceability
# compute(plan) SHOULD write `metadata$provenance$digest <- digest_plan(plan)`
```

---

## 6. Statistical Foundations

### 6.1 Variance Propagation Through Linear Maps

For a linear map `y = M x` where `M` is `[n_target × n_source]`:

#### 6.1.1 Independence Assumption

Assume `Cov(x) = diag(v)` (diagonal covariance).

```r
# Beta (location) mapping
beta_out = M %*% beta_in

# Variance propagation
var_out[i] = Σ_j M[i,j]^2 * var_in[j]
           = ((M ⊙ M) %*% var_in)[i]
```

**Implementation:**

```r
propagate_variance_independent <- function(M, beta, var) {
  # M: [n_target × n_source] sparse or dense matrix
  # beta: [n_source × n_subject × n_contrast]
  # var: [n_source × n_subject × n_contrast]

  dims <- dim(beta)
  n_target <- nrow(M)

  # Allocate output
  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  # Process each subject × contrast
  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      beta_out[, j, k] <- M %*% beta[, j, k]
      var_out[, j, k] <- (M^2) %*% var[, j, k]
    }
  }

  list(beta = beta_out, var = var_out)
}
```

#### 6.1.2 Full Covariance

When covariance structure `Σ` is known:

```
Var(y) = M Σ M^T
```

For row `i` of `M`:
```
Var(y[i]) = M[i,·] Σ M[i,·]^T
```

**Implementation:**

```r
propagate_variance_covariance <- function(M, beta, var, cov_provider) {
  # cov_provider: function(indices) -> Sigma_block

  dims <- dim(beta)
  n_target <- nrow(M)

  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      beta_out[, j, k] <- M %*% beta[, j, k]

      # Variance propagation per target row
      for (i in seq_len(n_target)) {
        idx <- which(M[i, ] != 0)
        w <- M[i, idx]

        # Get covariance block for source indices
        Sigma <- cov_provider(idx)

        # Compute variance for this target sample
        var_out[i, j, k] <- drop(t(w) %*% Sigma %*% w)
      }
    }
  }

  list(beta = beta_out, var = var_out)
}
```

### 6.2 Degrees of Freedom Aggregation

#### 6.2.1 Satterthwaite Approximation

When combining variances with different degrees of freedom:

```
df_combined ≈ (Σ a_i^2 s_i^2)^2 / Σ(a_i^4 s_i^4 / df_i)
```

where `a_i` are weights (e.g., `M[row, i]`) and `s_i^2` are variances.

**Implementation:**

```r
satterthwaite_df <- function(weights, variances, dfs) {
  # weights: [n] coefficients
  # variances: [n] unbiased variances
  # dfs: [n] degrees of freedom

  numerator <- sum(weights^2 * variances)^2
  denominator <- sum(weights^4 * variances^2 / dfs)

  df_combined <- numerator / denominator
  df_combined
}
```

**Application in Mapping:**

```r
# For each target row i
for (i in seq_len(n_target)) {
  idx <- which(M[i, ] != 0)
  w <- M[i, idx]
  s2 <- var[idx, j, k]
  df_source <- df[idx, j, k]

  df_out[i, j, k] <- satterthwaite_df(w, s2, df_source)
}
```

### 6.3 Test Statistic Derivations

#### 6.3.1 t-statistic

From `beta` and `var`:
```
t = beta / sqrt(var)
```

Requires `df` separately.

#### 6.3.2 z-score from t

Two-sided, sign-preserving conversion:
```
z = qnorm(1 - (2*pt(-abs(t), df))/2) * sign(t)
```

**Implementation (two-sided, sign-preserving):**

```r
t_to_z <- function(t, df) {
  # two-sided p-value
  p_two <- 2 * stats::pt(-abs(t), df)
  # two-sided z with sign of t
  z <- stats::qnorm(1 - p_two/2) * sign(t)
  z
}
```

#### 6.3.3 t from z

Requires known `df` and assumes z encodes a two-sided p-value with sign:
```
t = sign(z) * qt(pnorm(|z|), df)
```

#### 6.3.4 F-statistic

For multi-parameter tests. Not directly mappable; requires effect scale and covariance.

#### 6.3.5 p-value to z-score

```
z = qnorm(p / 2, lower.tail = FALSE) * sign(effect)
```

### 6.4 Meta-Analysis Combiners

#### 6.4.1 Fixed-Effects (Inverse-Variance Weighted)

```r
combine_fixed <- function(beta, var) {
  # beta: [n_samples × n_subjects × n_contrasts]
  # var: same dims

  dims <- dim(beta)
  n_samples <- dims[1]
  n_contrasts <- dims[3]

  # Weights: 1/variance
  w <- 1 / var
  W <- apply(w, c(1, 3), sum, na.rm = TRUE)

  # Weighted mean
  beta_pooled <- apply(beta * w, c(1, 3), sum, na.rm = TRUE) / W

  # Pooled variance
  var_pooled <- 1 / W

  # Reshape to [n_samples × 1 × n_contrasts]
  beta_out <- array(beta_pooled, dim = c(n_samples, 1, n_contrasts))
  var_out <- array(var_pooled, dim = c(n_samples, 1, n_contrasts))

  list(beta = beta_out, var = var_out)
}
```

#### 6.4.2 Stouffer's Method (z-scores)

For equal weights:
```
z_combined = Σ z_i / sqrt(n)
```

For weights `w`:
```
z_combined = (Σ w_i * z_i) / sqrt(Σ w_i^2)
```

**Implementation:**

```r
combine_stouffer <- function(z, weights = NULL) {
  # z: [n_samples × n_subjects × n_contrasts]

  dims <- dim(z)
  n_samples <- dims[1]
  n_subjects <- dims[2]
  n_contrasts <- dims[3]

  if (is.null(weights)) {
    # Equal weights
    z_combined <- apply(z, c(1, 3), sum, na.rm = TRUE) / sqrt(n_subjects)
  } else {
    # Weighted
    w <- weights
    numerator <- apply(z * w, c(1, 3), sum, na.rm = TRUE)
    denominator <- sqrt(apply(w^2, c(1, 3), sum, na.rm = TRUE))
    z_combined <- numerator / denominator
  }

  # Reshape
  array(z_combined, dim = c(n_samples, 1, n_contrasts))
}
```

#### 6.4.3 Fisher's Method

```r
combine_fisher <- function(z) {
  # Convert z to two-tailed p-values
  p <- 2 * pnorm(-abs(z))

  # Chi-square statistic
  chi2 <- -2 * apply(log(p), c(1, 3), sum, na.rm = TRUE)

  # Degrees of freedom
  n_subjects <- dim(z)[2]
  df <- 2 * n_subjects

  # Combined p-value
  p_combined <- pchisq(chi2, df, lower.tail = FALSE)

  # Convert back to z (preserving sign of mean effect)
  z_mean_sign <- sign(apply(z, c(1, 3), mean, na.rm = TRUE))
  z_combined <- qnorm(p_combined / 2, lower.tail = FALSE) * z_mean_sign

  array(z_combined, dim = c(dim(z)[1], 1, dim(z)[3]))
}
```

### 6.5 Statistical Invariants

**Enforced throughout pipeline:**

1. **Never average standard errors:**
   - Always work with variance (`var = se^2`)
   - Propagate variance through operations
   - Derive SE at the end: `se = sqrt(var)`

2. **Never linearly map test statistics:**
   - Do not map `t`, `z`, `F` directly
   - Map `beta` and `var`, then recompute test statistics

3. **Preserve effect scale:**
   - If only `t/z` available, use evidence combiners (Stouffer, Fisher)
   - Do not convert `t` to `beta` without explicit scale information

4. **Variance must be positive:**
   - Validate `var > 0` after all operations
   - NA values permitted but tracked

---

## 7. Storage & I/O

**Package Dependencies:**
- NIfTI I/O: Uses `neuroim2` package (local development version at `~/code/neuroim2`)
- CSV/TSV: Uses `data.table::fread()`
- Parquet: Uses `arrow::read_parquet()`
- HDF5: Uses `hdf5r` package

### 7.1 Storage Adapter Interface (Normative)

All adapters implement a single, minimal surface. Semantics (mapping and statistics) stay in gds; adapters provide bytes efficiently.

```r
# Minimal adapter API
register_adapter(name, detect, open, probe, read, close, ...)

# Required functions
detect(source) -> score in [0,1] or FALSE
open(source, mode = c("r","w"), ...) -> handle
probe(handle) -> list(
  assays = c("beta","var","t","df", ...),
  dims   = c(sample = I, subject = J, contrast = K),
  subjects = character(J),
  contrasts = character(K),
  space = <Space>,
  maps = list(),        # optional LinearMap objects (may be empty)
  metadata = list(schema_version = "0.1.0", units = list(...))
)
read(handle, assays, block = NULL) -> named list of arrays [I × J × K]
close(handle) -> NULL
```

### 7.2 Adapter Registration

```r
# Global adapter registry
.adapter_registry <- new.env(parent = emptyenv())

register_adapter <- function(name, detect, open, probe, read, close, ...) {
  .adapter_registry[[name]] <- list(
    name = name,
    detect = detect,      # function(source) -> score in [0,1] or FALSE
    open = open,          # function(source, options) -> handle
    probe = probe,        # function(handle) -> probe_result
    read = read,          # function(handle, assays, block) -> named list
    close = close,        # function(handle) -> NULL
    ...                   # Additional adapter-specific functions
  )

  invisible(NULL)
}

get_adapter <- function(name) {
  adapter <- .adapter_registry[[name]]
  if (is.null(adapter)) {
    stop("Adapter not found: ", name)
  }
  adapter
}

detect_adapter <- function(source, prefer = NULL) {
  scores <- sapply(ls(.adapter_registry), function(name) {
    adapter <- .adapter_registry[[name]]
    score <- adapter$detect(source)
    if (isFALSE(score)) 0 else score
  })

  if (!is.null(prefer) && prefer %in% names(scores) && scores[prefer] > 0) {
    return(prefer)
  }

  best <- names(scores)[which.max(scores)]
  if (scores[best] == 0) {
    stop("No adapter found for source: ", source)
  }

  best
}
```

### 7.3 Tabular Adapter (CSV/Parquet)

```r
# Adapter for CSV/Parquet ROI/parcel data
adapter_tabular <- function() {
  list(
    name = "tabular",

    detect = function(source) {
      if (is.character(source) && length(source) == 1) {
        ext <- tools::file_ext(source)
        if (ext %in% c("csv", "tsv", "parquet")) return(0.8)
      }
      FALSE
    },

    open = function(source, mode = "r", ...) {
      if (mode == "r") {
        if (tools::file_ext(source) == "parquet") {
          df <- arrow::read_parquet(source)
        } else {
          df <- data.table::fread(source, data.table = FALSE)
        }
      } else {
        df <- NULL  # Write mode
      }

      list(path = source, mode = mode, data = df)
    },

    probe = function(handle, effect_cols, subject_col, sample_col, contrast_col,
                    space = NULL, ...) {
      df <- handle$data

      # Detect columns if not specified
      # ... column detection logic ...

      # Pivot to wide format
      # beta: [sample × subject × contrast]
      # var: [sample × subject × contrast]

      samples <- unique(df[[sample_col]])
      subjects <- unique(df[[subject_col]])
      contrasts <- unique(df[[contrast_col]])

      list(
        assays = names(effect_cols),
        dims = c(sample = length(samples),
                subject = length(subjects),
                contrast = length(contrasts)),
        subjects = subjects,
        contrasts = contrasts,
        space = space %||% space_parcels(labels = samples),
        maps = list(),
        metadata = list(
          schema_version = "0.1.0",
          source_file = handle$path
        )
      )
    },

    read = function(handle, assays, block = NULL) {
      # Read specified assays from tabular data
      # Reshape to [sample × subject × contrast]
      # ... implementation ...
    },

    close = function(handle) {
      invisible(NULL)
    }
  )
}

# Register during package load
.onLoad <- function(libname, pkgname) {
  do.call(register_adapter, adapter_tabular())
  # ... register other adapters ...
}
```

### 7.4 NIfTI Adapter

**neuroim2 API Reference:**
- `read_header(file)` - Read NIfTI header, returns `NIFTIMetaInfo` with `@dims`, `@spacing`, `@origin`
- `read_vol(file, index=1)` - Read 3D volume, returns `NeuroVol` object
- `read_vec(file)` - Read 4D volume, returns `NeuroVec` object
- `sub_vector(vec, index)` - Extract single volume from 4D `NeuroVec`
- `space(vol)` - Get `NeuroSpace` object from volume
- `trans(space)` - Extract 4×4 affine transformation matrix
- `write_vol(vol, file)` - Write `NeuroVol` to NIfTI file

**Mask Convention:**
- Masks are separate NIfTI files with binary values: 1 = included, 0 = excluded
- NOT derived from NA values in data files
- Pass mask via `gds(source, mask = "path/to/mask.nii")` or `mask = NeuroVol` object

```r
adapter_nifti <- function() {
  list(
    name = "nifti",

    detect = function(source) {
      # source can be:
      # - Single .nii/.nii.gz file
      # - Vector of files
      # - Directory with structured layout

      if (is.character(source)) {
        files <- if (dir.exists(source[1])) {
          list.files(source[1], pattern = "\\.nii(\\.gz)?$", full.names = TRUE)
        } else {
          source
        }

        if (length(files) > 0 && all(file.exists(files))) {
          return(0.9)
        }
      }

      FALSE
    },

    open = function(source, mode = "r", layout = "bids", ...) {
      # Detect layout and organize files
      # layout: "bids", "fsl", "flat", "custom"

      if (dir.exists(source)) {
        files <- detect_nifti_layout(source, layout)
      } else {
        files <- source
      }

      list(
        files = files,
        mode = mode,
        layout = layout,
        handles = list()  # Lazy open
      )
    },

    probe = function(handle, stat_names = c("beta", "var", "t", "df"),
                    subjects = NULL, contrasts = NULL,
                    mask = NULL, ...) {
      # Read header using neuroim2 to get dimensions
      first_file <- handle$files[1]
      meta <- neuroim2::read_header(first_file)

      # Get dimensions - handle both 3D and 4D files
      dim_all <- meta@dims
      if (length(dim_all) == 3L) {
        dim_voxel <- dim_all
        n_contrasts <- 1L
      } else if (length(dim_all) == 4L) {
        dim_voxel <- dim_all[1:3]
        n_contrasts <- dim_all[4]
      } else {
        stop("Unsupported NIfTI dimensionality")
      }

      # Load first volume to get spatial information
      first_vol <- if (n_contrasts == 1L) {
        neuroim2::read_vol(first_file, index = 1)
      } else {
        # For 4D files, use sub_vector to extract first volume
        vec <- neuroim2::read_vec(first_file)
        neuroim2::sub_vector(vec, 1)
      }

      # Extract affine transformation from neuroim2 NeuroSpace
      nspace <- neuroim2::space(first_vol)
      affine <- neuroim2::trans(nspace)

      # Detect subjects and contrasts from file structure
      if (is.null(subjects)) {
        subjects <- extract_subjects_from_paths(handle$files)
      }
      if (is.null(contrasts)) {
        contrasts <- extract_contrasts_from_paths(handle$files)
      }

      # Load mask from separate NIfTI file or create full mask
      # Mask should be binary: 1 = included, 0 = excluded
      if (!is.null(mask)) {
        # Mask provided as file path or NeuroVol object
        mask_vol <- if (is.character(mask)) {
          neuroim2::read_vol(mask)
        } else if (inherits(mask, "NeuroVol")) {
          mask
        } else {
          stop("mask must be a file path or NeuroVol object")
        }

        # Convert to logical bitmap (1 = included, 0 = excluded)
        mask_bitmap <- as.array(mask_vol) > 0
      } else {
        # No mask provided: include all voxels
        mask_bitmap <- array(TRUE, dim = dim_voxel)
      }

      mask_idx <- which(as.vector(mask_bitmap))

      list(
        assays = stat_names,
        dims = c(sample = length(mask_idx),
                subject = length(subjects),
                contrast = length(contrasts)),
        subjects = subjects,
        contrasts = contrasts,
        space = space_voxel(dim_voxel, affine,
                           mask_bitmap = mask_bitmap,
                           mask_idx = mask_idx,
                           storage = "packed"),
        maps = list(),
        metadata = list(
          schema_version = "0.1.0",
          source_layout = handle$layout
        )
      )
    },

    read = function(handle, assays, block = NULL) {
      # Read NIfTI files for specified assays
      # Apply mask and reshape to [sample × subject × contrast]
      # block: list(sample = c(start, end))

      # ... implementation ...
    },

    close = function(handle) {
      # Close any open file handles
      invisible(NULL)
    }
  )
}
```

### 7.5 HDF5 Adapter (GDS Layout)

```r
adapter_h5_gds <- function() {
  list(
    name = "h5_gds",

    detect = function(source) {
      if (is.character(source) && length(source) == 1 && file.exists(source)) {
        if (tools::file_ext(source) %in% c("h5", "hdf5")) {
          h5 <- NULL
          on.exit({ if (!is.null(h5)) try(h5$close(), silent = TRUE) }, add = TRUE)
          h5 <- hdf5r::H5File$new(source, mode = "r")
          if (h5$exists("/gds")) return(1.0)
        }
      }
      FALSE
    },

    open = function(source, mode = "r", ...) {
      h5 <- hdf5r::H5File$new(source, mode = mode)
      list(h5 = h5, path = source, mode = mode)
    },

    probe = function(handle) {
      h5 <- handle$h5

      # Read axes
      subjects <- h5[["/gds/axes/subjects"]][]
      contrasts <- h5[["/gds/axes/contrasts"]][]

      # Read space
      space_type <- h5[["/gds/space/type"]][]
      space <- read_h5_space(h5, space_type)

      # List assays
      assays <- names(h5[["/gds/assays"]])

      # Read dims from first assay
      beta <- h5[["/gds/assays/beta"]]
      dims <- beta$dims

      list(
        assays = assays,
        dims = c(sample = dims[1], subject = dims[2], contrast = dims[3]),
        subjects = subjects,
        contrasts = contrasts,
        space = space,
        maps = list(),
        metadata = read_h5_metadata(h5)
      )
    },

    read = function(handle, assays, block = NULL) {
      h5 <- handle$h5

      # Read assays with optional block indexing
      result <- list()
      for (assay in assays) {
        dset <- h5[[paste0("/gds/assays/", assay)]]

        if (is.null(block)) {
          result[[assay]] <- dset[]
        } else {
          # Block indexing: block$sample, block$subject, block$contrast
          i <- block$sample %||% seq_len(dset$dims[1])
          j <- block$subject %||% seq_len(dset$dims[2])
          k <- block$contrast %||% seq_len(dset$dims[3])

          result[[assay]] <- dset[i, j, k]
        }
      }

      result
    },

    close = function(handle) {
      handle$h5$close()
      invisible(NULL)
    },

    # Write methods
    write_assay = function(handle, assay, x, block = NULL) {
      h5 <- handle$h5
      path <- paste0("/gds/assays/", assay)

      if (!h5$exists(path)) {
        # Create dataset
        h5$create_dataset(path, dtype = "float32",
                         dims = dim(x),
                         chunk_dims = c(min(10000, dim(x)[1]), 1, 1),
                         compression = "gzip")
      }

      dset <- h5[[path]]

      if (is.null(block)) {
        dset[] <- x
      } else {
        i <- block$sample %||% seq_len(dim(x)[1])
        j <- block$subject %||% seq_len(dim(x)[2])
        k <- block$contrast %||% seq_len(dim(x)[3])

        dset[i, j, k] <- x
      }
    }
  )
}
```

**HDF5 Layout Specification:**

```
/gds/
  /axes/
    /subjects            (string[S])
    /contrasts           (string[K])
  /space/
    /type                (string: "voxel" | "parcels" | "surface" | "basis")
    # If voxel:
    /dim                 (int[3])
    /affine              (double[4,4])
    /mask                (uint8[X,Y,Z]) or /mask_idx (int[V])
    /storage             (string: "dense" | "packed")
    # If parcels:
    /labels              (string[V])
    # If surface:
    /vertices            (double[N,3])
    /faces               (int[F,3])
    /hemi                (string)
    # If basis:
    /k                   (int)
    /basis_name          (string)
  /assays/
    /beta                (float32[V,S,K])
    /var                 (float32[V,S,K])
    /se                  (float32[V,S,K])
    /t                   (float32[V,S,K])
    /z                   (float32[V,S,K])
    /df                  (int[S,K] or [V,S,K])
    /meta/
      /beta.json         (JSON: role, units, etc.)
      /var.json          (JSON: variance_of, etc.)
  /contrast_info         (JSON)
  /provenance            (JSONL: append-only log)
  /metadata              (JSON)
  /version               (string: "gds-h5/1.0")
```

### 7.6 fmristore Integration

fmristore serves as an official storage backend via adapter.

**Path A (No Schema Change):**
- Use existing fmristore layouts
- Adapter composes subjects/contrasts across files
- Maps `/mask`, `/labels`, latent bases to GDS spaces

**Path B (GDS Shim Group - Recommended):**
- Add `/gds` group to fmristore files
- Perfect round-tripping
- Single file per cohort

**fmristore Adapter (Sketch):**

```r
adapter_fmristore <- function() {
  list(
    name = "fmristore",

    detect = function(source) {
      # Detect fmristore HDF5 files
      # Use fmristore::detect_h5_type()

      if (is.character(source) && file.exists(source)) {
        type <- fmristore::detect_h5_type(source)
        if (!is.null(type)) return(0.95)
      }
      FALSE
    },

    open = function(source, mode = "r", ...) {
      # Open using fmristore classes
      type <- fmristore::detect_h5_type(source)

      obj <- switch(type,
        "labeled_volume" = fmristore::H5NeuroVol(source),
        "parcellated" = fmristore::H5ParcellatedScan(source),
        "latent" = fmristore::LatentNeuroVec(source),
        stop("Unknown fmristore type: ", type)
      )

      list(obj = obj, type = type, path = source)
    },

    probe = function(handle, ...) {
      # Extract GDS-compatible info from fmristore object

      obj <- handle$obj
      type <- handle$type

      # Build space
      space <- switch(type,
        "labeled_volume" = {
          # Extract mask, affine from /header/*
          mask <- obj$mask()
          affine <- obj$affine()
          dim_voxel <- dim(mask)
          space_voxel(dim_voxel, affine,
                     mask_bitmap = mask,
                     storage = "packed")
        },
        "parcellated" = {
          clusters <- obj$cluster_map()
          labels <- unique(as.vector(clusters))
          space_parcels(labels = as.character(labels))
        },
        "latent" = {
          basis <- obj$basis()
          k <- ncol(basis)
          space_basis(k, basis_name = "latent", projector = t(basis))
        }
      )

      # ... extract subjects, contrasts, assays from file structure ...

      list(
        assays = c("beta", "var"),  # Detected from files
        dims = c(sample = ..., subject = ..., contrast = ...),
        subjects = ...,
        contrasts = ...,
        space = space,
        metadata = ...
      )
    },

    read = function(handle, assays, block = NULL) {
      # Use fmristore methods to read data
      # Reshape to GDS format
      # ...
    },

    close = function(handle) {
      # Close fmristore object
      invisible(NULL)
    }
  )
}
```

---

## 8. Spatial Operations

### 8.1 Mask Operations

#### 8.1.1 Union Mask

```r
compute_mask_union <- function(mask_list) {
  # mask_list: list of logical arrays [X,Y,Z]

  Reduce(`|`, mask_list)
}
```

#### 8.1.2 Intersection Mask

```r
compute_mask_intersection <- function(mask_list) {
  Reduce(`&`, mask_list)
}
```

#### 8.1.3 Threshold Mask

```r
compute_mask_threshold <- function(mask_list, threshold = 0.95) {
  # Count proportion of subjects with each voxel
  n_subjects <- length(mask_list)

  mask_sum <- Reduce(`+`, lapply(mask_list, as.integer))
  mask_prop <- mask_sum / n_subjects

  mask_prop >= threshold
}
```

#### 8.1.4 Pack/Unpack Assays

```r
pack_assays <- function(assays, mask_idx) {
  # assays: list of arrays [n_voxels × subject × contrast]
  # mask_idx: integer vector of indices to keep

  lapply(assays, function(A) {
    A[mask_idx, , , drop = FALSE]
  })
}

unpack_assays <- function(assays, mask_idx, full_n) {
  # assays: list of arrays [n_masked × subject × contrast]
  # mask_idx: integer vector of packed indices
  # full_n: total number of voxels

  lapply(assays, function(A) {
    dims <- dim(A)
    out <- array(NA_real_, dim = c(full_n, dims[2], dims[3]))
    out[mask_idx, , ] <- A
    out
  })
}
```

### 8.2 Alignment Transforms

#### 8.2.1 Apply Linear Alignment (Per-Subject)

```r
apply_alignment_linear <- function(M_by_subject, beta, var, uncertainty_rule) {
  # M_by_subject: named list of matrices [n_target × n_source]
  # beta, var: [n_source × n_subjects × n_contrasts]

  dims <- dim(beta)
  subjects <- dimnames(beta)[[2]]
  n_target <- nrow(M_by_subject[[1]])

  beta_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))
  var_out <- array(NA_real_, dim = c(n_target, dims[2], dims[3]))

  for (j in seq_along(subjects)) {
    subj <- subjects[j]
    M <- M_by_subject[[subj]]

    for (k in seq_len(dims[3])) {
      # Apply mapping
      beta_out[, j, k] <- M %*% beta[, j, k]

      # Propagate variance
      if (uncertainty_rule$mode == "independent") {
        var_out[, j, k] <- (M^2) %*% var[, j, k]
      } else if (uncertainty_rule$mode == "cov_provider") {
        var_out[, j, k] <- propagate_with_cov(M, var[, j, k],
                                              uncertainty_rule$cov_provider)
      }
    }
  }

  dimnames(beta_out) <- list(NULL, subjects, dimnames(beta)[[3]])
  dimnames(var_out) <- list(NULL, subjects, dimnames(var)[[3]])

  list(beta = beta_out, var = var_out)
}
```

#### 8.2.2 Orthogonal Family (Procrustes)

Used for latent basis alignment.

```r
# Example: align subject bases to consensus
align_orthogonal <- function(gds, R_by_subject, target_space) {
  # R_by_subject: named list of orthogonal matrices [k × k]

  family <- OrthogonalFamily(
    name = "procrustes_consensus",
    from_space = gds$space,
    to_space = target_space,
    matrices_by_subject = R_by_subject
  )

  # Register and align
  gds <- register_map(gds, family)
  plan <- align(gds, "procrustes_consensus")

  plan
}
```

#### 8.2.3 Optimal Transport (Wasserstein)

Mass-preserving alignment.

```r
align_ot <- function(gds, transport_plans, target_space) {
  # transport_plans: named list of row-stochastic matrices [k_target × k_source]

  family <- OTFamily(
    name = "wasserstein_alignment",
    from_space = gds$space,
    to_space = target_space,
    plans_by_subject = transport_plans
  )

  gds <- register_map(gds, family)
  plan <- align(gds, "wasserstein_alignment")

  plan
}
```

### 8.3 Space Transformations

#### 8.3.1 Voxel to Parcel

```r
voxel_to_parcel_map <- function(parcel_space, voxel_space, method = c("mean", "weighted_mean")) {
  method <- match.arg(method)

  # Build mapping matrix from membership
  membership <- parcel_space$membership
  n_parcels <- length(parcel_space$labels)
  n_voxels <- length(voxel_space$mask_idx)

  M <- Matrix::sparseMatrix(
    i = rep(seq_len(n_parcels), times = sapply(membership, length)),
    j = unlist(membership),
    x = if (method == "mean") {
      # Equal weights
      rep(1 / sapply(membership, length), times = sapply(membership, length))
    } else {
      # Could use distance-weighted, etc.
      1
    },
    dims = c(n_parcels, n_voxels)
  )

  map_linear(
    from_space = voxel_space,
    to_space = parcel_space,
    operator = M,
    uncertainty = UncertaintyRule("independent")
  )
}
```

#### 8.3.2 Basis to Voxel

```r
basis_to_voxel_map <- function(basis_space, voxel_space) {
  # Basis space has projector B [k × n_voxels]
  # Forward: y = B %*% x_voxel
  # Backward: x_voxel = t(B) %*% y

  B <- basis_space$projector

  map_linear(
    from_space = basis_space,
    to_space = voxel_space,
    operator = t(B),  # Transpose for backward projection
    uncertainty = UncertaintyRule("independent")
  )
}
```

### 8.4 Late-Binding Map Registration

```r
register_map <- function(x, family, overwrite = FALSE) {
  # x: GDS or plan
  # family: MapFamily object

  # Validate subjects match
  subjects_x <- if (inherits(x, "gds")) x$subjects else x$source$probe$subjects
  subjects_fam <- names(family$by_subject)

  if (!setequal(subjects_x, subjects_fam)) {
    stop("Subject mismatch between GDS and MapFamily")
  }

  # Validate dimensions
  # ... check that operator dims match source/target space sample counts ...

  # Add to space_graph
  if (inherits(x, "gds")) {
    if (is.null(x$metadata$space_graph)) {
      x$metadata$space_graph <- list(maps = list())
    }

    if (family$name %in% names(x$metadata$space_graph$maps) && !overwrite) {
      stop("MapFamily '", family$name, "' already registered. Use overwrite=TRUE.")
    }

    x$metadata$space_graph$maps[[family$name]] <- family

    # Add provenance
    x$metadata <- add_provenance_node(
      x$metadata,
      "register_map",
      list(family_name = family$name, type = family$type),
      inputs = list()
    )

  } else if (inherits(x, "gds_plan")) {
    # Store in plan metadata
    if (is.null(x$meta$space_graph)) {
      x$meta$space_graph <- list(maps = list())
    }

    x$meta$space_graph$maps[[family$name]] <- family
  }

  x
}
```

---

## 9. Implementation Modules

### 9.1 Package Structure

```
gdsfmri/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── gds-class.R
│   ├── gds-accessors.R
│   ├── gds-validators.R
│   ├── space-voxel.R
│   ├── space-parcels.R
│   ├── space-surface.R
│   ├── space-basis.R
│   ├── map-linear.R
│   ├── map-family.R
│   ├── uncertainty-rule.R
│   ├── plan-core.R
│   ├── plan-ops.R
│   ├── plan-optimizer.R
│   ├── plan-serialization.R
│   ├── verb-gds.R
│   ├── verb-subset.R
│   ├── verb-derive.R
│   ├── verb-align.R
│   ├── verb-mask.R
│   ├── verb-map-to.R
│   ├── verb-reduce.R
│   ├── verb-write-out.R
│   ├── verb-compute.R
│   ├── assay-registry.R
│   ├── derive-stats.R
│   ├── propagate-variance.R
│   ├── combine-fixed.R
│   ├── combine-stouffer.R
│   ├── combine-fisher.R
│   ├── adapter-interface.R
│   ├── adapter-tabular.R
│   ├── adapter-nifti.R
│   ├── adapter-h5-gds.R
│   ├── adapter-fmristore.R
│   ├── io-nifti.R
│   ├── io-h5.R
│   ├── mask-operations.R
│   ├── provenance.R
│   ├── metadata.R
│   ├── utils-misc.R
│   └── zzz.R (package hooks)
├── src/ (optional C++ for performance)
├── inst/
│   ├── schemas/
│   │   ├── h5-gds-layout.md
│   │   └── space-specification.md
│   └── extdata/ (test data)
├── tests/
│   └── testthat/
│       ├── test-gds-class.R
│       ├── test-spaces.R
│       ├── test-maps.R
│       ├── test-plan.R
│       ├── test-verbs.R
│       ├── test-derive.R
│       ├── test-variance-propagation.R
│       ├── test-adapters.R
│       ├── test-mask.R
│       ├── test-alignment.R
│       └── test-meta-analysis.R
├── vignettes/
│   ├── introduction.Rmd
│   ├── spaces.Rmd
│   ├── lazy-execution.Rmd
│   ├── alignment.Rmd
│   ├── fmristore-integration.Rmd
│   └── migration-from-group-data.Rmd
└── man/ (generated by roxygen2)
```

### 9.2 Key Implementation Files

#### 9.2.1 R/gds-class.R

```r
#' Create a new GDS object
#'
#' @param assays Named list of 3D arrays [sample × subject × contrast]
#' @param space Space object (space_voxel, space_parcels, etc.)
#' @param subjects Character vector of subject IDs
#' @param contrasts Character vector of contrast names
#' @param col_data Data frame of subject covariates (optional)
#' @param row_data Data frame of sample metadata (optional)
#' @param metadata List of metadata (optional)
#'
#' @return GDS object
#' @export
new_gds <- function(assays, space, subjects, contrasts,
                    col_data = NULL, row_data = NULL,
                    metadata = list()) {
  # Implementation as in §2.1
}

#' @export
print.gds <- function(x, ...) {
  cat("GDS object\n")
  cat("==========\n\n")

  dims <- dim(x$assays[[1]])
  cat(sprintf("Dimensions: %d samples × %d subjects × %d contrasts\n",
              dims[1], dims[2], dims[3]))

  cat(sprintf("Space: %s\n", x$space$type))
  cat(sprintf("Assays: %s\n", paste(names(x$assays), collapse = ", ")))

  if (nrow(x$col_data) > 0) {
    cat(sprintf("Subject covariates: %s\n",
                paste(names(x$col_data), collapse = ", ")))
  }

  invisible(x)
}

#' @export
summary.gds <- function(object, ...) {
  print(object, ...)

  cat("\nAssay summaries:\n")
  for (assay_name in names(object$assays)) {
    A <- object$assays[[assay_name]]
    cat(sprintf("  %s: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
                assay_name,
                mean(A, na.rm = TRUE),
                sd(A, na.rm = TRUE),
                min(A, na.rm = TRUE),
                max(A, na.rm = TRUE)))
  }

  invisible(object)
}
```

#### 9.2.2 R/verb-derive.R

```r
#' Derive additional assays
#'
#' @param x GDS, plan, or source
#' @param what Character vector of assays to derive
#' @param options List of derivation options
#'
#' @return Plan object
#' @export
derive <- function(x, what = c("var", "se", "t", "z"), options = list()) {
  plan <- as_plan(x)
  add_op(plan, op_derive(what, options))
}

# Internal: execute derivation
execute_derive <- function(assays, what, options) {
  for (target in what) {
    if (target %in% names(assays) && !isTRUE(options$overwrite)) {
      next  # Skip existing
    }

    assays[[target]] <- switch(target,
      var = derive_var(assays),
      se = derive_se(assays),
      t = derive_t(assays),
      z = derive_z(assays),
      stop("Unknown derivation: ", target)
    )
  }

  assays
}
```

#### 9.2.3 R/propagate-variance.R

```r
#' Propagate variance through linear mapping (independent assumption)
#'
#' @param M Mapping matrix [n_target × n_source]
#' @param beta Effect array [n_source × n_subjects × n_contrasts]
#' @param var Variance array (same dims as beta)
#'
#' @return List with beta_out and var_out
propagate_variance_independent <- function(M, beta, var) {
  # Implementation as in §6.1.1
}

#' Propagate variance with covariance provider
#'
#' @param M Mapping matrix
#' @param beta Effect array
#' @param var Variance array
#' @param cov_provider Function(indices) -> Sigma_block
#'
#' @return List with beta_out and var_out
propagate_variance_covariance <- function(M, beta, var, cov_provider) {
  # Implementation as in §6.1.2
}
```

### 9.3 Dependency Management

**Core Dependencies:**
- `Matrix`: Sparse matrices
- `digest`: Hashing for provenance
- `jsonlite`: JSON serialization
- `hdf5r`: HDF5 I/O (for HDF5 adapter)
- `RNifti` or `neuroim2`: NIfTI I/O

**Optional Dependencies:**
- `arrow`: Parquet support
- `data.table`: Fast CSV reading
- `yaml`: YAML serialization
- `future`: Parallel execution
- `fmristore`: Official storage backend

**No Bioconductor Required:**
- `DelayedArray` is NOT required
- Block streaming implemented directly

---

## 10. Test Specifications

### 10.1 Test Coverage Requirements

Minimum 90% code coverage across:
1. Core objects (GDS, Space, Map, Plan)
2. All verbs
3. Statistical derivations
4. Variance propagation
5. Adapters

### 10.2 Unit Tests

#### 10.2.1 test-gds-class.R

```r
test_that("GDS constructor validates dimensions", {
  beta <- array(rnorm(100 * 10 * 3), dim = c(100, 10, 3))
  var <- array(runif(100 * 10 * 3, 0.01, 1), dim = c(100, 10, 3))

  space <- space_parcels(labels = paste0("ROI_", 1:100))
  subjects <- paste0("sub-", sprintf("%02d", 1:10))
  contrasts <- c("A", "B", "C")

  gds <- new_gds(
    assays = list(beta = beta, var = var),
    space = space,
    subjects = subjects,
    contrasts = contrasts
  )

  expect_s3_class(gds, "gds")
  expect_equal(dim(assay(gds, "beta")), c(100, 10, 3))
})

test_that("GDS constructor enforces beta + var/se requirement", {
  beta <- array(rnorm(100 * 10 * 3), dim = c(100, 10, 3))

  expect_error(
    new_gds(
      assays = list(beta = beta),
      space = space_parcels(labels = paste0("ROI_", 1:100)),
      subjects = paste0("sub-", 1:10),
      contrasts = c("A", "B", "C")
    ),
    "must have either 'var' or 'se'"
  )
})
```

#### 10.2.2 test-variance-propagation.R

```r
test_that("Independent variance propagation is correct", {
  # Simple 2→1 mapping
  M <- matrix(c(0.6, 0.8), nrow = 1)
  beta <- matrix(c(10, 20), ncol = 1)
  var <- matrix(c(1, 4), ncol = 1)

  result <- propagate_variance_independent(M, beta, var)

  # Expected beta: 0.6*10 + 0.8*20 = 22
  expect_equal(result$beta[1,1], 22)

  # Expected var: 0.6^2*1 + 0.8^2*4 = 0.36 + 2.56 = 2.92
  expect_equal(result$var[1,1], 2.92)
})

test_that("Variance is always non-negative", {
  M <- matrix(runif(10 * 20), nrow = 10)
  beta <- array(rnorm(20 * 5 * 2), dim = c(20, 5, 2))
  var <- array(runif(20 * 5 * 2, 0.1, 2), dim = c(20, 5, 2))

  result <- propagate_variance_independent(M, beta, var)

  expect_true(all(result$var >= 0, na.rm = TRUE))
})
```

#### 10.2.3 test-derive.R

```r
test_that("Derive var from se", {
  se <- array(runif(100 * 10 * 3, 0.1, 2), dim = c(100, 10, 3))

  assays <- list(se = se)
  result <- execute_derive(assays, "var", list())

  expect_equal(result$var, se^2)
})

test_that("Derive t from beta and var", {
  beta <- array(rnorm(100 * 10 * 3), dim = c(100, 10, 3))
  var <- array(runif(100 * 10 * 3, 0.1, 2), dim = c(100, 10, 3))

  assays <- list(beta = beta, var = var)
  result <- execute_derive(assays, "t", list())

  expected_t <- beta / sqrt(var)
  expect_equal(result$t, expected_t)
})

test_that("Derive z from t and df", {
  t <- array(rnorm(100 * 10 * 3), dim = c(100, 10, 3))
  df <- array(30, dim = c(10, 3))  # Constant df across samples

  assays <- list(t = t, df = df)
  result <- execute_derive(assays, "z", list())

  # Check first element
  p_two <- 2 * pt(-abs(t[1, 1, 1]), df[1, 1])
  expected_z <- qnorm(1 - p_two/2) * sign(t[1, 1, 1])
  expect_equal(result$z[1, 1, 1], expected_z)
})
```

#### 10.2.4 test-plan.R

```r
test_that("Plan pushdown optimization works", {
  # Create mock source
  src <- gds_source("mock", list(data = mock_data))
  plan <- gds_plan(source = src)

  # Add two subset operations
  plan <- subset(plan, subject = c("sub-01", "sub-02"))
  plan <- subset(plan, contrast = c("A", "B"))

  # Optimize
  plan_opt <- optimize_plan(plan)

  # Should have merged subsets into source
  expect_null(plan_opt$source$probe$subset)
  expect_length(plan_opt$nodes, 0)  # All pushed down
})

test_that("Plan serialization round-trips", {
  src <- gds_source("tabular", "roi_stats.csv")
  plan <- gds_plan(source = src) %>%
    subset(subject = c("sub-01", "sub-02")) %>%
    derive(c("var", "t"))

  # Save and load
  tmp <- tempfile(fileext = ".json")
  save_plan(plan, tmp)
  plan2 <- load_plan(tmp)

  expect_equal(digest_plan(plan), digest_plan(plan2))
})
```

### 10.3 Integration Tests

#### 10.3.1 End-to-End Voxel Pipeline

```r
test_that("Voxel pipeline with alignment and mapping", {
  skip_if_not_installed("RNifti")

  # Create synthetic NIfTI data
  # ... setup ...

  plan <- gds(nifti_files) %>%
    subset(subject = paste0("sub-", 1:10)) %>%
    derive(c("var", "t")) %>%
    align(map = orthogonal_family) %>%
    mask(policy = MaskPolicy(rule = "threshold", threshold = 0.9)) %>%
    map_to(parcel_space, voxel_to_parcel_map) %>%
    reduce(method = "fixed")

  result <- compute(plan, sink = "memory")

  # Validate output
  expect_s3_class(result, "gds")
  expect_equal(result$space$type, "parcels")
  expect_equal(dim(assay(result, "beta"))[2], 1)  # Reduced across subjects
  expect_true("var" %in% names(result$assays))
})
```

#### 10.3.2 Stouffer Combination

```r
test_that("Stouffer combination of z-scores is correct", {
  # Create synthetic z-scores
  z <- array(rnorm(100 * 20 * 3, mean = 2), dim = c(100, 20, 3))

  assays <- list(z = z)
  result <- combine_stouffer(z, weights = NULL)

  # Check dimensions
  expect_equal(dim(result), c(100, 1, 3))

  # Manual calculation for first sample/contrast
  expected <- sum(z[1, , 1]) / sqrt(20)
  expect_equal(result[1, 1, 1], expected)
})
```

### 10.4 Edge Cases and Error Handling

```r
test_that("Handles NA values correctly", {
  beta <- array(rnorm(100 * 10 * 3), dim = c(100, 10, 3))
  beta[1:5, 1, 1] <- NA
  var <- array(runif(100 * 10 * 3, 0.1, 2), dim = c(100, 10, 3))

  # Should propagate NAs
  M <- matrix(runif(50 * 100), nrow = 50)
  result <- propagate_variance_independent(M, beta, var)

  expect_true(any(is.na(result$beta)))
})

test_that("Refuses to map t-statistics without effect scale", {
  t <- array(rnorm(100 * 10 * 3), dim = c(100, 10, 3))
  df <- array(30, dim = c(10, 3))

  gds_obj <- new_gds(
    assays = list(t = t, df = df),
    space = space_parcels(labels = paste0("ROI_", 1:100)),
    subjects = paste0("sub-", 1:10),
    contrasts = c("A", "B", "C")
  )

  M <- matrix(runif(50 * 100), nrow = 50)

  expect_error(
    map_to(gds_obj, target_space = space_parcels(labels = paste0("Cluster_", 1:50)),
           map = M),
    "Cannot linearly map.*without effect scale"
  )
})
```

---

## 11. Migration Strategy

### 11.1 Backward Compatibility with group_data

Keep existing `fmrireg::group_data` API intact while using GDS internally.

```r
# In fmrireg package:

#' @export
group_data <- function(...) UseMethod("group_data")

#' @export
group_data.default <- function(source, ...) {
  # Delegate to GDS
  plan <- gdsfmri::gds(source, ...)

  # Option to return plan or compute eagerly
  if (isTRUE(getOption("gds.lazy", FALSE))) {
    structure(plan, class = c("group_data_gds", "group_data", "gds_plan"))
  } else {
    result <- gdsfmri::compute(plan)
    structure(result, class = c("group_data_gds", "group_data", "gds"))
  }
}

# Keep existing group_data_from_* constructors
#' @export
group_data_from_tabular <- function(...) {
  plan <- gdsfmri::gds_from_tabular(...)

  if (isTRUE(getOption("gds.lazy", FALSE))) {
    structure(plan, class = c("group_data_gds", "group_data", "gds_plan"))
  } else {
    result <- gdsfmri::compute(plan)
    structure(result, class = c("group_data_gds", "group_data", "gds"))
  }
}

# Preserve existing methods
#' @export
n_subjects.group_data_gds <- function(x) {
  if (inherits(x, "gds_plan")) {
    length(x$source$probe$subjects)
  } else {
    length(x$subjects)
  }
}

#' @export
get_subjects.group_data_gds <- function(x) {
  if (inherits(x, "gds_plan")) {
    x$source$probe$subjects
  } else {
    x$subjects
  }
}

#' @export
get_covariates.group_data_gds <- function(x) {
  if (inherits(x, "gds_plan")) {
    x$source$probe$col_data
  } else {
    x$col_data
  }
}
```

### 11.2 Conversion Utilities

```r
#' Convert group_data to GDS
#'
#' @param x group_data object
#' @export
as_gds <- function(x) UseMethod("as_gds")

#' @export
as_gds.group_data_gds <- function(x) {
  # Already GDS-backed
  if (inherits(x, "gds_plan")) {
    gdsfmri::compute(x)
  } else {
    class(x) <- c("gds", "list")
    x
  }
}

#' @export
as_gds.group_data <- function(x) {
  # Convert legacy group_data to GDS
  # Extract assays, space, etc.
  # ...
}

#' Convert GDS to group_data
#'
#' @param x GDS object
#' @export
as_group_data <- function(x) {
  structure(x, class = c("group_data_gds", "group_data", "gds"))
}
```

### 11.3 Migration Guide for Users

**Step 1: Install gdsfmri**
```r
# install.packages("gdsfmri")  # When on CRAN
# remotes::install_github("user/gdsfmri")
```

**Step 2: Enable lazy mode (optional)**
```r
options(gds.lazy = TRUE)
```

**Step 3: Existing code continues to work**
```r
# Old code:
gd <- group_data_from_tabular("roi_stats.csv", ...)

# Still works! Now backed by GDS
```

**Step 4: Adopt new API gradually**
```r
# New code:
plan <- gds("roi_stats.csv") %>%
  subset(subject = my_subjects) %>%
  derive(c("var", "t")) %>%
  reduce(method = "fixed")

result <- compute(plan)
```

### 11.4 Breaking Changes (None Expected)

No breaking changes for existing `group_data` users. New GDS API is additive.

---

## 12. fmrireg Integration: Plugin Architecture

### 12.1 Design Principles

**The Bridge Pattern: fmrireg.gds as a Mounting Kit**

The purpose of `fmrireg.gds` is **not** to replace or duplicate anything in fmrireg. It exists purely as a **"bridge" layer** that lets two autonomous systems communicate cleanly without either depending on the other:

| Package | Primary Responsibility | Depends On |
|---------|------------------------|------------|
| **gdsfmri** | Data orchestration: GDS model, spaces, masks, alignment, lazy pipeline, I/O | — |
| **fmrireg** | Statistical computation: meta-analysis, evidence combiners, spatial FDR | — |
| **fmrireg.gds** | Translator/bridge: registers fmrireg's methods inside GDS pipeline | gdsfmri, fmrireg |

**Analogy:**
- `fmrireg` = high-performance engine (statistical methods)
- `gdsfmri` = car chassis (data structure, steering, wheels, fuel lines)
- `fmrireg.gds` = mounting kit that lets the engine bolt cleanly into the chassis

Without it, both are fine on their own. With it, you get a complete, driveable system.

**Clear Separation of Concerns:**

- **gdsfmri** owns:
  - Data model: `[sample × subject × contrast]` dimensional format
  - Spaces & maps: voxel/parcel/surface/basis with subject→group alignment
  - Uncertainty propagation across space maps (beta/var only; stats re-derived)
  - Lazy pipeline: Plan + block streaming + provenance + I/O adapters
  - Operator registries: pluggable reducers and post-hoc steps (generic infrastructure)

- **fmrireg** owns:
  - Meta-analysis & meta-regression engines (FE/REML/PM/DL, moderators, heterogeneity)
  - Evidence combiners (Stouffer/Lipták, Fisher, Lancaster) for t/z/F/χ² families
  - Multiple-comparison control (structure-adaptive spatial FDR, cluster-extent)
  - All statistical algorithms (single source of truth)

- **fmrireg.gds** owns:
  - **Only registration code** (no algorithms)
  - `.onLoad()` hook that calls `gdsfmri::register_reducer()` and `gdsfmri::register_posthoc()`
  - Thin wrappers that call `fmrireg::meta_random_core()`, `fmrireg::spatial_fdr_core()`, etc.
  - Optional convenience functions (e.g., `fmri_meta()` sugar)
  - **~100-300 lines of code total**

**Integration Strategy:** Plugin API, not object conversion

- fmrireg exposes stateless kernels that operate on blocks shaped `[subjects × samples_block]`
- GDS calls these kernels during `compute()` while streaming data
- fmrireg.gds makes fmrireg **discoverable and callable** inside the GDS execution engine
- No round-trip conversions between data models
- No algorithm duplication (all algorithms stay in fmrireg)

**Dependency Structure (one-way, acyclic):**

```
fmrireg.gds → Imports: fmrireg, gdsfmri
gdsfmri     → no dependency on either package (optional plugin)
fmrireg     → no dependency on gdsfmri or fmrireg.gds (autonomous)
```

**What happens without fmrireg.gds installed:**
- gdsfmri can read, align, and organize data but has no meta-analysis methods
- fmrireg can analyze data, but only after you reshape it into `group_data_*` format
- Both packages work independently

**What happens with fmrireg.gds installed:**
- When fmrireg.gds loads, it tells gdsfmri: "Here are reducers called `meta:pm`, `meta:reml`, `combine:stouffer`, `fdr:spatial`, etc. When a Plan asks for these, call my functions."
- Under the hood, those reducers just call `fmrireg::meta_random_core()` or `fmrireg::spatial_fdr_core()` - no data conversion, no copy, just streamed arrays
- Users write `reduce(method = "meta:pm")` and GDS seamlessly calls fmrireg

**Architectural Benefits:**
- **No circular dependencies** - clean acyclic graph
- **Optional integration** - install fmrireg.gds only if you need meta-analysis in GDS pipelines
- **Maintenance isolation** - statistical updates in fmrireg; data-model changes in gdsfmri
- **Open plugin system** - later you could have `nilearn.gds`, `afni.gds`, etc., all coexisting
- **Independent testing** - fmrireg tests pure statistical kernels; fmrireg.gds tests registration

**Visual Summary:**

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Code                               │
│  library(gdsfmri)                                              │
│  library(fmrireg.gds)  # Loads bridge, registers reducers     │
│                                                                 │
│  plan <- gds("data.nii.gz") %>%                                │
│    reduce(method = "meta:pm") %>%                              │
│    posthoc(method = "fdr:spatial")                             │
│  result <- compute(plan)                                       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        gdsfmri                                  │
│  • Data model [sample × subject × contrast]                    │
│  • Spaces, masks, alignment                                    │
│  • Lazy Plan + block streaming                                 │
│  • Reducer/posthoc registries (generic)                        │
│  • compute() executor                                          │
│                                                                 │
│  During compute():                                             │
│    1. Streams blocks [subjects × samples]                      │
│    2. Calls get_reducer("meta:pm")                             │
│    3. Executes registered kernel on block                      │
│    4. Writes results to sink                                   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ (looks up registered methods)
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      fmrireg.gds (bridge)                       │
│  • .onLoad() registers methods:                                │
│      register_reducer("meta:pm", wrapper_fn)                   │
│      register_posthoc("fdr:spatial", wrapper_fn)               │
│  • Wrapper functions immediately call fmrireg::*_core()        │
│  • Optional sugar: fmri_meta(plan, ...)                        │
│  • ~100-300 LOC total                                          │
│  • ZERO statistical algorithms                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ (delegates to)
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                         fmrireg                                 │
│  • meta_fixed_core(beta, var, X, df, opts)                     │
│  • meta_random_core(beta, var, X, df, opts)                    │
│    - DL, PM, REML τ² estimators                                │
│  • combine_stouffer_core(z, w, opts)                           │
│  • combine_fisher_core(p, opts)                                │
│  • spatial_fdr_core(z, graph, alpha, opts)                     │
│  • ALL STATISTICAL ALGORITHMS (single source of truth)         │
│  • No dependency on gdsfmri                                    │
└─────────────────────────────────────────────────────────────────┘

WITHOUT fmrireg.gds:
  gdsfmri: ✓ works (data orchestration)
  fmrireg: ✓ works (statistical analysis)
  Interop: ✗ manual reshaping required

WITH fmrireg.gds:
  gdsfmri: ✓ works (data orchestration)
  fmrireg: ✓ works (statistical analysis)
  Interop: ✓ seamless streaming via registered reducers
```

### 12.2 Reducer API (Per-Sample Meta-Analysis)

Each reducer operates on one contrast at a time over the subject axis.

**Registration Interface:**

```r
# In gdsfmri (core package)
.gds_reducers <- new.env(parent = emptyenv())

register_reducer <- function(name, fun, requires, provides,
                            options_schema = list(),
                            doc = NULL) {
  stopifnot(
    is.character(name), nzchar(name),
    is.function(fun),
    is.character(requires),
    is.character(provides)
  )

  .gds_reducers[[name]] <- list(
    name = name,
    fun = fun,
    requires = requires,     # c("beta","var") or c("z") or c("t","df")
    provides = provides,     # c("beta_g","var_g","tau2","Q","I2","p")
    options_schema = options_schema,
    doc = doc
  )

  invisible(name)
}

get_reducer <- function(name) {
  reducer <- .gds_reducers[[name]]
  if (is.null(reducer)) {
    stop("Reducer not found: ", name, "\n",
         "Available reducers: ", paste(ls(.gds_reducers), collapse = ", "))
  }
  reducer
}

list_reducers <- function() {
  names(.gds_reducers)
}
```

**Kernel Signature:**

```r
# Standard reducer kernel signature
# beta:  [subjects × samples]   (block; NULL if combiner-only)
# var:   [subjects × samples]   (block; NULL if combiner-only)
# X:     [subjects × p]          (moderators from col_data + formula)
# df:    [subjects × samples]    (degrees of freedom; NULL if not needed)
# opts:  named list              (method-specific options)
# return: named list of matrices [1 × samples] for each output assay

reducer_kernel <- function(beta, var, X = NULL, df = NULL, opts = list()) {
  # Returns: list(beta_g = ..., var_g = ..., tau2 = ..., Q = ..., I2 = ..., p = ...)
}
```

**Execution Flow:**

GDS's `reduce(method = ...)` resolves to a registered reducer and streams blocks:
1. Slices the current contrast
2. Passes a sample block (e.g., 100k voxels) across all subjects to the kernel
3. The kernel returns group summaries per sample; GDS writes them to the sink

No conversions. Pure functions. Easy to test and parallelize.

### 12.3 Post-Hoc API (Spatial Multiple Comparisons)

Post-hoc kernels operate on realized group-level maps (per contrast), using a spatial graph supplied by GDS from the active Space.

**Registration Interface:**

```r
# In gdsfmri (core package)
.gds_posthoc <- new.env(parent = emptyenv())

register_posthoc <- function(name, fun, requires = "z",
                            provides = c("q", "reject", "thresh"),
                            options_schema = list(),
                            doc = NULL) {
  stopifnot(
    is.character(name), nzchar(name),
    is.function(fun)
  )

  .gds_posthoc[[name]] <- list(
    name = name,
    fun = fun,
    requires = requires,
    provides = provides,
    options_schema = options_schema,
    doc = doc
  )

  invisible(name)
}

get_posthoc <- function(name) {
  ph <- .gds_posthoc[[name]]
  if (is.null(ph)) {
    stop("Post-hoc method not found: ", name)
  }
  ph
}
```

**Kernel Signature:**

```r
# Standard post-hoc kernel signature
# z:      [samples]           (group-level test statistic)
# graph:  adjacency object    (from space: vox26, mesh, parcel borders)
# alpha:  numeric             (family-wise error rate)
# opts:   named list          (method-specific options)
# return: named list with q-values, rejection mask, threshold

posthoc_kernel <- function(z, graph, alpha = 0.05, opts = list()) {
  # Returns: list(q = ..., reject = ..., thresh = ...)
}
```

**Adjacency Extraction:**

GDS provides spatial structure based on the active space:

```r
adjacency <- function(space, kind = c("auto", "vox6", "vox18", "vox26",
                                      "mesh", "parcel")) {
  kind <- match.arg(kind)

  if (kind == "auto") {
    kind <- switch(space$type,
      voxel = "vox26",
      surface = "mesh",
      parcels = "parcel",
      stop("No default adjacency for space type: ", space$type)
    )
  }

  switch(kind,
    vox6 = build_voxel_adjacency(space, connectivity = 6),
    vox18 = build_voxel_adjacency(space, connectivity = 18),
    vox26 = build_voxel_adjacency(space, connectivity = 26),
    mesh = build_mesh_adjacency(space),
    parcel = build_parcel_adjacency(space)
  )
}
```

### 12.4 fmrireg Kernel Implementations

In `fmrireg` package (no GDS dependency), expose stateless kernels:

**Meta-Analysis Kernels:**

```r
#' Fixed-effects meta-analysis kernel
#'
#' @param beta [subjects × samples] effect estimates
#' @param var [subjects × samples] variances
#' @param X [subjects × p] design matrix for moderators (NULL for intercept-only)
#' @param df [subjects × samples] degrees of freedom (NULL if not available)
#' @param opts named list (reserved for future options)
#' @return named list with beta_g, var_g, Q, I2, p
#' @export
meta_fixed_core <- function(beta, var, X = NULL, df = NULL, opts = list()) {
  n_subjects <- nrow(beta)
  n_samples <- ncol(beta)

  # Inverse-variance weights
  w <- 1 / var
  W <- colSums(w, na.rm = TRUE)

  # Weighted mean
  beta_g <- colSums(beta * w, na.rm = TRUE) / W

  # Pooled variance
  var_g <- 1 / W

  # Cochran's Q statistic
  Q <- colSums(w * (beta - rep(beta_g, each = n_subjects))^2, na.rm = TRUE)

  # I² heterogeneity
  df_Q <- n_subjects - 1
  I2 <- pmax(0, (Q - df_Q) / Q)

  # Two-tailed p-value
  se_g <- sqrt(var_g)
  z_g <- beta_g / se_g
  p <- 2 * pnorm(-abs(z_g))

  list(
    beta_g = matrix(beta_g, nrow = 1),
    var_g = matrix(var_g, nrow = 1),
    Q = matrix(Q, nrow = 1),
    I2 = matrix(I2, nrow = 1),
    p = matrix(p, nrow = 1)
  )
}

#' Random-effects meta-analysis kernel (DerSimonian-Laird, Paule-Mandel, REML)
#'
#' @param beta [subjects × samples] effect estimates
#' @param var [subjects × samples] variances
#' @param X [subjects × p] design matrix (NULL for intercept-only)
#' @param df [subjects × samples] degrees of freedom (NULL if not available)
#' @param opts named list; opts$tau2 = c("DL", "PM", "REML")
#' @return named list with beta_g, var_g, tau2, Q, I2, p
#' @export
meta_random_core <- function(beta, var, X = NULL, df = NULL,
                             opts = list(tau2 = "PM")) {
  tau2_method <- opts$tau2 %||% "PM"
  n_subjects <- nrow(beta)
  n_samples <- ncol(beta)

  # Compute τ² per sample
  tau2 <- numeric(n_samples)

  if (tau2_method == "DL") {
    # DerSimonian-Laird estimator (fast)
    w <- 1 / var
    W <- colSums(w, na.rm = TRUE)
    beta_fe <- colSums(beta * w, na.rm = TRUE) / W
    Q <- colSums(w * (beta - rep(beta_fe, each = n_subjects))^2, na.rm = TRUE)
    C <- W - colSums(w^2, na.rm = TRUE) / W
    tau2 <- pmax(0, (Q - (n_subjects - 1)) / C)

  } else if (tau2_method == "PM") {
    # Paule-Mandel estimator (bisection root-finding)
    for (i in seq_len(n_samples)) {
      tau2[i] <- estimate_tau2_pm(beta[, i], var[, i])
    }

  } else if (tau2_method == "REML") {
    # REML estimator (optimization)
    for (i in seq_len(n_samples)) {
      tau2[i] <- estimate_tau2_reml(beta[, i], var[, i])
    }
  }

  # Random-effects weights
  w_re <- 1 / (var + rep(tau2, each = n_subjects))
  W_re <- colSums(w_re, na.rm = TRUE)

  # Weighted mean
  beta_g <- colSums(beta * w_re, na.rm = TRUE) / W_re

  # Pooled variance
  var_g <- 1 / W_re

  # Heterogeneity statistics
  w_fe <- 1 / var
  W_fe <- colSums(w_fe, na.rm = TRUE)
  Q <- colSums(w_fe * (beta - rep(beta_g, each = n_subjects))^2, na.rm = TRUE)
  df_Q <- n_subjects - 1
  I2 <- pmax(0, (Q - df_Q) / Q)

  # Two-tailed p-value
  se_g <- sqrt(var_g)
  z_g <- beta_g / se_g
  p <- 2 * pnorm(-abs(z_g))

  list(
    beta_g = matrix(beta_g, nrow = 1),
    var_g = matrix(var_g, nrow = 1),
    tau2 = matrix(tau2, nrow = 1),
    Q = matrix(Q, nrow = 1),
    I2 = matrix(I2, nrow = 1),
    p = matrix(p, nrow = 1)
  )
}
```

**Evidence Combiner Kernels:**

```r
#' Stouffer's method (z-score combination)
#'
#' @param z [subjects × samples] z-scores
#' @param w [subjects × samples] weights (NULL for equal weights)
#' @param opts named list (reserved)
#' @return named list with z_g
#' @export
combine_stouffer_core <- function(z, w = NULL, opts = list()) {
  n_subjects <- nrow(z)
  n_samples <- ncol(z)

  if (is.null(w)) {
    # Equal weights
    z_g <- colSums(z, na.rm = TRUE) / sqrt(n_subjects)
  } else {
    # Weighted
    numerator <- colSums(z * w, na.rm = TRUE)
    denominator <- sqrt(colSums(w^2, na.rm = TRUE))
    z_g <- numerator / denominator
  }

  list(z_g = matrix(z_g, nrow = 1))
}

#' Fisher's method (p-value combination)
#'
#' @param p [subjects × samples] p-values
#' @param opts named list (reserved)
#' @return named list with chi2, df, p_g
#' @export
combine_fisher_core <- function(p, opts = list()) {
  n_subjects <- nrow(p)

  # Chi-square statistic
  chi2 <- -2 * colSums(log(p), na.rm = TRUE)
  df <- 2 * n_subjects

  # Combined p-value
  p_g <- pchisq(chi2, df, lower.tail = FALSE)

  list(
    chi2 = matrix(chi2, nrow = 1),
    df = matrix(df, nrow = 1),
    p_g = matrix(p_g, nrow = 1)
  )
}

#' Lancaster's method (weighted Fisher)
#'
#' @param p [subjects × samples] p-values
#' @param w [subjects × samples] weights
#' @param opts named list (reserved)
#' @return named list with chi2, df, p_g
#' @export
combine_lancaster_core <- function(p, w, opts = list()) {
  # Weighted log-sum
  chi2 <- -2 * colSums(w * log(p), na.rm = TRUE)

  # Effective degrees of freedom
  df <- 2 * colSums(w, na.rm = TRUE)^2 / colSums(w^2, na.rm = TRUE)

  # Combined p-value
  p_g <- pchisq(chi2, df, lower.tail = FALSE)

  list(
    chi2 = matrix(chi2, nrow = 1),
    df = matrix(df, nrow = 1),
    p_g = matrix(p_g, nrow = 1)
  )
}
```

**Spatial FDR Kernel:**

```r
#' Structure-adaptive spatial FDR
#'
#' @param z [samples] group-level z-scores
#' @param graph adjacency object from GDS
#' @param alpha numeric FDR level
#' @param opts named list; opts$method = c("BH", "BY", "adaptive")
#' @return named list with q-values, rejection mask, threshold
#' @export
spatial_fdr_core <- function(z, graph, alpha = 0.05,
                             opts = list(method = "adaptive")) {
  method <- opts$method %||% "adaptive"
  n_samples <- length(z)

  # Convert z to two-tailed p-values
  p <- 2 * pnorm(-abs(z))

  if (method == "adaptive") {
    # Structure-adaptive weighting based on local spatial autocorrelation
    weights <- compute_spatial_weights(z, graph)
    q <- weighted_bh(p, weights)
  } else if (method == "BH") {
    # Standard Benjamini-Hochberg
    q <- p.adjust(p, method = "BH")
  } else if (method == "BY") {
    # Benjamini-Yekutieli (dependency control)
    q <- p.adjust(p, method = "BY")
  }

  # Rejection mask
  reject <- q <= alpha

  # Discovery threshold (largest p among rejected)
  if (any(reject)) {
    thresh <- max(p[reject])
  } else {
    thresh <- NA_real_
  }

  list(
    q = q,
    reject = reject,
    thresh = thresh
  )
}

# Helper: compute spatial weights based on local autocorrelation
compute_spatial_weights <- function(z, graph) {
  # ... implementation of structure-adaptive weighting ...
  # Uses graph to compute local Moran's I or similar statistic
}
```

### 12.5 fmrireg.gds Bridge Package

**CRITICAL: This package contains ZERO statistical algorithms.**

The bridge package's sole purpose is to:
1. Call `gdsfmri::register_reducer()` and `gdsfmri::register_posthoc()` on load
2. Provide thin wrapper functions that immediately delegate to `fmrireg::*_core()` functions
3. Optionally provide convenience sugar (e.g., `fmri_meta()` that builds a Plan)

All statistical computation happens in `fmrireg`. All data orchestration happens in `gdsfmri`. This package is **only the wiring** between them.

**Package Structure:**

```
fmrireg.gds/
  DESCRIPTION        (Imports: fmrireg, gdsfmri)
  NAMESPACE
  R/
    zzz.R            (registration via .onLoad)
    sugar.R          (optional convenience functions)
  man/
  tests/
```

**DESCRIPTION:**

```
Package: fmrireg.gds
Title: Bridge Between fmrireg and gdsfmri
Version: 0.1.0
Authors@R: ...
Description: Registers fmrireg's meta-analytic and inferential kernels
    with gdsfmri's reducer and post-hoc APIs. Enables seamless streaming
    meta-analysis without object conversion.
Depends: R (>= 4.0.0)
Imports:
    fmrireg (>= 1.0.0),
    gdsfmri (>= 0.1.0)
Suggests:
    testthat
License: MIT + file LICENSE
```

**R/zzz.R (Registration):**

```r
.onLoad <- function(libname, pkgname) {
  # Version guards
  if (!requireNamespace("gdsfmri", quietly = TRUE) ||
      !requireNamespace("fmrireg", quietly = TRUE)) {
    return(invisible())
  }

  # Optional: check minimum versions
  if (utils::packageVersion("fmrireg") < "1.0.0") {
    warning("fmrireg.gds requires fmrireg >= 1.0.0")
    return(invisible())
  }

  # Register meta-analysis reducers
  try(gdsfmri::register_reducer(
    name = "meta:fixed",
    fun = function(beta, var, X, df, opts) {
      fmrireg::meta_fixed_core(beta, var, X, df, opts)
    },
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "Q", "I2", "p"),
    doc = "Fixed-effects inverse-variance weighted meta-analysis"
  ), silent = TRUE)

  try(gdsfmri::register_reducer(
    name = "meta:dl",
    fun = function(beta, var, X, df, opts) {
      fmrireg::meta_random_core(beta, var, X, df,
                                opts = c(opts, list(tau2 = "DL")))
    },
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "tau2", "Q", "I2", "p"),
    doc = "Random-effects meta-analysis (DerSimonian-Laird τ²)"
  ), silent = TRUE)

  try(gdsfmri::register_reducer(
    name = "meta:pm",
    fun = function(beta, var, X, df, opts) {
      fmrireg::meta_random_core(beta, var, X, df,
                                opts = c(opts, list(tau2 = "PM")))
    },
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "tau2", "Q", "I2", "p"),
    doc = "Random-effects meta-analysis (Paule-Mandel τ²)"
  ), silent = TRUE)

  try(gdsfmri::register_reducer(
    name = "meta:reml",
    fun = function(beta, var, X, df, opts) {
      fmrireg::meta_random_core(beta, var, X, df,
                                opts = c(opts, list(tau2 = "REML")))
    },
    requires = c("beta", "var"),
    provides = c("beta_g", "var_g", "tau2", "Q", "I2", "p"),
    doc = "Random-effects meta-analysis (REML τ²)"
  ), silent = TRUE)

  # Register evidence combiners
  try(gdsfmri::register_reducer(
    name = "combine:stouffer",
    fun = function(beta, var, X, df, opts) {
      # GDS passes z in opts$z_block when requires = "z"
      z_block <- opts$z_block
      weights <- opts$weights
      fmrireg::combine_stouffer_core(z_block, weights, opts)
    },
    requires = "z",
    provides = "z_g",
    doc = "Stouffer's method for z-score combination"
  ), silent = TRUE)

  try(gdsfmri::register_reducer(
    name = "combine:fisher",
    fun = function(beta, var, X, df, opts) {
      p_block <- opts$p_block
      fmrireg::combine_fisher_core(p_block, opts)
    },
    requires = "p",
    provides = c("chi2", "df", "p_g"),
    doc = "Fisher's method for p-value combination"
  ), silent = TRUE)

  try(gdsfmri::register_reducer(
    name = "combine:lancaster",
    fun = function(beta, var, X, df, opts) {
      p_block <- opts$p_block
      weights <- opts$weights
      fmrireg::combine_lancaster_core(p_block, weights, opts)
    },
    requires = "p",
    provides = c("chi2", "df", "p_g"),
    doc = "Lancaster's weighted Fisher method"
  ), silent = TRUE)

  # Register post-hoc methods
  try(gdsfmri::register_posthoc(
    name = "fdr:spatial",
    fun = function(z, graph, alpha, opts) {
      fmrireg::spatial_fdr_core(z, graph, alpha, opts)
    },
    requires = "z",
    provides = c("q", "reject", "thresh"),
    doc = "Structure-adaptive spatial FDR"
  ), silent = TRUE)

  invisible()
}
```

**R/sugar.R (Optional Convenience):**

```r
#' Convenience wrapper for fmri_meta on GDS plans
#'
#' @param plan gds_plan object
#' @param method meta-analysis method ("meta:pm", "meta:reml", "meta:fixed")
#' @param formula model formula for moderators (default: ~ 1)
#' @param alpha FDR level for spatial post-hoc (NULL to skip)
#' @param ... additional options passed to reducer
#' @return gds_plan with reduce() and optionally posthoc() added
#' @export
fmri_meta <- function(plan, method = "meta:pm", formula = ~ 1,
                     alpha = NULL, ...) {
  stopifnot(inherits(plan, "gds_plan"))

  # Add reduce operation
  plan <- gdsfmri::reduce(plan, method = method, formula = formula,
                         options = list(...))

  # Optionally add spatial FDR
  if (!is.null(alpha)) {
    plan <- gdsfmri::posthoc(plan, method = "fdr:spatial",
                            options = list(alpha = alpha))
  }

  plan
}
```

### 12.6 GDS Execution: Calling Registered Reducers

Inside `gdsfmri`, the `compute()` function calls registered reducers during block streaming:

**Reducer Execution Hook:**

```r
# Inside gdsfmri/R/compute.R

execute_reduce_node <- function(node, assays_block, col_data, space) {
  # node: op_reduce(method, formula, options)
  # assays_block: current block [samples × subjects × contrasts]
  # col_data: data.frame with subject-level covariates

  method <- node$method
  formula <- as.formula(node$formula)
  options <- node$options

  # Get reducer
  reducer <- get_reducer(method)

  # Build design matrix from formula + col_data
  if (identical(formula, ~ 1)) {
    X <- NULL
  } else {
    X <- model.matrix(formula, data = col_data)
  }

  # Extract required assays
  beta <- if ("beta" %in% reducer$requires) assays_block$beta else NULL
  var <- if ("var" %in% reducer$requires) assays_block$var else NULL
  df <- if ("df" %in% reducer$requires) assays_block$df else NULL

  # For evidence combiners (z/p-only), pass via options
  if ("z" %in% reducer$requires) {
    options$z_block <- assays_block$z
  }
  if ("p" %in% reducer$requires) {
    options$p_block <- assays_block$p
  }

  # Process each contrast separately
  n_contrasts <- dim(assays_block$beta)[3]
  results <- vector("list", n_contrasts)

  for (k in seq_len(n_contrasts)) {
    # Extract contrast slice: [samples × subjects]
    beta_k <- if (!is.null(beta)) t(beta[, , k]) else NULL  # [subjects × samples]
    var_k <- if (!is.null(var)) t(var[, , k]) else NULL
    df_k <- if (!is.null(df)) t(df[, , k]) else NULL

    # Call reducer kernel
    result_k <- reducer$fun(beta_k, var_k, X, df_k, options)

    # result_k: named list with [1 × samples] matrices
    results[[k]] <- result_k
  }

  # Combine results across contrasts: [samples × 1 × contrasts]
  combined <- list()
  for (assay_name in reducer$provides) {
    combined[[assay_name]] <- array(
      sapply(results, function(r) r[[assay_name]]),
      dim = c(nrow(results[[1]][[assay_name]]), 1, n_contrasts)
    )
  }

  combined
}
```

**Post-Hoc Execution Hook:**

```r
# Inside gdsfmri/R/compute.R

execute_posthoc_node <- function(node, assays_realized, space) {
  # node: op_posthoc(method, options)
  # assays_realized: realized GDS after reduce
  # space: gds_space object

  method <- node$method
  options <- node$options
  alpha <- options$alpha %||% 0.05

  # Get post-hoc method
  ph <- get_posthoc(method)

  # Build spatial graph
  graph <- adjacency(space, kind = options$adjacency %||% "auto")

  # Apply to each contrast
  n_contrasts <- dim(assays_realized$z)[3]
  results <- vector("list", n_contrasts)

  for (k in seq_len(n_contrasts)) {
    z_k <- assays_realized$z[, 1, k]  # [samples]

    # Call post-hoc kernel
    result_k <- ph$fun(z_k, graph, alpha, options)

    results[[k]] <- result_k
  }

  # Add post-hoc results as new assays
  for (assay_name in ph$provides) {
    assays_realized[[assay_name]] <- array(
      sapply(results, function(r) r[[assay_name]]),
      dim = c(length(results[[1]][[assay_name]]), 1, n_contrasts)
    )
  }

  assays_realized
}
```

### 12.7 User-Facing Workflow

**Example 1: Fixed-effects meta-analysis with spatial FDR**

```r
library(gdsfmri)
library(fmrireg.gds)  # Loads and registers reducers

plan <- gds(Sys.glob("sub-*/beta.nii.gz")) %>%
  subset(subject = paste0("sub-", sprintf("%02d", 1:30))) %>%
  derive(c("var", "t")) %>%
  align(map = mni_warp_family) %>%
  mask(MaskPolicy(scope = "group", rule = "threshold", threshold = 0.95)) %>%
  reduce(method = "meta:fixed") %>%
  posthoc(method = "fdr:spatial", options = list(alpha = 0.05)) %>%
  write_out("group_fixed_fdr.h5", format = "h5")

result <- compute(plan, sink = "memory")

# Access results
beta_group <- assay(result, "beta_g")  # [samples × 1 × contrasts]
q_values <- assay(result, "q")         # FDR q-values
reject_mask <- assay(result, "reject") # Significant voxels
```

**Example 2: Random-effects with Paule-Mandel estimator**

```r
plan <- gds("roi_stats.csv") %>%
  subset(contrast = c("A_main", "B_main", "A:B")) %>%
  derive(c("var", "t")) %>%
  reduce(method = "meta:pm") %>%  # Paule-Mandel τ²
  write_out("roi_meta_pm.csv", format = "csv")

result <- compute(plan)

# Access heterogeneity statistics
tau2 <- assay(result, "tau2")  # Between-study variance
I2 <- assay(result, "I2")      # Heterogeneity proportion
Q <- assay(result, "Q")        # Cochran's Q
```

**Example 3: Evidence combination (Stouffer) for z-scores only**

```r
# When you only have t/z statistics (no effect scale)
plan <- gds("sub-*/t_stat.nii.gz") %>%
  derive("z") %>%  # Convert t to z
  reduce(method = "combine:stouffer",
         options = list(weights = "equal")) %>%
  write_out("stouffer_combined.h5", format = "h5")

result <- compute(plan)
z_combined <- assay(result, "z_g")
```

**Example 4: Meta-regression with moderators**

```r
# col_data has subject-level covariates: age, group
plan <- gds(source) %>%
  derive(c("var", "t")) %>%
  align(map = mni_warp_family) %>%
  reduce(method = "meta:reml", formula = ~ group + age) %>%
  write_out("meta_regression.h5", format = "h5")

result <- compute(plan)

# Extract moderator effects
beta_intercept <- assay(result, "beta_g")[, 1, , 1]  # Intercept
beta_group <- assay(result, "beta_g")[, 1, , 2]      # Group effect
beta_age <- assay(result, "beta_g")[, 1, , 3]        # Age effect
```

### 12.8 Migration Checklist

**Phase A: Kernel Extraction (fmrireg) – 1-2 weeks**

- [ ] Extract `meta_fixed_core()` from existing `fmri_meta()` code
- [ ] Extract `meta_random_core(tau2 = "DL|PM|REML")` with consistent signature
- [ ] Implement `combine_stouffer_core()`, `combine_fisher_core()`, `combine_lancaster_core()`
- [ ] Extract `spatial_fdr_core()` from existing spatial FDR code
- [ ] Export kernels with roxygen2 documentation
- [ ] Add unit tests for kernels (pure function tests, no GDS dependency)

**Phase B: Bridge Package (fmrireg.gds) – 1 week**

- [ ] Create package skeleton with DESCRIPTION (Imports: fmrireg, gdsfmri)
- [ ] Implement `.onLoad()` registration in R/zzz.R
- [ ] Add optional sugar: `fmri_meta(plan, method, formula, alpha)`
- [ ] Add tests: verify registration succeeds, call reducers via GDS
- [ ] Document usage examples

**Phase C: GDS Registry & Executor (gdsfmri) – 1-2 weeks**

- [ ] Implement `register_reducer()`, `get_reducer()`, `list_reducers()`
- [ ] Implement `register_posthoc()`, `get_posthoc()`
- [ ] Add `execute_reduce_node()` hook in compute pipeline
- [ ] Add `execute_posthoc_node()` hook
- [ ] Implement `adjacency()` function for spatial graphs
- [ ] Add provenance tracking: record reducer name + package versions

**Phase D: Validation & Documentation – 2-3 weeks**

- [ ] Round-trip statistical parity tests on known datasets
- [ ] Benchmark streaming reducers vs legacy pipelines (speed & memory)
- [ ] Add vignette: "Using fmrireg with GDS"
- [ ] Add vignette: "Writing custom reducers"
- [ ] Document reducer/post-hoc API in gdsfmri
- [ ] CI matrix: test gdsfmri alone, with fmrireg, with fmrireg.gds

**Total Estimated Effort: 5-8 weeks**

### 12.9 Extension: F-Statistics and Future Test Families

**Supporting F-statistics (multi-parameter tests):**

F-statistics don't have an effect scale. Two approaches:

1. **If you have beta/var (preferred):**
   - Run meta-analysis on beta/var (effect scale preserved)
   - Derive F from pooled beta/var at group level

2. **If you only have F, df1, df2 (no effect scale):**
   - Use evidence combiners:
     - Convert F → p: `p = pf(F, df1, df2, lower.tail = FALSE)`
     - Apply `combine:fisher` or `combine:lancaster` on p
   - Or convert F → z and use `combine:stouffer`

**Reducer for F-to-p combiner:**

```r
# In fmrireg
combine_f_to_fisher_core <- function(F_stat, df1, df2, opts = list()) {
  # Convert F to p-values
  p <- pf(F_stat, df1, df2, lower.tail = FALSE)

  # Apply Fisher's method
  combine_fisher_core(p, opts)
}

# Register in fmrireg.gds
try(gdsfmri::register_reducer(
  name = "combine:f_fisher",
  fun = function(beta, var, X, df, opts) {
    F_block <- opts$F_block
    df1_block <- opts$df1_block
    df2_block <- opts$df2_block
    fmrireg::combine_f_to_fisher_core(F_block, df1_block, df2_block, opts)
  },
  requires = c("F", "df1", "df2"),
  provides = c("chi2", "df", "p_g"),
  doc = "Fisher's method for F-statistic combination"
), silent = TRUE)
```

**Future test families (Bayes factors, Wald/score tests):**

Add new reducers by:
1. Implementing kernel in fmrireg (or any package)
2. Registering via `register_reducer()` (declaring `requires` and `provides`)
3. GDS handles data movement; kernel handles statistics

### 12.10 Advantages of This Architecture

**1. No impedance mismatch:**
   - Data never bounces between `group_data` and `gds` formats
   - Streaming blocks go directly to statistical kernels

**2. Single source of truth:**
   - All algorithms live in fmrireg
   - No drift, no duplication
   - fmrireg can evolve independently

**3. Composability:**
   - Same reducers work on voxels, parcels, surfaces, basis
   - GDS handles alignment/masks/mapping; fmrireg handles statistics

**4. Testability:**
   - Kernels are pure functions: easy to unit test
   - No GDS dependency in fmrireg tests
   - Integration tests in fmrireg.gds

**5. Extensibility:**
   - Any package can register reducers/post-hoc methods
   - No need to modify gdsfmri core

**6. Performance:**
   - Block streaming with fused operations
   - No intermediate object creation
   - Parallelizable across contrasts and blocks

**7. Provenance:**
   - GDS records reducer name + fmrireg version
   - Reproducible pipelines with version locking

---

## End of Technical Specification

This completes the comprehensive technical specification for the `gdsfmri` R package. The document provides:

1. ✅ Complete API signatures for all 5 nouns and 8 verbs
2. ✅ Precise dimensional contracts for every operation
3. ✅ Mathematical formulas for statistical operations
4. ✅ Storage adapter interface and HDF5 layout
5. ✅ Implementation file structure with function skeletons
6. ✅ Comprehensive test specifications
7. ✅ Migration strategy for backward compatibility

**Next Steps for Implementation:**
1. Set up package skeleton with proper DESCRIPTION and NAMESPACE
2. Implement core objects (GDS, Space, Map) with validators
3. Implement lazy execution system (Plan, optimizer)
4. Implement verbs in order: gds() → subset() → derive() → compute()
5. Implement adapters (tabular first, then NIfTI, then HDF5)
6. Add statistical operations (variance propagation, meta-analysis)
7. Implement alignment and mask operations
8. Write comprehensive tests for each module
9. Add documentation and vignettes
10. Integration testing with real fMRI data
