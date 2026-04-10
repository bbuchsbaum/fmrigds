# fmrigds Compute Pipeline Analysis

**Date:** 2025-10-29
**Purpose:** Comprehensive trace of data flow from source to results in the fmrigds group-level analysis pipeline

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Data Flow Diagram](#data-flow-diagram)
3. [Pipeline Components](#pipeline-components)
4. [Step-by-Step Execution Trace](#step-by-step-execution-trace)
5. [Plan Structure Documentation](#plan-structure-documentation)
6. [Adapter Interface Specification](#adapter-interface-specification)
7. [Reducer Interface Specification](#reducer-interface-specification)
8. [Posthoc Interface Specification](#posthoc-interface-specification)
9. [Compute Result Structure](#compute-result-structure)
10. [Error Handling & Diagnostics](#error-handling--diagnostics)
11. [Common Failure Modes](#common-failure-modes)
12. [Debugging Recommendations](#debugging-recommendations)

---

## Executive Summary

The fmrigds pipeline implements a **lazy, plan-based computation model** for group-level fMRI analysis. Key characteristics:

- **Lazy Evaluation**: Operations build a plan without executing until `compute()` is called
- **Adapter Pattern**: Pluggable storage backends (H5, NIfTI, CSV/tabular, fmristore)
- **Registry Pattern**: Reducers and posthoc methods are registered functions
- **Immutable Plans**: Each verb returns a new plan with an added operation node
- **3D Array Structure**: All data flows as `[sample × subject × contrast]` arrays

**Critical Insight**: Assays are produced dynamically during compute execution. The `provides` field in reducer registrations determines what assays will be available, but this information is NOT stored in the plan metadata - you must inspect reducer registrations to know what assays a reducer produces.

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                         USER INTERFACE                               │
│  group_data(paths, format="h5") → fmri_meta() → fmri_ttest()       │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      1. PLAN CONSTRUCTION                            │
│                     (R/gds-verb.R, R/plan.R)                         │
├─────────────────────────────────────────────────────────────────────┤
│  gds(source, format) → gds_plan                                      │
│    ↓                                                                 │
│  1a. Detect adapter (or use explicit format)                        │
│  1b. adapter$open(source) → handle                                   │
│  1c. adapter$probe(handle, ...) → metadata                          │
│       Returns: assays, dims, subjects, contrasts, space, col_data   │
│  1d. adapter$close(handle)                                           │
│  1e. Create gds_source(adapter, source, probe_result)               │
│  1f. Create gds_plan(source, nodes=[], meta=list())                 │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                   2. PLAN MODIFICATION (LAZY)                        │
│              (R/verb-reduce.R, R/verb-posthoc.R)                     │
├─────────────────────────────────────────────────────────────────────┤
│  plan |> reduce(method, formula, weights) → plan                     │
│    Appends: op_reduce(method, weights, by, options, formula)        │
│                                                                      │
│  plan |> posthoc(method, options) → plan                            │
│    Appends: list(op="posthoc", method=method, options=options)      │
│                                                                      │
│  Each verb:                                                          │
│    - Calls as_plan(x) to ensure plan object                         │
│    - Calls add_op(plan, node) to append operation                   │
│    - Returns new plan (immutable pattern)                           │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      3. COMPUTE EXECUTION                            │
│                      (R/compute.R)                                   │
├─────────────────────────────────────────────────────────────────────┤
│  compute(plan) → gds object                                          │
│                                                                      │
│  3a. Optimize plan (plan-optimizer.R)                               │
│      - Fuse adjacent operations                                     │
│      - Push subsets to adapter read                                 │
│                                                                      │
│  3b. Open adapter and read all source assays                        │
│      handle <- adapter$open(source)                                 │
│      arrays <- adapter$read(handle, assays=all_assays, block=NULL) │
│      adapter$close(handle)                                           │
│      Result: Named list of 3D arrays [sample × subject × contrast]  │
│                                                                      │
│  3c. Apply plan nodes sequentially                                  │
│      For each node in plan$nodes:                                   │
│        - subset_axis: Slice arrays by dimension                     │
│        - derive: Add computed assays (z, p, t, se, var)             │
│        - align_to_group: Transform to common space                  │
│        - mask_policy: Apply masking rules                           │
│        - map: Transform spatial coordinates                         │
│        - reduce: Apply meta-analysis (SEE SECTION 3d)               │
│        - posthoc: Apply post-hoc corrections (SEE SECTION 3e)       │
│        - write: Queue export operation                              │
│                                                                      │
│  3d. Reduce Execution (R/reduce-exec.R)                             │
│      apply_reduce(node, arrays, weights, subjects, col_data)        │
│        ↓                                                             │
│      Get reducer: get_reducer(.normalize_reducer_name(method))     │
│        ↓                                                             │
│      Ensure required inputs: .ensure_required_arrays(arrays, req)   │
│        - Auto-derives z from {t,df} or {p,beta}                     │
│        - Auto-derives p from {t,df} or {z}                          │
│        ↓                                                             │
│      For each contrast k:                                           │
│        - Slice arrays[,,k] → subjects × samples matrices            │
│        - Build design matrix X from formula + col_data if needed    │
│        - Call reducer$fun(beta, var, X, z, p, df, df1, df2, opts)  │
│        - Collect outputs specified in reducer$provides              │
│        ↓                                                             │
│      Return new arrays with subject axis collapsed to "meta"        │
│                                                                      │
│  3e. Posthoc Execution (R/posthoc-exec.R)                           │
│      apply_posthoc(node, arrays)                                    │
│        ↓                                                             │
│      Get method: get_posthoc(method)                                │
│        ↓                                                             │
│      Ensure required inputs (e.g., derive p if needed)              │
│        ↓                                                             │
│      Call ph$fun(arrays, options) → list of new assays              │
│        ↓                                                             │
│      Merge new assays into arrays                                   │
│                                                                      │
│  3f. Build final GDS object                                         │
│      new_gds(                                                        │
│        assays = arrays,                                             │
│        space = space,                                               │
│        subjects = subjects,  # "meta" after reduce                  │
│        contrasts = contrasts,                                       │
│        col_data = col_data,                                         │
│        metadata = metadata                                          │
│      )                                                               │
│                                                                      │
│  3g. Handle sink operations                                         │
│      if (sink == "h5"): Write to HDF5                               │
│      if (write ops queued): Export to specified formats             │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     4. RESULT EXTRACTION                             │
│                      (R/gds-class.R)                                 │
├─────────────────────────────────────────────────────────────────────┤
│  gds <- compute(plan)                                                │
│                                                                      │
│  Accessors:                                                          │
│    assay(gds, "beta")     → 3D array [sample × 1 × contrast]        │
│    assay(gds, "se")       → 3D array                                │
│    assay(gds, "z_g")      → Group-level z-scores                    │
│    assay(gds, "p_g")      → Group-level p-values                    │
│    assay(gds, "Q")        → Heterogeneity Q statistic               │
│    assay(gds, "I2")       → I² heterogeneity index                  │
│    assay(gds, "tau2")     → Between-study variance (RE models)      │
│    assay(gds, "q")        → FDR-adjusted q-values (after posthoc)   │
│    assay(gds, "coef:X1")  → Regression coefficient for X1           │
│    assay(gds, "se_coef:X1") → SE of coefficient for X1              │
│                                                                      │
│  Other accessors:                                                    │
│    assays(gds)            → Named list of all arrays                │
│    subjects(gds)          → "meta" (or subject IDs if no reduce)    │
│    contrasts(gds)         → Character vector                        │
│    space(gds)             → Space descriptor                        │
│    col_data(gds)          → Subject-level covariates (if any)       │
│    metadata(gds)          → Provenance, software, units, etc.       │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Pipeline Components

### Core Architecture

The pipeline consists of six key components:

1. **Plans** (`gds_plan`): Lazy computation graphs
2. **Sources** (`gds_source`): Adapter bindings with probe metadata
3. **Adapters**: Storage backend interfaces (H5, NIfTI, tabular, fmristore)
4. **Reducers**: Meta-analysis kernels (registered functions)
5. **Posthoc Methods**: Post-hoc correction methods (registered functions)
6. **GDS Objects**: Realized group data with assays + metadata

### File Organization

```
fmrigds/R/
├── gds-verb.R           # gds() entry point, adapter dispatch
├── plan.R               # Plan constructors, operation nodes
├── compute.R            # Compute execution, node application
├── adapter-*.R          # Storage backend implementations
├── adapter-registry.R   # Adapter detection and lookup
├── verb-reduce.R        # reduce() verb
├── verb-posthoc.R       # posthoc() verb
├── reduce-exec.R        # Reducer execution engine
├── posthoc-exec.R       # Posthoc execution engine
├── reducer-registry.R   # Reducer registration system
├── posthoc-registry.R   # Posthoc registration system
├── reducers-core.R      # Built-in reducers (meta:fe, meta:re, etc.)
├── derive-stats.R       # Assay derivation (z, p, t, se, var)
└── gds-class.R          # GDS object constructors and accessors
```

---

## Step-by-Step Execution Trace

### Example: Simple Fixed-Effect Meta-Analysis

```r
# User code
df <- data.frame(
  sample = rep(c("ROI_1","ROI_2"), each = 4),
  subject = rep(c("s1","s2"), times = 4),
  contrast = rep(c("c1","c2"), each = 2, times = 2),
  beta = rnorm(8, mean = 0.5),
  var = runif(8, 0.01, 0.1)
)
write.csv(df, "/tmp/test.csv", row.names = FALSE)

result <- gds("/tmp/test.csv") |>
  reduce(method = "fixed", weights = "1/var") |>
  posthoc("fdr:bh") |>
  compute()

beta_g <- assay(result, "beta")
q_vals <- assay(result, "q")
```

### Execution Trace

#### Step 1: Plan Construction (`gds()`)

**File:** `R/gds-verb.R::gds()`

```r
# 1.1 Adapter detection
source <- "/tmp/test.csv"
format <- "auto"
adapter_name <- detect_adapter(source, prefer = NULL)
  # Checks all registered adapters via adapter$detect(source)
  # tabular adapter returns score > 0 for .csv files
  # Result: "tabular"

adapter <- get_adapter("tabular")
  # Retrieves from .adapter_registry environment

# 1.2 Open and probe
handle <- adapter$open(source)
  # tabular adapter: Reads CSV into memory, returns list(data=df)

probe_result <- adapter$probe(handle,
                               effect_cols = c(beta="beta", var="var"),
                               subject_col = "subject",
                               sample_col = "sample",
                               contrast_col = "contrast")
  # Returns:
  # list(
  #   assays = c("beta", "var"),
  #   dims = gds_dims(sample=2, subject=2, contrast=2),
  #   subjects = c("s1", "s2"),
  #   contrasts = c("c1", "c2"),
  #   space = space_sample_labels(c("ROI_1", "ROI_2")),
  #   metadata = list(),
  #   columns = list(effect_cols=..., subject_col=..., etc.)
  # )

adapter$close(handle)
  # tabular adapter: No-op

# 1.3 Build source and plan
src <- gds_source("tabular", "/tmp/test.csv", probe_result)
plan <- gds_plan(source = src, nodes = list(), meta = list())
plan$meta$subjects <- c("s1", "s2")
plan$metadata <- list(dims = probe_result$dims)

# Result: gds_plan with empty nodes
```

#### Step 2: Add Reduce Operation

**File:** `R/verb-reduce.R::reduce()`

```r
# Input plan from step 1
plan <- reduce(plan, method = "fixed", weights = "1/var")

# Inside reduce():
method <- "fixed"  # normalized to "meta:fe" later
weights <- "1/var"
by <- "contrast"
options <- list()
formula <- NULL

node <- op_reduce(method="fixed", weights="1/var", by="contrast",
                  options=list(), formula=NULL)
  # Creates: list(op="reduce", method="fixed", weights="1/var",
  #               by="contrast", options=list(), formula=NULL)

plan$nodes <- c(plan$nodes, list(node))

# Result: plan with 1 node
```

#### Step 3: Add Posthoc Operation

**File:** `R/verb-posthoc.R::posthoc()`

```r
plan <- posthoc(plan, method = "fdr:bh", options = list())

node <- list(op = "posthoc", method = "fdr:bh", options = list())
plan$nodes <- c(plan$nodes, list(node))

# Result: plan with 2 nodes
```

#### Step 4: Compute Execution

**File:** `R/compute.R::compute()`

```r
result <- compute(plan)

# 4.1 Optimize plan
plan <- optimize_plan(plan)
  # May fuse operations, but none in this simple case

# 4.2 Open adapter and read
adapter <- get_adapter("tabular")
handle <- adapter$open("/tmp/test.csv")

arrays <- adapter$read(
  handle,
  assays = c("beta", "var"),
  block = NULL,
  effect_cols = c(beta="beta", var="var"),
  subject_col = "subject",
  sample_col = "sample",
  contrast_col = "contrast"
)
# Returns:
# list(
#   beta = array(dim = c(2, 2, 2)),  # [sample × subject × contrast]
#   var  = array(dim = c(2, 2, 2))
# )

adapter$close(handle)

# 4.3 Apply plan nodes
node_result <- .apply_plan_nodes(arrays, plan, space, subjects, col_data=NULL)

# Node 1: Reduce (method="fixed")
# File: R/reduce-exec.R::apply_reduce()

name <- .normalize_reducer_name("fixed")  # "meta:fe"
reducer <- get_reducer("meta:fe")
  # Retrieves from .gds_reducers registry
  # reducer = list(
  #   name = "meta:fe",
  #   fun = function(beta, var, X, z, p, df, df1, df2, opts) {...},
  #   requires = c("beta", "var"),
  #   provides = c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2")
  # )

# Ensure required arrays
arrays <- .ensure_required_arrays(arrays, c("beta", "var"))
  # beta and var already present, no derivation needed

# Extract dimensions
dims <- c(2, 2, 2)  # sample, subject, contrast
n_samples <- 2; n_subject <- 2; n_contrast <- 2

# Prepare output arrays
out_arrays <- lapply(c("beta_g","var_g","se_g","z_g","p_g","Q","I2"),
                     function(.) array(NA_real_, dim = c(2, 1, 2)))
names(out_arrays) <- c("beta_g","var_g","se_g","z_g","p_g","Q","I2")

# For each contrast
for (k in 1:2) {
  # Slice and transpose to [subjects × samples]
  beta_mat <- t(arrays$beta[,,k])  # 2×2 matrix
  var_mat  <- t(arrays$var[,,k])   # 2×2 matrix

  # Call reducer function
  res <- core_meta_fe_kernel(beta_mat, var_mat, X=NULL, df=NULL, opts=list())
  # Returns: list(beta_g=vec, var_g=vec, se_g=vec, z_g=vec, p_g=vec, Q=vec, I2=vec)
  # Each vec has length = n_samples = 2

  # Store results
  out_arrays$beta_g[, 1, k] <- res$beta_g
  out_arrays$var_g[, 1, k]  <- res$var_g
  out_arrays$se_g[, 1, k]   <- res$se_g
  out_arrays$z_g[, 1, k]    <- res$z_g
  out_arrays$p_g[, 1, k]    <- res$p_g
  out_arrays$Q[, 1, k]      <- res$Q
  out_arrays$I2[, 1, k]     <- res$I2
}

# Merge into arrays
arrays[names(out_arrays)] <- out_arrays

# Legacy field mapping
arrays$beta <- arrays$beta_g
arrays$var  <- arrays$var_g
arrays$se   <- arrays$se_g
arrays$z    <- arrays$z_g
arrays$p    <- arrays$p_g

subjects <- "meta"  # Collapsed from c("s1","s2") to group level

# Node 2: Posthoc (method="fdr:bh")
# File: R/posthoc-exec.R::apply_posthoc()

ph <- get_posthoc("fdr:bh")
  # Retrieves from .gds_posthoc registry
  # ph = list(
  #   name = "fdr:bh",
  #   fun = .posthoc_fdr_template("BH"),
  #   requires = c("p"),
  #   provides = c("q")
  # )

# Ensure required arrays (p already exists from reduce)
arrays <- .ensure_required_arrays(arrays, c("p"))

# Call posthoc function
res <- ph$fun(arrays, options = list())
  # For each [j, k] in arrays$p:
  #   adj <- p.adjust(p[,j,k], method="BH")
  #   out[,j,k] <- adj
  # Returns: list(q = array(dim=c(2,1,2)))

# Merge into arrays
arrays$q <- res$q

# 4.4 Build final GDS
gds <- new_gds(
  assays = arrays,  # Now contains: beta, var, se, z, p, q, beta_g, var_g, se_g, z_g, p_g, Q, I2
  space = space_sample_labels(c("ROI_1", "ROI_2")),
  subjects = "meta",
  contrasts = c("c1", "c2"),
  col_data = NULL,
  metadata = probe_result$metadata
)

# Add provenance
for (node in plan$nodes) {
  params <- # extract params from node
  gds$metadata <- add_provenance_node(gds$metadata, node$op, params)
}

# Result: Realized gds object
```

#### Step 5: Assay Extraction

**File:** `R/gds-class.R::assay()`

```r
beta_g <- assay(result, "beta")
  # Returns: result$assays[["beta"]]
  # 3D array with dim = c(2, 1, 2)  [sample × subject × contrast]
  # subject dimension = 1 because reduce collapsed to "meta"

q_vals <- assay(result, "q")
  # Returns: result$assays[["q"]]
  # 3D array with dim = c(2, 1, 2)
```

---

## Plan Structure Documentation

### gds_plan Object

```r
plan <- list(
  source = <gds_source>,  # See below
  nodes = list(),         # List of operation nodes
  meta = list(            # Planner bookkeeping
    subjects = character(),          # Subject IDs
    adapter_columns = list(),        # Adapter-specific column mappings
    map_families = list(),           # Spatial mapping families
    col_data = data.frame(),         # Subject-level covariates
    temporal_policy = NULL,          # Temporal aggregation policy
    contrast_matrix = NULL,          # Contrast matrix for linear combinations
    contrast_names = NULL            # Names for contrasts
  ),
  metadata = list(        # Legacy compatibility
    dims = <gds_dims>     # Dimensions from probe
  )
)
class(plan) <- "gds_plan"
```

### gds_source Object

```r
source <- list(
  adapter = "h5",                # Adapter name
  source = "/path/to/data.h5",  # Source specification
  probe = list(                 # Probe result from adapter$probe()
    assays = c("beta", "var", "se", "t", "df"),
    dims = gds_dims(sample=1000, subject=20, contrast=5),
    subjects = c("sub01", "sub02", ...),
    contrasts = c("FaceVsPlace", "SceneVsObject", ...),
    space = <space_voxel or space_parcels or space_sample_labels>,
    maps = list(),              # Registered map families
    metadata = list(),          # Source-specific metadata
    columns = list(),           # Column mappings (for tabular adapter)
    col_data = data.frame()     # Subject-level covariates (if available)
  ),
  hash = "xxhash64-digest"      # Stable hash of adapter + source
)
class(source) <- "gds_source"
```

### Operation Nodes

All nodes are plain lists with `op` field:

```r
# Subset operation
list(op = "subset_axis", sample = 1:100, subject = NULL, contrast = "c1")

# Derive operation
list(op = "derive", what = c("z", "p"), options = list(overwrite = TRUE))

# Alignment operation
list(op = "align_to_group", family = <MapFamily>, family_name = "MNI152")

# Mask operation
list(op = "mask_policy", policy = <MaskPolicy>)

# Map operation
list(op = "map", target_space = <space>, map = <matrix>,
     uncertainty = <UncertaintyRule>, combine = "mean")

# Reduce operation
list(op = "reduce", method = "meta:fe", weights = "1/var",
     by = "contrast", options = list(), formula = NULL)

# Posthoc operation
list(op = "posthoc", method = "fdr:bh", options = list(alpha = 0.05))

# Write operation
list(op = "write", path = "/output/results.h5", format = "h5", options = list())
```

---

## Adapter Interface Specification

### Adapter Structure

Each adapter is a named list in `.adapter_registry`:

```r
adapter <- list(
  name = "h5",

  # Detection function: source -> score in [0,1] or FALSE
  detect = function(source) {
    if (!is.character(source) || !file.exists(source)) return(FALSE)
    ext <- tolower(tools::file_ext(source))
    if (ext %in% c("h5", "hdf5")) return(1.0)
    return(FALSE)
  },

  # Open function: source -> handle
  open = function(source, mode = "r", ...) {
    list(file = hdf5r::H5File$new(source, mode = mode), path = source)
  },

  # Probe function: (handle, ...) -> metadata
  probe = function(handle, ...) {
    # Read metadata from storage
    # Return:
    list(
      assays = c("beta", "var", "se", ...),
      dims = gds_dims(sample = N, subject = S, contrast = K),
      subjects = character(S),
      contrasts = character(K),
      space = <space object>,
      maps = list(),            # Optional: registered map families
      metadata = list(),        # Optional: provenance, software, etc.
      columns = list(),         # Adapter-specific (tabular only)
      col_data = data.frame()   # Optional: subject-level covariates
    )
  },

  # Read function: (handle, assays, block, ...) -> named list of arrays
  read = function(handle, assays, block = NULL, ...) {
    # Read requested assays from storage
    # block is a list(sample=idx, subject=idx, contrast=idx) for subsetting
    # Return: named list of 3D arrays [sample × subject × contrast]
    list(
      beta = array(dim = c(N, S, K)),
      var  = array(dim = c(N, S, K)),
      ...
    )
  },

  # Close function: handle -> NULL
  close = function(handle) {
    # Clean up resources
    invisible(NULL)
  }
)
```

### Registered Adapters

| Adapter | Format | detect() Priority | Notes |
|---------|--------|-------------------|-------|
| `h5` | HDF5 | 1.0 if has `/gds` group | Standard GDS format |
| `nifti` | NIfTI | 0.8 for .nii/.nii.gz | Per-subject volumes |
| `tabular` | CSV/TSV | 0.6 for .csv/.tsv | Long-form data |
| `fmristore` | fmristore | 1.0 if detectable | Custom binary format |

### Key Adapter Behaviors

**H5 Adapter:**
- Expects `/gds/assays/<name>` datasets
- Reads space from `/gds/space` group
- Supports packed voxel spaces with mask indices
- Can read map families from `/gds/alignments`
- Can read col_data from `/gds/axes/subjects_table`

**NIfTI Adapter:**
- Reads beta, SE, var from separate NIfTI files per subject
- Requires mask to define sample space
- Constructs space_voxel from NIfTI header
- No native col_data support (must pass separately)

**Tabular Adapter:**
- Reads CSV/TSV in long form: one row per sample × subject × contrast
- Requires column mappings: `effect_cols`, `subject_col`, `sample_col`, `contrast_col`
- Constructs space_sample_labels from unique sample IDs
- Can extract col_data from duplicate subject rows

---

## Reducer Interface Specification

### Reducer Structure

Reducers are registered functions in `.gds_reducers`:

```r
reducer <- list(
  name = "meta:fe",

  # Reducer function signature
  fun = function(beta, var, X, z, p, df, df1, df2, opts) {
    # beta, var, z, p: [subjects × samples] matrices (transposed and sliced)
    # X: [subjects × p] design matrix (optional, for meta-regression)
    # df, df1, df2: degrees of freedom (various shapes)
    # opts: named list of options

    # Return: named list of results
    list(
      beta_g = numeric(n_samples),
      var_g  = numeric(n_samples),
      se_g   = numeric(n_samples),
      z_g    = numeric(n_samples),
      p_g    = numeric(n_samples),
      Q      = numeric(n_samples),
      I2     = numeric(n_samples)
    )
  },

  requires = c("beta", "var"),  # Assays that must be present
  provides = c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2"),  # Assays produced
  options_schema = list()  # Optional schema for opts validation
)
```

### Built-in Reducers

| Name | Description | Requires | Provides |
|------|-------------|----------|----------|
| `meta:fe` | Fixed-effect meta-analysis | `beta`, `var` | `beta_g`, `var_g`, `se_g`, `z_g`, `p_g`, `Q`, `I2` |
| `meta:re` | Random-effects (DL) | `beta`, `var` | `beta_g`, `var_g`, `se_g`, `z_g`, `p_g`, `tau2`, `Q`, `I2` |
| `meta:fe_reg` | Fixed-effect meta-regression | `beta`, `var`, `X` | `coef`, `se_coef`, `Q`, `df_res` |
| `meta:re_reg` | Random-effects meta-regression | `beta`, `var`, `X` | `coef`, `se_coef`, `tau2`, `Q`, `df_res` |
| `combine:stouffer` | Stouffer's Z combination | `z` | `z_g`, `p_g` |
| `combine:fisher` | Fisher's combination | `p` | `p_g`, `chi2`, `df` |
| `combine:lancaster` | Lancaster combination | `p` | `p_g`, `chi2`, `df` |
| `ols:voxelwise` | Voxelwise OLS regression | `beta`, `X` | `coef`, `se_coef`, `t_coef`, `p_coef`, `sigma2`, `df_res` |

### Reducer Execution Details

**Input Preparation:**

1. Arrays are sliced by contrast: `arrays[,,k]` → 2D matrix
2. Matrices are transposed: `[sample × subject]` → `[subject × sample]`
3. Design matrix X is built from formula + col_data if needed
4. Required assays are auto-derived if missing (z from t+df, p from z, etc.)

**Output Handling:**

- **Scalar outputs** (e.g., `Q`, `I2`, `tau2`): Single value per sample, replicated if needed
- **Vector outputs** (e.g., `beta_g`, `var_g`): One value per sample
- **Matrix outputs** (e.g., `coef`, `se_coef`): One value per (parameter × sample)
  - Matrix outputs are unpacked into param-suffixed assays: `coef:X1`, `coef:X2`, `se_coef:X1`, etc.
  - This allows flexible access to individual parameter estimates

**Subject Axis Collapse:**

After reduce, `subjects` changes from `c("s1","s2",...)` to `"meta"`, and the subject dimension becomes 1.

---

## Posthoc Interface Specification

### Posthoc Structure

Posthoc methods are registered in `.gds_posthoc`:

```r
posthoc_method <- list(
  name = "fdr:bh",

  # Posthoc function: (arrays, opts) -> named list of new assays
  fun = function(arrays, opts) {
    # arrays: Named list of 3D arrays [sample × subject × contrast]
    # opts: Named list of options

    # Process p-values (or derive from z/t)
    p <- arrays$p

    # Apply correction
    q <- array(NA_real_, dim = dim(p))
    for (j in seq_len(dim(p)[2])) {
      for (k in seq_len(dim(p)[3])) {
        q[,j,k] <- p.adjust(p[,j,k], method = "BH")
      }
    }

    # Return new assays
    list(q = q)
  },

  requires = c("p"),  # Assays needed
  provides = c("q")   # Assays produced
)
```

### Built-in Posthoc Methods

| Name | Description | Requires | Provides |
|------|-------------|----------|----------|
| `fdr:bh` | Benjamini-Hochberg FDR | `p` | `q` |
| `fdr:by` | Benjamini-Yekutieli FDR | `p` | `q` |

### Extensibility

Users can register custom posthoc methods:

```r
register_posthoc(
  name = "bonferroni",
  fun = function(arrays, opts) {
    p <- arrays$p
    dims <- dim(p)
    n_tests <- dims[1] * dims[3]  # samples × contrasts
    p_adj <- pmin(p * n_tests, 1)
    list(p_bonf = p_adj)
  },
  requires = c("p"),
  provides = c("p_bonf")
)
```

---

## Compute Result Structure

### GDS Object

```r
gds <- list(
  assays = list(
    beta = array(dim = c(N, 1, K)),  # After reduce: subject dim = 1
    var  = array(dim = c(N, 1, K)),
    se   = array(dim = c(N, 1, K)),
    z    = array(dim = c(N, 1, K)),
    p    = array(dim = c(N, 1, K)),
    q    = array(dim = c(N, 1, K)),  # If posthoc applied
    # Meta-analysis specific:
    beta_g = array(dim = c(N, 1, K)),
    var_g  = array(dim = c(N, 1, K)),
    se_g   = array(dim = c(N, 1, K)),
    z_g    = array(dim = c(N, 1, K)),
    p_g    = array(dim = c(N, 1, K)),
    Q      = array(dim = c(N, 1, K)),
    I2     = array(dim = c(N, 1, K)),
    tau2   = array(dim = c(N, 1, K)),  # If random-effects
    # Meta-regression with p parameters:
    `coef:X1`    = array(dim = c(N, 1, K)),
    `se_coef:X1` = array(dim = c(N, 1, K)),
    `t_coef:X1`  = array(dim = c(N, 1, K)),
    `p_coef:X1`  = array(dim = c(N, 1, K)),
    `coef:X2`    = array(dim = c(N, 1, K)),
    ...
  ),

  space = <space object>,  # space_voxel, space_parcels, or space_sample_labels

  subjects = "meta",  # Character vector, "meta" after reduce

  contrasts = c("c1", "c2", ...),  # Character vector

  col_data = data.frame(  # Subject-level covariates (NULL after reduce)
    # Keyed by subject, rownames = subjects
  ),

  row_data = data.frame(  # Sample-level metadata
    # One row per sample
  ),

  metadata = list(
    schema_version = "0.1.0",

    units = list(
      beta = "beta coefficient",
      var = "variance",
      ...
    ),

    provenance = list(
      graph = list(  # Provenance nodes
        list(op = "reduce", params = list(...), timestamp = ..., hash = ...),
        list(op = "posthoc", params = list(...), timestamp = ..., hash = ...)
      ),
      log = c(  # Human-readable log
        "[2025-10-29 12:00:00] reduce(method=meta:fe, weights=1/var)",
        "[2025-10-29 12:00:01] posthoc(method=fdr:bh)"
      ),
      digest = "xxhash64-digest-of-plan"
    ),

    software = list(
      package = "gdsfmri",
      version = "0.1.0",
      R_version = "4.3.0"
    ),

    alignment = NULL,  # Alignment metadata if align_to_group was used

    map_families = list(),  # Registered spatial mappings

    mask_info = NULL,  # Mask metadata if masking was applied

    contrast_info = NULL,  # Contrast metadata

    design_mats = list(  # Design matrices from meta-regression
      list(
        method = "meta:fe_reg",
        formula = "~ age + sex",
        columns = c("(Intercept)", "age", "sex"),
        hash = "..."
      )
    ),

    attachments = list(  # Large auxiliary data
      "reduce/ols:voxelwise/contrast=c1" = list(
        type = "cov_tri",
        terms = c("X1", "X2", "X3"),
        pack = "upper",
        cov_tri = matrix(...)  # Packed upper triangle of covariance
      )
    ),

    notes = NULL,  # User notes

    created = <POSIXct timestamp>
  )
)

class(gds) <- c("gds", "group_data")
```

### Assay Naming Conventions

| Pattern | Example | Meaning |
|---------|---------|---------|
| Plain | `beta`, `var`, `se` | Standard effect size arrays |
| `_g` suffix | `beta_g`, `z_g`, `p_g` | Group-level (meta-analysis) results |
| `coef:` prefix | `coef:age`, `coef:(Intercept)` | Regression coefficient for named parameter |
| `se_coef:` prefix | `se_coef:age` | Standard error of coefficient |
| `t_coef:` prefix | `t_coef:age` | t-statistic for coefficient |
| `p_coef:` prefix | `p_coef:age` | p-value for coefficient |
| Other | `Q`, `I2`, `tau2`, `chi2`, `df`, `df_res` | Meta-analysis statistics |

---

## Error Handling & Diagnostics

### Common Error Messages

#### 1. "No adapter detected for source"

**Cause:** No registered adapter's `detect()` function returned a score > 0.

**Location:** `R/adapter-registry.R::detect_adapter()`

**Debug:**
```r
# Check registered adapters
ls(fmrigds:::.adapter_registry)

# Manually test detection
source <- "/path/to/data.csv"
adapters <- ls(fmrigds:::.adapter_registry)
for (name in adapters) {
  adapter <- fmrigds:::get_adapter(name)
  score <- adapter$detect(source)
  cat(name, ": ", score, "\n")
}
```

#### 2. "Evidence reducer requires z or p"

**Cause:** Reducer requires specific assays that are not present and cannot be auto-derived.

**Location:** `R/reduce-exec.R::.reduce_evidence()`

**Explanation:**
- Stouffer/Fisher combiners need z or p values
- Auto-derivation attempted via `.ensure_required_arrays()`:
  - `z` derived from `{t, df}` or `{p, beta}`
  - `p` derived from `{t, df}` or `{z}`
- If neither path succeeds, error is raised

**Debug:**
```r
# Check available assays in probe result
plan$source$probe$assays

# Check what reducer requires
reducer <- fmrigds:::get_reducer("combine:stouffer")
reducer$requires  # c("z")

# Manually derive if possible
arrays <- adapter$read(handle, assays = c("beta", "se", "t", "df"), ...)
arrays$z <- fmrigds:::derive_z(arrays)
```

#### 3. "assay() returns NULL"

**Cause:** Requested assay name is not in `gds$assays`.

**Location:** `R/gds-class.R::assay.gds()`

**Debug:**
```r
# List all available assays
names(gds$assays)

# Check what reducer provides
reducer <- fmrigds:::get_reducer("meta:fe")
reducer$provides
# c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2")

# Assay naming after reduce
# Legacy fields are mapped: beta <- beta_g, var <- var_g, etc.
```

**Common Issue:** Requesting `assay(gds, "beta")` when reducer produced `beta_g`. Solution: Check `reducer$provides` or request `beta_g` directly (though legacy mapping should handle this).

#### 4. "meta:fe_reg requires X (subjects x p) in options$X"

**Cause:** Meta-regression reducer requires design matrix X, but none was provided.

**Location:** `R/reducers-core.R::register_core_reducers()`

**Solution:**
```r
# Option 1: Use formula (auto-builds X)
plan <- reduce(plan, method = "meta:fe_reg", formula = ~ age + sex)
# X is built from formula + col_data during reduce execution

# Option 2: Pass X explicitly
X <- model.matrix(~ age + sex, data = col_data)
plan <- reduce(plan, method = "meta:fe_reg", options = list(X = X))
```

#### 5. "Unknown post-hoc method: fdr:custom"

**Cause:** Requested posthoc method is not registered.

**Location:** `R/posthoc-exec.R::apply_posthoc()`

**Debug:**
```r
# List registered posthoc methods
fmrigds:::list_posthoc()

# Register custom method
register_posthoc(
  name = "fdr:custom",
  fun = function(arrays, opts) { ... },
  requires = c("p"),
  provides = c("q")
)
```

---

## Common Failure Modes

### 1. Missing Assays After Compute

**Symptom:** `assay(gds, "tau2")` returns `NULL` after fixed-effect meta-analysis.

**Cause:** The `meta:fe` reducer does not produce `tau2` (only `meta:re` does).

**Solution:** Check `reducer$provides` to see what assays are produced:
```r
fmrigds:::get_reducer("meta:fe")$provides
# c("beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2")

fmrigds:::get_reducer("meta:re")$provides
# c("beta_g", "var_g", "se_g", "z_g", "p_g", "tau2", "Q", "I2")
```

### 2. Dimension Mismatch After Reduce

**Symptom:** `dim(assay(gds, "beta"))` is `c(1000, 1, 5)` but expected `c(1000, 20, 5)`.

**Cause:** After `reduce()`, the subject dimension collapses from S subjects to 1 ("meta").

**Expected Behavior:** This is intentional. Reduce operations collapse across subjects.

**Solution:** To access subject-level data, compute the plan *before* reduce:
```r
# Subject-level
subject_gds <- gds(source) |> compute()
dim(assay(subject_gds, "beta"))  # c(1000, 20, 5)

# Group-level
group_gds <- gds(source) |> reduce("meta:fe") |> compute()
dim(assay(group_gds, "beta"))  # c(1000, 1, 5)
```

### 3. Regression Coefficients Not Accessible

**Symptom:** After `reduce(method = "meta:fe_reg", formula = ~ age + sex)`, cannot find coefficient assays.

**Cause:** Coefficient assays use param-suffixed naming: `coef:age`, `se_coef:age`, etc.

**Solution:**
```r
# List all assays
names(gds$assays)
# Look for: "coef:(Intercept)", "coef:age", "coef:sex",
#           "se_coef:(Intercept)", "se_coef:age", "se_coef:sex", etc.

# Extract specific coefficient
beta_age <- assay(gds, "coef:age")
se_age   <- assay(gds, "se_coef:age")
t_age    <- assay(gds, "t_coef:age")
p_age    <- assay(gds, "p_coef:age")
```

### 4. FDR Correction Not Applied

**Symptom:** `assay(gds, "q")` returns `NULL` even though `posthoc("fdr:bh")` was called.

**Cause:** `posthoc()` only adds a node to the plan; must call `compute()` to execute.

**Solution:**
```r
# Incorrect (posthoc not computed)
plan <- gds(source) |> reduce("meta:fe") |> posthoc("fdr:bh")
gds <- compute(plan)  # Wait, this is correct...

# Debug: Check if posthoc was actually executed
gds$metadata$provenance$log
# Should contain: "[timestamp] posthoc(method=fdr:bh)"

# Check if p values exist (required for FDR)
"p" %in% names(gds$assays)  # Must be TRUE

# Check registered posthoc methods
fmrigds:::list_posthoc()
# Must contain "fdr:bh"
```

### 5. Col_data Mismatch After Reduce

**Symptom:** `col_data(gds)` is `NULL` after reduce, even though it was present before.

**Cause:** After reduce, subjects collapse to "meta", making subject-level covariates meaningless.

**Expected Behavior:** `col_data` is automatically dropped when `subjects == "meta"`.

**Location:** `R/compute.R::compute()` lines 86-96

**Workaround:** Store col_data separately before reduce:
```r
pre_reduce <- gds(source) |> compute()
col_data_saved <- col_data(pre_reduce)

post_reduce <- gds(source) |> reduce("meta:fe") |> compute()
# col_data(post_reduce) is NULL, but col_data_saved still has it
```

---

## Debugging Recommendations

### Enabling Verbose Output

Currently, fmrigds does not have a built-in verbose mode, but you can add tracing:

```r
# Trace adapter operations
trace(fmrigds:::get_adapter("h5")$probe, tracer = function() {
  cat("Probing H5 file...\n")
})

# Trace reducer calls
trace(fmrigds:::apply_reduce, tracer = function() {
  cat("Applying reduce with method:", node$method, "\n")
})

# Untrace when done
untrace(fmrigds:::get_adapter("h5")$probe)
untrace(fmrigds:::apply_reduce)
```

### Inspecting Plans Before Compute

```r
plan <- gds(source) |> reduce("meta:fe") |> posthoc("fdr:bh")

# Examine plan structure
str(plan, max.level = 2)

# List operation nodes
lapply(plan$nodes, function(node) {
  c(op = node$op, method = node$method %||% NA)
})

# Check source metadata
plan$source$probe$assays
plan$source$probe$dims
plan$source$probe$subjects
```

### Debugging Reducer Execution

```r
# List all registered reducers
fmrigds:::list_reducers()

# Inspect a specific reducer
reducer <- fmrigds:::get_reducer("meta:fe")
str(reducer)

# Check what it requires and provides
reducer$requires
reducer$provides

# Manually test reducer kernel
beta <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 subjects × 10 samples
var  <- matrix(runif(100, 0.01, 0.1), nrow = 10, ncol = 10)

result <- fmrigds:::core_meta_fe_kernel(beta, var, X = NULL, df = NULL, opts = list())
names(result)  # Should match reducer$provides
```

### Debugging Assay Derivation

```r
# Check derivation logic
arrays <- list(
  beta = array(rnorm(100), dim = c(10, 5, 2)),
  se   = array(runif(100, 0.1, 0.5), dim = c(10, 5, 2))
)

# Derive var from se
arrays$var <- fmrigds:::derive_var(arrays)

# Derive t from beta and var
arrays$t <- fmrigds:::derive_t(arrays)

# Derive p from t (but need df!)
arrays$df <- array(20, dim = c(10, 5, 2))
arrays$p <- fmrigds:::derive_p(arrays)
```

### Logging Assay Names During Compute

Add temporary logging to understand assay flow:

```r
# Before compute
plan <- gds(source) |> reduce("meta:fe")

# Add tracer
trace(fmrigds:::apply_reduce, exit = function() {
  cat("Assays after reduce:", paste(names(arrays), collapse = ", "), "\n")
})

# Compute
gds <- compute(plan)

# Check final assays
cat("Final assays:", paste(names(gds$assays), collapse = ", "), "\n")

# Cleanup
untrace(fmrigds:::apply_reduce)
```

### Interactive Plan Execution

For deep debugging, manually step through compute:

```r
plan <- gds(source) |> reduce("meta:fe") |> posthoc("fdr:bh")

# Open adapter
adapter <- fmrigds:::get_adapter(plan$source$adapter)
handle <- adapter$open(plan$source$source)

# Read arrays
arrays <- adapter$read(handle,
                       assays = plan$source$probe$assays,
                       block = NULL)
adapter$close(handle)

# Apply nodes manually
for (i in seq_along(plan$nodes)) {
  node <- plan$nodes[[i]]
  cat("Node", i, ":", node$op, "\n")

  if (node$op == "reduce") {
    result <- fmrigds:::apply_reduce(node, arrays, node$weights,
                                     plan$source$probe$subjects,
                                     plan$meta$col_data)
    arrays <- result$arrays
    cat("  Assays after reduce:", paste(names(arrays), collapse = ", "), "\n")
  } else if (node$op == "posthoc") {
    result <- fmrigds:::apply_posthoc(node, arrays)
    arrays <- result$arrays
    cat("  Assays after posthoc:", paste(names(arrays), collapse = ", "), "\n")
  }
}

# Inspect final arrays
str(arrays, max.level = 1)
```

### Recommended Debugging Workflow

1. **Start Simple**: Test with minimal CSV data (2 samples, 2 subjects, 1 contrast)
2. **Check Adapter**: Verify `probe()` returns expected metadata
3. **Inspect Plan**: Examine `plan$nodes` before compute
4. **Trace Reducers**: Use `trace()` to log reducer inputs/outputs
5. **Check Provides**: Always verify `reducer$provides` matches expected assays
6. **Manual Stepping**: For complex issues, step through compute manually
7. **Validate Dimensions**: After each operation, check array dimensions
8. **Log Assay Names**: Track which assays exist at each pipeline stage

---

## Recommendations for Improvements

### 1. Add Verbose Mode

**Proposal:** Add `verbose = TRUE` parameter to `compute()` that logs:
- Which adapter is used
- Dimensions of arrays read from adapter
- Each operation node being applied
- Assays produced by each operation
- Final assay inventory

**Implementation:**
```r
compute <- function(x, verbose = FALSE, ...) {
  if (verbose) cat("Opening adapter:", plan$source$adapter, "\n")
  # ...
  if (verbose) cat("Read assays:", paste(names(arrays), collapse=", "), "\n")
  # ...
  for (node in plan$nodes) {
    if (verbose) cat("Applying:", node$op, "\n")
    # ...
    if (verbose && node$op == "reduce") {
      cat("  Produced:", paste(names(result$arrays), collapse=", "), "\n")
    }
  }
}
```

### 2. Assay Discovery API

**Proposal:** Add function to predict what assays will be available after compute, without actually computing:

```r
predict_assays <- function(plan) {
  # Start with source assays
  assays <- plan$source$probe$assays

  # Track through operations
  for (node in plan$nodes) {
    if (node$op == "derive") {
      assays <- c(assays, node$what)
    } else if (node$op == "reduce") {
      reducer <- get_reducer(.normalize_reducer_name(node$method))
      if (!is.null(reducer)) {
        assays <- c(assays, reducer$provides)
      }
    } else if (node$op == "posthoc") {
      ph <- get_posthoc(node$method)
      if (!is.null(ph)) {
        assays <- c(assays, ph$provides)
      }
    }
  }

  unique(assays)
}
```

**Usage:**
```r
plan <- gds(source) |> reduce("meta:fe") |> posthoc("fdr:bh")
predict_assays(plan)
# c("beta", "var", "beta_g", "var_g", "se_g", "z_g", "p_g", "Q", "I2", "q")
```

### 3. Plan Validation

**Proposal:** Add `validate_plan()` to check for common errors before compute:

```r
validate_plan <- function(plan, strict = FALSE) {
  errors <- character()
  warnings <- character()

  # Check reducer requirements
  for (node in plan$nodes) {
    if (node$op == "reduce") {
      reducer <- get_reducer(.normalize_reducer_name(node$method))
      if (is.null(reducer)) {
        errors <- c(errors, paste("Unknown reducer:", node$method))
      } else {
        available <- plan$source$probe$assays
        missing <- setdiff(reducer$requires, available)
        if (length(missing) > 0) {
          # Check if can be auto-derived
          for (assay in missing) {
            if (assay == "z" && all(c("t","df") %in% available)) next
            if (assay == "p" && "z" %in% available) next
            errors <- c(errors, paste("Reducer", node$method, "requires", assay, "but not available"))
          }
        }
      }
    }
  }

  if (length(errors) > 0) {
    stop(paste("Plan validation failed:\n", paste("  -", errors, collapse="\n")), call. = FALSE)
  }

  if (length(warnings) > 0 && strict) {
    warning(paste("Plan validation warnings:\n", paste("  -", warnings, collapse="\n")), call. = FALSE)
  }

  invisible(TRUE)
}
```

### 4. Assay Documentation

**Proposal:** Add `?assays` help page documenting standard assay naming conventions:

- Effect sizes: `beta`, `se`, `var`, `t`, `z`, `p`
- Group-level: `beta_g`, `var_g`, `se_g`, `z_g`, `p_g`
- Meta-analysis: `Q`, `I2`, `tau2`, `df_res`
- Regression: `coef:<param>`, `se_coef:<param>`, `t_coef:<param>`, `p_coef:<param>`
- Post-hoc: `q` (FDR-adjusted p-values)
- Evidence: `chi2`, `df`, `Z`

### 5. Helpful Error Messages for Missing Assays

**Current:** `assay(gds, "tau2")` returns `NULL` silently.

**Proposed:** Check if assay is a common typo or expected from different reducer:

```r
assay.gds <- function(x, name = "beta", ...) {
  if (!name %in% names(x$assays)) {
    # Check if it's a known assay from different reducer
    suggestions <- character()
    if (name == "tau2") {
      suggestions <- "tau2 is only produced by random-effects reducers (meta:re, meta:re_reg)"
    } else if (grepl("^coef:", name)) {
      suggestions <- "Regression coefficients require meta-regression reducers (meta:fe_reg, meta:re_reg)"
    }

    msg <- paste0("Assay '", name, "' not found in GDS object.\n",
                  "Available assays: ", paste(names(x$assays), collapse = ", "))
    if (length(suggestions) > 0) {
      msg <- paste0(msg, "\n\nHint: ", suggestions)
    }

    warning(msg, call. = FALSE)
    return(NULL)
  }
  x$assays[[name]]
}
```

### 6. Plan Serialization

**Current:** Plans can be serialized but no helper functions exist.

**Proposal:** Add `save_plan()` and `load_plan()`:

```r
save_plan <- function(plan, path) {
  saveRDS(plan, path)
}

load_plan <- function(path) {
  plan <- readRDS(path)
  if (!inherits(plan, "gds_plan")) {
    stop("Loaded object is not a gds_plan", call. = FALSE)
  }
  plan
}
```

This would allow saving expensive probe results and reusing plans across sessions.

---

## Conclusion

The fmrigds pipeline implements a sophisticated lazy evaluation system for group-level fMRI analysis. Key takeaways:

1. **Plans are immutable**: Each verb returns a new plan
2. **Execution is delayed**: Nothing happens until `compute()` is called
3. **Assays are dynamic**: The set of available assays depends on which reducers were applied
4. **Reducers are modular**: New meta-analysis methods can be registered
5. **Errors are informative**: Most failures happen at compute time with clear messages

The most common debugging challenge is understanding what assays will be available after compute. The `provides` field in reducer registrations is the source of truth, but it's not exposed in plan metadata. The recommended improvements would make this more transparent.

For users transitioning from legacy `fmri_meta()`, the key mental model shift is:

```r
# Legacy: All-in-one execution
result <- fmri_meta(group_data, formula = ~ 1, method = "fe")

# fmrigds: Lazy plan + compute
result <- group_data |>           # Build plan
  fmrigds::as_plan() |>           #   (implicit via reduce)
  fmrigds::reduce("meta:fe") |>   # Add reduce op
  fmrigds::compute()              # Execute plan → GDS

# Extract via assay()
beta <- assay(result, "beta")
```

The pipeline's flexibility comes from separating plan construction (what to do) from execution (doing it), enabling optimization, serialization, and extensibility.
