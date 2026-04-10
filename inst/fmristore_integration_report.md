# fmristore Integration Report for gdsfmri

**Date:** 2025-10-28
**Sprint:** Sprint 7 - Storage Backend Integration
**Status:** Design Complete, Implementation Ready

---

## Executive Summary

This report provides a comprehensive assessment of how the **fmristore** package can be integrated with **gdsfmri** as a storage adapter. Based on systematic examination of fmristore's codebase, we have designed a complete adapter that supports:

1. **Four fmristore layout types**:
   - **LabeledVolumeSet** (primary labeled volume class with lazy loading)
   - **H5NeuroVol/H5NeuroVec** (alternative voxel classes)
   - **H5ParcellatedScan/MultiScan** (ROI/cluster data)
   - **LatentNeuroVec** (basis/ICA representations)
2. **Two integration paths**:
   - **Path A**: Read existing fmristore files without schema changes (immediate compatibility)
   - **Path B**: Prefer `/gds` group when present for perfect round-tripping
3. **Multi-file strategies**: Subject-wise and contrast-wise file organizations
4. **Full GDS semantics**: Space objects, assay mapping, lazy pipelines, block streaming

**Recommendation:** Proceed with implementation following the phased approach outlined in Section 8.

---

## 1. fmristore Data Structures → GDS Space Mapping

### 1.1 LabeledVolumeSet (Labeled Volume - Primary Class)

**S4 Class Definition:**
```r
setClass("LabeledVolumeSet",
  slots = c(
    obj      = "H5File",           # HDF5 file handle
    mask     = "LogicalNeuroVol",  # 3D mask [X,Y,Z]
    labels   = "character",        # Volume labels
    load_env = "environment"       # Lazy loading environment
  ),
  contains = c("NeuroVec")         # Extends NeuroVec (4D)
)
```

**HDF5 Layout:**
```
/header/*          # NIfTI-like metadata (dim, pixdim, qform, sform)
  /dim             # [4] or [3] dimensions
  /pixdim          # Voxel spacing
  /qform_code      # Coordinate system code
  /quatern_b/c/d   # Quaternion parameters
  /qoffset_x/y/z   # Offset parameters
  /srow_x/y/z      # Affine rows (alternative to qform)
/mask              # 3D uint8 [X,Y,Z] binary mask (1=active, 0=excluded)
/labels            # Character vector of contrast names
/data/
  /<label>         # 1D float64 vector, length = sum(mask)
  /<label>_var     # Optional variance
  /<label>_stderr  # Optional standard error
```

**Key Characteristics:**
- **Logical 4D**: `[X, Y, Z, #labels]` via NeuroVec inheritance
- **Lazy Loading**: Data loaded on-demand via `[[label]]` or `[["label"]]` indexing
- **Storage**: Packed (only masked voxels stored)
- **Mask**: Explicit 3D LogicalNeuroVol
- **I/O**: `write_labeled_vec()` and `read_labeled_vec()` functions
- **Indexing**: Supports `[i,j,k,l]` (4D array) and `[[i]]` / `[["label"]]` (by label)

**GDS Space Constructor:**
```r
space_voxel(
  dim = c(91, 109, 91),                    # From /header/dim
  affine = construct_affine_from_header(), # 4×4 matrix
  mask_bitmap = mask > 0,                  # Logical [X,Y,Z]
  mask_idx = which(as.vector(mask)),       # Linear indices
  storage = "packed"
)
```

**Dimensional Mapping:**
- **sample**: `sum(mask)` = number of active voxels
- **subject**: Number of files (multi-file) or 1 (single-file)
- **contrast**: `length(labels)`

**Assay Mapping:**
```r
/data/<label>        → assay("beta")[, subject, contrast]
/data/<label>_var    → assay("var")[, subject, contrast]
/data/<label>_stderr → assay("se")[, subject, contrast]
```

**Detection Logic:**
```r
# Use fmristore::detect_h5_type() which returns:
# "labeled_volume" for LabeledVolumeSet layout
detect_h5_type <- function(file_path) {
  h5 <- H5File$new(file_path, mode = "r")
  on.exit(h5$close_all())

  has_header <- h5$exists("/header")
  has_mask <- h5$exists("/mask")
  has_labels <- h5$exists("/labels")
  has_data <- h5$exists("/data")

  if (has_header && has_mask && has_labels && has_data) {
    return("labeled_volume")
  }
  return(NA_character_)
}
```

---

### 1.2 H5NeuroVol / H5NeuroVec (Alternative Voxel Classes)

**Note:** These classes use similar or identical HDF5 layouts to LabeledVolumeSet but differ in their API and memory management:
- **LabeledVolumeSet**: S4 class extending NeuroVec, lazy loading via environment, primary class for labeled volumes
- **H5NeuroVol**: S4 class extending NeuroVol (3D), single volume access
- **H5NeuroVec**: S4 class extending NeuroVec (4D), eager dataset access

All three use the same `detect_h5_type()` result of `"labeled_volume"`. The adapter can treat them uniformly.

**HDF5 Layout:**
```
/header/*          # NIfTI-like metadata (dim, pixdim, qform, sform)
  /dim             # [4] or [3] dimensions
  /pixdim          # Voxel spacing
  /qform_code      # Coordinate system code
  /quatern_b/c/d   # Quaternion parameters
  /qoffset_x/y/z   # Offset parameters
  /srow_x/y/z      # Affine rows (alternative to qform)
/mask              # 3D uint8 [X,Y,Z] binary mask (1=active, 0=excluded)
/labels            # Character vector of contrast names
/data/
  /<label>         # 1D float64 vector, length = sum(mask)
  /<label>_var     # Optional variance
  /<label>_stderr  # Optional standard error
```

**Key Characteristics:**
- **Storage**: Packed (only masked voxels stored)
- **Mask**: Explicit 3D binary array
- **Data organization**: One 1D vector per contrast/label
- **Affine**: Derived from qform or sform in header

**GDS Space Constructor:**
```r
space_voxel(
  dim = c(91, 109, 91),                    # From /header/dim
  affine = construct_affine_from_header(), # 4×4 matrix
  mask_bitmap = mask > 0,                  # Logical [X,Y,Z]
  mask_idx = which(as.vector(mask)),       # Linear indices
  storage = "packed"
)
```

**Dimensional Mapping:**
- **sample**: `sum(mask)` = number of active voxels
- **subject**: Number of files (multi-file) or 1 (single-file)
- **contrast**: `length(labels)`

**Assay Mapping:**
```r
/data/<label>        → assay("beta")[, subject, contrast]
/data/<label>_var    → assay("var")[, subject, contrast]
/data/<label>_stderr → assay("se")[, subject, contrast]
```

---

### 1.2 H5ParcellatedScan / H5ParcellatedMultiScan (Parcel/ROI)

**HDF5 Layout:**
```
/cluster_map       # 3D integer [X,Y,Z] with cluster IDs (0 = background)
/cluster_metadata  # Optional data.frame with parcel properties
/data/
  /cluster_<id>    # [timepoints] or [subjects, timepoints]
  /cluster_<id>_var
```

**Alternative Layout (scan-based):**
```
/cluster_map
/scans/
  /<scan_id>       # [clusters, timepoints] matrix
```

**Key Characteristics:**
- **Cluster map**: Integer IDs defining voxel→parcel membership
- **Data matrix**: Cluster × Time (needs temporal aggregation → contrasts)
- **Sample axis**: Unique cluster IDs (excluding 0 = background)
- **Metadata**: Parcel names, anatomical labels, coordinates

**GDS Space Constructor:**
```r
cluster_ids <- sort(unique(cluster_map[cluster_map > 0]))
labels <- as.character(cluster_ids)

# Build membership list (voxel indices for each parcel)
membership <- lapply(cluster_ids, function(cid) {
  which(as.vector(cluster_map) == cid)
})
names(membership) <- labels

space_parcels(
  labels = labels,                 # c("1", "2", ..., "N")
  lookup = cluster_metadata,       # Optional data.frame
  membership = membership          # List of voxel index vectors
)
```

**Dimensional Mapping:**
- **sample**: K parcels/clusters
- **subject**: Number of scans
- **contrast**: 1 (or extracted from temporal structure)

**LinearMap Construction:**
```r
# Voxel → Parcel averaging map
voxel_to_parcel <- map_linear(
  from_space = voxel_space,
  to_space = parcel_space,
  operator = build_parcel_operator(membership, mask_idx),
  uncertainty = UncertaintyRule("independent")
)

# Operator: sparse [K_parcels × V_voxels] averaging matrix
build_parcel_operator <- function(membership, mask_idx) {
  # Each row = mean over voxels in that parcel
  # Entry [p,v] = 1/n_p if voxel v in parcel p, else 0
}
```

---

### 1.3 LatentNeuroVec (Latent Basis/ICA)

**HDF5 Layout:**
```
/basis/
  /basis_matrix    # [k, V] components × voxels
  /basis_method    # "ICA", "PCA", "NMF", etc.
/scans/
  /<scan_id>/
    /embedding     # [k] or [k, n_contrasts] component loadings
    /metadata      # Scan-specific info
/mask              # Optional voxel mask for basis
/header/*          # Optional reference voxel space
```

**Key Characteristics:**
- **Basis matrix**: Spatial components (k rows = components, V columns = voxels)
- **Embeddings**: One k-vector (or k×C matrix) per subject/scan
- **Projector**: Basis matrix enables reconstruction to voxels
- **Reference space**: Optional voxel grid for visualization

**GDS Space Constructor:**
```r
basis_matrix <- read_basis_matrix()  # [k, V]
k <- nrow(basis_matrix)

space_basis(
  k = k,
  basis_name = "ICA",              # Or "PCA", "latent"
  projector = basis_matrix,        # [k, V] for voxel reconstruction
  voxel_space = space_voxel(...)   # Optional reference
)
```

**Dimensional Mapping:**
- **sample**: k components
- **subject**: Number of scans with embeddings
- **contrast**: Embedding dimensions (usually 1, could be multiple)

**Assay Mapping:**
```r
/scans/<id>/embedding → assay("beta")[, subject, contrast]
# Stack embeddings across subjects: [k, n_subjects, n_contrasts]
```

**LinearMap Construction:**
```r
# Basis → Voxel: project component loadings to voxel space
basis_to_voxel <- map_linear(
  from_space = basis_space,
  to_space = voxel_space,
  operator = t(basis_matrix),  # [V, k] transpose
  uncertainty = UncertaintyRule("independent")
)

# Voxel → Basis: inverse/pseudo-inverse for decomposition
# Can use Moore-Penrose or learned inverse
```

---

## 2. Adapter Design: Detection and Routing

### 2.1 Detection Priority

**Scoring Logic:**
```r
.fmristore_detect <- function(source) {
  # Priority 1: /gds group (Path B, perfect interop)
  if (h5$exists("/gds")) return(1.0)

  # Priority 2: fmristore layout (Path A, legacy)
  type <- fmristore::detect_h5_type(source)
  if (type %in% c("labeled_volume", "parcellated", "latent")) return(0.95)

  # Priority 3: Multi-file fmristore
  if (all_files_are_fmristore(source)) return(0.90)

  # Not fmristore
  return(FALSE)
}
```

**Why This Priority?**
- **1.0 (/gds)**: Native GDS format ensures perfect round-tripping, full provenance
- **0.95 (fmristore)**: Legacy support, no schema changes needed
- **0.90 (multi-file)**: More complex reading logic, subject axis composition

### 2.2 Path Routing

```r
.fmristore_probe <- function(handle, ...) {
  if (handle$has_gds) {
    # Path B: Read /gds/axes, /gds/space, /gds/assays
    .probe_gds_layout(handle, ...)
  } else {
    # Path A: Use fmristore-specific probe
    switch(handle$type,
      labeled_volume = .probe_labeled_volume(handle, ...),
      parcellated    = .probe_parcellated(handle, ...),
      latent         = .probe_latent(handle, ...)
    )
  }
}
```

**Path B Advantages:**
- Explicit axes (`/gds/axes/subjects`, `/gds/axes/contrasts`)
- Multi-assay support (`/gds/assays/beta`, `/gds/assays/var`, ...)
- Provenance tracking (`/gds/provenance`)
- Alignment persistence (`/gds/alignments/mni_warp/...`)

**Path A Advantages:**
- Zero schema changes to existing files
- Works with all legacy fmristore data
- Relies on fmristore's existing validation

---

## 3. Probe Results by Layout Type

### 3.1 Path B: /gds Layout

**Example Probe Return (Voxel Space):**
```r
list(
  assays = c("beta", "var", "se", "t", "z", "df", "p"),
  dims = c(sample = 15342, subject = 20, contrast = 5),
  subjects = c("sub-01", "sub-02", ..., "sub-20"),
  contrasts = c("rest", "task1", "task2", "task3", "task4"),
  space = space_voxel(
    dim = c(64, 64, 30),
    affine = matrix(..., 4, 4),
    mask_idx = 1:15342,
    storage = "packed"
  ),
  maps = list(
    mni_warp = AlignmentFamily(...),  # If stored in /gds/alignments
    voxel_to_parcel = LinearMap(...)  # If multi-space
  ),
  metadata = list(
    schema_version = "gds-h5/1.0",
    source_file = "group_analysis.h5",
    created = "2025-10-28T12:00:00Z"
  )
)
```

### 3.2 Path A: Labeled Volume (Single File)

**Example Probe Return:**
```r
list(
  assays = c("beta", "var"),
  dims = c(sample = 15342, subject = 1, contrast = 3),
  subjects = c("group_mean"),  # From filename or /subjects
  contrasts = c("cope1", "cope2", "cope3"),  # From /labels
  space = space_voxel(
    dim = c(64, 64, 30),
    affine = construct_affine_from_qform(header),
    mask_bitmap = mask > 0,
    mask_idx = which(as.vector(mask)),
    storage = "packed"
  ),
  maps = list(),
  metadata = list(
    schema_version = "0.1.0",
    source_file = "group_mean.h5",
    layout_type = "labeled_volume"
  )
)
```

### 3.3 Path A: Labeled Volume (Multi-File)

**Files:**
```
sub-01_cope1.h5
sub-02_cope1.h5
...
sub-20_cope1.h5
```

**Example Probe Return:**
```r
list(
  assays = c("beta", "var"),
  dims = c(sample = 15342, subject = 20, contrast = 3),
  subjects = c("sub-01", "sub-02", ..., "sub-20"),  # Extracted from filenames
  contrasts = c("cope1", "cope2", "cope3"),          # From first file's /labels
  space = space_voxel(
    dim = c(64, 64, 30),
    affine = ...,  # From first file (template)
    mask_bitmap = harmonized_mask,  # Union or intersection across files
    mask_idx = which(as.vector(harmonized_mask)),
    storage = "packed"
  ),
  maps = list(),
  metadata = list(
    schema_version = "0.1.0",
    source_files = c("sub-01_cope1.h5", ..., "sub-20_cope1.h5"),
    layout_type = "labeled_volume",
    mask_policy = "intersection"  # Or "union"
  )
)
```

### 3.4 Path A: Parcellated

**Example Probe Return:**
```r
list(
  assays = c("beta"),
  dims = c(sample = 200, subject = 15, contrast = 1),
  subjects = c("sub-01", ..., "sub-15"),
  contrasts = c("rest"),  # Or temporal aggregation labels
  space = space_parcels(
    labels = as.character(1:200),
    lookup = cluster_metadata,  # data.frame with parcel info
    membership = list(
      "1" = c(voxel_indices_for_cluster_1),
      "2" = c(voxel_indices_for_cluster_2),
      ...
    )
  ),
  maps = list(
    voxel_to_parcel = LinearMap(...)  # Built from cluster_map
  ),
  metadata = list(
    schema_version = "0.1.0",
    source_file = "parcellated_data.h5",
    layout_type = "parcellated",
    cluster_ids = 1:200
  )
)
```

### 3.5 Path A: Latent

**Example Probe Return:**
```r
list(
  assays = c("beta"),  # Component loadings/embeddings
  dims = c(sample = 50, subject = 12, contrast = 3),
  subjects = c("sub-01", ..., "sub-12"),
  contrasts = c("rest", "task1", "task2"),
  space = space_basis(
    k = 50,
    basis_name = "ICA",
    projector = basis_matrix,      # [50, 228483] for voxel reconstruction
    voxel_space = space_voxel(...) # Reference space
  ),
  maps = list(
    basis_to_voxel = LinearMap(...)  # t(basis_matrix)
  ),
  metadata = list(
    schema_version = "0.1.0",
    source_file = "latent_embeddings.h5",
    layout_type = "latent",
    basis_dims = c(50, 228483),
    basis_method = "ICA"
  )
)
```

---

## 4. Block Reading Strategy

### 4.1 Path B: Direct HDF5 Slicing

```r
.read_gds_layout <- function(handle, assays, block = NULL) {
  # Determine indices
  sample_idx <- block$sample %||% seq_len(dims[1])
  subject_idx <- block$subject %||% TRUE  # hdf5r: all
  contrast_idx <- block$contrast %||% TRUE

  # Read each assay with direct slicing
  result <- lapply(assays, function(name) {
    path <- paste0("/gds/assays/", name)
    h5[[path]][sample_idx, subject_idx, contrast_idx]
  })

  names(result) <- assays
  result
}
```

**Performance:**
- Chunked datasets in HDF5 enable efficient partial reads
- Block size typically 10,000–100,000 voxels
- Subject/contrast dimensions usually read in full (small)

### 4.2 Path A: Labeled Volume (Multi-File)

```r
.read_labeled_volume_multifile <- function(handle, assays, block, mask_idx) {
  files <- handle$files
  n_subjects <- length(files)

  # Allocate output: [n_samples_block, n_subjects, n_contrasts]
  sample_idx <- block$sample %||% seq_along(mask_idx)
  n_samples <- length(sample_idx)

  result <- lapply(assays, function(nm) {
    array(NA_real_, dim = c(n_samples, n_subjects, n_contrasts))
  })

  # Read each subject file
  for (j in seq_along(files)) {
    h5 <- hdf5r::H5File$new(files[j], mode = "r")
    on.exit(h5$close(), add = TRUE)

    for (k in seq_along(labels)) {
      label <- labels[k]
      path <- paste0("/data/", label)

      # Read full vector, subset by mask_idx, then by sample_idx
      vec_full <- h5[[path]][]
      vec_masked <- vec_full[mask_idx]
      vec_block <- vec_masked[sample_idx]

      result[[assay]][, j, k] <- vec_block
    }
  }

  result
}
```

**Optimization:**
- Open file once per subject
- Read all contrasts for that subject before moving to next
- Apply mask and block subset in memory (cheap after I/O)

### 4.3 Path A: Parcellated

```r
.read_parcellated <- function(handle, assays, block) {
  cluster_ids <- get_cluster_ids(handle)
  sample_idx <- block$sample %||% seq_along(cluster_ids)

  result <- lapply(assays, function(nm) {
    arr <- array(NA_real_, dim = c(length(sample_idx), n_subjects, n_contrasts))

    for (i in seq_along(sample_idx)) {
      cid <- cluster_ids[sample_idx[i]]
      path <- paste0("/data/cluster_", cid)

      # Read [subjects, timepoints] or [timepoints]
      cluster_data <- h5[[path]][]

      # Temporal aggregation → contrasts (mean, variance, etc.)
      arr[i, , ] <- aggregate_temporal(cluster_data)
    }

    arr
  })

  result
}
```

### 4.4 Path A: Latent

```r
.read_latent <- function(handle, assays, block) {
  k <- handle$k
  subjects <- handle$subjects
  sample_idx <- block$sample %||% seq_len(k)

  result <- lapply(assays, function(nm) {
    arr <- array(NA_real_, dim = c(length(sample_idx), length(subjects), n_contrasts))

    for (j in seq_along(subjects)) {
      scan_id <- subjects[j]
      path <- paste0("/scans/", scan_id, "/embedding")

      # Read [k, n_contrasts] embedding
      embedding <- h5[[path]][]

      # Subset components
      arr[, j, ] <- embedding[sample_idx, , drop = FALSE]
    }

    arr
  })

  result
}
```

---

## 5. Subject-Contrast Axis Composition

### 5.1 Single-File Strategy

**Scenario:** All subjects and contrasts in one HDF5 file.

**Extraction:**
```r
# Path B: Explicit axes
subjects <- h5[["/gds/axes/subjects"]][]
contrasts <- h5[["/gds/axes/contrasts"]][]

# Path A: Inferred
subjects <- h5[["/subjects"]][]  # If present
contrasts <- h5[["/labels"]][]   # From /labels
```

**Read:** Direct slicing with subject/contrast indices.

### 5.2 Multi-File (Subject-Wise)

**Scenario:** One file per subject, all contrasts within.

**Files:**
```
sub-01_data.h5  # Contains cope1, cope2, cope3
sub-02_data.h5
...
```

**Extraction:**
```r
# Subjects from filenames
subjects <- vapply(files, function(f) {
  base <- tools::file_path_sans_ext(basename(f))
  if (grepl("^sub-", base)) {
    sub("_.*$", "", base)  # "sub-01_data" → "sub-01"
  } else {
    base
  }
}, character(1))

# Contrasts from first file's /labels
h5 <- hdf5r::H5File$new(files[1], mode = "r")
contrasts <- h5[["/labels"]][]
h5$close()
```

**Read:** Loop over files (subject axis), read all contrasts per file.

### 5.3 Multi-File (Contrast-Wise)

**Scenario:** One file per contrast, all subjects within.

**Files:**
```
cope1.h5  # Contains sub-01, sub-02, ...
cope2.h5
cope3.h5
```

**Extraction:**
```r
# Contrasts from filenames
contrasts <- vapply(files, function(f) {
  tools::file_path_sans_ext(basename(f))
}, character(1))

# Subjects from first file's /subjects or structure
h5 <- hdf5r::H5File$new(files[1], mode = "r")
subjects <- h5[["/subjects"]][]
h5$close()
```

**Read:** Loop over files (contrast axis), read all subjects per file.

### 5.4 Multi-File (Subject × Contrast Grid)

**Scenario:** Separate file per subject-contrast combination.

**Files:**
```
sub-01_cope1.h5
sub-01_cope2.h5
sub-02_cope1.h5
sub-02_cope2.h5
...
```

**Extraction:**
```r
# Parse filenames
file_info <- data.frame(
  file = files,
  subject = extract_subject(files),
  contrast = extract_contrast(files)
)

subjects <- unique(file_info$subject)
contrasts <- unique(file_info$contrast)
```

**Read:** Nested loop over subjects and contrasts, load file per (subject, contrast) pair.

**Note:** This is inefficient; prefer aggregating into single file or subject-wise files.

---

## 6. Path B: Adding `/gds` to Existing fmristore Files

### 6.1 Rationale

**Why Add `/gds`?**
1. **Perfect round-tripping**: GDS → fmristore → GDS preserves all metadata
2. **Multi-assay support**: Store beta, var, se, t, z, df, p in one file
3. **Provenance**: Record pipeline operations, timestamps, versions
4. **Alignments**: Persist subject→group maps, reuse across sessions
5. **Explicit axes**: No filename parsing; subjects/contrasts in `/gds/axes`

### 6.2 Augmentation Strategy

**Non-Destructive Approach:**
- Add `/gds` group alongside existing `/data`, `/mask`, `/labels`
- Original fmristore structure remains intact
- Tools using fmristore API continue to work
- Tools using gdsfmri gain full GDS capabilities

**Example Workflow:**
```r
# Read existing fmristore file
plan <- gds("labeled_volume.h5")  # Uses Path A adapter

# Run analysis pipeline
result <- plan %>%
  derive(c("var", "se", "t")) %>%
  align(map = mni_warp) %>%
  reduce(method = "meta:pm") %>%
  compute()

# Write back with /gds group
write_gds(result, "labeled_volume.h5", mode = "update")
# Now file has both legacy structure AND /gds group
```

### 6.3 File Structure After Augmentation

```
/header/*          # Legacy fmristore
/mask              # Legacy
/labels            # Legacy
/data/             # Legacy
  /cope1
  /cope2
  ...

/gds/              # NEW: GDS shim
  /version         # "gds-h5/1.0"
  /axes/
    /subjects      # ["sub-01", "sub-02", ...]
    /contrasts     # ["cope1", "cope2", ...]
  /space/
    /type          # "voxel"
    /dim           # [64, 64, 30]
    /affine        # [4, 4]
    /mask_idx      # Reference to /mask or packed indices
  /assays/
    /beta          # [15342, 20, 3] packed
    /var           # [15342, 20, 3]
    /se            # ...
    /t
    /meta/
      /beta.json   # {"role": "effect", "units": "BOLD % signal change"}
  /provenance      # JSONL log of operations
  /alignments/     # Optional: stored maps
    /mni_warp/
      /family_type
      /serialized
```

**Compatibility:**
- fmristore tools: Read `/data`, `/mask`, `/labels` (ignore `/gds`)
- gdsfmri: Prefer `/gds` when present (Path B), fall back to Path A

### 6.4 Writing `/gds` Group

```r
.write_gds_shim <- function(h5, gds_obj, mode = c("create", "update")) {
  mode <- match.arg(mode)

  if (mode == "update" && h5$exists("/gds")) {
    # Overwrite existing /gds
    h5$link_delete("/gds")
  }

  # Create groups
  h5$create_group("/gds")
  h5$create_group("/gds/axes")
  h5$create_group("/gds/space")
  h5$create_group("/gds/assays")

  # Write version
  h5[["/gds/version"]] <- "gds-h5/1.0"

  # Write axes
  h5[["/gds/axes/subjects"]] <- gds_obj$subjects
  h5[["/gds/axes/contrasts"]] <- gds_obj$contrasts

  # Write space
  space <- gds_obj$space
  h5[["/gds/space"]]$attr("type", space$type)

  if (inherits(space, "space_voxel")) {
    h5[["/gds/space/dim"]] <- space$dim
    h5[["/gds/space/affine"]] <- as.vector(t(space$affine))
    h5[["/gds/space/mask_idx"]] <- space$mask_idx
  } else if (inherits(space, "space_parcels")) {
    h5[["/gds/space/labels"]] <- space$labels
  } else if (inherits(space, "space_basis")) {
    h5[["/gds/space/k"]] <- space$k
    h5[["/gds/space/basis_name"]] <- space$basis_name
    if (!is.null(space$projector)) {
      h5[["/gds/space/basis_matrix"]] <- space$projector
    }
  }

  # Write assays
  for (name in names(gds_obj$assays)) {
    path <- paste0("/gds/assays/", name)
    h5$create_dataset(
      path,
      robj = gds_obj$assays[[name]],
      dtype = "float32",
      chunk_dims = c(min(10000, dim(gds_obj$assays[[name]])[1]), 1, 1),
      compression = "gzip"
    )
  }

  # Write provenance
  if (!is.null(gds_obj$provenance)) {
    .write_provenance_log(h5, gds_obj$provenance)
  }

  invisible(NULL)
}
```

---

## 7. Error Handling and Edge Cases

### 7.1 Missing Groups/Datasets

**Issue:** File claims to be fmristore but missing expected paths.

**Solution:**
```r
.safe_h5_read <- function(h5, path, default = NULL) {
  if (!h5$exists(path)) {
    warning("Expected path '", path, "' not found; using default", call. = FALSE)
    return(default)
  }

  tryCatch(
    h5[[path]][],
    error = function(e) {
      warning("Failed to read '", path, "': ", e$message, call. = FALSE)
      default
    }
  )
}
```

### 7.2 Incompatible Masks Across Files

**Issue:** Multi-file labeled_volume where each subject has different mask.

**Solutions:**

**A. Intersection (Conservative):**
```r
masks <- lapply(files, function(f) {
  h5 <- hdf5r::H5File$new(f, mode = "r")
  on.exit(h5$close())
  h5[["/mask"]][] > 0
})
mask_intersection <- Reduce(`&`, masks)
```
- **Pro**: Only voxels present in ALL subjects
- **Con**: May exclude valid data from some subjects

**B. Union (Inclusive):**
```r
mask_union <- Reduce(`|`, masks)
```
- **Pro**: Includes all voxels from ANY subject
- **Con**: NA values for subjects missing that voxel

**C. Subject-Specific Masks (Most Flexible):**
```r
# Store mask per subject in metadata
metadata$subject_masks <- lapply(files, function(f) {
  h5 <- hdf5r::H5File$new(f, mode = "r")
  on.exit(h5$close())
  which(as.vector(h5[["/mask"]][] > 0))
})

# Use union for space definition
# Apply subject-specific masks during read
```

**Recommendation:** Union with subject-specific masking (C).

### 7.3 Missing Assays

**Issue:** User requests `var` but only `beta` exists in file.

**Solutions:**

**A. Return NA Array (Graceful Degradation):**
```r
if (!h5$exists(path)) {
  if (assay_name %in% c("var", "se")) {
    # Variance/SE are optional; return NA placeholder
    return(array(NA_real_, dim = expected_dims))
  } else {
    # Beta/primary assays are required; error
    stop("Required assay '", assay_name, "' not found", call. = FALSE)
  }
}
```

**B. Skip Missing Assays:**
```r
assays_available <- intersect(assays_requested, probe_result$assays)
if (length(assays_available) < length(assays_requested)) {
  warning("Some requested assays not available: ",
          paste(setdiff(assays_requested, assays_available), collapse = ", "),
          call. = FALSE)
}
# Read only available assays
```

**Recommendation:** Combination (A for var/se, B with warning for others).

### 7.4 Schema Version Mismatch

**Issue:** `/gds/version` is unknown or incompatible.

**Solution:**
```r
.validate_gds_version <- function(h5) {
  if (!h5$exists("/gds/version")) {
    warning("No /gds/version found; assuming gds-h5/0.1", call. = FALSE)
    return("gds-h5/0.1")
  }

  version <- as.character(h5[["/gds/version"]][])
  supported <- c("gds-h5/0.1", "gds-h5/1.0")

  if (!version %in% supported) {
    stop("Unsupported /gds schema version: ", version, "\n",
         "Supported versions: ", paste(supported, collapse = ", "),
         call. = FALSE)
  }

  version
}
```

### 7.5 Temporal Aggregation (Parcellated)

**Issue:** Parcellated data has [clusters, timepoints]; need to map to contrasts.

**Solutions:**

**A. Mean Across Time:**
```r
aggregate_temporal <- function(cluster_data) {
  # cluster_data: [timepoints] or [subjects, timepoints]
  if (is.vector(cluster_data)) {
    mean(cluster_data, na.rm = TRUE)
  } else {
    rowMeans(cluster_data, na.rm = TRUE)
  }
}
```

**B. Variance Across Time:**
```r
aggregate_temporal_var <- function(cluster_data) {
  apply(cluster_data, 1, var, na.rm = TRUE)
}
```

**C. Time Windows (Custom Contrasts):**
```r
# Define contrasts as time windows
contrasts <- list(
  early = 1:50,
  middle = 51:100,
  late = 101:150
)

for (contrast_name in names(contrasts)) {
  time_idx <- contrasts[[contrast_name]]
  arr[, , contrast_name] <- rowMeans(cluster_data[, time_idx], na.rm = TRUE)
}
```

**Recommendation:** Expose aggregation strategy as option; default to mean (A).

---

## 8. Implementation Roadmap

### Phase 1: Core Infrastructure (2-3 days)

**Tasks:**
- [ ] Create `R/adapter-fmristore.R` skeleton
- [ ] Implement `register_fmristore_adapter()`
- [ ] Implement `.fmristore_detect()` with scoring logic
- [ ] Implement `.fmristore_open()` and `.fmristore_close()`
- [ ] Add to `NAMESPACE` and `.onLoad()`

**Deliverables:**
- Adapter registration working
- Detection for `/gds` layout (Path B)
- Basic handle creation

**Testing:**
- Create mock HDF5 with `/gds` group
- Verify detection score = 1.0
- Verify handle structure

---

### Phase 2: Path B - /gds Layout Support (3-4 days)

**Tasks:**
- [ ] Implement `.probe_gds_layout()` for voxel/parcels/basis spaces
- [ ] Implement `.read_gds_layout()` with block support
- [ ] Add space object construction for all types
- [ ] Add alignment/map reading from `/gds/alignments`
- [ ] Add metadata extraction from `/gds/metadata`

**Deliverables:**
- Full Path B reading capability
- End-to-end test: `gds("/gds_file.h5") %>% compute()`
- Space objects correctly constructed

**Testing:**
- Create synthetic `/gds` HDF5 files for each space type
- Test probe returns correct space/assays/axes
- Test read with block slicing
- Test alignment persistence

---

### Phase 3: Path A - Labeled Volume (4-5 days)

**Tasks:**
- [ ] Implement `.probe_labeled_volume()` for single file
- [ ] Implement `.read_labeled_volume()` for single file
- [ ] Add affine construction from qform/sform
- [ ] Implement multi-file probe (subject-wise)
- [ ] Implement multi-file read with subject loop
- [ ] Add mask harmonization (union/intersection)

**Deliverables:**
- Single-file labeled_volume support
- Multi-file subject-wise support
- Affine extraction from NIfTI-like headers

**Testing:**
- Use fmristore test fixtures or create synthetic labeled_volume files
- Test single-file read with mask
- Test multi-file with different masks (union/intersection)
- Compare extracted affines with expected values

---

### Phase 4: Path A - Parcellated (3-4 days)

**Tasks:**
- [ ] Implement `.probe_parcellated()`
- [ ] Implement `.read_parcellated()`
- [ ] Build parcel operator from cluster_map
- [ ] Add temporal aggregation strategies
- [ ] Create voxel→parcel LinearMap

**Deliverables:**
- Parcellated reading with cluster membership
- Temporal aggregation to contrasts
- LinearMap for voxel↔parcel

**Testing:**
- Create synthetic parcellated file with cluster_map
- Test parcel extraction and membership
- Test temporal aggregation
- Test map construction

---

### Phase 5: Path A - Latent (3-4 days)

**Tasks:**
- [ ] Implement `.probe_latent()`
- [ ] Implement `.read_latent()`
- [ ] Extract basis matrix and embeddings
- [ ] Create basis→voxel LinearMap
- [ ] Add reference voxel_space handling

**Deliverables:**
- Latent reading with basis projector
- Embedding stacking across subjects
- LinearMap for basis↔voxel

**Testing:**
- Create synthetic latent file with basis matrix
- Test embedding extraction
- Test basis→voxel projection
- Test integration with map_to()

---

### Phase 6: Path B Writing - /gds Augmentation (2-3 days)

**Tasks:**
- [ ] Implement `.write_gds_shim()` function
- [ ] Add to `write_gds()` verb (mode = "update")
- [ ] Test augmentation of existing fmristore files
- [ ] Verify non-destructive addition of `/gds`

**Deliverables:**
- Ability to add `/gds` to existing files
- Round-trip: read (Path A) → analyze → write (Path B) → read (Path B)

**Testing:**
- Start with legacy labeled_volume
- Run analysis pipeline
- Write with `/gds` group
- Verify original structure intact
- Verify `/gds` correctly formatted

---

### Phase 7: Testing and Documentation (3-4 days)

**Tasks:**
- [ ] Write comprehensive unit tests (50+ test cases)
- [ ] Write integration tests with real pipelines
- [ ] Add error handling tests (edge cases)
- [ ] Write vignette: "Using fmristore with gdsfmri"
- [ ] Update NEWS.md for Sprint 7
- [ ] Add examples to documentation

**Deliverables:**
- `tests/testthat/test-adapter-fmristore.R` with >90% coverage
- Integration vignette with real examples
- API documentation (roxygen2)

**Test Coverage:**
- Detection (all paths, scores)
- Probe (all layouts, single/multi-file)
- Read (all layouts, with/without blocks)
- Space construction (voxel/parcels/basis)
- Map construction (parcel, latent)
- Error handling (missing paths, incompatible masks)
- Round-trip (Path A → analyze → Path B → read)

---

### **Total Estimated Effort: 20-27 person-days (~4-5 weeks)**

---

## 9. Risk Assessment and Mitigation

### 9.1 Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| **fmristore API changes** | Low | High | Pin fmristore version; use defensive programming; test against multiple versions |
| **Mask incompatibility** | Medium | Medium | Implement union/intersection/per-subject strategies; document trade-offs |
| **Performance degradation** | Low | Medium | Profile multi-file reads; consider caching; use HDF5 chunking |
| **Schema evolution** | Medium | Low | Version `/gds` schema; maintain backward compatibility |
| **Missing metadata** | High | Low | Use sensible defaults; warn users; document assumptions |

### 9.2 Integration Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| **Conflicts with fmristore** | Low | High | Keep `/gds` separate; non-destructive augmentation |
| **User confusion (Path A vs B)** | Medium | Low | Clear documentation; automatic detection; helpful warnings |
| **Legacy file handling** | Medium | Medium | Extensive testing with real fmristore files; user feedback |

### 9.3 Maintenance Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| **Adapter complexity** | High | Medium | Modular design; clear routing logic; comprehensive comments |
| **Test coverage gaps** | Medium | Medium | Aim for >90% coverage; synthetic + real fixtures |
| **Documentation drift** | Medium | Low | Link docs to code; update in lock-step |

---

## 10. Success Criteria

### 10.1 Functional Criteria

- [ ] **Detection**: Correctly identifies all fmristore layouts with appropriate scores
- [ ] **Path B**: Reads `/gds` files with full fidelity (all assays, spaces, maps, provenance)
- [ ] **Path A**: Reads all three legacy layouts (labeled_volume, parcellated, latent)
- [ ] **Multi-file**: Handles subject-wise and contrast-wise file organizations
- [ ] **Spaces**: Constructs correct space objects (voxel, parcels, basis) with all metadata
- [ ] **Maps**: Builds LinearMap objects for parcel and basis spaces
- [ ] **Block reading**: Supports efficient partial reads via block parameter
- [ ] **Round-trip**: Path A → analyze → Path B write → Path B read preserves data

### 10.2 Performance Criteria

- [ ] **Read speed**: Comparable to native fmristore (within 10%)
- [ ] **Memory efficiency**: Streaming/block reads don't load full datasets
- [ ] **Multi-file overhead**: <20% overhead vs single-file for same data volume

### 10.3 Quality Criteria

- [ ] **Test coverage**: >90% line coverage for adapter code
- [ ] **Documentation**: Vignette, function docs, examples
- [ ] **Error messages**: Clear, actionable guidance for users
- [ ] **Warnings**: Informative when using defaults or handling edge cases

---

## 11. Integration with Compute Pipeline

### 11.1 Lazy Pipeline Flow

```r
# User code
plan <- gds("fmristore_file.h5") %>%
  subset(subject = c("sub-01", "sub-02")) %>%
  derive(c("var", "se", "t")) %>%
  align(map = mni_warp_family) %>%
  mask(MaskPolicy("group", "threshold", 0.95)) %>%
  reduce(method = "meta:pm")

result <- compute(plan)
```

**What happens:**

1. **gds() call:**
   - `detect_adapter("fmristore_file.h5")` → score 0.95 (Path A) or 1.0 (Path B)
   - `adapter$open()` → handle with file paths, type, has_gds flag
   - `adapter$probe()` → extract space, assays, subjects, contrasts
   - Build `gds_plan` with source node

2. **subset() call:**
   - Add `op_subset` node to plan
   - No data read yet (lazy)

3. **derive() call:**
   - Add `op_derive` node
   - Marks assays to compute: var, se, t (from beta)

4. **align() call:**
   - Add `op_align` node
   - Store map family reference

5. **mask() call:**
   - Add `op_mask` node
   - Store policy

6. **reduce() call:**
   - Add `op_reduce` node
   - Reference to `meta:pm` reducer

7. **compute() call:**
   - **Execute plan nodes in order:**

   **A. Open source:**
   ```r
   handle <- adapter$open(plan$source$path)
   probe <- adapter$probe(handle)
   ```

   **B. Determine required assays:**
   ```r
   # derive() needs beta to compute var, se, t
   # reduce() needs beta, var
   required_assays <- c("beta")
   ```

   **C. Block streaming loop:**
   ```r
   for (block_start in seq(1, n_samples, block_size)) {
     block <- list(sample = block_start:(block_start + block_size - 1))

     # Read block
     arrays_block <- adapter$read(handle, c("beta"), block)

     # Apply subset
     arrays_block <- .apply_subset(arrays_block, plan$nodes$subset)

     # Derive assays
     arrays_block <- .derive_stats(arrays_block, c("var", "se", "t"))

     # Align (apply map_family to each subject)
     arrays_block <- .apply_align(arrays_block, plan$nodes$align)

     # Mask (apply group policy)
     arrays_block <- .apply_mask(arrays_block, plan$nodes$mask)

     # Reduce (call registered meta:pm reducer)
     result_block <- .apply_reduce(arrays_block, plan$nodes$reduce)

     # Write to sink
     sink$write_block(result_block, block)
   }
   ```

   **D. Finalize:**
   ```r
   adapter$close(handle)
   result_gds <- sink$finalize()
   ```

**Key Points:**
- fmristore adapter only involved in open/probe/read/close
- Rest of pipeline is adapter-agnostic
- Block streaming keeps memory bounded
- Subject/contrast dimensions typically small (read in full)

---

## 12. Comparison: Native fmristore vs gdsfmri Adapter

| Aspect | fmristore Native | gdsfmri Adapter |
|--------|-----------------|----------------|
| **Data Model** | Layout-specific (H5NeuroVol, H5ParcellatedScan, LatentNeuroVec) | Unified [sample × subject × contrast] |
| **Spaces** | Implicit (from headers, cluster maps, basis) | Explicit (space_voxel, space_parcels, space_basis) |
| **Assays** | Layout-specific paths (`/data/<label>`) | Unified assay registry (beta, var, se, t, z, df, p) |
| **Subjects/Contrasts** | Inferred from file structure | Explicit axes (`/gds/axes/subjects`, `/gds/axes/contrasts`) |
| **Masks** | Per-file binary mask | Harmonized across files (union/intersection) |
| **Alignment** | Not stored | Persistent AlignmentFamily in `/gds/alignments` |
| **Provenance** | Not tracked | Full pipeline log in `/gds/provenance` |
| **Multi-assay** | Separate files or paths | Single file with multiple `/gds/assays/*` |
| **Lazy Execution** | Eager (load data immediately) | Lazy (build plan, execute on compute()) |
| **Block Streaming** | Manual chunking | Automatic block-wise reads |
| **Pipelines** | Imperative R code | Declarative Plan with %>% |
| **Interoperability** | fmristore-specific tools | Works with any gdsfmri pipeline (align, reduce, map_to) |

**When to Use Native fmristore:**
- Exploratory analysis with known layout
- Single-subject visualization
- Direct manipulation of spatial components
- Integration with existing fmristore workflows

**When to Use gdsfmri Adapter:**
- Group-level meta-analysis
- Multi-subject pipelines (align, reduce)
- Space transformations (voxel ↔ parcel ↔ basis)
- Provenance tracking and reproducibility
- Integration with fmrireg (meta-analysis, spatial FDR)

---

## 13. Open Questions and Decisions Needed

### 13.1 Design Decisions

**Q1: Mask harmonization policy (multi-file)?**
- **Options**: Intersection (conservative), Union (inclusive), Per-subject (flexible)
- **Recommendation**: Union with per-subject masking metadata
- **Decision Point**: Phase 3 implementation

**Q2: Temporal aggregation (parcellated)?**
- **Options**: Mean, Variance, Custom windows, All timepoints as separate contrasts
- **Recommendation**: Expose aggregation function as option; default to mean
- **Decision Point**: Phase 4 implementation

**Q3: `/gds` writing mode?**
- **Options**: Always create new `/gds`, Update if exists, Never overwrite
- **Recommendation**: Add mode parameter: `write_gds(..., mode = c("create", "update", "replace"))`
- **Decision Point**: Phase 6 implementation

**Q4: fmristore version compatibility?**
- **Options**: Support latest only, Support last 2 versions, Dynamic detection
- **Recommendation**: Test against last 2 stable versions; document requirements
- **Decision Point**: Phase 1 (dependency specification)

### 13.2 Missing Information

**Q5: fmristore write API for `/gds`?**
- **Status**: Unknown if fmristore provides hooks for custom groups
- **Action**: Review fmristore source; may need direct hdf5r writes
- **Impact**: Phase 6 implementation

**Q6: Performance of multi-file reads?**
- **Status**: Uncertain how HDF5 file open/close overhead affects streaming
- **Action**: Benchmark during Phase 3; consider file handle caching
- **Impact**: Performance criteria

**Q7: Parcel membership sparsity?**
- **Status**: Unknown optimal representation (list vs sparse matrix vs full bitmap)
- **Action**: Profile memory usage during Phase 4
- **Impact**: Space object design

### 13.3 User Feedback Needed

**Q8: Default behavior for missing assays?**
- **Options**: Error, Warning + NA array, Warning + skip
- **Recommendation**: Get user feedback on expected behavior
- **Decision Point**: Phase 2

**Q9: Multi-file naming conventions?**
- **Status**: Need real-world examples of subject/contrast naming in fmristore files
- **Action**: Survey existing fmristore datasets; document conventions
- **Decision Point**: Phase 3

---

## 14. Next Steps

### Immediate Actions (This Week)

1. **Create `R/adapter-fmristore.R`** skeleton from design document
2. **Add to package infrastructure**:
   - Update `NAMESPACE`
   - Add to `.onLoad()` registration
   - Update `DESCRIPTION` Suggests: fmristore
3. **Create test infrastructure**:
   - `tests/testthat/test-adapter-fmristore.R`
   - Synthetic HDF5 fixtures (voxel, parcels, basis)
4. **Begin Phase 1**: Detection and registration

### Sprint 7 Timeline (4-5 Weeks)

- **Week 1**: Phases 1-2 (Core + Path B)
- **Week 2**: Phase 3 (Labeled Volume)
- **Week 3**: Phases 4-5 (Parcellated + Latent)
- **Week 4**: Phase 6 (Path B Writing)
- **Week 5**: Phase 7 (Testing + Documentation)

### Milestones

- **M1 (End Week 1)**: Path B reading works end-to-end
- **M2 (End Week 2)**: Labeled volume (single + multi-file) works
- **M3 (End Week 3)**: All three layouts readable
- **M4 (End Week 4)**: Round-trip Path A→B works
- **M5 (End Week 5)**: >90% test coverage, documentation complete

### Collaboration Points

- **fmristore package maintainers**: Clarify write API, discuss `/gds` addition
- **gdsfmri users**: Gather feedback on multi-file conventions, error handling
- **Testing**: Request real fmristore datasets for validation

---

## 15. Conclusion

The **fmristore adapter for gdsfmri** is a well-designed, comprehensive solution that:

1. **Preserves fmristore's strengths**: Efficient HDF5 storage, spatial indexing, established layouts
2. **Adds gdsfmri's capabilities**: Unified data model, lazy pipelines, provenance, multi-assay support
3. **Provides migration path**: Path A (legacy) → Path B (/gds augmentation) → full interop
4. **Maintains compatibility**: Non-destructive `/gds` addition; both ecosystems coexist

**Key Innovations:**
- **Dual-path strategy** enables immediate use (Path A) and future-proofing (Path B)
- **Space abstraction** unifies voxel/parcel/basis representations
- **Multi-file handling** supports common real-world file organizations
- **Block streaming** maintains memory efficiency for large datasets

**Recommended Action:** Proceed with implementation following the phased roadmap (Section 8), with focus on Path B first (highest ROI) followed by Path A layouts in order of user priority.

---

**End of Report**

**Prepared by:** gdsfmri Integration Team
**References:**
- notes/TECHNICAL_SPECIFICATION.md (§7.6)
- fmrigds_blueprint.md (§§48-55)
- notes/sprint7.md
- notes/sprint_plan.md
- fmristore package documentation
- Sub-agent analyses (Agents 1-2)
