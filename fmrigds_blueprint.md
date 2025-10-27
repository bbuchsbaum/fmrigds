# fmrigds Blueprint: A Format‑Agnostic Group Data Set (GDS) for fMRI

This blueprint defines a compact, reusable, lazy, and format‑agnostic representation for first‑level fMRI results suitable for group/meta analysis. It unifies dense voxel maps, parcels/clusters, latent bases, and surfaces under one object model while integrating cleanly with an existing `group_data` API and enabling a small standalone package (e.g., `gdsfmri`).

The core of the design is a single object, a Group Data Set (GDS), whose data model guarantees predictable dimensions, explicit spaces, and derivable statistics. Adapters normalize inputs from common formats (NIfTI, HDF5, tabular, in‑memory fits) into a stable schema. The system is big‑data friendly by construction via delayed/on‑disk arrays.

---

## 1. Vision and Goals

- One object for first‑level results across spaces (voxels, parcels, bases, surfaces, etc.).
- Fixed, predictable dimensionality: `[sample × subject × contrast]` to simplify downstream code.
- Lazy/on‑disk arrays to scale to whole‑brain voxels and large cohorts.
- Clean migration path from current `group_data_*` constructors to a future standalone package.
- Format‑agnostic adapters that normalize diverse inputs into the canonical schema.
- Space‑aware transforms and I/O, without forcing all data into voxel space.

---

## 2. Core Object: GDS

The `gds` object carries:

- Assays (aka statistical layers): `beta`, `se` or `var`, optionally `t`, `z`, `df`, `n_eff`, etc.
- Three canonical axes in fixed order: `sample` (locations in a space), `subject` (first‑level fits), `contrast` (within‑subject contrast names/ids).
- A `space` descriptor that defines what a “sample” means and how to index/map it.
- Covariates and metadata: `col_data` (subject covariates), `row_data` (sample metadata), `metadata` (schema_version, units, provenance, software, alignment, etc.).

Data model guarantees:

- If `{beta + se|var}` are present, `t`, `z`, etc. are derivable lazily.
- If `{t + df}` are present but `beta` is not, combination methods (e.g., Stouffer) remain valid; optional converters can be used when scale info exists.

---

## 3. Minimal Schema (Stable “Wire Format”)

Required members:

- `assays`: named list of 3‑D arrays (or delayed/on‑disk arrays) with shape `[sample × subject × contrast]`.
  - Required combinations: `{beta, se}` or `{beta, var}` or `{t, df}`.
  - Optional: `z`, `stderr_beta` (alias of `se`), `s2` (alias of `var`), `n_eff`, etc.
- `space`: object describing the `sample` axis (see Spaces).
- `subjects`: character vector of length `dim(assays[[1]])[2]`.
- `contrasts`: character vector of length `dim(assays[[1]])[3]`.

Optional but recommended:

- `col_data`: data.frame keyed by subjects (covariates like group, age, etc.).
- `row_data`: data.frame keyed by samples (ROI names, cluster extent, vertex coords, basis meta, etc.).
- `metadata`: list with `schema_version`, `provenance`, `units`, `mask_info`, `design_mats`, etc.

Invariants:

- All assay arrays share identical dims and dimnames.
- If `beta` exists, you must have either `se` or `var`.
- If `t` exists, you should have `df`.
- Units recorded per assay in `metadata$units` (e.g., `% BOLD`, a.u.).

---

## 4. Spaces: One Abstraction, Multiple Kinds

A `space` object explains what a `sample` is and how to index or map it across spaces. Examples (S3 constructors):

- `space_voxels(affine, dim, mask_idx = NULL, template_id = NULL)`
- Alias: `space_voxel()` in the technical spec mirrors this constructor.
- `space_surface(vertices, faces, hemi, template_id = NULL)`
- `space_parcels(labels, lookup, membership = NULL)`
- `space_basis(B, basis_name, voxel_space = NULL)`

Details:

- Voxel space: carries `dim`, `affine` (or qform/sform), and optional sparse `mask_idx` selecting samples.
- Parcel/cluster space: `labels` for samples, optional `membership` mapping to voxel indices (sparse list or `[samples × voxels]` sparse matrix) for aggregation.
- Basis space: holds `(voxels × k)` projector `B` (possibly sparse) and optional reference voxel space for projection.
- Surface space: mesh geometry (vertices, faces) and hemisphere information.

Mapping between spaces uses linear operators (sparse matrices/callbacks). Examples:

- Voxel → parcels: mean/weighted mean by membership.
- Voxel → basis: `Bᵀ y` or `y = B c` depending on direction.
- Voxel ↔ surface: barycentric/ribbon mapping (pluggable; not all methods need to ship initially).

Mapping is optional. If first‑level results are ROI‑level, the space is `space_parcels()`. Methods that require voxels can request a mapping or error gracefully.

---

## 5. Lazy, On‑Disk Data

- Support base arrays or delayed arrays:
  - Bioconductor `DelayedArray` + `HDF5Array` for R work well initially.
  - Thin backend interface allows alternates (TileDB, Arrow) later.
- Adapters from NIfTI/IMG are adapter‑based and yield on‑disk assays—no large in‑RAM copies required.

---

## 6. Canonical API (S3)

- Constructor:
  - `new_gds(assays, space, subjects, contrasts, col_data = NULL, row_data = NULL, metadata = list())`
- Accessors:
  - `assay(gds, name = "beta")` → array `[sample × subject × contrast]`
  - `assays(gds)`, `space(gds)`, `subjects(gds)`, `contrasts(gds)`
  - `col_data(gds)`, `row_data(gds)`, `metadata(gds)`
- Verbs:
  - `subset_gds(gds, sample = NULL, subject = NULL, contrast = NULL)`
  - `derive(gds, what = c("var","se","t","z"), overwrite = FALSE)`
  - `as_long(gds, stats = c("beta","se"), drop_na = TRUE)`
  - `map_to(gds, target_space, op = c("mean","matmul","sum"), ...)` (space transforms)
- Meta‑analysis/helpers:
  - `weights(gds, kind = c("fixed","random"), clamp = TRUE)`
  - `combine_fixed(gds)` (group maps in same space)
  - `combine_stouffer(gds)` (if `z` or `t+df` only)

Dimensional order is fixed: `[sample × subject × contrast]` for predictable vectorization.

---

## 6A. Lazy vs Eager API

- Lazy‑first (recommended): all verbs return a `gds_plan`; `compute()` materializes to a `gds_result` (see §§25–27).
  - Accepts `gds_source` or `gds_plan` on input; returns an augmented `gds_plan`.
  - Example: `res <- compute(map_to(derive(subset_gds(src, subject=ids), c("var","t")), parcels, map, UncertaintyRule("independent")))`.
- Eager wrappers (compatibility/convenience): call the lazy verb and immediately `compute()`.
  - `subset_gds_eager(x, ...)`, `derive_eager(x, ...)`, `map_to_eager(x, ...)` mirror lazy signatures and accept `compute()` options.
- Back‑compat shims: existing `group_data_from_*()` can return either a realized `gds` or a `gds_plan` based on an option (e.g., `options(gds.lazy = TRUE)`), with `as_gds()` and `as_group_data()` bridging during migration.
- Provenance: eager wrappers embed the plan digest and compute options into `metadata$provenance` of the realized object for traceability.

---

## 7. Adapters (“Drivers”) for Input Formats

All adapters create a `gds` without eager loading:

- NIfTI: `gds_from_nifti(beta, se|var, t, df, mask, subjects, contrasts, space_voxels(...))`
  - Validates dims/alignment, stores assays as delayed/on‑disk arrays.
- HDF5: `gds_from_h5(h5_paths, stat_names = c("beta","se","t","df"), ...)`
  - Uses a small, documented HDF5 layout or reads an existing one.
- CSV/Parquet (ROI/basis/cluster summaries): `gds_from_tabular(df_or_path, effect_cols, subject_col, sample_col, contrast_col, space = space_parcels(...))`
- In‑memory first‑level fits: `gds_from_fmrilm(list_of_fits, mask, space_voxels(...))`

Compatibility layer: keep existing `group_data_from_*()` and have them call new adapters and return a `gds`. Provide `as_gds.group_data_*` and `as_group_data.gds` to support migration.

---

## 8. Standardized Derivations

- `{beta, se} → var = se^2`.
- `{beta, var} → se = sqrt(var)`.
- `{t, df}` with known design scale (e.g., `metadata$scale_per_contrast`) allows explicit back‑conversion to effect scale (optional, not automatic).
- Helpers: `to_z(t, df)` and `to_t(z)`; record transforms in `metadata$provenance`.

---

## 9. Space‑Aware Outputs

- `to_nifti(gds, stat, contrast, subject = NULL)` — if space is voxel or can map to voxel.
- `to_surface(gds, ...)` — if space is surface.
- `to_tabular(gds, ...)` — long tibble for parcels/basis outputs.

---

## 10. Package Layout (new package: `gdsfmri`)

R/

- `gds-class.R` — constructor, validators, accessors, print/summary.
- `gds-space.R` — `space_*` constructors, mapping interface.
- `gds-derive.R` — `derive()`, `to_z()`, `to_t()`, `weights()`.
- `gds-subset.R` — `subset_gds()`, `as_long()`.
- `gds-meta.R` — `combine_fixed()`, (`combine_random()` later).
- `adapters-nifti.R` — `gds_from_nifti()`.
- `adapters-h5.R` — `gds_from_h5()`.
- `adapters-tabular.R` — `gds_from_tabular()`.
- `adapters-fmrilm.R` — `gds_from_fmrilm()`.
- `io-nifti.R` — `to_nifti()`, header utilities.
- `utils-validate.R` — invariants and schema checks.

inst/

- `schemas/hdf5.md` — tiny spec documenting the HDF5 layout.
- `schemas/space.md` — space object spec + mapping API.

tests/testthat/

- `test-gds-class.R`, `test-adapters.R`, `test-derive.R`, `test-meta.R`.

---

## 11. R Skeletons (Abbreviated Signatures)

- Class and constructor: `new_gds(assays, space, subjects, contrasts, col_data = NULL, row_data = NULL, metadata = list())` with invariant checks and minimal stat sanity.
- Spaces: `space_voxels()`, `space_parcels()`, `space_basis()`, `space_surface()`; mapping API accepts sparse matrices or callbacks.
- Subsetting and tidy view: `subset_gds()`, `as_long()` preserving dimensionality and dimnames.
- Derivations: `derive()` computes `var`, `se`, `t`, `z` as available; provenance recorded in metadata.
- Adapters: thin adapters that validate alignment, assemble assays, and populate spaces/axes.

---

## 12. Migration Path from Current `group_data`

- Keep existing generics (`group_data()`, print/summary methods, `n_subjects`, `get_subjects`, `get_covariates`).
- Internally, `group_data_from_*()` builds a `gds` and returns it (or class `c("group_data_gds","group_data")`).
- Provide `as_gds.group_data_*()` to retrieve embedded `gds` and `as_group_data.gds` to shim back, ensuring downstream compatibility during transition.

---

## 13. HDF5/Zarr Mini‑Spec (Pragmatic Layout)

Groups/paths:

- `/assays/beta` (float32, chunked) `[sample × subject × contrast]`
- `/assays/se`, `/assays/var`, `/assays/t`, `/assays/df` (as present)
- `/space/type` (attr: `"voxel"|"parcels"|"basis"|"surface"`)
- `/space/voxel/dim` (int, len=3)
- `/space/voxel/affine` (float64, 4×4)
- `/space/parcels/labels` (string, len = sample)
- `/subjects` (string), `/contrasts` (string)
- `/row_data/...` (columns), `/col_data/...` (columns)
- `/metadata/json` (utf8)

Notes:

- If physical order differs, record logical order `[sample, subject, contrast]` as an attribute.
- Chunk along `sample` for voxel/basis; along `subject` for ROI/cluster; compression via zstd/gzip.
- Same logical layout applies if targeting Zarr/TileDB later.

---

## 14. Why This Works Across Voxels, Clusters, Bases, Surfaces

- Voxels: `space_voxels` + on‑disk assays for scale; mapping optional.
- ROIs/clusters: `space_parcels(labels=roi_names)` + tabular adapter; small dense arrays.
- Latent bases: `space_basis(k, projector=B_or_fun)`; first‑level outputs already in component space; `map_to()` projects to voxels/surfaces as needed.
- Surfaces: `space_surface(vertices, faces, hemi)`; assays dimension equals number of vertices; exporters produce GIFTI/curv, etc.

All higher‑level analyses use the same object to subset/derive stats, compute weights, combine (fixed/random), and emit results in the native or mapped space.

---

## 15. Future‑Proof Extras (Optional)

- `contrast_info`: per‑contrast scale (units) and design row to enable disciplined `t`→effect conversions.
- Mask semantics: distinguish “masked out” vs “missing” via NA policy and binary sample mask in `row_data`.
- QC utilities: variance histograms, missingness maps, leverage diagnostics—identical across spaces.
- Typed IDs with `vctrs::rcrd` (`SubjectId`, `ContrastId`) for stricter joins.

---

## 16. What to Ship Immediately

1. Define `gds` (constructor + validators + accessors).
2. Implement `gds_from_tabular()` fully (including wide→long normalization where needed).
3. Implement a minimal `gds_from_nifti()` (can read into RAM first; swap to `HDF5Array` later).
4. Provide `derive()` + `combine_fixed()` (weighted mean with `1/var`).
5. Add `as_long()` for tidy plotting/debugging.
6. Add a compat wrapper in `fmrireg::group_data()` that returns a `gds` under the hood.

---

## 17. Example Usage

ROI table → gds

```
gd <- gds_from_tabular(
  data = "roi_stats.csv",
  effect_cols = c(beta = "beta", se = "se"),
  subject_col = "sub",
  sample_col = "roi",
  contrast_col = "contrast",
  space = space_parcels(labels = NULL)
)

gd <- derive(gd, c("var","t"))   # fills var and t
fixed <- combine_fixed(gd)          # group-level maps in same space
```

NIfTI pair → gds (voxel space)

```
gd_nii <- gds_from_nifti(
  beta = B, se = S,
  subjects = sub_ids,
  contrasts = c("Faces>Places"),
  voxel_dim = c(91,109,91),
  affine = diag(4)
)
```

---

## 18. Alignment With Current Code

- Existing `group_data_*` constructors are the adapters in this blueprint.
- Print/summary map directly to `gds` methods.
- `extract_csv_data()` functionality decomposes into `as_long(gds)`, `subset_gds()`, and `assay(gds, "beta")` slices.
- Keep `group_data` as a thin shim on top of `gds` to avoid breaking user code; new packages can depend directly on `gdsfmri`.

---

## 19. Non‑Goals and Scope Notes

- Random‑effects meta (`combine_random`) and advanced space mappings (full voxel↔surface) are scoped for later; the API is designed to accept them without breaking changes.
- Back‑conversion from `t` to effect size requires explicit scale info and is never automatic.

---

## 20. Glossary

- Assay: A 3‑D statistical layer (e.g., `beta`, `se`, `t`).
- Sample: A location index within a defined space (voxel, vertex, ROI, basis component, etc.).
- Space: Object that defines what samples mean and how to map them.
- Adapter: A reader/normalizer that produces a `gds` from an external format.

---

## 21. Uncertainty Propagation Through Spatial Mappings

Statistical core for linear spatial mappings `y = M x`:

- Location mapping (per subject × contrast): `beta' = M beta`.
- Variance mapping: `Var(beta') = M Σ M^T`, where `Σ = Cov(beta)`.
- Independence default (diagonal `Σ = diag(v)`): `Var(beta')_i = sum_j (M_ij^2 * v_j)`.
- Weighted parcel mean (row weights `w` that sum to 1): `beta' = w^T beta`, `Var(beta') = w^T Σ w`.
- Never average standard errors. Always propagate via variance/covariance.
- Do not map `t`/`z` directly. Recompute after mapping: `t' = beta'/SE'` with `SE' = sqrt(Var(beta'))`.
- If only `z` (or `t+df` without effect scale) are available, use explicit combiners (not spatial maps). Example (Stouffer, independence): `z' = sum_i w_i z_i / sqrt(sum_i w_i^2)`.
- Optional df aggregation (Satterthwaite/Welch) for approximate `df'` when combining variances: `df' ≈ ((sum_i a_i^2 s_i^2)^2) / (sum_i (a_i^4 s_i^4 / df_i))`, with `a_i = M[row,i]` and `s_i^2` unbiased variances.

---

## 22. Assay Registry and Uncertainty Rules

Make mapping assay‑aware via two layers:

- Assay registry defines semantics and units:
  - `register_assay(role = c("location","variance","stdev","z","t","df","n_eff"), name, units = NULL)`
  - Examples: `beta → location`, `var → variance`, `se → stdev` (derived), `t,z → distributional` (never linearly mapped).

- UncertaintyRule controls propagation:
  - `UncertaintyRule(mode = c("independent","cov_provider","kernel","none"), df_rule = c("satterthwaite","none"), cov_provider = NULL, kernel = NULL)`
  - Modes:
    - `independent`: treat `Σ` as diagonal (fast path).
    - `cov_provider`: function supplying `Σ` blocks for source indices per target row.
    - `kernel`: build `Σ` from a spatial correlation model (e.g., exp(−dist/ℓ)).
    - `none`: skip variance propagation (allowed but discouraged).

`map_to` inspects registered assays to choose behavior:

- If location present: map with `M`.
- If variance or stdev present: propagate with rule; derive complementary forms as needed.
- If only `t/df` or `z`: use explicit combiners (e.g., Stouffer) or refuse unless confirmed.

---

## 23. Mapping API: Minimal Surface

Signature (plan‑first API shown; eager wrappers can materialize):

```
map_to(x, target_space,
       map = LinearMap(...),
       uncertainty = UncertaintyRule("independent", df_rule = "satterthwaite"),
       combine = NULL)
```

Independent propagation implementation sketch:

```
.propagate_independent <- function(M, beta, var) {
  beta_out <- M %*% beta
  var_out  <- (M^2) %*% var
  list(beta = beta_out, var = var_out)
}
```

With a covariance provider per target row `r`:

```
# idx <- which(M[r,] != 0); w <- M[r, idx]
# Sigma_block <- cov_provider(idx)
# var_r <- drop(t(w) %*% Sigma_block %*% w)
```

This isolates linear algebra/geometry in `LinearMap` and keeps statistical propagation explicit and extensible.

---

## 24. Provenance: Auditable and Reproducible

Two complementary mechanisms:

- Provenance DAG (always on): every verb adds a node with `op_name`, package version, normalized params, timestamp, seed (if any), input content hashes, and parent pointers. Stored at `metadata$provenance$graph` with a human log at `metadata$provenance$log`.
- Declarative plan: verbs return a plan; the operation history is the object. Serialize/rehydrate via `save_plan()`/`load_plan()`; `digest(plan)` yields a stable hash; `explain(plan)` prints a fused, readable tree.

Utilities:

- `explain(plan)`, `digest(plan)`, `save_plan(plan, "analysis.json")`, `load_plan("analysis.json")`.

---

## 25. Fully Lazy, Declarative Processing Graph

Make every verb build a plan; only `compute()` touches data. Benefits: pushdown optimization, streaming, reproducibility.

Objects:

- `gds_source` — adapter probe/open (no data loaded).
- `gds_plan` — DAG of ops; root is `gds_source`.
- `gds_result` — realized `gds` (in‑memory or on‑disk arrays).

Core verbs (all return plans):

- `subset_gds(plan_or_source, sample = NULL, subject = NULL, contrast = NULL)`
- `derive(plan, what = c("var","se","t","z"), options = list())`
- `map_to(plan, target_space, map, uncertainty = UncertaintyRule("independent"), combine = NULL)`
- `mutate_assay(plan, name, fn)`
- `reduce_subjects(plan, method = c("fixed","random"), weights = c("1/var","n_eff"), by = "contrast")`
- `write_out(plan, path, format = c("h5","nifti","csv","parquet"))`

Execution:

- `compute(plan, sink = c("memory","HDF5"), path = NULL, block = list(sample = 100000), scheduler = c("sequential","multicore","future"), cache = TRUE) -> gds_result`

Planner rewrites and fusion (initial rules):

- Push `subset(subject, contrast)` to the source.
- Push `subset(sample)` before `map_to` when source indices are given.
- Coalesce chains of `derive` into a single pass.
- Fuse `t→z→combine` into one combiner kernel when `beta/var` are absent.

Executor streams blocks from the adapter, applies fused kernels, and writes to the sink.

---

## 26. Minimal Op Types and Plan Serialization

Internal op types:

- `op_read` (adapter root)
- `op_subset_axis` (any axis)
- `op_derive` (var, se, t, z)
- `op_map` (LinearMap + UncertaintyRule or combiner)
- `op_reduce` (meta‑analysis across subjects)
- `op_write` (materialize/export)

Serialized form (JSON example):

```
{
  "version": "0.1.0",
  "nodes": [
    {"id":"n0","op":"read","adapter":"nifti","source":[...],"hash":"..."},
    {"id":"n1","op":"subset_axis","parent":"n0","subject":["sub-01","sub-02"]},
    {"id":"n2","op":"derive","parent":"n1","what":["var","t"]},
    {"id":"n3","op":"map","parent":"n2","space":"parcels","map":"mean",
     "uncertainty":{"mode":"independent","df_rule":"satterthwaite"}}
  ]
}
```

`compute()` on any machine with registered adapters reproduces results deterministically (avoid nondeterministic kernels).

---

## 27. R Skeletons: Plan Core and Execution (Sketch)

Plan core and verbs:

```
gds_plan <- function(source, nodes = list(), meta = list()) {
  structure(list(source = source, nodes = nodes, meta = meta), class = "gds_plan")
}

add_op <- function(plan, node) {
  stopifnot(inherits(plan, "gds_plan"))
  plan$nodes <- append(plan$nodes, list(node))
  plan
}

op_subset_axis <- function(sample=NULL, subject=NULL, contrast=NULL) {
  list(op="subset_axis", sample=sample, subject=subject, contrast=contrast)
}
op_derive <- function(what, options=list()) {
  list(op="derive", what=as.character(what), options=options)
}
op_map <- function(space, map, uncertainty, combine=NULL) {
  list(op="map", space=space, map=map, uncertainty=uncertainty, combine=combine)
}

subset_gds <- function(x, sample=NULL, subject=NULL, contrast=NULL) {
  plan <- if (inherits(x, "gds_plan")) x else gds_plan(source = x)
  add_op(plan, op_subset_axis(sample, subject, contrast))
}

derive <- function(x, what, options=list()) {
  plan <- if (inherits(x, "gds_plan")) x else gds_plan(source = x)
  add_op(plan, op_derive(what, options))
}

map_to <- function(x, target_space, map, uncertainty = UncertaintyRule("independent"), combine=NULL) {
  plan <- if (inherits(x, "gds_plan")) x else gds_plan(source = x)
  add_op(plan, op_map(space = target_space, map = map, uncertainty = uncertainty, combine = combine))
}
```

Planner rewrite (seeded; extend as needed):

```
optimize_plan <- function(plan) {
  nodes <- plan$nodes
  # Implement deterministic rewrites: push subsets, coalesce derives, etc.
  plan$nodes <- nodes
  plan
}
```

Executor (block streaming; pseudocode outline):

```
compute <- function(plan, sink = c("memory","HDF5"), path=NULL, block=list(sample=1e5), cache=TRUE) {
  sink <- match.arg(sink)
  plan <- optimize_plan(plan)
  cur <- open_stream(plan$source, plan)   # adapter cursor with pushed‑down filters
  out <- execute_stream(cur, plan$nodes, sink, path, block)
  out
}
```

---

## 28. Package Layout Additions

- `R/gds-map.R` — LinearMap, UncertaintyRule, assay registry, mapping kernels.
- `R/gds-plan.R` — plan object, ops, optimizer, compute/executor stubs.
- `R/gds-provenance.R` — provenance DAG, explain/digest/save_plan/load_plan.
- `tests/testthat/test-plan.R` — plan construction, rewrite, and serialization.
- `tests/testthat/test-mapping.R` — independent propagation and Stouffer checks.

---

## 29. Practical Next Steps (Updated)

1. Add assay registry + `UncertaintyRule` and wire into `map_to()`; start with `independent` mode.
2. Wrap readers as `gds_source` and have adapters return a plan by default (`gds_from_*()` returns a `gds_plan`).
3. Implement a minimal planner (pushdown subset; coalesce derives) and a block executor that streams source assays → derives → map.
4. Record provenance nodes on each verb; implement `explain(plan)`.
5. Ship Stouffer as a named combiner for `t/z`‑only pipelines with clear docs on assumptions.

---

## 30. Subject→Group Alignment Transforms

Sometimes a shared group space requires per‑subject alignment before analysis:

- Orthogonal/Procrustes: rotate latent bases to a consensus basis.
- Affine/diffeomorphic (3D): volumetric spatial registration with resampling.
- Optimal Transport (OT/Wasserstein): mass‑preserving reweighting toward a barycenter/template.

Model per‑subject transforms as subject‑aware map families attached to the space graph.

Map family abstraction:

```
MapFamily <- function(name, from_space, to_space,
                      type = c("linear","orthogonal","affine3d","deform3d","ot"),
                      maps = NULL,         # list named by subject id: matrix or function
                      traits = list(linear = TRUE, isometry = FALSE, mass_preserving = FALSE),
                      uncertainty = UncertaintyRule("independent", df_rule = "satterthwaite")) {
  type <- match.arg(type)
  structure(list(name = name, from = from_space, to = to_space,
                 type = type, maps = maps, traits = traits,
                 uncertainty = uncertainty),
            class = "gds_map_family")
}
```

- For subject `s`, `maps[[s]]` holds that subject’s operator:
  - Orthogonal/linear: sparse/dense matrix `M_s`.
  - Affine3d/deform3d: resampling function for 3D arrays.
  - OT: transport plan `P_s` (row‑stochastic; mass‑preserving).

Plan verb:

```
align_to_group <- function(plan, family) {
  add_op(plan, list(op = "align_to_group", family = family))
}
```

Execution semantics (per subject):

- Linear/orthogonal/OT families, for each subject `s` and contrast `k`:
  - `beta'(·,s,k) = M_s beta(·,s,k)`
  - `Var(beta') = M_s Σ_s M_s^T`
  - Independence shortcut (`Σ_s = diag(v_s)`): `Var(beta') = (M_s ⊙ M_s) v_s` (elementwise square then matrix–vector multiply).
- Affine3d/deform3d (volumes): resample 3D fields per assay.
  - Resample variance maps (not SE) using the same interpolator; recompute `t/z` afterward.
  - If local covariance is available, use the general `M_s Σ_s M_s^T` on a local stencil.

Traits:

- `isometry=TRUE` (orthogonal) indicates Euclidean norm preservation (QC); variance still changes unless `Σ ∝ I`.
- `mass_preserving=TRUE` (OT) guarantees row sums = 1; use variance rule with `M_s = P_s`.

Integration:

- Add families to the space graph, e.g., `sg$maps <- c(sg$maps, list(family_orth, family_ot))`.
- Planner pushes `align_to_group()` as early as possible (right after subject/contrast subsetting) so downstream `map_to()` and `reduce_subjects()` operate in group space.

---

## 31. Masks and Packed Storage for Voxel Spaces

Distinguish two notions:

- Storage mask: samples physically stored in assays (e.g., in‑brain voxels only).
- Analysis mask: samples used by the analysis step (subset or QC‑derived).

Extend voxel space to encode masks and storage mode:

```
space_voxels <- function(dim, affine,
                         mask_bitmap = NULL,      # logical array [dim], optional
                         mask_idx    = NULL,      # integer vector of linear indices
                         storage = c("dense","packed"),   # packed = only mask voxels stored
                         template_id = NULL) {
  storage <- match.arg(storage)
  if (!is.null(mask_bitmap) && is.null(mask_idx)) mask_idx <- which(as.vector(mask_bitmap))
  structure(list(type = "voxel", dim = dim, affine = affine,
                 mask_bitmap = mask_bitmap, mask_idx = mask_idx,
                 storage = storage, template_id = template_id),
            class = c("space_voxels","gds_space"))
}
```

- `dense`: assays are full `[X×Y×Z]` along sample.
- `packed`: assays store only `length(mask_idx)` samples. The executor expands/remaps as needed.
- Adapters should advertise storage; the executor may repack eagerly into a preferred sink (e.g., packed HDF5).

Mask policy API:

```
MaskPolicy <- function(scope = c("subject","group"),
                       rule  = c("intersection","union","threshold","custom"),
                       threshold = 0.95, custom = NULL) {
  structure(list(scope=match.arg(scope), rule=match.arg(rule),
                 threshold=threshold, custom=custom), class="gds_mask_policy")
}
apply_mask_policy <- function(plan, policy) {
  add_op(plan, list(op="mask_policy", policy=policy))
}
```

- Examples: `scope="group", rule="intersection"` (classic group mask), `rule="threshold", threshold=0.9` (≥90% subjects).
- Planner pushes mask computation after `align_to_group()`, then packs to minimize I/O.

Mapping with masks: linear maps use only source indices in the active mask; uncertainty propagation uses the corresponding rows/columns of `M` and variances.

---

## 32. Plan Integration and Optimizer Rules

New plan nodes:

- `op_align_to_group(family)` — subject‑aware per‑space alignment.
- `op_mask_policy(policy)` — derive analysis mask; may repack storage.

Rewrite rules (optimizer):

1. Push `subset(subject, contrast)` before `align_to_group`.
2. Push `align_to_group` before `mask_policy` (derive group mask in group space).
3. Push `mask_policy` before `derive/map/reduce` (smaller working set).
4. Fuse multiple `derive` steps; recompute any `t/z` after mapping from `beta/var`.

Streaming order per block: read → align (per subject) → mask/pack → derive → (optional) map_to → reduce/write.

---

## 33. Implementation Sketches: Alignment and Masks

Map families (helpers):

```
OrthogonalFamily <- function(name, from_space, to_space, mats_by_subject, uncertainty = UncertaintyRule("independent")) {
  MapFamily(name, from_space, to_space, type="orthogonal",
            maps = mats_by_subject,
            traits = list(linear=TRUE, isometry=TRUE, mass_preserving=FALSE),
            uncertainty = uncertainty)
}

OTFamily <- function(name, from_space, to_space, plans_by_subject, uncertainty = UncertaintyRule("independent")) {
  MapFamily(name, from_space, to_space, type="ot",
            maps = plans_by_subject,
            traits = list(linear=TRUE, isometry=FALSE, mass_preserving=TRUE),
            uncertainty = uncertainty)
}
```

Executor kernel (per block, per subject; independence rule):

```
.apply_family_linear <- function(M, beta_block, var_block, rule) {
  beta_out <- M %*% beta_block                # [target x ncols]
  var_out  <- (M^2) %*% var_block             # propagate variance only
  list(beta = beta_out, var = var_out)
}
```

Mask packing utilities:

```
pack_assays <- function(assays, mask_idx) {
  lapply(assays, function(A) {
    A[mask_idx, , , drop = FALSE]
  })
}

unpack_assays <- function(assays, mask_idx, full_n) {
  lapply(assays, function(A) {
    d <- dim(A); out <- array(NA_real_, dim = c(full_n, d[2], d[3]))
    out[mask_idx, , ] <- A
    out
  })
}
```

---

## 34. Practical Guidance and Defaults

- Default alignment: none (identity). Apply only if a `MapFamily` is provided and `align_to_group()` is in the plan.
- Default mask policy: if adapters provide per‑subject masks, auto‑compute a group intersection mask and pack immediately (configurable).
- Uncertainty propagation:
  - For orthogonal/linear/OT families: use variance rule described above.
  - For 3D resampling (affine/deform): resample variance, never SE; recompute `t/z` after mapping.
- t/z‑only inputs: refuse alignment‑as‑linear unless `beta/var` are present or reconstructable; otherwise use explicit combiners (e.g., Stouffer) only when alignment semantics make sense for z/t (rare).

---

## 35. Why This Solves the Two Cases

1. General subject→group transforms: handled via subject‑aware `MapFamily` nodes (orthogonal, affine, deformable, OT) that integrate with `UncertaintyRule` and optimizer pushdown.
2. 3D mask handling: encoded in `space_voxels` with storage vs analysis masks and packed storage. Adapters can emit packed assays; the plan computes a group mask and repacks in one pass; mappings respect the active mask.

---

## 36. Checklist to Wire Up Now

- Add `MapFamily` and `align_to_group()` (start with `OrthogonalFamily` and `OTFamily`).
- Extend executor to apply per‑subject linear maps and propagate variance (`independent` rule initially).
- Upgrade `space_voxels()` with storage + `mask_idx` semantics; add `pack_assays()` utilities.
- Add `MaskPolicy` + `apply_mask_policy()`; push it early in the optimizer.
- Document “don’t map t/z; recompute from beta/var”; provide a Stouffer combiner as an explicit, separate verb.

---

## 37. Minimal Surface: 5 Nouns + 8 Verbs

This section distills the full design into a compact surface that covers multi‑space data, masks, alignment to a consensus, uncertainty propagation, factorial contrasts, adapters, and full laziness/provenance without heaviness.

---

## 38. The 5 Nouns

1) GDS — a realized group dataset

- Assays: named 3‑D arrays `[sample × subject × contrast]`. Canonical: `beta` (location), `var` (variance). Optional: `se` (derived), `t`, `z`, `df`, `n_eff`.
- Axes: `sample`, `subject`, `contrast` (with dimnames).
- Space: a Space object defining what `sample` means.
- Metadata: `col_data` (per‑subject covariates), `row_data` (per‑sample), `contrast_info`, `provenance`, `units`, `schema_version`.

2) Space — what the sample axis is

- `space_voxel(dim, affine, mask_idx = NULL, storage = c("dense","packed"))`  [alias: `space_voxels()`]
- `space_parcels(labels)`
- `space_surface(vertices, faces, hemi)`
- `space_basis(k, projector = NULL)` (optional mapping to voxels/other spaces)

Masking is first‑class in Space: voxel spaces can be packed (store only `mask_idx`) or dense.

3) Map — linear mapping between spaces (global or per‑subject)

- `map_linear(from, to, operator, by_subject = NULL, traits = list(orthogonal = FALSE, mass_preserving = FALSE), uncertainty = UncertaintyRule("independent"))`
- `operator` can be a sparse/dense matrix or a function; `by_subject` (named list) supports subject‑specific transforms (orthogonal/Procrustes, OT/Wasserstein, affine/deform resampling via function).

4) UncertaintyRule — how to propagate variance through maps

- `UncertaintyRule(mode = c("independent","cov_provider","kernel"), df_rule = c("satterthwaite","none"), cov_provider = NULL, kernel = NULL)`
- Default `independent` uses `(M^2) %*% var`. If covariance is available, plug it via `cov_provider`.

5) Plan — a lazy pipeline (the “recipe”)

- Immutable DAG of steps applied to a source; serializable; auditable provenance.
- `compute(plan, ...)` realizes a GDS (in memory or on disk).

---

## 39. The 8 Verbs (lazy; only compute() touches data)

```
gds()        # front door: create a Plan using the Adapter API
subset()     # select samples / subjects / contrasts
derive()     # fill {var, se, t, z} from what you have (never map t/z directly)
align()      # apply subject-specific Map (orthogonal/OT/affine) into group space
mask()       # derive/apply group mask policy; pack/unpack storage
map_to()     # change spaces using a Map (assay-aware uncertainty propagation)
reduce()     # meta-analysis across subjects (fixed/random; or z-combiners like Stouffer)
write_out()  # export to NIfTI/HDF5/CSV/Parquet (still lazy; executes at compute)
compute()    # execute optimized pipeline; returns a realized GDS
```

Keep the rules simple and safe:

- Only location (`beta`) and uncertainty (`var/se`) are mapped linearly.
- `t` and `z` are recomputed after mapping (or combined explicitly in `reduce(method = "stouffer")`).
- Masks live in Space; `mask()` picks/derives the analysis mask (intersection/threshold/custom) and can pack storage.

---

## 40. Core Invariants

- All assays share identical dimensions and dimnames `[sample × subject × contrast]`.
- If `beta` exists, you must have `var` or `se` (`se` is `sqrt(var)`).
- If only `t/z` are provided, generic spatial mapping is not allowed; use `reduce(method = "stouffer")` or provide `beta/var`.
- Spaces are explicit; mapping between spaces is linear via Map.
- Provenance is structural: every verb adds a node to the plan; plans are serializable; `compute()` stores a full audit trail in `metadata$provenance`.

---

## 41. Minimal Adapter API (public, extensible)

```
register_adapter(
  name,
  detect(source) -> score in [0,1] or FALSE,
  open(source, options=list()) -> handle,
  probe(handle) -> list(
      assays = c("beta","var","t","df", ...),
      dims   = c(sample=I, subject=J, contrast=K),
      subjects, contrasts,
      space,                         # a Space object
      covariates = NULL, row_data = NULL,
      maps = list(),                 # optional Space-to-Space Map(s)
      metadata = list(units = list(...), schema_version = "0.1.0")
  ),
  read(handle, assays, block = NULL) -> named list of arrays [I × J × K],
  close(handle)
)

plan <- gds(source, format = "auto", prefer = NULL)   # uses registered adapters; returns a Plan
```

---

## 42. Factorial / Multilevel Designs

Keep the contrast axis flat; attach optional formal structure:

```
contrast_info(
  contrasts,                      # must match the axis
  design_spec = ~ A * B,          # formula or list
  factors     = list(A = c("a1","a2"), B = c("b1","b2")),
  C           = NULL,             # contrast matrix, if desired
  type        = c("main","main","interaction"),
  depends_on  = list("A:B" = c("A_main","B_main")),
  scale       = NULL,             # per-contrast units/scales
  notes       = NULL
)
```

Downstream modeling can consult this without expanding the API.

---

## 43. Uncertainty Through Mappings (recap)

For a linear map `y = M x` on the sample axis:

- `beta' = M beta`
- Independence: `Var(beta') = (M ⊙ M) Var(beta)`
- With covariance: `Var(beta') = M Σ M^T` via `cov_provider`
- Recompute `SE' = sqrt(Var(beta'))` and then `t' = beta'/SE'` if requested.

This is assay‑aware inside `map_to()` via the assay registry.

---

## 44. Masks: Storage vs Analysis (recap)

- Voxel Space stores `mask_idx` (active samples) and `storage = "packed"|"dense"`.
- `mask(policy = MaskPolicy(scope = "group", rule = "intersection" | "threshold", threshold = 0.9))` computes the analysis mask, repacking assays if beneficial.
- All mappings respect the active mask automatically.

---

## 45. End‑to‑End Examples

1) Voxel NIfTI → align to group → group mask → parcels → fixed‑effects

```
plan <- gds(Sys.glob("sub-*/cope1.nii.gz")) %>%
  subset(subject = paste0("sub-", sprintf("%02d", 1:30))) %>%
  derive(c("var","t")) %>%                              # ensures var present; t recomputable
  align(map = map_linear("subject_voxel","group_voxel",
                         by_subject = ortho_mats,
                         traits = list(orthogonal = TRUE))) %>%
  mask(policy = MaskPolicy(scope = "group", rule = "threshold", threshold = 0.95)) %>%
  map_to(target_space = "parcels", map = voxel_to_parcel_map,
         uncertainty = UncertaintyRule("independent")) %>%
  reduce(method = "fixed", weights = "1/var") %>%
  write_out(path = "results/cope1_parcels.h5", format = "h5")

g <- compute(plan)  # realized GDS in parcel space
```

2) ROI CSV with factorial contrasts → Stouffer on z‑scores

```
plan <- gds("roi_stats.csv") %>%
  subset(contrast = c("A_main","B_main","A:B")) %>%
  derive("z") %>%                                # if only t/df present; or pass-through if z exists
  reduce(method = "stouffer", weights = "equal") %>%
  write_out("group_roi.csv", format = "csv")

g <- compute(plan)
```

---

## 46. Tiny R Skeleton (exported signatures)

```
# nouns
new_gds(assays, space, subjects, contrasts, col_data=NULL, row_data=NULL, metadata=list())
space_voxel(dim, affine, mask_idx=NULL, storage=c("dense","packed"))
space_parcels(labels)
map_linear(from, to, operator, by_subject=NULL, traits=list(), uncertainty=UncertaintyRule("independent"))
UncertaintyRule <- function(mode=c("independent","cov_provider","kernel"), df_rule=c("satterthwaite","none"),
                            cov_provider=NULL, kernel=NULL)

# adapters
register_adapter(...)
gds <- function(source, format="auto", prefer=NULL)  # returns Plan

# verbs (all return Plan)
subset <- function(plan, sample=NULL, subject=NULL, contrast=NULL)
derive <- function(plan, what=c("var","se","t","z"))
align  <- function(plan, map)                         # subject-aware linear maps welcome
mask   <- function(plan, policy = MaskPolicy(scope="group", rule="intersection", threshold=0.9))
map_to <- function(plan, target_space, map, uncertainty=UncertaintyRule("independent"))
reduce <- function(plan, method=c("fixed","random","stouffer"), weights=c("1/var","equal","n_eff"))
write_out <- function(plan, path, format=c("h5","nifti","csv","parquet"))
compute <- function(plan, sink=c("memory","HDF5"), path=NULL, block=list(sample=1e5))

# extras
contrast_info <- function(...)
provenance    <- function(plan)
explain       <- function(plan)
save_plan     <- function(plan, file)
load_plan     <- function(file)
```

---

## 47. Packaging and Migration

- Package name: `gdsfmri` (or similar). Keep `fmrireg::group_data()` as a shim that calls `gds()`.
- Ship two adapters (NIfTI, CSV) and one HDF5 writer. Add HDF5/Zarr mini‑spec in `inst/schemas/`.
- Start with `UncertaintyRule("independent")`; the interface already supports covariance later.
- Keep documentation focused on the 5 nouns + 8 verbs; everything else hangs off these.

---

## 48. Integration With fmristore (Roles and Boundaries)

- Roles remain crisp and non-overlapping:
  - gds owns analysis semantics, assay awareness, lazy plans/DAG, provenance, mapping + uncertainty propagation, and meta-analysis.
  - fmristore is the HDF5 storage engine and an official Storage Adapter. It provides efficient, typed, chunked I/O for voxel/parcel/latent data.
- Outcome: gds becomes format-agnostic with a thin fmristore adapter; fmristore stays the canonical HDF5 hub rather than an analysis layer.

---

## 49. What fmristore Already Provides (Mapped to gds Nouns)

- Voxels (with mask):
  - `write_labeled_vec()` layout with NIfTI-like header under `/header/*`, a 3D mask under `/mask`, and 1-D data vectors per label under `/data/<label>` in mask order.
  - Maps cleanly to `space_voxel(mask_idx=...)` with per-contrast slabs; array-like access via `H5NeuroVol/H5NeuroVec` (`as_h5()` conversions).

- Parcels / clusters:
  - `H5ParcellatedScan`/`H5ParcellatedMultiScan` store cluster maps, read from `/cluster_map`, and handle cluster-time data efficiently.
  - Natural `space_parcels` with an accompanying voxel→parcel `LinearMap`; `cluster_metadata()` can populate `row_data`.

- Latent bases:
  - `LatentNeuroVec` stores basis/loadings: spatial basis under `/basis/basis_matrix`, temporal embeddings under `/scans/.../embedding`.
  - Textbook `space_basis` with a linear map back to voxels.

- HDF5 navigation and detection:
  - Centralized constants (`H5_PATHS`, `H5_DSETS`, `H5_ATTRS`) and `detect_h5_type()` provide robust probing and layout traversal.

Conclusion: fmristore already covers Spaces (voxels with mask, parcels, latent bases) with fast, chunked I/O—ideal backend material.

---

## 50. Two Integration Paths

Path A — Zero schema changes (ship a FmriStoreAdapter now)

- Keep current fmristore layouts; compose at the gds layer.
- Axes:
  - Subjects: vector of file paths (one per subject) + `subjects` vector.
  - Contrasts: if `LabeledVolumeSet`, use `/labels` as names and `/data/<label>` as vectors; if `H5NeuroVec`, treat dim 4 as contrasts; if parcellated, the sample axis is cluster IDs; if latent, component index.
  - Samples: for voxels, flatten `/mask` order; for parcels, cluster IDs; for latent, component index.
- Assays:
  - Expose what’s present: `beta`, `var`, `t`, `z`, `df`. May be split across files by assay or bundled as label sets.
- Space and maps:
  - Build `space_voxel`, `space_parcels`, or `space_basis` from HDF5 headers/paths; return `LinearMap` objects so gds can propagate uncertainty at compute time.
- Probe:
  - Use `fmristore::detect_h5_type()` within the adapter’s `probe()` to dispatch to the right reader.

Pros: no changes to fmristore; ingest immediately. Cons: multi-assay/multi-subject within a single HDF5 is less tidy; axes/provenance are implicit.

Path B — Tiny “GDS shim” group inside fmristore (recommended)

Add a small optional group so fmristore becomes an official GDS backend with perfect round-tripping. This coexists with existing fmristore groups.

```
/gds
  /axes
    /subjects            (string)                # length S
    /contrasts           (string)                # length K
  /space
    /type                (string: "voxel" | "parcels" | "latent")
    /mask                (uint8 [X,Y,Z])        # if voxel space
    /mask_idx            (int [V])              # optional cached linear indices (mask order)
    /affine              (double [4,4])         # from header or s/qform
    /cluster_map         (int [X,Y,Z])          # if parcels; or vector in mask-order
    /basis/basis_matrix  (double [k, V])        # if latent
  /assays
    /beta                (float [V,S,K])        # samples × subjects × contrasts
    /var                 (float [V,S,K])        # optional; or /se
    /df                  (int   [S,K])          # optional
    /meta/beta.json      (JSON; kind="beta", units, ...)
    /meta/var.json       (JSON; variance_of="beta", ...)
  /contrast_info         (JSON)                 # factorial structure, dependencies
  /provenance            (JSONL)                # append-only plan log
  /version               (string)               # e.g., "gds-h5/0.1"
```

Minimal fmristore helpers to add:

- `write_gds_assays(file, space, subjects, contrasts, assays = list(beta=…, var=…, df=…))`
- `read_gds_assays(file, which = c("beta","var","df"))`
- `contains_gds(file)` / `validate_gds(file)`
- (Optional) provenance append helper to `/gds/provenance` (JSONL), which `gds::compute()` can call.

Pros: single canonical file per study/cohort split, perfect interop, clean provenance, multi-assay support, zero ambiguity. Cons: small initial work.

Note: This layout aligns with the HDF5 mini-spec (§13) by namespacing under `/gds` to coexist with fmristore’s existing paths.

---

## 51. Minimal StorageAdapter API (What gds Expects)

Pseudocode interface in R that enables lazy, block-wise I/O while keeping semantics in gds:

```
StorageAdapter <- list(
  # capabilities
  probe   = function(path) -> TRUE/FALSE,
  open    = function(paths, mode="r", ...) -> handle,
  close   = function(handle),

  # axes & space
  subjects    = function(handle) -> character(S),
  contrasts   = function(handle) -> character(K),
  space       = function(handle) -> Space,             # voxel | parcels | latent
  linear_maps = function(handle) -> list(LinearMap),   # e.g., vox->parcels, latent->vox

  # assays (lazy blocks)
  list_assays = function(handle) -> c("beta","var","t","df", ...),
  read_assay  = function(handle, assay, i=NULL, s=NULL, k=NULL) -> array,
  write_assay = function(handle, assay, x, i=NULL, s=NULL, k=NULL),

  # metadata & provenance
  metadata    = function(handle) -> list(),
  append_prov = function(handle, entries) -> invisible(TRUE)
)
```

Why thin: fmristore already offers type detection (`detect_h5_type()`), space ingredients (`/mask`, header dims/affine, `/cluster_map`, latent basis), and chunked array classes. gds retains responsibility for statistically correct mapping and meta-analysis.

---

## 52. Practical Recipes

Path A (no schema change):

- Per-subject files with `write_labeled_vec()`:
  - `beta`: one file per subject, `labels =` contrasts, `/data/<contrast>` = mask-ordered β.
  - `var`: same layout (or parallel file set), same mask and labels.
- `FmriStoreAdapter` composes subjects×contrasts across files and exposes `assay("beta")`, `assay("var")`, `space_voxel(mask=…)`, and any parcel/latent maps for later `map_to()`.
- Parcellated first-level summaries: point to `H5Parcellated*`, expose `space_parcels` and an optional vox→parcel `LinearMap`.

Path B (tiny shim):

- Add the `/gds` group as in §50 and provide `write_gds_assays()`/`read_gds_assays()`; implement `contains_gds()`/`validate_gds()`.
- The gds adapter prefers `/gds/assays/*` when present and falls back to Path A heuristics otherwise.

---

## 53. Why This Avoids Overlap

- gds: 5 nouns, 8 verbs, lazy DAG, provenance, assay semantics, correct uncertainty propagation.
- fmristore: HDF5 engine with performant data structures and optional `/gds/*` conventions any ecosystem can write (FSL/AFNI/Python).

---

## 54. Nice-to-Have Tweaks in fmristore (Optional)

1) Cache sample order: write `/gds/space/mask_idx` (int) to avoid recomputing `which(mask)` across stacks.
2) Expose transforms to group space: per-subject `/gds/align/<subject_id>` transform (orthogonal/Procrustes/Wasserstein) for consensus alignment metadata.
3) Assay semantics metadata: JSON next to each dataset (e.g., `/gds/assays/var` with `{ "variance_of": "beta" }`) so gds can re-derive `t/z` post-mapping.
4) Provenance append: helper that appends JSONL to `/gds/provenance` (function, arguments, package version, timestamp) when `gds::compute()` runs.

---

## 55. Suggested Next Steps (Concrete)

- Week 1: implement read-only `FmriStoreAdapter` (Path A) using `detect_h5_type()`, `/labels`, `/mask`, `/cluster_map`, and `LatentNeuroVec` to populate Space + axes; expose `assay("beta")`, `assay("var")` from today’s files.
- Week 2: add the `/gds` writer/reader in fmristore (Path B), plus `contains_gds()` / `validate_gds()`.
- Week 3: prefer `/gds/assays/*` in the adapter when present, fall back to Path A otherwise.
- Docs: “Using fmristore as a GDS backend” vignette in both projects (voxel, parcel, latent examples).

---

## 56. Introduction & Design Overview (Executive Summary)

GDS (Group Data Set) is a small, composable standard for representing first‑level fMRI results so they flow cleanly into group‑level and meta‑analytic workflows—no matter the upstream toolchain or spatial representation (voxels, parcels, surfaces, latent bases). It is format‑agnostic, space‑aware, and statistically principled, giving you one consistent shape, a tiny set of verbs to describe intent, and a lazy Plan that compiles intentions into efficient, reproducible computation—without hard dependencies on Bioconductor stacks.

Purpose (why GDS exists):

- One lingua franca for first‑level outputs across SPM/FSL/AFNI/Python/R.
- One shape that works across spaces (dense voxels, ROI/cluster, surfaces, latent bases).
- One set of uncertainty rules: map effects/variances correctly; never average SE; re‑derive test statistics after transforms.
- One lazy pipeline with pushdown optimization and full provenance for large datasets.

Design in one picture: five nouns describe the world; eight verbs describe intent; `compute()` materializes. See §§37–46 for the minimal surface.

---

## 57. Statistics Families (including F)

GDS separates assays from statistical roles so new families can be added without expanding the surface:

- Location & Uncertainty: `beta` (effect), `var` (variance). `se` is derived from `var`.
- Test‑stat families (extensible): `t` with `df`, `z`; `F` with `df1`, `df2`; `chi2` with `df`; future families can register.
- `derive()` handles family → p/z (e.g., `pf(F, df1, df2)` → `p` → `z`).

Mapping policy:

- If `beta/var` exist → map linearly, then re‑derive `se/t/z/F` as needed.
- If only test statistics exist (e.g., just `F`) → do not linearly map; use explicit evidence combiners in `map_to(combine=…)` (within‑subject across samples) or `reduce(method=…)` (across subjects).

This “role + family” approach future‑proofs the API (e.g., Bayes factors, Wald tests) by registering family semantics and p/z derivations where applicable.

---

## 58. Storage & Dependencies (No Bioconductor Required)

GDS is not tied to `DelayedArray/HDF5Array`. Adapters talk directly to light backends:

- `hdf5r` for HDF5 datasets, `RNifti`/`neuroim2` for NIfTI, Arrow/Parquet or TileDB, and CSV/TSV for tabular; fmristore via a thin adapter (§§48–55).
- Block streaming is built into `compute()` so large arrays are processed in chunks irrespective of backend. No delayed‑array framework is required.
- Minimal Storage Adapter interface exposes axis labels, Space (with `mask_idx` and storage mode for voxels), and chunked `read_assay()`/`write_assay()`; GDS handles mapping, uncertainty, and statistics on top.

---

## 59. Masks & Alignment (Recap)

- Masks live in Space (voxel). Files can arrive with a storage mask; `mask()` derives the analysis mask by policy (e.g., group intersection ≥ 95%) and may repack to packed storage for speed (§§31, 44).
- Alignment is a subject‑aware `Map` used by `align()` to bring subjects into a consensus space (orthogonal/Procrustes, affine/deform, OT). After alignment, you can map to parcels/surfaces and `reduce()` across subjects (§§30–36).

---

## 60. Provenance & Reproducibility (Recap)

Every verb adds a node to the Plan’s provenance graph with operation name, canonicalized parameters, package version, timestamp, and input/operator hashes. Plans are serializable (JSON/YAML), explainable (`explain()`), and re‑runnable anywhere with the same adapters; `compute()` writes realized GDS and can append a provenance log to the output store (§§24–27).

---

## 61. Interoperation with fmristore (Recap)

Use fmristore as an official storage backend via a small adapter today (no schema change), or add an optional `/gds` group (beta/var arrays, axes, space, and provenance) to make it a first‑class GDS store. GDS stays in charge of semantics/computation; fmristore stays in charge of fast bytes (§§48–55).

---

## 62. Late‑Binding Alignment Transforms (Register After Stats)

Goal: you have first‑level stats (beta/var) in subject space. Later you compute/receive alignment transforms (orthogonal rotations, affine/deformable warps, or OT/Wasserstein plans) to a consensus group space. You want to attach these transforms to an existing GDS so future `align()`/`map_to()` can use them lazily and reproducibly.

Minimal, clean API (two verbs + one noun):

Noun — MapFamily (subject‑aware transform set):

```
MapFamily <- function(name, from_space, to_space,
                      type = c("orthogonal","linear","affine3d","deform3d","ot"),
                      by_subject,                 # named list: subject_id -> operator descriptor
                      traits = list(orthogonal = FALSE, mass_preserving = FALSE),
                      uncertainty = UncertaintyRule("independent")) {
  structure(list(name=name, from=from_space, to=to_space, type=match.arg(type),
                 by_subject=by_subject, traits=traits, uncertainty=uncertainty),
            class="gds_map_family")
}
```

Operator descriptors (serializable, no closures):

- orthogonal/linear: `list(kind="matrix", dim=c(Nt,Ns), storage="dense|csr", data=...)`
- affine3d: `list(kind="affine", mat=matrix(…,4,4), interp="trilinear|bspline")`
- deform3d: `list(kind="warp", path="warp_sub01.h5", layout="vx,vy,vz", interp="…")`
- ot: `list(kind="transport", storage="csr", rowsums=1, data=...)`

Verb — register_map(): attach one or more families to a GDS or Plan (immutably; returns modified object).

```
register_map <- function(x, family, overwrite = FALSE) {
  # 1) validate subjects and dimensions against x$space
  # 2) add to x$space_graph$maps[[family$name]]
  # 3) append a provenance node ("register_map", family$name, digest)
  # return modified x (GDS or Plan)
}
```

Verb — align(): use a registered family to transform each subject into the target space (lazy plan node).

```
align <- function(x, family, name = NULL) {
  fam <- if (inherits(family, "gds_map_family")) family else x$space_graph$maps[[family]]
  add_op(as_plan(x), list(op="align_to_group", family=fam$name))
}
```

Semantics & validation:

- Subjects: `names(by_subject)` must match `subjects(x)` (reordering allowed).
- Spaces: `family$from_space` equals current Space; `family$to_space` is the post‑align group Space.
- Operators: matrix dims `[N_target × N_source]` must match sample sizes; voxel affine/warp compatibility validated; OT rows sum to 1.
- Uncertainty: `family$uncertainty` controls propagation (`independent` → `(M ⊙ M) var`; `cov_provider`/`kernel` if available). Registration is lazy; assays are untouched until compute.

Persistence (optional, HDF5 via fmristore):

```
/gds/alignments/<family_name>/
  attrs:  type=..., from=..., to=...
  /subjects            (string[S])
  /traits              (JSON)
  /uncertainty         (JSON)
  # one of:
  /matrices/  /affine/  /warp/  /transport/   # subject-specific operators
```

Usage scenarios:

1) Orthogonal basis alignment computed later

```
G <- load_gds("firstlevel.h5")
R <- list("sub-01"=list(kind="matrix", dim=c(k,k), storage="dense", data=R1),
          "sub-02"=list(kind="matrix", dim=c(k,k), storage="dense", data=R2))
fam <- MapFamily("consensus_basis_v1", from_space="basis_k", to_space="basis_consensus",
                 type="orthogonal", by_subject=R, traits=list(orthogonal=TRUE))
G2 <- register_map(G, fam)
plan <- align(G2, family="consensus_basis_v1") %>% reduce(method="fixed")
res  <- compute(plan)
```

2) Deformable warps (FNIRT/ANTs) arriving after stats

```
G <- gds("beta_var_subject_space.h5")
warps <- read_warp_index("warplist.csv")
by_sub <- lapply(names(warps), function(s) list(kind="warp", path=warps[[s]], layout="vx,vy,vz", interp="trilinear"));
names(by_sub) <- names(warps)
fam <- MapFamily("mni_fnirt_2025", from_space="native_voxel", to_space="mni_voxel", type="deform3d", by_subject=by_sub)
plan <- register_map(G, fam) %>%
        align("mni_fnirt_2025") %>%
        mask(policy = MaskPolicy(scope="group", rule="threshold", threshold=0.95)) %>%
        map_to(target_space="parcels", map=voxel_to_parcel_map) %>%
        reduce(method="fixed")
res <- compute(plan)
```

Notes on stats & correctness:

- Never map `t/z/F` directly. Align/map `beta/var` (with the registered family), then re‑derive `se/t/z/F` after transforms.
- If only test statistics exist (no effect scale), use evidence combiners (`reduce(method="stouffer"|"fisher")`) rather than spatial maps.
- For `independent` propagation use `(M ⊙ M) var`; provide covariance via `cov_provider` when available.

Why this solves the timing issue: transforms can be created anytime and registered later; the object simply gains a named family in its `space_graph`. Alignment stays lazy and auditable; persistence is optional and backward‑compatible. No heavy dependencies are introduced; works with `hdf5r`, `RNifti/neuroim2`, or fmristore.
