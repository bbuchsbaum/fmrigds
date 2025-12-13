# GDS HDF5 Layout (`gds-h5`)

```
/gds
  /axes
    /subjects    (string dataset, len = #subjects)
    /contrasts   (string dataset, len = #contrasts)
  /space
    /type        (string attribute: "voxel" | "parcels" | "basis")
    /voxel/mask_idx (int dataset)        # optional, mask-ordered indices
    /voxel/dim      (int dataset, len=3)
    /voxel/affine   (double dataset, 4×4)
    /parcels/labels (string dataset)
  /assays
    /beta       (float dataset [sample × subject × contrast])
    /var        (... optional)
    /se         (... optional)
    /t, /z, /df (optional)
  /metadata
    /json       (string dataset, JSON blob)
  /provenance
    /log        (string dataset or JSONL)
  /alignments
    /<family>/serialized (string dataset; ASCII-serialised MapFamily)
```

- All assay datasets share identical dimensions and chunking along the sample axis.
- Logical order is sample × subject × contrast; store as attribute if physical order differs.
- Additional assays may be added under `/gds/assays/<name>` with matching dimensions.
- Alignment families are stored as R-serialised (`serialize(..., ascii = TRUE)`) objects; adapters may expose them via `plan$meta$map_families`.
- Writers may also emit tabular/CSV (`write_out(..., format = "csv")`), Parquet (optional `arrow`), or NIfTI volumes (requires `RNifti`) when exporting realised results. These formats are derived via `compute()` after data normalises into the `gds` structure.
- Streaming: `compute(block = list(sample = ...))` forwards block indices to adapters; HDF5 adapters honour `[sample × subject × contrast]` slicing, enabling chunked processing downstream.
