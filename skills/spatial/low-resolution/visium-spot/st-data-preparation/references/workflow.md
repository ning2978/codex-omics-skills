# Visium Spot Preparation Workflow

Use this reference when implementing or reviewing a Visium spot data-preparation script.

## Inputs

- Raw 10x Visium sample folders with expression matrix and `spatial/` assets, or
- Per-sample Seurat spatial objects that must be converted to Python/Scanpy-compatible `.h5ad`.

Before writing code, identify:

- sample IDs and display names
- metadata columns needed downstream
- raw input root
- standardized h5ad output root
- intended QC thresholds
- whether images and scalefactors must be preserved

## Raw 10x Route

1. Read the filtered feature-barcode matrix.
2. Read `spatial/tissue_positions*`, scalefactors, and available images.
3. Match expression barcodes to spatial barcodes before saving coordinates.
4. Store full-resolution coordinates in `obsm["spatial"]` as `[x, y]`.
5. Store 10x-style spatial metadata under `uns["spatial"][library_id]`.
6. Preserve raw integer counts in `layers["raw_counts"]`.
7. Add sample metadata before concatenating or saving.
8. Compute QC metrics and apply the selected thresholds.
9. Save one `.h5ad` per sample or one merged object only if downstream code expects that shape.
10. Validate the saved objects by reopening them from disk.

## Seurat Spatial Route

Prefer a transparent intermediate when converting Seurat spatial objects:

- `matrix.mtx.gz`
- `normalized_data.mtx.gz` when available
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `metadata.csv`
- `spatial/`
- `source_manifest.json`

Then build `.h5ad` from that intermediate with Python. This makes barcode order, feature order, spatial coordinates, and raw counts auditable.

## Validation Checklist

- `n_obs` equals the expected spot count.
- expression matrix rows align with `obs_names`.
- `layers["raw_counts"]` is integer-like and full-gene when needed for deconvolution.
- `obsm["spatial"]` exists and has two columns.
- `uns["spatial"][library_id]` includes scalefactors and any images required by plotting code.
- sample metadata columns are present and not silently inferred from stale filenames.
- QC metrics and thresholds are recorded.
- saved `.h5ad` files reopen with `anndata.read_h5ad(..., backed="r")`.

## Common Failure Modes

- Seurat coordinate columns named `x`/`y` may not match pixel-axis semantics; validate orientation against tissue images.
- Barcode suffixes such as `-1` can be dropped or duplicated during R/Python conversion.
- HVG-only objects are not safe deconvolution references unless a full raw-count layer is preserved.
- Concatenating samples before validating per-sample `uns["spatial"]` can hide broken image metadata.
