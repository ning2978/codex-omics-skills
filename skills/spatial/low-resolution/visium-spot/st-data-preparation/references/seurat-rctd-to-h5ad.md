# Seurat Spatial or RCTD Objects to h5ad

Use this reference when ST data already exists as per-sample Seurat spatial/RCTD objects and needs Python-compatible `.h5ad` files.

## Preferred Conversion Shape

Use a two-stage export:

1. R exports each sample to a transparent 10x-like directory.
2. Python reads the directory and writes `.h5ad`.

The intermediate directory should contain:

- `matrix.mtx.gz` for raw counts
- `normalized_data.mtx.gz` when normalized expression must be preserved
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `metadata.csv`
- `spatial/tissue_positions.csv` or equivalent
- `spatial/scalefactors_json.json` and images when available
- `source_manifest.json`

## Coordinate Guidance

Check Seurat spatial coordinates carefully. Some Seurat/Visium accessors return columns whose names do not match downstream pixel-axis expectations.

Before saving final h5ad:

- render a small coordinate overlay against the image
- verify x increases left-to-right
- verify y increases top-to-bottom in pixel coordinates
- confirm `obsm["spatial"]` is `[x, y]`
- confirm the library ID used in `uns["spatial"][library_id]`

## R Export Rules

- Preserve barcode order exactly as used in the expression matrix.
- Preserve feature order exactly as used in the expression matrix.
- Export raw counts from the count slot/layer, not scaled data.
- Export normalized data separately when needed.
- Write metadata with one row per barcode.
- Write a manifest with source object path, assay/layer names, dimensions, and timestamp.

## Python Import Rules

- Read the sparse matrix without densifying.
- Join metadata by barcode and assert no missing or duplicated rows.
- Create `layers["raw_counts"]` from the raw-count matrix.
- Put normalized/log data in `X` only when that is the intended project convention.
- Store coordinates in `obsm["spatial"]`.
- Store scalefactors/images in `uns["spatial"][library_id]` when available.

## Validation

- Matrix shape equals metadata and feature table dimensions.
- Raw counts are non-negative and integer-like.
- Spot coordinates are finite and in expected pixel range.
- A saved `.h5ad` can be reopened.
- A quick spatial plot or overlay is nonblank and correctly oriented.
