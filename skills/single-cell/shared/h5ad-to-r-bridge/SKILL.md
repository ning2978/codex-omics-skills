---
name: h5ad-to-r-bridge
description: Use when converting AnnData h5ad files into R, Seurat, spacexr, RCTD, or Matrix Market inputs while preserving sparse raw counts, barcode order, metadata, labels, reductions, or spatial coordinates.
---

# H5AD To R Bridge

## Purpose

Use this skill for Python/AnnData `.h5ad` -> R workflows. The default route is a transparent sparse Matrix Market export, not dense CSV and not an opaque direct RDS conversion.

Use with `algorithm-rctd` when the exported data will feed `spacexr::Reference`, `spacexr::SpatialRNA`, or a pilot/full RCTD run.

## Core Contract

- Never export expression with `adata.X.todense()` for large data.
- RCTD counts must come from an explicit raw-count layer, not normalized/log `adata.X`.
- Intermediate counts are always genes x observations in `counts.mtx.gz`.
- `barcodes.tsv.gz` order must match matrix columns.
- `features.tsv.gz` order must match matrix rows.
- `metadata.csv` must contain `barcode`.
- ST exports must write `spatial/tissue_positions.csv` with `x = pxl_col_in_fullres` and `y = pxl_row_in_fullres`.
- If `obs["pxl_col_in_fullres"]` and `obs["pxl_row_in_fullres"]` exist, cross-check them against `obsm["spatial"]`.

## Scripts

All bundled scripts live under this skill's `scripts/` directory.

### Export h5ad to sparse R input

```bash
python \
  <skill_dir>/scripts/export_h5ad_sparse_for_r.py \
  --input-h5ad path/to/input.h5ad \
  --output-dir path/to/export_dir \
  --mode st \
  --counts-layer raw_counts
```

For single-cell references, use `--mode sc --counts-layer counts --label-key <celltype_column>`.

### Build Seurat object

```bash
Rscript <skill_dir>/scripts/read_sparse_export_to_seurat.R \
  --input-dir path/to/export_dir \
  --output-rds path/to/object.rds
```

Seurat spatial images are optional and not part of the default RCTD path.

### Build RCTD objects

```bash
Rscript <skill_dir>/scripts/read_sparse_export_to_rctd.R \
  --input-dir path/to/reference_export \
  --type reference \
  --label-key cell2location_label \
  --output-rds path/to/reference.rds

Rscript <skill_dir>/scripts/read_sparse_export_to_rctd.R \
  --input-dir path/to/spatial_export \
  --type spatial \
  --output-rds path/to/spatial_rna.rds
```

### Run an RCTD pilot

```bash
Rscript <skill_dir>/scripts/run_rctd_from_sparse_exports.R \
  --reference-dir path/to/reference_export \
  --spatial-dir path/to/spatial_export \
  --output-dir path/to/rctd_pilot \
  --label-key cell2location_label \
  --max-cores 4 \
  --test-mode TRUE
```

## Project Defaults To Define

Before exporting, define:

- ST h5ad path or directory
- ST counts layer
- reference h5ad path
- reference counts layer
- reference label key
- export root

## Validation Checklist

- Re-open the Matrix Market export in R and confirm matrix dimensions.
- Confirm all counts are nonnegative and integer-like.
- Confirm metadata row order matches barcodes.
- For ST, confirm coordinate rownames match matrix colnames.
- For RCTD reference, sanitize invalid cell-type names and write a reversible mapping CSV.
- For pilot RCTD, require an RDS object, weights CSV if produced, and a run manifest.
