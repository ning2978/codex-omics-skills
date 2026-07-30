---
name: visiumhd-qc-analysis
description: Use when preparing, QCing, threshold-tuning, or exporting 10x Visium HD binned spatial transcriptomics data, especially Matrix Market plus spatial parquet outputs that need per-bin QC tables, spatial threshold comparison plots, or QC-filtered `.h5ad` files.
---

# Visium HD QC Analysis

## Purpose

Use this skill for 10x Visium HD binned outputs, especially `outs/binned_outputs/square_008um/` directories containing `filtered_feature_bc_matrix/` and `spatial/tissue_positions.parquet`.

This skill complements `spatial/low-resolution/visium-spot/st-data-preparation`: use that skill for spot-level Visium contracts, and this one for Visium HD-specific QC threshold comparison and h5ad export.

## Workflow

1. Confirm the Visium HD input structure:
   - `filtered_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}`
   - `spatial/{tissue_positions.parquet,scalefactors_json.json,tissue_hires_image.png,tissue_lowres_image.png}`
2. Confirm sample metadata and output naming.
3. Compute or reuse per-bin QC metrics:
   - `total_counts`
   - `n_genes_by_counts`
   - `mt_counts`
   - `pct_counts_mt`
   - spatial columns from `tissue_positions.parquet`
4. Before final filtering, compare candidate `min_genes` thresholds with spatial plots.
5. Export QC-filtered `.h5ad` only after the user chooses the threshold.
6. Validate saved h5ad files before claiming they are ready.

## Recommended Defaults

For 8um Visium HD bins, start with threshold exploration instead of guessing.
Use the target repository's adapter when one exists; the command below shows
the required interface rather than a bundled script:

```bash
python /path/to/project/compare_visiumhd_qc_thresholds.py \
  --thresholds 20 50 80 100 \
  --max-pct-mt 20 \
  --max-counts 50000
```

Formal export after threshold selection:

```bash
python /path/to/project/export_visiumhd_qced_h5ad.py \
  --min-genes 100 \
  --max-pct-mt 20 \
  --max-counts 50000
```

Use the current project's validated Python environment. Use conda/pip environments only.

## Output Contract

Threshold comparison should produce:

- summary CSV with one row per sample per threshold
- one overview count plot
- one multi-panel spatial pass-bin plot per sample

QC-filtered h5ad should preserve:

- raw counts in `adata.X`
- raw counts copy in `adata.layers["raw_counts"]`
- full-resolution spatial coordinates in `adata.obsm["spatial"]` as `[x, y]`
- sample-keyed spatial metadata in `adata.uns["spatial"][sample_id]`
- QC parameters in `adata.uns["*_qc"]` or an equivalent explicit key

## Validation Checklist

After threshold comparison:

- summary rows equal `n_samples * n_thresholds`
- totals match the expected filtering logic
- all sample spatial plots exist and use `x = pxl_col_in_fullres`, `y = pxl_row_in_fullres`, with y-axis inverted

After h5ad export:

- each h5ad can be opened with `anndata.read_h5ad(..., backed="r")`
- shapes match the summary CSV
- `layers["raw_counts"]`, `obsm["spatial"]`, `uns["spatial"]`, and QC metadata exist
- total retained bins match the selected threshold comparison total

## Boundaries

- Do not reuse a published or old project threshold without checking spatial pass-bin plots.
- Do not treat Visium spot QC thresholds as Visium HD defaults.
- Keep binned output resolution explicit, for example `square_008um`.
