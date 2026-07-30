# Visium Spot h5ad Data Contract

Use this contract when deciding whether a prepared Visium spot `.h5ad` is ready for downstream deconvolution, visualization, or neighborhood analysis.

## Required Structure

- `adata.X`: normalized/log data or raw counts, depending on project convention. State the convention explicitly.
- `adata.layers["raw_counts"]`: pre-normalization integer-like counts when downstream model fitting needs counts.
- `adata.obs_names`: spot barcodes, unique within the object.
- `adata.var_names`: gene identifiers or symbols, unique enough for the downstream workflow.
- `adata.obsm["spatial"]`: full-resolution spot coordinates in `[x, y]` order.
- `adata.uns["spatial"][library_id]`: 10x-style nested spatial metadata with scalefactors and available images.
- `adata.obs`: sample metadata, QC metrics, and any condition/timepoint/batch labels needed downstream.

## Recommended QC Fields

- `total_counts`
- `n_genes_by_counts`
- `pct_counts_mt`
- sample/library identifier
- tissue/region annotation if available
- selected QC-threshold metadata in `uns`

## Readiness Checks

- Reopen the file from disk and inspect shape, layers, `obsm`, and `uns`.
- Confirm spatial coordinates render over the tissue image with the expected orientation.
- Confirm raw-count layer is not HVG-restricted if the file is intended for cell2location or RCTD.
- Confirm sample metadata values are sourced from the current project mapping, not inferred from obsolete filenames.
- Confirm the output directory and naming convention are stable before handing off to model-fitting skills.
