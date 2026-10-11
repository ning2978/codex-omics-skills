---
name: codex-imc-neighborhood-analysis
description: Use when preparing CODEX, IMC, or mIHC cell-level protein tables as AnnData and running cellular neighborhood, RCN/CN, Squidpy, or SOAPy-style spatial niche analysis.
---

# CODEX / IMC Neighborhood Analysis

## Purpose

Use this skill for cell-level spatial proteomics where each row is a segmented cell with marker intensities, centroid coordinates, image/FOV id, and a phenotype/cluster label.

It distills the course `codex.py`, `pro.RCN.py`, and `SOAPy.pro.py` patterns into a reusable workflow.

## Environment

Use the project's validated Python environment. Pair with
`omics-analysis-runtime-preflight` to verify the actual interpreter and required
imports. Invoke any project runner with that interpreter; no environment wrapper
or dataset-specific runner is bundled here.

The supported default stack is AnnData/Scanpy + scimap + Squidpy. SOAPy is optional and should use a separate environment if a concrete SOAPy script is being run, because `SOAPy-st` pins older Scanpy/AnnData/scikit-image/scikit-learn versions that conflict with the verified `codex-imc` stack.

## Input Contract

Require a table or h5ad with:

- marker intensity columns for `adata.X`;
- centroid columns, e.g. `Centroid X`, `Centroid Y`;
- image/FOV column, e.g. `IMAGE` or `fov`;
- phenotype/cluster column, e.g. `cluster`, `cell_type`, or `cluster_ID_*`;
- optional existing RCN/CN assignment table.

AnnData contract:

- `adata.X`: numeric marker matrix;
- `adata.obs[cluster_key]`: categorical phenotype label;
- `adata.obs[fov_key]`: FOV/image boundary;
- `adata.obsm["spatial"]`: numeric `[x, y]`;
- optional `adata.obsm["X_umap"]` from upstream clustering.

## Workflow

1. Convert source CSV/export to AnnData with marker columns in `X`.
2. Preserve centroid and FOV columns in `obs`; write `obsm["spatial"]`.
3. If an upstream neighbor label table exists, align by stable cell id, not row position unless that is the documented key.
4. Build spatial graph per FOV:
   - Squidpy: `sq.gr.spatial_neighbors(..., coord_type="generic", radius=<cutoff>)`;
   - SOAPy: `sp.pp.make_network(..., sample_key=fov_key, method="radius", cutoff=<cutoff>)`.
5. Run neighborhood enrichment, CN/RCN clustering, or niche clustering from local phenotype composition.
6. Export h5ad, spatial PNGs, composition tables, and parameters.

## Guardrails

- Never connect cells across different FOVs/images.
- Keep marker intensity preprocessing separate from spatial neighborhood analysis.
- Record physical units for radius. Many CODEX/IMC exports use pixels, not microns.
- If row alignment is uncertain, stop and require a cell-id join key.

## Output Contract

- `codex_or_imc_prepared.h5ad`;
- `neighborhood_assignment.csv` or `cn_assignment.csv`;
- `neighborhood_composition_by_phenotype.csv`;
- `neighborhood_fraction_by_fov.csv`;
- spatial plots by FOV and neighborhood/CN label;
- `run_params.json` with marker columns, coordinate columns, FOV key, cluster key, method, radius/KNN, and package versions when available.

## Validation

- `adata.n_obs` equals source cell table rows after any documented filtering.
- `obsm["spatial"]` is numeric, shape `(n_obs, 2)`, and nonzero.
- FOV labels are present for all cells.
- Spatial graph has no cross-FOV edges.
- Fractions sum to about 1 within each FOV.
- PNGs are nonblank and use the same coordinate orientation as the source images.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
