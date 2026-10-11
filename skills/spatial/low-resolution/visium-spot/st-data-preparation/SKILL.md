---
name: st-data-preparation
description: Use this skill when preparing raw 10x Visium spatial transcriptomics inputs or Seurat/RCTD spatial objects into standardized `.h5ad` files, including unified R/Seurat-to-Python conversion with sc/ST mode detection, especially when checking barcode matching, spatial metadata structure, raw-count preservation, QC assumptions, Matrix Market exports, or a project's canonical ST preparation workflow.
---

# Skill: st-data-preparation

## Purpose

Use this skill when the task is to prepare ST inputs from raw 10x Visium outputs, or from already prepared Seurat spatial objects, into standardized `.h5ad` files that remain suitable for downstream deconvolution, plotting, neighborhood analysis, or integration.

This skill is for raw 10x Visium intake, Seurat/RCTD ST object conversion, QC, barcode matching, spatial metadata preservation, and validation of the generated `.h5ad` files.

Route first through `spatial/shared/workflow-st-analysis`; use this skill for concrete Visium spot implementation details.

## Use This Skill When

- the task starts from raw 10x Visium outputs
- the task starts from per-sample Seurat spatial RDS objects and needs Python/Scanpy `.h5ad`
- you need a transparent `barcodes.tsv.gz` / `features.tsv.gz` / `matrix.mtx.gz` intermediate before `.h5ad`
- you need to validate barcode matching, spatial metadata, QC thresholds, or saved `.h5ad` structure
- you want the resulting files to be ready for downstream cell2location use

Do not use this skill when:

- the user already has valid standardized ST `.h5ad` files and wants model fitting
- the task is mainly about cell2location parameter choice or result interpretation
- the task is high-resolution Visium HD, Xenium, CODEX, or IMC

## Quick Start

When this skill is used:

1. Identify the raw input root, output root, sample metadata source, and intended `.h5ad` naming convention.
2. Confirm the Python/R environment from the current project before running conversion or QC scripts.
3. Preserve raw integer counts and spatial metadata before normalization or feature selection.
4. Validate barcode matching, coordinate orientation, `obsm['spatial']`, and `uns['spatial']` before declaring files ready.
5. Record the chosen QC thresholds and output contract in the project work log if one exists.

## Core Data Contract

Generated ST `.h5ad` files should preserve this downstream contract:

- `adata.layers['raw_counts']` as the pre-normalization integer count matrix for downstream cell2location use
- `adata.X` as normalized and log-transformed data only if that matches the chosen repository convention
- `adata.obsm['spatial']` as full-resolution spot coordinates
- `adata.uns['spatial'][GSM_ID]` with nested `scalefactors` and available `images`
- project-level sample annotations such as `adata.obs['sample']` and, when available, fields like `donor`, `condition`, `timepoint`, `slide`, or `batch`
- QC metrics such as `total_counts`, `n_genes_by_counts`, and `pct_counts_mt`

Important:

- For 10x Visium projects, keep implementation aligned with the actual saved data structure, not only with generic Scanpy examples.
- Most tools expect `adata.uns['spatial'][library_id]`, not a flat `adata.uns['spatial']` payload.

## Working Procedure

For raw 10x Visium intake, follow this order unless the user asks for a narrower check:

1. Confirm the task really starts from raw 10x Visium data.
2. Verify the environment and use the `scanpy` wrapper when running checks or scripts.
3. Use the project's sample-mapping source of truth if one exists.
4. Read raw counts from `filtered_feature_bc_matrix.h5`.
5. Read spatial files from each sample `spatial/` directory.
6. Match expression barcodes to spatial barcodes before saving coordinates.
7. Compute QC metrics and apply project QC thresholds unless the user requests alternatives.
8. Save `.h5ad` files into the project-standard output directory.
9. Validate the saved structure before declaring the output ready for cell2location.

For Seurat/RCTD spatial objects that already exist, prefer the transparent two-stage route:

1. Read per-sample Seurat RDS files.
2. Export each sample to a 10x-like Matrix Market directory with:
   `matrix.mtx.gz`, `normalized_data.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`, `metadata.csv`, `source_manifest.json`, and `spatial/`.
3. Build `.h5ad` from that intermediate directory with the scanpy Python environment.
4. Validate shape, spot count, integer-like counts, `obsm['spatial']`, `uns['spatial']`, tissue images, and any annotation joins.

For generic Seurat RDS conversion, prefer a converter that exposes these mode rules:

- `--mode auto` detects sc vs ST from Seurat image slots and writes mode-scoped output.
- ST outputs go to a spatial output folder; sc outputs go to a separate single-cell output folder.
- Use `--mode st` only when the object should be treated as spatial despite auto-detection ambiguity.
- In ST mode, validate `adata.obsm['spatial']`, `adata.uns['spatial'][library_id]`, images/scalefactors, and the spatial coordinate column order.

## QC Defaults

Start with conservative exploratory thresholds, not universal constants. Common Visium spot checks include:

- `min_genes`
- `max_counts`
- `max_pct_mt`
- tissue/edge spot status where available

If the tissue type, chemistry, spot density, or study design changes, call out that thresholds need review.

## Related Skills

- Use `spatial/shared/workflow-deconvolution-run` + algorithm skill after standardized ST `.h5ad` files are ready.
- Use `spatial/low-resolution/visium-spot/algorithm-cell2location` for cell2location deconvolution.
- Use `spatial/low-resolution/visium-spot/algorithm-rctd` for RCTD deconvolution.

## What to Read Next

Read these references only when needed:

- For the detailed workflow, validation checklist, and common failure modes:
  `references/workflow.md`
- For the expected output structure:
  `references/data-contract.md`
- For the verified Seurat/RCTD RDS -> 10x-like Matrix Market -> `.h5ad` route:
  `references/seurat-rctd-to-h5ad.md`
  - Must read the Seurat `GetTissueCoordinates()` coordinate-orientation guidance before writing or reviewing any Seurat -> h5ad converter.

## Boundaries

- Prefer current project docs and current outputs over old notebooks and archived scripts.
- If you mention implementation history, say clearly whether a file is a reference implementation or the recommended path.
- Do not present archived scripts as the default execution entrypoint unless the user explicitly asks for legacy reproduction.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
