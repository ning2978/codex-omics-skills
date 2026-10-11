---
name: visiumhd-spatial-domain-neighborhood
description: Use when analyzing 10x Visium HD or other sub-spot/binned spatial data for spatial domains, BANKSY-style spatial clustering, Squidpy neighborhood enrichment, differential neighborhood enrichment, ligand-receptor communication, or niche profiles from cell-type/bin labels.
---

# Visium HD Spatial Domain + Neighborhood

## Purpose

Use this skill after Visium HD bins have been loaded into an AnnData object with `obsm["spatial"]`. It covers the reusable patterns found in the "single-cell spatial symphony" HD scripts:

- basic HD Scanpy clustering and spatial plotting;
- BANKSY-like spatial domain clustering;
- cell-type/bin neighborhood enrichment and between-sample differential enrichment;
- Squidpy ligand-receptor analysis constrained by a spatial graph;
- simple neighborhood-composition niche profiles.

For raw binned Matrix Market plus `tissue_positions.parquet` QC/export, use `visiumhd-qc-analysis` first. For mapping existing CellCharter labels onto H&E, use `visiumhd-cellcharter-visualization`.

When assigning or reviewing biological domain/niche names, read [the joint evidence review](../../../../single-cell/python/sc-marker-annotation/references/annotation-evidence-review.md). Combine differential ranking, sourced canonical/competing programs, biological-unit support and appropriate QC with anatomy/composition. Neighborhood or communication computation on supplied labels does not itself request new annotation.

## Input Contract

Require:

- one `.h5ad` per sample or a combined `.h5ad`;
- raw or normalized counts with `var_names` made unique;
- `obsm["spatial"]` as numeric `[x, y]`, preferably full-resolution pixel coordinates for 10x HD;
- `obs[cluster_key]` for neighborhood/communication work, usually cell type, domain, or niche label;
- `obs[sample_key]` or `obs[library_key]` before combining samples.

Before analysis, assert coordinates are numeric and not strings. If using 10x H&E plotting, confirm `uns["spatial"][sample_id]["scalefactors"]` exists and do not hard-code a library id.

Declare whether observations are bins or documented segmented-cell profiles. Bin size alone does not establish single-cell purity. For biological review, locate explicit unintegrated gene expression, sample/unit design, applicable QC and available histology/composition; missing checks remain limitations.

## Workflow Router

1. **Domain discovery**: run Scanpy PCA/neighbors/Leiden or BANKSY. Use this when the output is a spatial domain label such as `molniche`, `banksy_domain`, or `leiden`.
2. **Neighborhood enrichment**: build a Squidpy graph with `coord_type="generic"` and either `radius` or `n_neighs`, then run `sq.gr.nhood_enrichment`.
3. **Differential enrichment**: compute each sample/condition separately or use a method that preserves library boundaries. Never let cross-sample edges enter a spatial graph.
4. **Ligand-receptor**: run `sq.gr.ligrec` only after cell-type labels are stable. Filter source/target groups explicitly for the biological question.
5. **Neighborhood niche profile**: convert each bin/cell's local label composition into a feature matrix, cluster that matrix, and interpret clusters as microenvironment states, not new cell types.

For biological naming, rank markers from all available expression-filtered genes in original unintegrated expression, preserving complete statistics plus a candidate view. Do not test scaled/HVG-only, BANKSY-augmented or neighborhood-composition features as gene-expression markers. Check measured canonical/competing programs beyond DE cutoffs, anatomy/composition, unit recurrence and QC. Retain mixed/interface domain hypotheses for coherent multi-lineage bins; domain/niche labels stay separate from cell identities and states.

## Parameter Defaults

- HD bin QC: start from `min_genes` chosen by `visiumhd-qc-analysis`; do not reuse course defaults blindly.
- Squidpy radius: explore a small grid such as `50, 100, 150, 200` in the same coordinate unit as `obsm["spatial"]`.
- KNN graph: use `n_neighs=20` or `30` for dense HD bins when a physical radius is not interpretable.
- BANKSY: treat `lambda=0.2`, `k_geom=15`, `pca_dims=20`, `resolution=0.5` as exploratory starting points, not fixed project defaults.

## Output Contract

Write outputs under an explicit run directory:

- `run_params.json` or equivalent parameter record;
- domain/niche labels in a new `.h5ad` copy, not by overwriting the source object;
- neighborhood enrichment matrices as CSV plus heatmap PNG;
- differential enrichment table with condition labels and sample/library scope;
- ligand-receptor table and filtered plot for named source/target groups;
- spatial plots with y-axis/image orientation checked.
- for biological naming/review, an evidence table with program sources/availability, positive/competing support, biological-unit support, QC, spatial/composition observations, alternatives and decision/status.

## Validation

- `obsm["spatial"]` exists, numeric, shape `(n_obs, 2)`, and is not all zeros.
- Spatial graph has no cross-sample edges when multiple libraries are present.
- Enrichment matrices contain the same label categories as `obs[cluster_key]`.
- Group-level differential results use declared independent biological units; ROIs/libraries from one specimen retain their technical or nested structure rather than becoming extra independent replicates.
- Biological names reflect bin/segmented-cell semantics and joint evidence or an explicit limitation; a mixed bin is not automatically contamination or a pure cell identity.
- PNGs are nonblank and labels are readable.

## Source Mapping

This skill distills these Baidu Netdisk script patterns from `/01_生信培训/单细胞交响乐/2025单细胞交响乐空间转录组`:

- lessons 11-12: `scanpy.spatial.HD.py`, `tissue_positions.py`, `banksy.HD.py`;
- lessons 14-18: `scniche.hd.py`, `coor.hd.py`, `coor.hd.diff.py`, `squidpy_cci.py`, `spatial.nmf.py`.

The course scripts are reference patterns only; many contain hard-coded paths, example labels, or plotting assumptions that must be replaced with project-specific parameters.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
