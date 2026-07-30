---
name: workflow-st-analysis
description: Use when starting, auditing, or staging a spatial transcriptomics project and the current data stage, stable input-output contracts, environment preflight, or handoff to preparation, deconvolution, visualization, or platform-specific analysis must be determined.
---

# Spatial Transcriptomics Workflow

## Role

Route spatial work after `skills-taxonomy` has identified the modality and resolution. Use this skill when the project stage is unclear or spans multiple stages. When a named algorithm already has ready inputs, use this skill only for a brief contract/handoff check; the algorithm skill owns method execution and completion.

Keep project paths and sample conventions in repository-local adapters; keep this workflow cross-project.

## Triage

1. Identify the platform and resolution: Visium spot, Visium HD/bin, Xenium/cell, or spatial proteomics.
2. Locate the current authoritative object: raw vendor output, QCed h5ad, reference h5ad, fitted model, or result tables.
3. Record stable sample, batch, condition, coordinate, count-layer, and annotation keys.
4. Confirm the expected output boundary and execution environment.
5. Hand off to the narrowest implementation and algorithm skill.

## Stage Router

| Current state | Next action | Handoff |
|---|---|---|
| Raw 10x Visium spot output or Seurat spatial object | Prepare and validate standardized h5ad | `st-data-preparation` |
| Raw Visium HD binned output | Tune bin QC and export h5ad | `visiumhd-qc-analysis` |
| Xenium XOA 2.0+ segmentation parquet | Inspect cells and polygon boundaries | `xenium-segmentation-new` |
| Per-ROI Xenium h5ad | Filter cell and segmentation QC | `xenium-roi-h5ad-qc-segmentation` |
| QCed segmented-cell Xenium expression subset | Recluster retained cells from counts | Pair with `sc-subcluster-analysis` |
| QCed spot-level ST plus reference | Launch deconvolution | `workflow-deconvolution-run` plus `algorithm-cell2location` or `algorithm-rctd` |
| Existing deconvolution object | Validate and summarize results | Algorithm-specific visualization skill |
| Existing spatial labels or embedding | Run domain, neighborhood, communication, or annotation analysis | Matching platform or `algorithm-*` skill |

## Stable Contracts

Require the following before a handoff:

- raw integer counts preserved in an explicit matrix or layer;
- coordinates stored in `obsm["spatial"]` or a documented equivalent with units and orientation;
- sample/library boundaries that prevent cross-sample spatial edges;
- stable sample, condition, batch, and annotation keys;
- input, output, log, and parameter-record paths;
- a machine-readable run manifest for long or multi-sample analyses.

For scRNA-seq or segmented-cell Xenium expression reclustering after any cell
filter, require retained `layers["counts"]` and recompute HVG, PCA, the
selected batch correction, neighbors, UMAP, and Leiden. Existing `X_umap`
coordinates may be retained only as explicitly named historical
projection/overlay fields; their presence never skips the new UMAP.

Do not apply this expression-UMAP gate to scimap RCN/CN, target-cell NMF
neighborhood clustering, or CellCharter niche clustering. Those methods
cluster local phenotype composition or aggregated spatial features rather
than segmented-cell gene-expression profiles.

For Visium-like deconvolution, keep this boundary explicit:

`raw spatial input -> QCed spatial h5ad -> reference h5ad -> deconvolution result -> summary tables and figures`

## Environment Preflight

- Prefer a repository-provided wrapper or locked environment.
- Verify imports and accelerator availability before a heavy run.
- Apply `environment-installation-policy` when an environment must be created or changed.
- Keep notebook exploration separate from the reproducible CLI or script used for the final run.
- Do not choose scVI for general Xenium expression refinement. Use it only
  when the user explicitly requests CellCharter expression latent mode and
  record that choice in the run manifest.

## Completion Gate

This gate proves that a stage produced a valid handoff artifact. It does not replace method-specific completion checks in an algorithm skill.

Do not declare a stage complete until:

- the expected artifact can be reopened;
- required AnnData keys and shapes match the recorded contract;
- an expression recluster contains counts, current two-dimensional `X_umap`,
  neighbor graph, cluster key, and true
  `recomputed_hvg/pca/neighbors/umap/leiden` provenance flags;
- multi-sample spatial graphs have zero cross-sample edges when applicable;
- parameters, package versions, and sample scope are recoverable;
- logs or summaries contain a clear successful completion marker.

## Boundaries

- Do not duplicate algorithm parameters here.
- Do not encode one tissue, cohort, or repository layout.
- Do not keep routing inside this skill after a specialized skill has taken ownership.
