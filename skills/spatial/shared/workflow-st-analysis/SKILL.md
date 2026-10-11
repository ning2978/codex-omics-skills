---
name: workflow-st-analysis
description: Stage and route spatial omics by platform, resolution, current artifact, data contract, and next handoff. Use when a Xenium, Visium/Visium HD, deconvolution, or spatial-analysis stage is unclear or spans steps.
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
| Requested new Xenium masks/transcript assignment | Select method and representative pilot | `xenium-resegmentation` |
| Per-ROI Xenium h5ad | Filter cell and segmentation QC | `xenium-roi-h5ad-qc-segmentation` |
| QCed segmented-cell Xenium expression subset | Recluster retained cells from counts | Pair with `sc-subcluster-analysis` |
| Existing segmented-cell Xenium expression clusters | Review markers and candidate identities without new clustering | `sc-marker-annotation` plus platform/project contracts |
| Requested segmented-cell subcluster–annotation–cleanup loop | Round review, concrete human confirmation, optional mixed whole-cluster proposal with target-cell loss | `sc-subcluster-analysis` plus `sc-marker-annotation`; spots/bins remain mixture workflows |
| Existing Xenium transfer predictions | Audit measured expression, program evidence, sample support and QC | `xenium-scrna-label-transfer` |
| Existing Visium spot niches or HD bin domains | Review anatomy, programs, composition and QC | `visium-spot-niche-annotation` or `visiumhd-spatial-domain-neighborhood` |
| QCed spot-level ST plus reference | Launch deconvolution | `workflow-deconvolution-run` plus `algorithm-cell2location` or `algorithm-rctd` |
| Existing deconvolution object | Validate and summarize results | Algorithm-specific visualization skill |
| Trusted segmented-cell phenotypes and geometry | RCN or target phenotype-composition NMF | `xenium-scimap-neighborhood-analysis` |
| Fixed anchors and candidate labels | Conditional proximity/cooccurrence/enrichment null | `spatial-null-model` |
| Explicit measured expression and spatial geometry | Panel/complex audit then communication model | `algorithm-commot` |
| Existing spatial labels or embedding | Domain or annotation analysis | Matching platform/algorithm skill |

## Stable Contracts

Require the following before a handoff:

- expression stages: explicit measured matrix source/meaning; preserve raw integer counts for counts-dependent work; geometry/composition stages need no expression matrix;
- coordinates stored in `obsm["spatial"]` or a documented equivalent with units and orientation;
- sample/library boundaries that prevent cross-sample spatial edges;
- unique cell keys and ROI nesting in genuine biological units; only the condition/batch/label fields needed by the requested stage, with label evidence grade;
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

## Biological review handoff

For new biological naming or annotation review, read [the joint evidence review](../../../single-cell/python/sc-marker-annotation/references/annotation-evidence-review.md). Handoff combines full available-gene differential statistics and a candidate view with independently checked canonical/competing programs, biological-unit support, and applicable QC. Record program sources and measured availability; missing targeted genes are not negative evidence.

Segmented Xenium profiles support cell-identity interpretation with segmentation alternatives. Visium spots/HD bins retain mixture and anatomy-first interpretation unless documented as segmented-cell profiles. Domain, RCN and CellCharter niche assignments stay separate from identities; transferred labels remain predictions. Existing-label plotting or marker export does not automatically request annotation, removal, or reclustering.

## Environment Preflight

- Prefer a repository-provided wrapper or locked environment.
- Verify imports and accelerator availability before a heavy run.
- Apply `environment-installation-policy` when an environment must be created or changed.
- Keep notebook exploration separate from the reproducible CLI or script used for the final run.
- General Xenium expression clustering may use scVI only on explicit user
  request, after exact counts/data/GPU preflight and a real scVI CUDA pilot
  pass, with CPU fallback disabled. Record batch key, model/training
  parameters, trainer device, peak CUDA allocation, training history, saved
  model, and complete expression-recluster provenance. CellCharter latent
  mode remains a compatible explicit-request alias, not the only allowed use.

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
- a biological review records program availability, unit support, applicable
  QC, material alternatives and decision/status, or explicit evidence limits;
- logs or summaries contain a clear successful completion marker.

## Boundaries

- Do not duplicate algorithm parameters here.
- Do not encode one tissue, cohort, or repository layout.
- Do not keep routing inside this skill after a specialized skill has taken ownership.

## Independent stage handoffs

Read [stage-contract.md](references/stage-contract.md) when linking stages or preparing portable tools. Each stage can run alone with its own minimum inputs; artifacts include upstream provenance and validation scope. A geometry PASS, supplied evidence audit, model smoke test and biological validation are separate states.

Current package commands, repository paths, cohort codes and registry roles belong in a repository-local adapter. Generic bundled tools do not import project packages. Resegmentation creates new cell/transcript-assignment versions and invalidates downstream counts/QC/annotation/spatial products for that new cell set while preserving old outputs. Fixed-selection spatial nulls do not correct same-data selection. Group inference requires genuine biological units, not cells or nested ROIs.

## Segmented-cell Xenium expression delivery

Only for a complete segmented-cell Xenium expression reduction, clustering,
annotation or subclustering handoff, read [the shared expression delivery contract](../../../single-cell/workflow/workflow-sc-analysis/references/delivery-contract.md).
Keep full stage artifacts/validation internal; `delivery/` selects final UMAP
PNGs, marker DotPlot PNGs, cluster review CSV and README. Standalone plots,
marker exports or review tasks keep their requested scope. This boundary does
not apply to spatial geometry, RCN/NMF/niches, Visium spots/HD bins or other
spatial-analysis outputs.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
