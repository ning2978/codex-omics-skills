---
name: workflow-sc-analysis
description: Route conventional scRNA-seq or snRNA-seq from count intake and QC through integration, annotation, sample-level comparison and reproducible delivery. Use for multi-stage requests; subset reruns route to sc-subcluster-analysis.
---

# Conventional single-cell workflow

Use Python for the AnnData main workflow and R for statistics/specialist methods. Read [references/data-contract.md](references/data-contract.md) before acting. This entry owns stage routing; leaves own methods. Read selected leaves completely, even if discovery-disabled. Paths below resolve from the skills root (`../../..` from this directory).

## Stage router

| Stage / leaf | Inputs | Completion evidence |
|---|---|---|
| Intake/QC: `single-cell/python/sc-preprocessing-qc` | Counts, platform/reference, capture metadata | Unique keys, immutable counts, design table, per-library calling/QC/doublets, per-cell exclusion ledger |
| Representation: `single-cell/python/sc-integration-clustering` | Retained counts and design | Correct log/HVG/PCA inputs, unintegrated baseline, integration assessment, named graph, resolution/seed candidates, selected Leiden/UMAP |
| Annotation: `single-cell/python/sc-marker-annotation` | All-gene log expression, clusters, sample/unit design and applicable QC | Complete markers/candidate view, sourced positive/competing programs with availability, unit support/QC, versioned identity/subtype/state and explicit unknown/ambiguous |
| Comparison: `single-cell/python/sc-sample-level-comparison` | Raw counts, annotations, biological-unit design | Estimable model or descriptive-only status, pseudobulk DEG, proportions, effects/FDR, configured enrichment/scores |
| Subset/filter rerun: `single-cell/python/sc-subcluster-analysis` | Retained cell keys and counts | Recomputed normalization/HVG/PCA/integration/graph/Leiden/UMAP, all-gene markers, removal evidence |
| Reviewed subcluster–annotation–cleanup loop: `single-cell/python/sc-subcluster-analysis` | Full counts, target lineage, sourced programs and human review | Frozen per-round actions/exact cells, partial or full confirmation, counts rerun, cumulative ledger and separate technical/scientific stop states |
| R conversion: `single-cell/shared/h5ad-to-r-bridge` | Declared counts and requested metadata/reductions | Sparse count/key-order/feature-map/reduction parity |
| Deconvolution reference | Counts, reviewed label key/sample map | Reference manifest with genes, annotation version/confidence, units, exclusions/limitations; hand off to `spatial/shared/workflow-deconvolution-run` and algorithm |

Run only requested stages; supplied validated artifacts are prerequisites, not newly recomputed stages. Missing biological units block group inference, not exploration. Distinguish library/sample/donor/condition before correction/testing; `sample_id` does not automatically justify Harmony.

For annotation or review, the marker leaf owns [the joint evidence review](../../python/sc-marker-annotation/references/annotation-evidence-review.md). Handoff includes matrix/filter provenance, canonical/competing program sources and availability, biological-unit support, applicable QC, and material alternatives. A shared canonical gene need not rank differentially; a single-sample group needs a support limitation rather than automatic rejection. Review-only work may finish with evidence tables; marker export or plotting alone does not authorize new labels or cell removal.

## Optional question-driven handoffs

- Ligand/receptor predictions: `single-cell/python/sc-cell-communication`.
- Velocity/topology/fate: `single-cell/python/sc-velocity-paga-cellrank-analysis`.
- Monocle2: `single-cell/r/sc-r-monocle2-analysis`.
- GRN/motif/AUCell: `single-cell/python/sc-pyscenic-regulon-analysis`.
- TCR/BCR: `single-cell/python/sc-multiomics-tcr-vdj-analysis`.
- CNV only when requested and assay-suitable: declare reference and method; expression-derived CNV alone does not prove malignancy.
- Heavy/failing runs: `_meta/omics-analysis-runtime-preflight`; formal manuscript claims additionally use `publication-grade-omics-analysis`.

Targeted spatial platforms route via `skills-taxonomy` to their own QC/gene-availability rules. Filtering plus a new embedding/reclustering requires retained-count recomputation; sidecar relabeling or historical overlays cannot complete it. Pure UMAP layout changes use the unchanged-graph exception and preserve clusters.

Retain internally the full-count annotated H5AD, sample/design table, exclusion ledger, complete statistics, annotation evidence/version, parameters/software/seeds/sources and machine-readable validation. Optional scaled HVG objects declare limited scope. Expression plots explicitly select log expression; publication DotPlots use `render-scrna-dotplot` with native Seurat/Python parity. Source H5AD, counts, registry and prior results remain protected from implicit mutation.

## Final user handoff

Read [the shared delivery contract](references/delivery-contract.md) for complete
expression-analysis handoffs. Required internal completion artifacts above remain
mandatory; select only final UMAP PNGs, marker DotPlot PNGs, cluster review CSV and
README into `delivery/`. Narrow plotting/export/review requests keep their own
scope. Packaging does not rerun analysis; final replies link the delivery content.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
