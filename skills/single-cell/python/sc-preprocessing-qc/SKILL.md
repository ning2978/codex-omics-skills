---
name: sc-preprocessing-qc
description: Receive whole-transcriptome scRNA-seq or snRNA-seq counts and perform capture-library-aware cell calling, QC, doublet assessment and optional ambient-RNA correction. Use before normalization; targeted spatial assays need platform QC.
---

# Count intake and QC

Read [../../workflow/workflow-sc-analysis/references/data-contract.md](../../workflow/workflow-sc-analysis/references/data-contract.md). Write a derivative; preserve sources and exact retained/excluded cell keys.

## Intake and calling

Inspect species/reference, assay/chemistry, count origin, feature types and filtered/unfiltered status. Preserve stable gene IDs/symbols and capture-aware keys; separate gene expression from antibody/VDJ features. Build library/sample/donor/batch/condition/pair metadata; missing design fields limit later inference.

For unfiltered droplets, inspect barcode ranks and assess empty droplets with an ambient-aware method such as DropletUtils `emptyDrops`, recording FDR, lower-count limit, high-count retention rule and barcodes. Without unfiltered droplets report upstream calling and `cell_calling=not_reassessed`; do not claim empty-droplet validation. Non-droplet platforms need their own rules.

## QC per physical library

- Derive counts, detected genes and mitochondrial fraction from original counts and species/reference annotations; optionally ribosomal/hemoglobin metrics. Plot distributions before/after per capture library and biological class when available.
- Select explainable platform/tissue-aware absolute or robust outlier thresholds, recording metric/direction/value/rationale per library. No universal mitochondrial or minimum-gene cutoff. High counts alone do not prove doublets; rare lineages can have different distributions.
- For snRNA-seq record intronic counting/nuclear chemistry, reduced mitochondrial/cytoplasmic signal, debris and complexity. Do not inherit scRNA mitochondrial thresholds. Stress/cell cycle are not default exclusions.
- Export a ledger for **every original cell**: key/library, metrics, rule flags, doublet decision, final inclusion and all reasons. Include calling decisions where available and attrition by sample/condition.

## Doublets

For droplet scRNA/snRNA default to Scrublet per **physical capture library**, on unnormalized raw counts after declared basic cell/low-quality screening. Do not split multiplexed captures by donor or use integrated/log/scaled input. Direct Scrublet or Scanpy's wrapper is acceptable in the checked environment.

Record loading/chemistry-derived expected rate, simulation parameters, version/seed, score/threshold, detected rate and uncertainty per library. Review competing-lineage coexpression, complexity, demultiplexing and sample distributions together; homotypic doublets may escape. Tiny/unsuitable libraries can be `SKIP` with reasons, not fabricated passing scores. Keep prediction, score and final exclusion separate; examine borderline-threshold sensitivity.

## Optional ambient RNA

Trigger from background/soup evidence. SoupX typically needs droplets/clusters; CellBender needs unfiltered counts; DecontX has cell-count/cluster assumptions. Check method requirements, record contamination and parameter sensitivity, and review low-RNA lineages. Preserve raw `counts`, store correction separately and declare normalization/model input. Fractional correction must not be rounded into purported raw counts.

## Handoff

Deliver retained full-count H5AD, input/sample manifest, QC plots/tables, all-cell ledger and doublet/background evidence including skipped steps. Verify finite nonnegative counts, unique keys, exact retained set and count immutability; integer-like counts for UMI models. Pass means traceable decisions, not absence of all contamination. Hand off to `sc-integration-clustering`; deletion invalidates existing downstream analysis.

## Validated count intake in an available downstream checkout

Use the platform QC skill for targeted per-ROI Xenium counts. Verify raw-count
provenance, unique keys and finite nonnegative integer values before QC; numerical
resemblance alone does not establish raw counts. This handoff does not replace
whole-transcriptome cell calling, doublet or mitochondrial QC. Require a complete
cell decision ledger, including cells in any excluded ROI.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
