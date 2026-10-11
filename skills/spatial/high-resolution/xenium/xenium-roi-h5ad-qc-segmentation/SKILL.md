---
name: xenium-roi-h5ad-qc-segmentation
description: QC-filter existing segmented-cell Xenium per-ROI H5AD counts with a complete cell decision ledger and strict cell/gene merge checks. Does not rerun segmentation.
---

# Xenium ROI count and segmentation QC

Use existing segmented cells. Boundary stain is an accepted method; morphology supports QC decisions and does not establish biological identity. Choose real-data thresholds from distributions and tissue context; this entrypoint does not choose them scientifically.

## Input and current entrypoint

Use a project QC runner implementing the contract below. This package provides
the QC specification, not a project-specific `xenium_downstream` CLI. Configure
inputs, counts provenance, thresholds, segmentation methods, metadata and new
output directories explicitly; verify the runner against these checks.

Require one H5AD per ROI, unique stable `obs_names` across all inputs, unique genes, `obs[cell_area,nucleus_area,segmentation_method]` and two-dimensional `obsm["spatial"]` in microns. Default counts are `layers["counts"]`; recompute metrics from that matrix even when metrics already exist. Reject nonfinite, negative or noninteger counts. For raw `X`, explicitly specify `--counts-source X --counts-origin origin.json`; the JSON requires `kind=raw_counts`, `path`, `sha256`, and `matrix` (`X` or `layers/<name>`). Verify the source hash and key-aligned count values. Numerical resemblance alone is insufficient.

Before filtering/merging reject duplicate cell keys and different gene sets; equal gene sets in different order align to the first input. Preserve identifiers; never rename barcodes or intersect genes. Metadata must not overwrite an existing identifier.

Use an existing validated environment. Thresholds are dataset-specific. Optional
ROI minimum-cell exclusion, maximum nucleus count and per-tier area thresholds
are explicit decisions. Dropped ROIs remain in the all-cell ledger and receive
empty per-ROI outputs. Global percentiles use all input ROIs before filtering.
Existing output artifacts are protected from overwrite.

## Outputs and acceptance

- Per-ROI H5AD and optional `<sample>_filtered.h5ad`, preserving all genes and retained counts.
- `all_cell_qc_ledger.csv.gz`: every original `obs_name`, original metadata, recalculated metrics, every pass/fail rule, `retained`, and all `exclusion_reasons`.
- `n_genes_threshold_summary.csv` with TOTAL, `filter_retention_per_roi.csv`, `segmentation_method_breakdown.csv`, `global_area_thresholds.json`, `run_manifest.json` with sources/hashes/parameters.

Reopen outputs, match exact retained keys/count values, check merged row count equals the per-ROI sum and gene order matches the first input. Include excluded ROIs in attrition. Hand off counts-bearing full-gene H5AD to `sc-integration-clustering`; select integration and final resolution explicitly. For boundaries in existing parquet use `xenium-segmentation-new`; rerunning a segmentation model is a separate task.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
