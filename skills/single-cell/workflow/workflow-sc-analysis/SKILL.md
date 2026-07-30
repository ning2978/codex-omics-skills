---
name: workflow-sc-analysis
description: Use this skill for single-cell RNA-seq workflow staging and contracts, including preprocessing, reference preparation, conversion, downstream handoff, and routing post-filtering subset reruns with new UMAP/Leiden/markers to sc-subcluster-analysis.
---

# Skill: workflow-sc-analysis

## Purpose

Use this skill for stage-driven scRNA analysis workflows, independent of project.

## Stage Router

1. raw counts intake
2. QC + normalization + annotation
3. reference-ready export for deconvolution
4. optional Seurat/Scanpy conversion
5. filtered subset rerun with new UMAP/Leiden/markers

When the user asks to filter/remove contamination or clean a subset and then rerun UMAP, Leiden, marker discovery, reclustering, re-embedding, 重新降维, 重新聚类, or 新的 UMAP, route to `sc-subcluster-analysis`. The rerun must start from retained-cell counts and recompute HVG/PCA/batch correction/neighbors/UMAP/Leiden; a sidecar-only filter table or projection plot does not satisfy this workflow stage.

## Stable Contracts

- preserve raw counts for downstream model-based methods
- keep annotation keys explicit and stable
- export machine-readable manifest with sample and label mapping

## Handoff Targets

- if goal is spatial deconvolution, hand off to:
  - `workflow-deconvolution-run`
  - plus selected algorithm skill
- if goal is subclustering or post-filtering rerun with new UMAP/Leiden/markers, hand off to:
  - `sc-subcluster-analysis`

This workflow skill is only a stage router and contract checker for subclustering requests. Do not treat sidecar-only filtering, relabeling, or global-coordinate overlay plots as completion of a requested subset rerun.
