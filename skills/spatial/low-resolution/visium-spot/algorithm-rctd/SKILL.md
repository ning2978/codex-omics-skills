---
name: algorithm-rctd
description: Use this skill for RCTD-specific deconvolution workflow, input requirements, run patterns, and output validation across projects.
---

# Skill: algorithm-rctd

## Purpose

Use this skill when the deconvolution algorithm is RCTD, regardless of project.

Typical scope:

- RCTD input contract checks (reference + spatial + platform assumptions)
- running RCTD in Seurat/R pipelines and exporting reusable outputs
- converting RCTD outputs for downstream Python/ST workflows

## Use With Workflow Skills

- Pair with `workflow-deconvolution-run` for execution orchestration.
- Pair with `spatial/shared/workflow-st-analysis` when ST preprocessing or conversion is needed.

## Minimum Input Contract

- spatial counts with barcode consistency
- matched reference cell states/labels
- explicit metadata keys for downstream grouping

## Output Contract (Recommended)

- per-sample deconvolution table (spot x cell type)
- serialized run object (RDS or equivalent)
- exported standardized matrix/table for downstream Python analysis
- run manifest with sample ID, reference version, and key parameters

## Operational Guidance

- use explicit output roots and stable naming
- keep conversion scripts deterministic and log-rich
- avoid hidden notebook state as the only execution path
