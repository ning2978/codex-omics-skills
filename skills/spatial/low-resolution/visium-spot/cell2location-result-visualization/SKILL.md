---
name: cell2location-result-visualization
description: Use this skill when summarizing, plotting, or regenerating figures from existing cell2location `spatial.h5ad` outputs, especially for abundance summaries, QC tables, grouped comparisons, spatial panels, and server-side non-interactive figure generation.
---

# Skill: cell2location-result-visualization

## Status

Algorithm-specific visualization skill in the new taxonomy.

Recommended routing:

- stage routing by `spatial/shared/workflow-st-analysis`
- algorithm routing by `spatial/low-resolution/visium-spot/algorithm-cell2location`

## Purpose

Use this skill after cell2location outputs already exist and the goal is to regenerate summaries or plots from `spatial.h5ad` files.

This is a good reusable skill candidate because plotting after deconvolution usually converges on a small set of stable inputs and outputs:

- `spatial.h5ad`-like result objects
- abundance matrices
- sample metadata
- exported summary tables and figures

It is also the right place to absorb reusable plotting patterns from the author ST scripts without inheriting their hard-coded Rmd structure.

## Use This Skill When

- the user wants to regenerate summary tables or plots from existing deconvolution results
- the task is about QC summaries, grouped comparisons, abundance trends, or optional spatial panels
- the plotting should run non-interactively on a server

## Do Not Use This Skill When

- the task is to rerun cell2location training
- the task starts from raw Visium data and still needs ST preparation
- the task is mainly about choosing model parameters rather than visualizing outputs

## Quick Start

1. Confirm deconvolution outputs already exist.
2. Confirm the abundance payload and coordinate payload names.
3. Use a reproducible environment wrapper if the project provides one.
4. Prefer non-interactive runs on the server.

## Stable Inputs

- Generic inputs:
  - sample-level `spatial.h5ad` or equivalent result objects
  - an abundance payload in `obsm`
  - coordinates for spatial maps
- Common payloads:
  - `obsm['q05_cell_abundance_w_sf']` or another explicit abundance matrix
  - `obsm['spatial']` for spatial maps
  - sample/condition metadata in `obs`

## Stable Outputs

- Generic outputs:
  - figures root
  - abundance summary tables
  - QC summary tables
  - grouped comparison plots
- Common tables:
  - `cell_abundance_summary.csv`
  - `fold_change_summary.csv`
  - `qc_metrics.csv`
- Common plots:
  - composition summaries
  - cell-type heatmaps
  - abundance boxplots
  - time/group trend plots
  - fold-change heatmaps
  - QC metric plots

## Packaging Guidance

This plotting layer is suitable for a reusable skill because:

- input structure is usually stable once deconvolution finishes
- outputs are easy to standardize across projects
- scripts can usually be made runnable non-interactively with CLI flags

For direct server testing, start with summary plots and QC plots. Spatial panel generation is more expensive and should be opt-in.

## Recommended Run Pattern

Use the validated project Python interpreter and a plotting runner that consumes
the result contract above. Supply its results directory explicitly. No
project-specific renderer or environment wrapper is bundled. Start with summary
and QC panels; generate per-sample spatial maps when requested.

## Author Workflow References

When you want plotting to stay close to the paper-source ST scripts, read:

- `references/author-inspired-panels.md`

That file extracts the reusable parts of the author workflow:

- marker-panel organization
- grouped spatial plot conventions
- deconvolution summary display patterns

## Related Skills

- Use `spatial/shared/workflow-deconvolution-run` + `spatial/low-resolution/visium-spot/algorithm-cell2location` when the missing piece is generating `spatial.h5ad`.

## Boundaries

- Use this skill for result regeneration, not for rerunning cell2location training.
- Do not treat old notebook-only plotting fragments as the default interface.
- If the plotting logic needs richer CLI controls later, extend the script rather than reintroducing interactive prompts.
- Treat author Rmd files as design references, not direct execution entrypoints.
- Keep biological groupings reusable; avoid hard-coding disease- or project-only comparison labels into the generic skill definition.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
