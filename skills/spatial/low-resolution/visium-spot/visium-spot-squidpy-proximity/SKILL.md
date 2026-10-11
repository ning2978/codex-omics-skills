---
name: visium-spot-squidpy-proximity
description: Use when designing, reviewing, running, or packaging 10x Visium spot-level Squidpy spatial proximity analyses from AnnData objects, especially anchor-vs-neighbor tests using existing spatial graphs, deconvolution abundance matrices such as RCTD/cell2location, topology hops, physical distance rings, sample-aware backgrounds, threshold previews, and filtered heatmaps.
---

# Skill: Visium Spot Squidpy Proximity

## Purpose

Use this skill for **spot-level** spatial proximity questions in 10x Visium-like data, for example:

- Which deconvolved cell types are enriched near a target anchor spot population?
- Which RCTD/cell2location abundances are highest within 1-hop, 2-hop, or physical rings?
- How should low-support timepoints or low-background fold changes be filtered before plotting?

This skill is algorithm-oriented and cross-project. For CellCharter niche interpretation, combine it with `spatial/shared/algorithm-cellcharter` and `spatial/low-resolution/visium-spot/visium-spot-niche-annotation`.

## Interpretation Boundary

State this boundary in outputs and final answers:

> This is spot-level deconvolution proximity in Visium spots, not true single-cell physical distance.

Do not claim direct cell-cell contact unless the input is single-cell spatial data with appropriate segmentation.

## Start Here

1. Identify the analysis object and keys:
   - AnnData `.h5ad`
   - abundance matrix, commonly `obsm['celltype_abundance']` or cell2location abundance
   - abundance column names, often in `uns[...]` if `obsm` is not a DataFrame
   - anchor variable or cell type abundance column
   - optional niche/region column in `obs`
   - sample/donor/condition/timepoint metadata
   - existing Squidpy graph: `obsp['spatial_connectivities']`, `obsp['spatial_distances']`
   - coordinates for physical rings: `obsm['spatial']`
2. If the anchor threshold is not already fixed, run a threshold preview first.
3. Before full analysis, validate that the graph does not connect spots across samples.
4. Keep output tables metadata-rich: preserve `sample`, `donor`, `condition`, `timepoint`, and analysis labels.

## Canonical Workflow

1. **Threshold preview**
   - Compute per-sample anchor abundance distribution.
   - Report candidate top fractions such as Top5%, Top10%, Top20%.
   - Report how many anchors remain inside the target region/niche per sample.
   - Stop here if the user requested script-first review or threshold review.

2. **Anchor definition**
   - Define anchors per sample, not globally, unless the user explicitly requests a global cutoff.
   - For deconvolution abundance anchors, prefer continuous abundance top fractions over hard binary labels.
   - Restrict anchors to the target region/niche only when the biology question is region-specific.
   - Write `anchor_spot_summary_by_sample.csv` and flag low anchor counts.

3. **Neighbor definitions**
   - Topology mode: 1-hop direct graph neighbors; 2-hop neighbors of 1-hop, excluding anchors and 1-hop spots.
   - Physical mode: Euclidean rings from anchors using `obsm['spatial']`, for example `0-600`, `600-1200`.
   - Always exclude anchor spots from neighbor sets.
   - For multi-sample objects, compute neighbor sets within each sample.

4. **Background definitions**
   - Internal/region-specific analysis: background is same sample, same target region, non-anchor, excluding current neighbor set.
   - Boundary analysis: background is same sample or same sample+condition, non-anchor, excluding current neighbor set.
   - Avoid overlapping neighbor and background groups.

5. **Cell-type statistics**
   - For each cell type / abundance column calculate:
     - neighbor mean abundance
     - background mean abundance
     - fold change and log2 fold change
     - Mann-Whitney U p-value
     - BH-adjusted p-value
     - anchor, neighbor, and background spot counts
   - Keep raw counts of support in every table; these often explain surprising heatmap colors.

6. **Plots and interpretation**
   - Use abundance heatmaps for "what is most present near anchors".
   - Use log2FC heatmaps for "what is enriched near anchors versus background".
   - Treat log2FC cautiously when background abundance is near zero.
   - Mask or filter low-support summaries before ranking plots.

## Default Plot Filters

For timepoint/grouped heatmaps, apply support filters unless the user explicitly asks for raw diagnostics:

- `min_samples_with_neighbors >= 2`
- `min_total_neighbor_spots >= 20`
- for log2FC plots, `min_background_mean_abundance >= 0.001`

Show filtered or missing heatmap cells with a soft low-attention color, not a saturated color. Keep a separate filter summary CSV with failure reasons.

## Implementation Template

Use project-local scripts as templates only after parameterizing:

- input h5ad path
- abundance key
- anchor column or anchor abundance
- sample/group metadata keys
- graph key and coordinate key
- output directory
- threshold and support filters

Run post-hoc Scanpy/Squidpy inspection in the current project's validated Python environment, with writable cache directories when needed.

## Validation Checklist

Before calling a result complete:

- `py_compile` passes.
- `--help` works for CLI scripts.
- Required AnnData keys exist.
- Target anchor column exists.
- Target region/niche has spots, if used.
- Spatial graph exists and has shape `n_obs x n_obs`.
- Cross-sample graph edges are absent or explicitly handled.
- Threshold preview or anchor summary was reviewed.
- Output tables preserve sample metadata.
- Plots were regenerated from the intended filtered/unfiltered tables.

## More Detail

Read `references/spot-proximity-contract.md` when implementing a new script, reviewing output tables, or debugging disagreement between abundance and log2FC plots.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
