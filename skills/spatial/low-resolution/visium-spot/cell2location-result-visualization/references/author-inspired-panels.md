# Author-Style Deconvolution Result Panels

Use this reference when a user wants cell2location plots to stay close to a paper or upstream workflow while still using project-local, scriptable code.

## Reusable Panel Types

- per-cell-type spatial abundance panels
- grouped abundance summaries by sample/condition/timepoint
- composition barplots or stacked area plots
- fold-change heatmaps against a chosen baseline
- QC metric panels for model outputs and abundance payloads
- marker or lineage panels that combine biological grouping with spatial abundance

## Implementation Rules

- Treat paper scripts as design references, not default execution entrypoints.
- Parameterize sample groups, color palettes, marker lists, and output roots.
- Keep plot generation non-interactive.
- Write the tables used for each figure.
- Record the abundance key, coordinate key, and any filtering rules in the output metadata or run summary.

## Validation

- Each expected figure exists and is non-empty.
- Each plotted table exists and has the expected grouping columns.
- Spatial panels use the intended coordinate orientation.
- Fold-change plots handle near-zero background values explicitly.
