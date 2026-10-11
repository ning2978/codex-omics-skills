# Shared DotPlot style

`style_defaults.json` is the shared source of defaults. Both Python entrypoints
pass resolved style settings to the same R renderer via the grouping contract.
Explicit `output.height_inches`, `output.width_inches`, `output.dpi` and `style`
values override defaults; width still follows the existing auto-expansion rule.
Height defaults separately per figure to max(3.2, 2.2 + 0.30 × group count).

White background, no grid, thin left and bottom axes; a 0.5-inch pastel strip,
horizontal bold black labels (newlines and explicit short labels supported).
Gene labels: 60°, 9 pt. Groups/strip/legends: 11 pt. Title: 14 pt. PNG: 240 DPI.
Native Seurat sizes and Percent Expressed legend are retained. Blue-white-red
expression colors use [-2.5, 2.5], with display-only squishing outside the range.
Actual settings, program colors/short labels, dimensions and scaling mode are
recorded in run/layout contracts and marker audits. Single-group expression is
unscaled native log1p mean, with missing z-scores. Record this in the layout/numeric contracts and
accompanying documentation; retain the Expression (unscaled) legend.

## Figure presentation

Do not add explanatory captions, method descriptions, footnotes or provenance
below figures by default. Put explanations in accompanying documentation and
machine-readable result records. Preserve axis/tick/panel labels, program
strips, legends and color-bar labels. Titles follow the task. Inspect exported
PNGs for overlap, cropping and excessive margins before delivery.

Historical redraws use a separate output directory with explicit source and
candidate status and retain original statistics and labels. They do not bypass
current annotation gates or promote candidate labels.

For long or single-gene programs, supply multiline `short_label` or an explicit
abbreviation and inspect the PNG for overlap and cropping; increase width or
strip height when necessary. Long gene names may also require an explicit
figure height to leave enough room for the data rows at 60 degrees. No silent truncation or automatic label merging.
