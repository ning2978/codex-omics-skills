---
name: render-scrna-dotplot
description: Render and validate publication-style, PNG-only scRNA-seq marker DotPlots with native Seurat DotPlot statistics, its native Percent Expressed size legend, a fixed blue-white-red expression scale, and a patchwork marker-program strip. Use for counts-H5AD plus registry-declared unique sidecar-key visualization, configurable marker programs/groupings, marker audits, and Python-Seurat numeric parity. Do not use to create or promote annotations, recluster cells, or mutate source files.
---

# Render scRNA-seq DotPlot

Create a derivative-only figure package from complete `layers["counts"]`, an
explicit unique sidecar key mapped to AnnData `obs_name`, and declarative
grouping and marker contracts.

Set `SKILLS_ROOT` to the installed skills root (the directory containing this
skill and `_shared/`) before running the commands below.

## Workflow

1. Read every applicable project `AGENTS.md`. Before inspecting labels, read the
   annotation guide and registry, resolve exactly one use case per semantic role, and stop if its
   semantic contract conflicts.
2. Read [references/input-contract.md](references/input-contract.md) completely
   before preparing a run YAML. Use the registry-declared scRNA join exactly.
3. Use a new output directory. Never overwrite source data, sidecars,
   annotations, registries, guides, or prior result packages.
4. Run from the project root with a Python environment containing h5py, NumPy,
   pandas, Pillow, PyYAML, and SciPy:

   ```bash
   python ${SKILLS_ROOT}/render-scrna-dotplot/scripts/render_scrna_dotplot.py render --config /absolute/path/run.yaml
   ```

5. Accept the result only after `PASS`, then verify it independently:

   ```bash
   python ${SKILLS_ROOT}/render-scrna-dotplot/scripts/render_scrna_dotplot.py validate --output /absolute/path/result --verify-inputs
   ```

## Fixed semantics

- Require an explicit unique `sidecar_key -> obs_name` exact-set join. A
  registry-declared sample-prefixed sidecar column named `cell_id` is allowed
  only as that sidecar key; never join it to a non-unique AnnData obs column or
  use row order.
- Scan complete CSR `layers["counts"]`; compute linear CP10K, group mean
  `avg.exp`, `log1p(avg.exp)`, per-gene sample z-score (`ddof=1`), clipped
  display values, and `pct.exp` independently in Python.
- For each grouping, set Seurat `Idents`, call exactly
  `DotPlot(object, features=markers, assay="RNA")`, and export unmodified
  `plot$data` before the display scale is changed.
- Preserve native Percent Expressed point sizes and legend. Override only the
  final colour scale with `#2166AC -> #FFFFFF -> #B2182B`, midpoint `0`, and
  limits `[-2.5, 2.5]`.
- Add a legend-free patchwork program strip and recompute strip lengths and
  positions after filtering.
- Audit every marker slot, duplicate, unmeasured gene, threshold decision,
  role, program, short label, color, and strip span. Never treat an unmeasured
  marker as zero.
- Generate only opaque-white PNG. Render x-axis gene labels at 60 degrees with
  `hjust=1`, `vjust=1`, and 9 pt text. Treat the configured width (default 20
  inches) as a minimum and auto-expand by default to
  `max(configured_width, 6.0 + marker_count * max(0.42, 0.18 + 0.025 *
  longest_gene_length))`. If auto-expansion is disabled, fail when the
  configured width is below that lower bound. Default height per grouping is `max(3.2, 2.2 + 0.30 * n_groups)` inches and
  default resolution is 240 DPI.
- Fail closed on invalid joins or counts, R/Seurat version mismatch, numeric
  parity failure, invalid artifacts, changed protected inputs, or manifest
  drift. Reject legacy midpoint `#F7F7F7`; migrate to `#FFFFFF`.

Rendering never relabels, backfills, formalizes, promotes, clusters, or mutates
cells. Keep native top-marker labels, fine/state evidence, transfer evidence,
claim classes, and cross-modality labels separate.

## Bundled resources

- `scripts/render_scrna_dotplot.py`: `render`, `validate`, and `self-test` CLI.
- [Shared R renderer](../_shared/dotplot/seurat_dotplot_renderer.R): native Seurat renderer and program strip.
- `references/input-contract.md`: complete run and marker YAML contracts.

## Shared implementation and source validation

Both DotPlot entrypoints execute `../_shared/dotplot/engine.py` in independent
module namespaces and use one shared R renderer. Keep CLI and YAML entrypoints.
The renderer compares the complete requested track against its marked JSON
summary in the guide, ignoring only `updated_at`. Missing summaries fail closed.
It verifies labels against registered source rows and applicability gates.
See the input contract for per-grouping source/role selection.

## Default visual style and single-group semantics

Use the shared [style contract](../_shared/dotplot/style.md) for defaults,
configuration overrides, and the historical visual reference. White background,
no grid, thin left/bottom axes; a fixed 0.5-inch pastel strip with horizontal,
bold black labels; gene labels 60°/9 pt; groups, strip and legends 11 pt; title
14 pt; PNG at 240 DPI. Explicit dimensions and style settings take precedence.
Record actual style, dimensions, labels and colors in result contracts.

A single group is supported: native Seurat disables cross-group scaling and
returns log1p mean expression. Keep z-score missing, retain unmodified native
values, and record the unscaled-expression explanation in result records and documentation. Do not add a bottom caption. Constant multi-group genes
retain native missing scaled values; parity requires identical missing positions.
The display scale squishes values outside [-2.5, 2.5] without changing statistics.
Use task-specific titles/subtitles reflecting lineage, grouping and annotation
status; never default to PCa, v2 or formal=false. Historical candidate redraws
must name their source/status, retain all labels and statistics, and use an
independent output directory; current registry gates remain required.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
