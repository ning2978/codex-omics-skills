# Biomedical Skill Migration Map

Use this map only when an old skill name appears in a prompt, script, or document. New work must use the canonical taxonomy.

## Retired Active Skills

| Retired name | Canonical replacement | Reason |
|---|---|---|
| `spatial-project-analysis` | `workflow-st-analysis` | Duplicate spatial project/stage router |
| `cell2location-params-official` | `algorithm-cell2location` | Parameters belong to the algorithm skill |
| `xenium-scimap-rcn` | `xenium-scimap-neighborhood-analysis` | RCN and target-NMF share the same spatial-count contract |
| `xenium-scimap-nmf-neighborhood` | `xenium-scimap-neighborhood-analysis` | Folded into the unified scimap neighborhood skill |
| `st-deconvolution-cell2location` | `workflow-deconvolution-run` + `algorithm-cell2location` | Separate execution lifecycle from method logic |
| `st-cellcharter-niche-analysis` | `algorithm-cellcharter` + `visium-spot-niche-annotation` | Separate clustering method from biological annotation |
| `wound-st-analysis` | `workflow-st-analysis` plus the relevant implementation skill | Project-specific adapter retired |

## Stable Renames and Handoffs

- Old generic `st-data-preparation` calls remain `st-data-preparation`, routed through `workflow-st-analysis` when stage discovery is needed.
- Cell2location plotting remains `cell2location-result-visualization`; training and parameters use `algorithm-cell2location`.
- Project-local paths and cohort labels belong in repository adapters, not global skills.

## Forward Convention

- method-specific logic -> `algorithm-*`
- execution/stage logic -> `workflow-*`
- single-cell workflows -> `single-cell/`
- spot-level spatial -> `spatial/low-resolution/visium-spot/`
- bin/cell-level spatial -> `spatial/high-resolution/<platform>/`
- cross-platform spatial methods -> `spatial/shared/`
- retired or historical adapters -> `~/.codex/skills-archive/`

Do not create a second skill solely for a parameter preset, runtime wrapper, sensitivity branch, tissue, cohort, or one-off plot when the owning canonical skill can express it as a mode or reference.

## 43-skill publication refresh

The 34-skill release expands with nine canonical additions:

- Single-cell: `sc-preprocessing-qc`, `sc-integration-clustering`,
  `sc-marker-annotation`, `sc-sample-level-comparison`, `sc-cell-communication`.
- Spatial: `spatial-null-model`, `xenium-resegmentation`.
- Visualization: `render-scrna-dotplot`, `render-xenium-dotplot`.

Existing canonical names remain stable. The two DotPlot entrypoints now resolve
one shared engine and R renderer under `_shared/dotplot`; their CLI names remain
unchanged. Spatial executable helpers require `_shared/spatial`. Use the installer
to retain dependency paths. Project-specific scripts and historical cohort
benchmarks are not distributed; keep those in the consuming project's adapter.

Typed cohort predicates replace disease-code inference in the published DotPlot
engine. Legacy prose-only scopes require an explicit configuration mapping and
fail closed when absent. Annotation source, exact-key, priority and gate checks
remain mandatory. Single-group explanatory text moves from the figure caption to
result records without changing native statistics or the expression legend.

## 44-skill environment setup addition

Added `_meta/omics-environment-setup` for isolated CPU Scanpy and Seurat setup.
No existing skill is retired or renamed. Installation policy and runtime-preflight
retain their separate responsibilities.
