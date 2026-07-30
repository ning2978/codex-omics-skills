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
