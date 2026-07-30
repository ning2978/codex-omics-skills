---
name: xenium-roi-dapi-overlays
description: Use this skill when generating, auditing, rerunning, or validating Xenium ROI split-tile DAPI overlays from morphology_focus crops, especially per-ROI morphology_focus DAPI crop PNGs colored by scVI cell_type_draft centroid scatter points, ROI manifests, sample_roi selections, and tmux/background runs. This skill is for centroid dot overlays only; if the task asks for cell_boundaries.parquet, nucleus_boundaries.parquet, irregular cell polygons, segmentation outlines, or polygon-filled cell overlays, use or combine with xenium-segmentation-new instead.
---

# Xenium ROI DAPI Overlays

## Purpose

Render per-ROI Xenium `morphology_focus/ch0000_dapi.ome.tif` crops with cell labels drawn as colored centroid scatter points. This skill captures the generic workflow for ROI split-tile figures and manifest validation.

Use this for requests such as:

- "Xenium ROI DAPI overlay"
- "ROI 拆片图"
- "morphology_focus overlay"
- "scVI cell type 圆点图"
- "audit ROI DAPI overlay manifest"
- "run ROI DAPI overlays in tmux"

Do not present these outputs as cell segmentation polygon overlays. The renderer does not read `cell_boundaries.parquet` or `nucleus_boundaries.parquet`; it plots cell centroids from `pixel_x_int` / `pixel_y_int`. For irregular cell outlines, filled cells, segmentation-method coloring, or boundary QC, switch to `xenium-segmentation-new`.

## Input Contract

Required inputs:

- Per-sample or merged cell/morphology feature tables with centroid coordinates and labels.
- Xenium DAPI OME-TIFFs under each sample's `morphology_focus/` directory:
  - `ch0000_dapi.ome.tif`
- Required table columns:
  - `sample_id`
  - `sample_roi`
  - `roi_id`
  - `cell_type_draft`
  - `pixel_x_int`
  - `pixel_y_int`
  - `morphology_in_bounds`

Use the current project's validated Python environment.

## Optional repository adapters

Use project scripts directly when they already follow this contract; do not
copy project paths or cohort defaults into this global skill.

- `xenium_morphology_focus_roi_dapi_overlays.py`
  - Lazily crops `ch0000_dapi.ome.tif`.
  - Uses per-sample morphology parquet rows.
  - Filters to `morphology_in_bounds`.
  - Computes an ROI bbox from centroid pixel coordinates plus `--margin-um`.
  - Downsamples by `--target-max-dim`.
  - Draws colored centroid scatter points by `--label-col`, default `cell_type_draft`.
  - Writes one PNG per `sample_roi` plus `roi_dapi_overlay_manifest.csv`.

- `run_xenium_roi_dapi_overlays_tmux.sh`
  - Starts, checks, tails, or stops the full renderer in a tmux session.
  - Uses `PYTHON`, `ANALYSIS_DIR`, `OUT_DIR`, `TMUX_SESSION`, and `EXTRA_ARGS` environment variables when provided.

## Recommended Workflow

First run a small smoke test with representative small and large ROIs. Keep the output in a temporary or explicit test folder unless intentionally overwriting the main outputs.

```bash
python \
  /path/to/project/xenium_morphology_focus_roi_dapi_overlays.py \
  --sample-rois <sample_roi_1> <sample_roi_2> \
  --out-dir <output_dir>/roi_dapi_overlays_smoke \
  --force
```

Then run the full job in tmux:

```bash
bash /path/to/project/run_xenium_roi_dapi_overlays_tmux.sh start
bash /path/to/project/run_xenium_roi_dapi_overlays_tmux.sh status
bash /path/to/project/run_xenium_roi_dapi_overlays_tmux.sh tail
```

To pass debug or selection options through the tmux launcher:

```bash
EXTRA_ARGS="--sample-rois <sample_roi_1> <sample_roi_2>" \
  bash /path/to/project/run_xenium_roi_dapi_overlays_tmux.sh start
```

Example final output directory:

```text
<output_dir>/roi_dapi_overlays/
```

Example manifest:

```text
<output_dir>/roi_dapi_overlays/roi_dapi_overlay_manifest.csv
```

## Validation

Expected full-run checks:

- Manifest row count matches the selected ROI count.
- `sample_id` row counts match the ROI manifest.
- `n_cells` totals match filtered in-bounds cells used for plotting.
- There is one non-empty PNG per successfully rendered ROI.
- Manifest statuses are expected to be `written` on forced renders or `exists` on non-forced reruns; any `missing_dapi`, `unknown_sample`, or `no_in_bounds_cells` row needs investigation.
- PNG paths listed in `png_path` exist and are non-empty.

Fast local validation snippet:

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

out_dir = Path("<output_dir>/roi_dapi_overlays")
manifest = pd.read_csv(out_dir / "roi_dapi_overlay_manifest.csv", dtype={"sample_id": str, "sample_roi": str, "roi_id": str})
pngs = sorted(out_dir.glob("*_dapi_scvi_cell_type_overlay.png"))

print("rows", len(manifest))
print("sample_counts", manifest["sample_id"].astype(str).value_counts().sort_index().to_dict())
print("cells", int(manifest["n_cells"].sum()))
print("pngs", len(pngs))
print("empty_pngs", sum(p.stat().st_size == 0 for p in pngs))
print("statuses", manifest["status"].value_counts().to_dict())
print("missing_png_paths", sum(not Path(p).exists() or Path(p).stat().st_size == 0 for p in manifest["png_path"]))
PY
```

## Interpretation Guardrails

- These figures are qualitative ROI context plots: DAPI grayscale background plus scVI label centroids.
- Point size, alpha, and downsampling affect visual density; they do not change the underlying cell table.
- Projected/working label keys should not be treated as final annotation evidence by themselves. If the user asks to audit or revise biological labels, combine with label-transfer or annotation-review workflows.
- When the user asks whether the plot is a segmentation plot, answer explicitly: it is a centroid scatter overlay, not a polygon segmentation overlay.
