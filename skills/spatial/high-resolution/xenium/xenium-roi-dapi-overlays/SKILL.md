---
name: xenium-roi-dapi-overlays
description: Use this skill when generating, auditing, rerunning, or validating Xenium ROI split-tile DAPI overlays from morphology_focus crops, especially per-ROI morphology_focus DAPI crop PNGs colored by a configured cell-label field, ROI manifests, sample_roi selections, and tmux/background runs. This skill is for centroid dot overlays only; if the task asks for cell_boundaries.parquet, nucleus_boundaries.parquet, irregular cell polygons, segmentation outlines, or polygon-filled cell overlays, use or combine with xenium-segmentation-new instead.
---

# Xenium ROI DAPI Overlays

## Purpose

Render per-ROI Xenium `morphology_focus/ch0000_dapi.ome.tif` crops with cell labels drawn as colored centroid scatter points. This skill captures the generic workflow for ROI split-tile figures and manifest validation.

Use this for requests such as:

- "Xenium ROI DAPI overlay"
- "ROI 拆片图"
- "morphology_focus overlay"
- "cell type 圆点图"
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

## Rendering workflow

Use a validated project renderer or implement this contract with tifffile,
NumPy and matplotlib. This skill does not bundle a project-specific renderer or
tmux launcher. Configure the image paths, label/coordinate columns, ROI selection,
pixel size, crop margin, maximum image dimension and new output directory.

1. Validate the pixel-to-physical coordinate transform and DAPI channel.
2. Lazily crop each image to the in-bounds centroids' bounding box plus the
   declared physical margin, then downsample without changing coordinate alignment.
3. Draw grayscale DAPI and centroid points using one explicit label-to-hex mapping
   across ROIs. Keep working/projected labels distinct from formal identities.
4. Write one PNG per ROI and a manifest containing sample/ROI, exact cell count,
   source image, transform, crop/downsampling, palette, PNG path and status.
5. Smoke-test representative small and large ROIs, visually inspect alignment,
   then run the full selection. Use a project background launcher if needed;
   protect existing artifacts unless overwrite is explicitly requested.

Example output: `<output_dir>/roi_dapi_overlays/`, with
`roi_dapi_overlay_manifest.csv` and `*_dapi_cell_type_overlay.png` files.

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
pngs = sorted(out_dir.glob("*_dapi_cell_type_overlay.png"))

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

- These figures are qualitative ROI context plots: DAPI grayscale background plus configured label centroids.
- Point size, alpha, and downsampling affect visual density; they do not change the underlying cell table.
- Projected/working label keys should not be treated as final annotation evidence by themselves. If the user asks to audit or revise biological labels, combine with label-transfer or annotation-review workflows.
- When the user asks whether the plot is a segmentation plot, answer explicitly: it is a centroid scatter overlay, not a polygon segmentation overlay.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
