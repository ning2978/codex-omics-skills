---
name: xenium-segmentation-new
description: Use this skill when working with new-version 10x Xenium (XOA 2.0+) cell segmentation outputs, especially the parquet trio cells.parquet, cell_boundaries.parquet, nucleus_boundaries.parquet, including loading polygons in long-format efficiently, segmentation QC (boundary-stain vs interior-stain vs nucleus-expansion fallback, multi-nucleated cells), and rendering polygon overlays for cropped regions.
---

# Skill: xenium-segmentation-new

## Purpose

Handle 10x Xenium **new-version (XOA 2.0+)** cell segmentation outputs, where boundaries are stored as **parquet long-tables** (replacing the legacy `cell_boundaries.csv.gz` / `nucleus_boundaries.csv.gz`) and `cells.parquet` carries multi-stain segmentation provenance.

Use when:

- inputs include `output-XETG*/cell_boundaries.parquet` and/or `nucleus_boundaries.parquet`
- `cells.parquet` has `segmentation_method` and `nucleus_count` columns
- you need cell/nucleus polygons for QC, visualization, niche overlays, or downstream rasterization (masks)

Do NOT use this skill for:

- legacy CSV.gz boundaries from XOA ≤ 1.6 (different grouping rules; convert first)
- transcripts.parquet analyses (separate concern)
- cell × gene matrix QC (covered by `cell_feature_matrix.h5` workflows)

## New-Version Output Contract

Each `output-XETG{instrument}__{slide}__Region_*__{date}__{time}/` folder MUST contain:

| File | Schema |
| --- | --- |
| `cells.parquet` | `cell_id` (str), `x_centroid`, `y_centroid` (µm, double), `transcript_counts`, `control_probe_counts`, `genomic_control_counts`, `control_codeword_counts`, `unassigned_codeword_counts`, `deprecated_codeword_counts`, `total_counts` (int64), `cell_area`, `nucleus_area` (µm², double; `nucleus_area` is NaN for anucleate cells), `nucleus_count` (int64; 0 = anucleate, ≥2 = multi-nucleated), `segmentation_method` (str) |
| `cell_boundaries.parquet` | long-format vertices: `cell_id` (str), `vertex_x`, `vertex_y` (float32, µm), `label_id` (int64; >1 means multiple sub-polygons per cell, rare) |
| `nucleus_boundaries.parquet` | same schema as above; one polygon per nucleus, multiple per cell when `nucleus_count > 1` |

**Coordinate system**: microns, top-left origin (image convention). Always invert the y-axis when plotting (`ax.set_ylim(ymax, ymin)` or `ax.invert_yaxis()`).

**Polygon resolution**: 10x writes a **fixed ~25-vertex** smoothed ring per cell. Vertex count is NOT a polygon-complexity signal in new XOA; treat as constant-resolution.

**`segmentation_method` values** (typical multi-stain Xenium 5K Prime / Multimodal Cell Segmentation kit):

- `Segmented by interior stain (18S)` — primary path (~70 % of cells in our prostate runs)
- `Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)` — multi-protein boundary path (~29 %)
- `Segmented by nucleus expansion of 5.0µm` — fallback when no usable stain (rare, ~1.6 %; flag for QC)

## Memory Discipline

`cell_boundaries.parquet` for one Xenium slide is typically **20–30 M rows / 23 row groups / ~80 MB**. Strategies:

1. **Project columns**: read only `cell_id`, `vertex_x`, `vertex_y`, `label_id`.
2. **Filter early**: pass a `cell_ids` set or a spatial `bbox` to drop vertices before grouping.
3. **Stream by row group**: use `pyarrow.parquet.ParquetFile.read_row_group(rg)`, not `read_table()`, for slide-wide passes.
4. **Vertex-level bbox is approximate**: for exact polygon clipping, convert to shapely / GeoPandas downstream.
5. Avoid converting all polygons to shapely up front; keep `dict[cell_id -> ndarray (N,2) float32]` until shapely is required.

## Optional implementation modules

When the target repository provides them, use a reusable IO module and two
CLIs that satisfy this interface:

- `xenium_segmentation_io.py`
  - `find_xenium_output_dir(sample_root)` → unique `output-XETG*` path
  - `load_cells_table(cells_parquet, columns=...)` with per-column fallback for resilient reads
  - `load_polygons(boundary_parquet, cell_ids=None, bbox=None)` → `dict[cell_id -> (N,2) float32]`
  - `iter_polygons(boundary_parquet, chunk_row_groups=1)` → streaming generator (cell_id boundary spanning row groups handled)
  - `polygons_to_geodataframe(polys)` (lazy `geopandas` / `shapely` import)
  - `bbox_filter_cells(cells_df, xmin, xmax, ymin, ymax)` for centroid filtering

- `xenium_segmentation_qc.py` — per-sample CLI
  - segmentation_method breakdown (count + fraction)
  - anucleate / mono / multi-nucleated counts (from `nucleus_count`)
  - cell_area, nucleus_area, nucleus/cell area ratio percentiles
  - optional: vertex-count summary per polygon (sanity check that 10x didn't ship truncated rings)
  - 6-panel histogram PNG

- `xenium_segmentation_plot.py` — per-sample CLI
  - bbox-cropped or full-FOV rendering
  - color polygons by any cells.parquet column (numeric → `viridis` w/ optional log; categorical → `tab20`)
  - overlay nucleus boundaries as thin black outline
  - non-interactive Agg backend, deterministic PNG names

## Typical Invocations

Activate the project's spatial env first:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate sc_spatial
```

Per-sample QC:

```bash
python /path/to/project/xenium_segmentation_qc.py --samples <sample_a> <sample_b>
```

Cropped polygon render (300 × 300 µm), color by transcript counts on log scale, with nucleus overlay:

```bash
python /path/to/project/xenium_segmentation_plot.py \
  --sample <sample_id> \
  --bbox 5400 5700 9600 9900 \
  --color-by total_counts --log-color \
  --overlay-nucleus
```

Same crop, color by segmentation_method (categorical legend):

```bash
python /path/to/project/xenium_segmentation_plot.py \
  --sample <sample_id> \
  --bbox 5400 5700 9600 9900 \
  --color-by segmentation_method
```

## QC Interpretation

| Signal | What to flag |
| --- | --- |
| `frac_nucleus_expansion > 0.1` | weak interior + boundary staining; downstream celltype assignments will lean on nucleus territory only |
| `frac_anucleate > 0.05` | possible over-segmentation of cytoplasmic territories without DAPI support; consider filtering for downstream phenotyping |
| `frac_multinucleated > 0.05` | tissue clumps / undersegmentation; check `nucleus_count` per region |
| `nuc/cell area ratio` peak ≈ 1.0 | thin or absent cytoplasm; filter or label as "near-nucleus only" |
| `cell_area p99 > 5 × p50` | a few outsized polygons (large cells, syncytia, or merge artifacts) |

## Output Contract for Downstream Use

- **Polygons → AnnData uns**: `adata.uns["xenium_polygons"] = {"cell": dict[cell_id->(N,2)], "nucleus": dict[cell_id->(N,2)]}` is acceptable for in-memory passes; do NOT serialize to h5ad (slow). Persist to parquet/GeoParquet instead.
- **Polygons → GeoJSON / GeoParquet**: use `polygons_to_geodataframe(...).to_parquet(...)` or `.to_file(..., driver="GeoJSON")` for QuPath / napari handoff.
- **Masks**: rasterize with `rasterio.features.rasterize` over a target image grid; persist as compressed `uint32` TIFF keyed by integer cell index (build `cell_id ↔ int` table separately).

## Validation Checklist (before downstream)

- `len(load_cells_table(...)) == cells.parquet num_rows`
- every `segmentation_method` value matches one of the three documented categories (else: log-warn, do not silently bucket)
- `cells_with_polygon = set(load_polygons(cell_b_pq).keys())` is a subset of `cells.parquet["cell_id"]` (Xenium guarantees this; assert)
- the difference `cells.parquet["cell_id"] - cells_with_polygon` is empty for `cell_boundaries.parquet`; for nuclei, it equals the anucleate set
- `bbox`-cropped renders show cells filling the box and y-axis is inverted (anatomy not mirrored)

## Out of Scope

- transcripts.parquet → polygon assignment (already done by Xenium ranger; revisit only for custom resegmentation)
- alignment to H&E / IF (separate skill — needs the per-Xenium transform matrix)
- legacy CSV.gz boundary loading (write a one-shot converter; don't mix in the new path)
