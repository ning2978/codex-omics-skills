---
name: visiumhd-cellcharter-visualization
description: Use when plotting, summarizing, or regenerating 10x Visium HD CellCharter niche/cluster outputs, especially mapping CellCharter labels from h5ad back onto H&E images, creating spatial overview plots, bin-to-HE coordinate CSVs, per-sample niche count/fraction tables, and validating nonblank high-DPI PNG outputs.
---

# Skill: Visium HD CellCharter Visualization

## Purpose

Use this workflow skill with `algorithm-cellcharter` for post-hoc visualization of 10x Visium HD CellCharter results. It assumes niches/clusters already exist in an `.h5ad`; it does not train scVI, run CellCharter clustering, or change labels.

## Input Contract

Require an AnnData object with:

- `obs[sample_key]`: sample or library id for each bin.
- `obs[cluster_key]`: CellCharter niche/cluster labels.
- `obsm["spatial"]`: full-resolution pixel coordinates as `[pxl_col_in_fullres, pxl_row_in_fullres]`.
- `uns["spatial"][sample_id]["images"]["hires"|"lowres"]`: embedded 10x H&E image arrays.
- `uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"|"tissue_lowres_scalef"]`: scale from fullres coordinates to the selected image.

If the HE images are not embedded, add them to `uns["spatial"]` from Space Ranger `spatial/tissue_hires_image.png`, `tissue_lowres_image.png`, and `scalefactors_json.json` before plotting.

## Environment

Prefer the CellCharter environment for scripts that read `.h5ad` with the project stack:

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mpl \
python
```

Use server-side plotting only: set `matplotlib.use("Agg")`, write PNG/CSV files, and validate with PIL rather than visual inspection alone.

## Plotting Defaults

- For a compact HE handoff, generate one all-sample overview PNG plus CSV tables; do not generate per-sample or per-niche image grids unless explicitly requested.
- For Visium HD bins over HE, use high DPI and small-but-visible points. A good default is `dpi=360`, `spot_size=0.9`, `alpha=0.82`.
- Preserve all cluster labels as strings and sort numeric labels numerically (`0, 1, 2, ... 10`) for palettes and tables.
- Use `tab20` for up to 20 niches; switch to a continuous categorical fallback only when there are more than 20 labels.
- Always export:
  - a bin-to-HE coordinate table with original fullres coordinates, selected image scale, and scaled HE coordinates;
  - a niche count table by sample;
  - a parameters JSON recording input path, keys, image key, DPI, spot size, and generated outputs.

## Workflow

1. Confirm the h5ad contains `sample_key`, `cluster_key`, `obsm["spatial"]`, and sample-keyed `uns["spatial"]` images/scalefactors.
2. Build a plotting table with `obs_name`, `sample_id`, `niche`, `fullres_col`, `fullres_row`, `he_col`, and `he_row`.
3. Draw the selected HE image per sample, then scatter bins at `he_col/he_row` colored by niche.
4. For overview plots, use a 2x2 grid for four samples; keep the legend outside the panels.
5. When changing output mode, remove stale plot subdirectories first so old per-niche images do not look current.
6. Never modify the source h5ad during visualization unless the user explicitly asks to embed missing images.

## Bundled Script

Use the bundled generic mapper when the repository does not already have a project-specific plotting script:

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mpl \
python \
<skill_dir>/scripts/map_visiumhd_cellcharter_niches_to_he.py \
  --input-h5ad path/to/cellcharter_result.h5ad \
  --outdir path/to/plots/he_niche_mapping \
  --sample-key sample_id \
  --cluster-key cluster_cellcharter \
  --output-prefix visiumhd_cellcharter \
  --overwrite
```

In repositories with existing scripts, reuse local code if it already follows the same input contract and validation rules.

## Validation

After editing plotting scripts, run:

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mpl \
python -m py_compile <script.py>

NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mpl \
python <script.py> --help
```

After plotting, verify:

- expected PNG/CSV/JSON count and no stale unwanted directories;
- coordinate table row count equals `adata.n_obs`;
- observed niche labels match `obs[cluster_key]`;
- output PNG opens with PIL, has nonzero file size, expected dimensions, and non-flat RGB extrema.

## References

Read `references/visiumhd-he-coordinate-contract.md` when implementing or debugging coordinate mapping from fullres Visium HD positions to hires/lowres H&E images.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
