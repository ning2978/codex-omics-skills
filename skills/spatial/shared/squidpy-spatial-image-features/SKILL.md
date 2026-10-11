---
name: squidpy-spatial-image-features
description: Use when working with Squidpy spatial image containers, Visium or multiplex image metadata in AnnData, spatial graph construction, image cropping, segmentation layers, image feature extraction, Moran spatial autocorrelation, co-occurrence, neighborhood enrichment, centrality, Ripley statistics, or Squidpy ligand-receptor plots.
---

# Skill: Squidpy Spatial Image Features

## Purpose

Use this skill for Squidpy-based spatial image handling and general spatial statistics that are not already covered by a more specific CellCharter, Xenium, Visium HD, or proximity skill.

## Environment

- Use a Scanpy/Squidpy environment; keep Napari/GUI dependencies separate on servers.
- Use `conda` and `pip` only:

```bash
conda create -n squidpy_spatial python=3.10
conda activate squidpy_spatial
pip install scanpy squidpy scikit-image
```

Only install `napari`, `pyqt`, or `napari-spatialdata` when the user explicitly needs an interactive GUI.

## AnnData Spatial Contract

For spatial plotting, confirm:

- `adata.obsm['spatial']` exists and is in the coordinate system expected by images.
- `adata.uns['spatial'][library_id]['images']` contains image layers such as `hires` or `segmentation`.
- `adata.uns['spatial'][library_id]['scalefactors']` includes `tissue_hires_scalef` and `spot_diameter_fullres` when using Visium-style images.
- multi-library objects have a `library_id`/`library_key` column and do not mix image keys across samples.

## ImageContainer Workflow

Use `sq.im.ImageContainer` when manipulating image layers:

```python
img = sq.im.ImageContainer(adata.uns["spatial"][library_id]["images"]["hires"], library_id=library_id, layer="image")
img.add_img(adata.uns["spatial"][library_id]["images"]["segmentation"], library_id=library_id, layer="segmentation")
img["segmentation"].attrs["segmentation"] = True
```

Typical operations:

- `img.show(layer="image")`
- `img.show("image", segmentation_layer="segmentation", segmentation_alpha=0.5)`
- `img.crop_center(...)` or `img.crop_corner(...)`
- `adata_crop = crop.subset(adata)`
- `sq.im.ImageContainer.concat(imgs)` for multiple libraries

## Image Feature Extraction

Use image features when histology or multiplex image morphology is part of the analysis:

```python
sq.im.calculate_image_features(
    adata,
    img,
    features="summary",
    key_added="summary_features",
    mask_circle=True,
)
```

Common feature families:

- `summary`: channel quantiles/intensity summaries.
- `histogram`: binned channel intensities.
- `texture`: contrast and related texture features.
- `segmentation`: label count, area, centroid, and per-segment features.

For segmentation:

```python
sq.im.process(crop, layer="image", method="smooth", sigma=4)
sq.im.segment(crop, layer="image_smooth", method="watershed", channel=0)
sq.im.calculate_image_features(
    adata_crop,
    crop,
    layer="image",
    features="segmentation",
    key_added="segmentation_features",
    features_kwargs={"segmentation": {"label_layer": "segmented_watershed", "props": ["label", "area", "centroid"]}},
    mask_circle=True,
)
```

## Spatial Statistics

Build spatial neighbors with an explicit coordinate type:

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```

Use sample-aware construction when working across multiple libraries. For high-resolution or cell-level coordinates, inspect graph edge lengths before interpreting results.

Common analyses:

- `sq.gr.centrality_scores(adata, cluster_key=...)`
- `sq.gr.co_occurrence(adata_subsample, cluster_key=...)`
- `sq.gr.nhood_enrichment(adata, cluster_key=...)`
- `sq.gr.ripley(adata, cluster_key=..., mode="L")`
- `sq.gr.spatial_autocorr(adata, mode="moran", n_perms=...)`
- `sq.gr.ligrec(adata, cluster_key=..., clusters=...)`

Subsample large cell-level data for expensive permutation or co-occurrence runs, and record the subsampling fraction.

## Output Checks

- image layers align visually with coordinates.
- cropped AnnData contains only spots/cells inside the crop.
- image feature tables are in `adata.obsm[key_added]` with named columns.
- graph statistics were run on a valid `spatial_connectivities` graph.
- plots use consistent `library_id`, image key, and coordinate basis.

## Boundaries

- For anchor-vs-neighbor deconvolution questions, use `visium-spot-squidpy-proximity`.
- For Xenium-specific segmentation parquet or ROI QC, use the Xenium skills.
- Do not run Napari on a headless server unless the user specifically asks for GUI setup.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
