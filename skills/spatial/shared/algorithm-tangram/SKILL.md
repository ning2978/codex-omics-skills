---
name: algorithm-tangram
description: Use when mapping single-cell RNA-seq references onto spatial transcriptomics with Tangram, including marker-gene selection, `tg.pp_adatas`, cell/cluster mapping modes, annotation projection, gene projection validation, segmentation-aware constrained mapping, and Tangram output checks.
---

# Skill: algorithm-tangram

## Purpose

Use this skill for Tangram-based scRNA-to-ST alignment and projection. Pair it with a project preparation skill when raw Visium, Xenium, or other ST inputs still need normalization, QC, or metadata cleanup.

## Environment

- Prefer a dedicated Tangram environment because `tangram-sc`, `scanpy`, `squidpy`, `torch`, and `scikit-image` versions can conflict with CellCharter/cell2location environments.
- Use `conda` and `pip` only. Example:

```bash
conda create -n tangram_spatial python=3.10
conda activate tangram_spatial
pip install tangram-sc scanpy squidpy scikit-image
```

- Choose GPU only after a small CPU smoke test proves that shared genes and AnnData keys are valid.

## Inputs

Confirm these before mapping:

- `adata_sc`: scRNA AnnData with raw or log-normalized expression and a trusted annotation column.
- `adata_st`: spatial AnnData with `obsm['spatial']`, spatial images if plotting is needed, and expression in the same gene namespace.
- Shared genes: use symbols or stable IDs consistently. Do not silently mix Ensembl IDs and symbols.
- Annotation key: e.g. `cell_type`, `cell_subclass`, or `cluster`.
- Optional image segmentation / cell-count prior for constrained deconvolution.

## Canonical Workflow

1. **Select training genes**
   - Prefer marker genes from `sc.tl.rank_genes_groups(adata_sc, groupby=annotation_key)`.
   - Use the top N per annotation, deduplicate, then intersect with spatial genes.
   - Avoid using all genes by default; it is slower and can dilute cell-type signal.

2. **Prepare AnnData objects**

```python
tg.pp_adatas(adata_sc, adata_st, genes=markers)
```

3. **Map cells or clusters**
   - Use `mode='cells'` for cell-level mapping and annotation projection.
   - Use `mode='clusters'` with `cluster_label=annotation_key` when the reference is large or the question is cell-type-level.
   - Use `mode='constrained'` only when a defensible target count or density prior exists.

4. **Project annotations**

```python
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_st,
    mode="cells",
    density_prior="rna_count_based",
    num_epochs=500,
    device="cpu",
)
tg.project_cell_annotations(ad_map, adata_st, annotation=annotation_key)
```

5. **Project and validate genes**

```python
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
df_gene_score = tg.compare_spatial_geneexp(ad_ge, adata_st, adata_sc)
```

Review training scores and predicted-vs-measured spatial expression before interpreting projected cell types.

## Segmentation-Aware Mapping

Use constrained mapping only when image-derived cell counts are available per spot/bin:

```python
sq.im.calculate_image_features(
    adata_st,
    img,
    layer="image",
    features="segmentation",
    key_added="image_features",
    features_kwargs={
        "segmentation": {
            "label_layer": "segmented_watershed",
            "props": ["label", "centroid"],
            "channels": [1, 2],
        }
    },
    mask_circle=True,
)
adata_st.obs["cell_count"] = adata_st.obsm["image_features"]["segmentation_label"]
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_st,
    mode="constrained",
    target_count=adata_st.obs["cell_count"].sum(),
    density_prior=adata_st.obs["cell_count"].to_numpy() / adata_st.obs["cell_count"].sum(),
)
```

If Tangram fails when `density_prior` is a NumPy array, check the installed `tangram-sc` version before patching package internals. Prefer upgrading or pinning a known-compatible version over editing site-packages in place.

## Output Checks

Before reporting results:

- `ad_map` exists and contains mapping probabilities.
- `adata_st.obsm` has projected annotation scores or counts.
- projected annotations sum to plausible spot/bin totals.
- `ad_ge` predicted expression was compared against measured ST expression.
- at least one diagnostic plot was generated for reference annotation, projected annotation, and selected marker genes.
- constrained runs record target count, density prior source, and segmentation parameters.

## Boundaries

- Tangram projection is not direct experimental co-localization; treat it as model-based mapping.
- Do not use Tangram as a replacement for cell2location/RCTD when the task specifically needs abundance uncertainty or robust deconvolution benchmarks.
- Do not copy course paths or demo marker lists into production scripts.
