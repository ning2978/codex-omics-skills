---
name: algorithm-commot
description: Use when running or reviewing COMMOT spatial communication analysis in Python, including ligand-receptor database filtering, `ct.tl.spatial_communication`, sender/receiver pathway summaries, communication direction, cluster communication networks, communication-dependent genes, and output validation.
---

# Skill: algorithm-commot

## Purpose

Use this skill for COMMOT ligand-receptor communication on spatial transcriptomics data. It is best for spatially constrained signaling questions where spot/cell coordinates are part of the model.

## Environment

- Prefer a dedicated COMMOT environment.
- Use `conda` and `pip` only:

```bash
conda create -n commot_spatial python=3.10
conda activate commot_spatial
pip install commot scanpy squidpy
conda install -n commot_spatial -c conda-forge pygraphviz
```

- Some downstream DEG functions may need R packages through `rpy2`; document those separately and do not hide R dependencies inside a "pure Python" run.

## Inputs

Required:

- AnnData with expression, `obsm['spatial']`, and spot/cell metadata.
- raw counts preserved in `layers['counts']` or `raw` if DEG steps are planned.
- clustering column, commonly `leiden` or a curated cell-type label.
- species and ligand-receptor database choice, e.g. CellChat human/mouse secreted signaling.

Check that the coordinate scale matches the communication distance threshold. `dis_thr=500` means different things for pixels, microns, and array coordinates.

## Canonical Workflow

1. **Prepare expression**

```python
adata.raw = adata
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
```

Keep a full-expression copy for COMMOT and a reduced copy for clustering if needed.

2. **Cluster if no trusted labels exist**

```python
sc.pp.highly_variable_genes(adata_work)
sc.tl.pca(adata_work)
sc.pp.neighbors(adata_work, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_work)
sc.tl.leiden(adata_work, resolution=0.4)
adata.obs["leiden"] = adata_work.obs["leiden"]
```

3. **Build and filter LR database**

```python
df_lr = ct.pp.ligand_receptor_database(
    species="mouse",
    signaling_type="Secreted Signaling",
    database="CellChat",
)
df_lr = ct.pp.filter_lr_database(df_lr, adata, min_cell_pct=0.05)
```

4. **Run spatial communication**

```python
ct.tl.spatial_communication(
    adata,
    database_name="cellchat",
    df_ligrec=df_lr,
    dis_thr=500,
    heteromeric=True,
    pathway_sum=True,
)
```

5. **Summarize a pathway**

```python
ct.tl.communication_direction(adata, database_name="cellchat", pathway_name="PSAP", k=5)
ct.pl.plot_cell_communication(
    adata,
    database_name="cellchat",
    pathway_name="PSAP",
    plot_method="grid",
    summary="sender",
    background="image",
    clustering="leiden",
)
```

6. **Cluster-level network**

```python
ct.tl.cluster_communication(
    adata,
    database_name="cellchat",
    pathway_name="PSAP",
    clustering="leiden",
    n_permutations=100,
)
ct.tl.cluster_position(adata, clustering="leiden")
```

## Optional DEG / Impact Workflow

Use only after the core communication result is stable:

- `ct.tl.communication_deg_detection`
- `ct.tl.communication_deg_clustering`
- `ct.pl.plot_communication_dependent_genes`
- `ct.tl.communication_impact`
- `ct.pl.plot_communication_impact`

Record whether DEG uses receiver or sender summaries and preserve the genes used as background.

## Output Checks

- `adata.uns` has COMMOT result keys for the chosen database/pathway.
- `adata.obsm` has sender/receiver summary matrices, e.g. pathway-specific received/sent signal.
- cluster communication outputs use the same clustering column as spatial plots.
- network plots were regenerated after `pygraphviz` installation if needed.
- distance threshold, species, database, and `min_cell_pct` are written into run metadata.

## Boundaries

- COMMOT reports modeled communication potential, not direct proof of signaling.
- Compare results against expression maps for ligand and receptor genes.
- Use stLearn CCI or Squidpy `ligrec` as complementary methods, not interchangeable drop-in replacements.
