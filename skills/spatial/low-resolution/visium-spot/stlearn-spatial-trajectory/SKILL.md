---
name: stlearn-spatial-trajectory
description: Use when running or reviewing stLearn workflows for 10x Visium spatial trajectory, morphology-adjusted PCA, SME normalization, louvain clustering, pseudotime root selection, pseudotimespace trajectories, transition marker detection, stLearn ligand-receptor CCI, or stLearn CCI maps.
---

# Skill: stLearn Spatial Trajectory

## Purpose

Use this skill for stLearn-based Visium workflows, especially when histology-informed morphology features, spatial trajectory inference, or stLearn ligand-receptor CCI are required.

## Environment

stLearn is version-sensitive. Prefer a dedicated environment:

```bash
conda create -n stlearn_spatial python=3.9
conda activate stlearn_spatial
pip install stlearn scanpy
```

If CCI GO enrichment is required, document the R path and required Bioconductor packages separately.

## Input Contract

Confirm:

- 10x Visium-style directory exists with matrix and `spatial/`.
- high-resolution image is available when using morphology or CCI plots over image.
- `library_id` and `count_file` are correct.
- raw counts are preserved before normalization.

```python
data = st.Read10X(path=visium_dir, library_id=library_id, count_file="filtered_feature_bc_matrix.h5")
data.layers["raw_count"] = data.X.copy()
data.var_names_make_unique()
```

## Trajectory Workflow

1. **Preprocess**

```python
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.raw = data
st.pp.scale(data)
```

2. **Morphology-adjusted embedding**

```python
st.em.run_pca(data, n_comps=50, random_state=0)
st.pp.tiling(data, out_path="tiles", crop_size=40)
st.pp.extract_feature(data)
st.spatial.morphology.adjust(data, use_data="X_pca", radius=50, method="mean")
st.pp.neighbors(data, n_neighbors=25, use_rep="X_pca_morphology", random_state=0)
st.tl.clustering.louvain(data, random_state=0)
```

3. **Choose a root**

Pick a biologically defensible root cluster from histology, markers, or cluster maps. Do not blindly reuse demo cluster numbers.

```python
data.uns["iroot"] = st.spatial.trajectory.set_root(data, use_label="louvain", cluster=root_cluster, use_raw=True)
st.spatial.trajectory.pseudotime(data, eps=50, use_rep="X_pca", use_label="louvain")
```

4. **Run global trajectory**

```python
st.spatial.trajectory.pseudotimespace_global(data, use_label="louvain", list_clusters=target_clusters)
st.pl.cluster_plot(data, use_label="louvain", show_trajectories=True, list_clusters=target_clusters)
st.pl.trajectory.tree_plot(data)
```

5. **Transition markers**

```python
st.spatial.trajectory.detect_transition_markers_clades(data, clade=clade_id, use_raw_count=False, cutoff_spearman=0.4)
```

Remove ribosomal genes from ranked transition markers before interpretation when they dominate:

```python
data.uns[key] = data.uns[key][~data.uns[key]["gene"].str.upper().str.startswith(("RPL", "RPS"))]
```

## stLearn CCI Workflow

Use this only after spot cell-type labels or mixtures are available:

```python
data.obs["cell_type"] = labels.astype("category")
data.uns["cell_type"] = spot_mixtures
lrs = st.tl.cci.load_lrs(["connectomeDB2020_lit"], species="human")
st.tl.cci.run(data, lrs, min_spots=20, distance=None, n_pairs=10000, n_cpus=4)
st.tl.cci.adj_pvals(data, correct_axis="spot", pval_adj_cutoff=0.05, adj_method="fdr_bh")
st.tl.cci.run_cci(data, "cell_type", min_spots=3, spot_mixtures=True, cell_prop_cutoff=0.2, sig_spots=True, n_perms=1000)
```

Start with low permutations for smoke tests only; use higher values for final analysis.

## Output Checks

- `data.layers['raw_count']` exists.
- morphology embedding exists when clustering used `X_pca_morphology`.
- root cluster choice and target trajectory clusters are recorded.
- pseudotime is visible on both non-spatial and spatial plots.
- transition marker tables exist and ribosomal dominance was checked.
- CCI results include `data.uns['lr_summary']`, adjusted p-values, LR maps, and cell-type interaction maps when requested.

## Boundaries

- stLearn trajectory is spatial-pseudotime inference, not lineage proof.
- For RNA velocity or fate probabilities, use `sc-velocity-paga-cellrank-analysis`.
- For COMMOT-specific communication, use `algorithm-commot`.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
