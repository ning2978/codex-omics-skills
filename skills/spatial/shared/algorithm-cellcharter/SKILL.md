---
name: algorithm-cellcharter
description: Use when running, reviewing, fixing, or packaging CellCharter spatial niche workflows for Visium-like spatial transcriptomics, including CellCharter API choices, spatial graph construction, neighborhood aggregation, Cluster/ClusterAutoK settings, scVI/RCTD feature inputs, and CellCharter output validation.
---

# Skill: algorithm-cellcharter

## Purpose

Use this skill for CellCharter-specific decisions. Pair it with a project or workflow skill for dataset paths and sample metadata.

## Environment

- Prefer a dedicated CellCharter environment, not the default scanpy environment:
  `NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mpl python`
- If a command only inspects finished `.h5ad` or CSV outputs and does not import `cellcharter`, the scanpy environment may be enough.
- Always run at least `py_compile` and `--help` for edited CellCharter scripts in the CellCharter environment.

## API Pattern

Use this current API shape when available:

```python
sq.gr.spatial_neighbors(..., library_key="sample", spatial_key="spatial")
cc.gr.remove_long_links(adata, percentile=99.0)
cc.gr.aggregate_neighbors(
    adata,
    n_layers=n_layers,
    use_rep=feature_key,
    sample_key="sample",
    out_key=cellcharter_rep_key,
)
clusterer = cc.tl.Cluster(n_clusters=k, random_state=seed, trainer_params=...)
clusterer.fit(adata, use_rep=cellcharter_rep_key)
labels = clusterer.predict(adata, use_rep=cellcharter_rep_key)
```

For automatic K, use `cc.tl.ClusterAutoK(...)`, then save both long-form and summarized stability outputs when available.

Avoid older or unverified call shapes such as `cc.tl.cluster.auto_k`, `cc.tl.nhood_annotation`, and `cc.pl.nhood_heatmap`.

## Input Representations

- Expression/scVI CellCharter: raw counts -> normalized/log expression -> HVGs -> scVI latent -> `X_cellcharter`.
- RCTD cell-type CellCharter: RCTD weights or abundance columns -> row-normalized proportions -> `clr_zscore` -> `X_celltype` -> `X_cellcharter_celltype`.
- Hybrid CellCharter: `[zscore(X_scVI) * scvi_weight, transform(RCTD weights) * rctd_weight]` -> `X_scVI_RCTD_hybrid` -> `X_cellcharter_hybrid`.

Use cell-type composition workflows when the question is about mixtures of annotated cell types forming niches. Use expression/scVI workflows when the question is about transcriptomic spatial domains.

Select the expression/scVI route only when the user explicitly requests
CellCharter latent mode. Do not introduce scVI as a default Xenium expression
refinement or label-transfer step. CellCharter clustering operates on
aggregated spatial features and is not subject to the expression-reclustering
UMAP gate.

## Spatial Graph Defaults

- Multi-sample Visium must be sample-aware: set `library_key`/`sample_key` so spots from different slides are never connected.
- For image pixel coordinates, prefer `coord_type="generic"`, `delaunay=True`, then remove long links around the 99th percentile.
- For canonical array/grid coordinates, use `coord_type="grid"` and set `n_neighs` explicitly.
- Log graph edge counts and cross-sample edge checks after graph construction.

## Output Checks

After a run, confirm:

- result `.h5ad` exists and contains the cluster key in `obs`
- `obsm[cellcharter_rep_key]` has expected dimensions
- `spatial_connectivities` exists and has zero cross-sample edges
- cluster labels are ordered categoricals or stable strings
- `run_parameters.json` and `run_summary.csv` record K, transform, graph, accelerator, and input paths
- spatial plots and niche composition tables were written

## References

Read `references/repository-api-notes.md` when modifying or reviewing CellCharter scripts.
