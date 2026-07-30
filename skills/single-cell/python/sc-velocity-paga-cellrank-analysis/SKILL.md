---
name: sc-velocity-paga-cellrank-analysis
description: Use when running or reviewing Python single-cell trajectory workflows with scVelo, PAGA, or CellRank, including loom-to-AnnData alignment, spliced/unspliced layer checks, velocity graph construction, latent time, PAGA topology, CellRank fate probabilities, and trajectory output validation.
---

# Skill: scVelo PAGA CellRank Analysis

## Purpose

Use this skill for RNA velocity and graph-based fate analysis in Python. It covers scVelo, Scanpy PAGA, and CellRank workflows from prepared AnnData/loom inputs.

## Environment

Use a dedicated environment because scVelo/CellRank dependency versions can conflict with modern Scanpy stacks:

```bash
conda create -n scvelo_analysis python=3.10
conda activate scvelo_analysis
pip install scanpy scvelo cellrank loompy
```

If the source scripts request an older Python, pin only when required by reproducibility or package compatibility. Use `conda` and `pip` only.

## Inputs

Confirm:

- AnnData `.h5ad` has filtered cells and trusted metadata/cluster labels.
- Loom or AnnData velocity object contains `spliced`, `unspliced`, and ideally `ambiguous` layers.
- Cell barcodes are aligned between expression and velocity objects.
- Gene IDs/symbols are in the same namespace.
- If using a reduced "mini AnnData", it preserves `obs_names`, `var_names`, layers, `obs` annotations, and embedding coordinates.

## Loom / AnnData Alignment

Never assume barcode formats match. Normalize barcodes explicitly and audit overlap:

```python
adata = sc.read_h5ad("annotated.h5ad")
ldata = scv.read("sample.loom", cache=True)

shared = adata.obs_names.intersection(ldata.obs_names)
if len(shared) == 0:
    raise ValueError("No shared cells between AnnData and loom; inspect barcode formatting.")

adata = adata[shared].copy()
ldata = ldata[shared, adata.var_names.intersection(ldata.var_names)].copy()
```

For multiple samples, process each sample's loom with sample-aware barcode prefixes before concatenating.

## scVelo Workflow

```python
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata)
scv.tl.velocity_embedding(adata, basis="umap")
scv.tl.latent_time(adata)
```

Use `mode="dynamical"` only after `scv.tl.recover_dynamics` succeeds on a tractable subset or final gene set.

## PAGA Workflow

Use PAGA when the question is topology among clusters/states:

```python
sc.tl.paga(adata, groups="cell_type")
sc.pl.paga(adata, threshold=0.03)
sc.tl.draw_graph(adata, init_pos="paga")
```

If PAGA is combined with velocity, compare PAGA edges with velocity stream direction instead of treating all graph edges as directional.

## CellRank Workflow

Use CellRank when terminal states, lineage drivers, or fate probabilities are requested:

```python
vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = 0.8 * vk + 0.2 * ck

est = cr.estimators.GPCCA(kernel)
est.fit(cluster_key="cell_type")
est.predict_terminal_states()
est.compute_fate_probabilities()
```

Record terminal-state selection rules and avoid manually forcing states without biological justification.

## Output Checks

- layers `spliced` and `unspliced` are present and nonempty.
- cell/gene overlap report is saved.
- `velocity_graph` exists and matches `n_obs`.
- velocity plots use the intended embedding (`X_umap`, `X_draw_graph_fa`, etc.).
- PAGA group key exists and has enough cells per group.
- CellRank fate probabilities are present and terminal states were recorded.
- final `.h5ad` preserves velocity layers and run parameters.

## Boundaries

- RNA velocity needs unspliced/spliced data; do not run it on ordinary expression-only h5ad.
- PAGA alone is topology, not velocity.
- CellRank fate probabilities are model outputs and should be validated against marker genes and known biology.
