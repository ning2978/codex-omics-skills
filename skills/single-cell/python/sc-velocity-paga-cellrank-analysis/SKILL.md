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
ldata = sc.read_loom("sample.loom")
# First apply an audited capture-aware key mapping to both objects.
# sc_contracts is the shared workflow's resolved scripts/sc_contracts.py module.
from sc_contracts import merge_velocity_layers, validate_velocity
adata_velocity = merge_velocity_layers(adata, ldata, allow_gene_subset=False)
```

For multiple captures use library-aware keys before concatenation. Require unique stable gene IDs, not duplicate symbols. The helper reorders both axes and actually copies velocity layers. Missing analyzed cells fail; use an explicit chosen subset/attrition ledger rather than silent intersection. Gene restriction needs an explicit derivative and omitted-gene report. Preserve original expression/counts object: velocity normalization modifies a working derivative.

## scVelo Workflow

```python
validate_velocity(adata_velocity, "stochastic")
scv.pp.filter_and_normalize(adata_velocity, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_velocity, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_velocity, mode="stochastic")
scv.tl.velocity_graph(adata_velocity)
scv.tl.velocity_embedding(adata_velocity, basis="umap")
scv.tl.velocity_pseudotime(adata_velocity)  # not dynamical latent time
```

For **dynamical** mode run `scv.tl.recover_dynamics(adata_velocity)` after preparation, inspect fitted genes/convergence and then `validate_velocity(..., "dynamical")`, `velocity(mode="dynamical")`, `velocity_graph`, and `latent_time`. Never call `latent_time` in a stochastic-only run. Record fitted gene count/failed fits, kinetic assumptions, preprocessing/graph source, root/terminal rules and exact scVelo/CellRank versions (API checks validated with 0.3.3/2.0.7). PAGA/CellRank examples below assume `adata` is the validated velocity derivative with an appropriate graph; PAGA alone can use an expression graph without velocity layers.

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

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
