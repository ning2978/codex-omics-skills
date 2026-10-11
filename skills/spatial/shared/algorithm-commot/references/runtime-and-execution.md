---
name: algorithm-commot
description: Run, review, diagnose, or optimize COMMOT spatial communication analysis in Python, including ligand-receptor filtering, `ct.tl.spatial_communication`, sender/receiver summaries, communication direction, cluster networks, communication-dependent genes, CPU/GPU capability checks, safe multi-process execution, reproducibility, and output validation.
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

## Runtime Diagnosis

Inspect the exact interpreter and installed source before promising acceleration:

1. Record the COMMOT, NumPy, SciPy, POT, and threadpool versions.
2. Inspect `ct.tl.spatial_communication` and downstream function signatures for
   `n_jobs`, device, or backend arguments.
3. Inspect the active process with `%CPU`, `NLWP`, per-thread CPU, affinity, and
   elapsed time. A wide CPU affinity does not mean the code is parallel.
4. Inspect `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS`, then
   confirm active math-library limits with `threadpoolctl`.
5. Check GPU process PIDs, not only global GPU utilization; another user's job
   can make a CPU-only COMMOT run appear GPU-active.

For COMMOT 0.0.3, `spatial_communication` has no `n_jobs`, device, or backend
parameter. Its collective optimal-transport implementation uses Python,
NumPy, and SciPy sparse operations and serial loops. Increasing BLAS thread
variables therefore does not make the main solver multi-core. POT reporting
only `NumpyBackend` is direct evidence that this path is CPU-only.

Do not infer GPU support from an NVIDIA device, a working driver, a
CUDA-enabled framework in another environment, or POT's general ability to
support other backends. Treat a GPU port or alternative implementation as a
new method requiring matrix-level and biological revalidation, not a runtime
flag for COMMOT 0.0.3.

## Safe CPU Parallelism

Parallelize independent samples, libraries, or ROIs at the process level when
all of these conditions hold:

- each task builds its own spatial graph and has zero permitted cross-boundary
  edges;
- each task writes unique cache and metadata paths;
- the parent reassembles outputs in a deterministic declared order;
- resume logic accepts only a complete result-plus-metadata pair;
- peak memory has been measured on a representative large boundary.

Do not split ligand-receptor pairs inside one collective optimal-transport
model: COMMOT jointly couples ligand and receptor distributions, so those
subproblems are not generally independent.

Before starting Python workers, set a fixed `PYTHONHASHSEED`. COMMOT 0.0.3
constructs ligand and receptor orders from Python sets, so different process
hash seeds can change ordering and floating-point results. Also set math
libraries to one thread per worker before spawned workers import NumPy/SciPy:

```bash
PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python run_commot.py --workers 4
```

Prefer a spawn-based process pool and avoid nested process/thread parallelism.
Use 4 independent workers as the normal starting point. When the available CPU
allocation explicitly permits it and the measured memory gate passes,
process-level parallelism may scale to at most 8 independent workers. Eight is
a hard cap, not a target: automatically downshift below 4 when CPU affinity,
independent task count, or measured memory requires it. Bound workers by both
independent task count and memory:

`workers <= floor(0.7 * available_memory / measured_peak_worker_memory)`

Scaling above four workers does not relax any reproducibility requirement:
use spawn-isolated workers, retain one math thread per worker and a fixed
`PYTHONHASHSEED`, and require the deterministic serial/parallel sparse-output
probe to pass before the larger pool is used.

Use the largest-ROI pilot when the distance matrix is dense because memory can
grow approximately quadratically with cells per boundary. Record worker count,
threads per worker, hash seed, task partition, peak RSS, and elapsed time in
the run manifest.

Run the bundled deterministic synthetic probe in the exact COMMOT environment
before enabling process parallelism:

```bash
python <skill-dir>/scripts/runtime_probe.py \
  --cells 250 --lr-pairs 3 --tasks 4 --workers 4 \
  --threads-per-worker 1 --python-hash-seed 0 \
  --output /path/to/commot_runtime_probe.json
```

Require PASS, distinct parallel worker PIDs, one active math thread per worker,
and exact sparse-output hashes. Treat the measured speedup as host- and
workload-specific. The synthetic probe does not prove full-data memory safety.

## Inputs

Required:

- AnnData with expression, `obsm['spatial']`, and spot/cell metadata.
- raw counts preserved in `layers['counts']` or `raw` if DEG steps are planned.
- trusted grouping label only when grouped summaries are requested; no clustering is required for a cell-level model.
- species and ligand-receptor database choice, e.g. CellChat human/mouse secreted signaling.

Check that the coordinate scale matches the communication distance threshold. `dis_thr=500` means different things for pixels, microns, and array coordinates.

## Canonical Workflow

1. **Use explicit prepared expression**

Read input-and-panel.md. Load a prepared per-ROI working H5AD and its audited LR table; preserve the source counts and matrix contract. No implicit normalization, species choice or expression imputation is allowed. Geometry and threshold use the same declared unit.

2. **Choose summary labels only if needed**

Cell-level COMMOT can run without clusters. Use existing trusted labels for grouped summaries; missing labels do not automatically request new expression clusters. A separately requested expression recluster follows the counts-first platform contract.

3. **Build and filter LR database**

```python
df_lr = ct.pp.ligand_receptor_database(
    species=declared_species,
    signaling_type="Secreted Signaling",
    database=declared_database,
)
df_lr = ct.pp.filter_lr_database(df_lr, adata, min_cell_pct=0.05)
```

4. **Run spatial communication**

```python
ct.tl.spatial_communication(
    adata,
    database_name="cellchat",
    df_ligrec=df_lr,
    dis_thr=declared_distance_threshold,
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
- multi-sample or multi-ROI runs have no cross-boundary edges unless explicitly modeled;
- parallel and serial pilots agree exactly after fixing the hash seed, or agree
  within a predeclared numerical tolerance when exact agreement is impossible;
- runtime metadata records the exact interpreter, COMMOT/POT backend, worker
  count, per-worker thread limit, hash seed, task partition, peak memory, and
  resume status.

## Completion Gate

Do not declare an optimized run ready until:

- a serial baseline is valid;
- the runtime probe passes in the exact environment;
- a representative real-data pilot reproduces the serial matrices and summary
  tables;
- the worker count passes the largest-boundary memory estimate;
- interrupted-task and resume behavior has been tested without duplicate or
  partial cache acceptance;
- full-run artifacts satisfy the ordinary COMMOT output checks above.

## Boundaries

- COMMOT reports modeled communication potential, not direct proof of signaling.
- Compare results against expression maps for ligand and receptor genes.
- Use stLearn CCI or Squidpy `ligrec` as complementary methods, not interchangeable drop-in replacements.
- Do not mutate an active run's script when its manifest signatures include a
  script hash; test in isolated outputs and change the live launcher only at a
  safe restart boundary.
