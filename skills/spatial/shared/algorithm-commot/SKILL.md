---
name: algorithm-commot
description: Prepare, run or audit COMMOT spatial communication with explicit measured expression, species-matched LR/complex panel coverage, physical distance, ROI isolation and verified runtime/parallel execution.
---

# COMMOT spatial communication

COMMOT models spatially constrained communication potential; it does not prove signaling. Read [input/panel and method handoff](references/input-and-panel.md) for preparation. Before model/runtime work, read [execution and runtime diagnosis](references/runtime-and-execution.md), including the existing `scripts/runtime_probe.py`. Do not load model/runtime detail for a panel-only audit.

Declare an existing expression matrix and its counts/normalized/log1p meaning, species, database/version/source, unique gene/cell keys, ROI boundaries, biological units and distance units. Never use scaled/integrated/HVG-only expression as full measured LR evidence. Preserve counts separately when needed. Do not recluster solely because communication summaries have no trusted label; cell-level modeling can precede biological grouping.

Audit the complete measured Xenium panel for both LR members and every required complex subunit, then audit detection within each modeled ROI. Unmeasured genes are unavailable evidence; no transfer-imputed expression is used. Keep measured availability, expression eligibility and inferred communication distinct.

```bash
python <skill-dir>/scripts/prepare_inputs.py --config <json> --out-dir <new-output>
```

This portable tool emits mapped, explicit-expression ROI H5ADs, eligible LR tables and panel/subunit audits; its status is `PREPARED_MODEL_NOT_RUN`. Inspect them before calling `ct.tl.spatial_communication`. Fit independent ROIs without cross-ROI edges and preserve sender/receiver denominators and unit support. Other methods use the method-selection reference; results are not interchangeable.

Use an existing verified environment. COMMOT 0.0.3 has no device/n_jobs solver flag; inspect the exact version before promising GPU/multicore acceleration. The existing deterministic serial/parallel runtime probe verifies a synthetic process path, not real-data memory safety. Full runs require a representative real-data memory/output pilot. Do not split jointly coupled LR pairs as independent models.

Export full LR/pathway matrices and summaries, distances/units, filtering and complex availability, parameters/environment, manifests, validation state and PNGs. Resume requires matching complete artifacts. Replication/group inference needs genuine biological units; nested ROIs/cells are not independent replicates. Keep missing/unavailable LR and failed ROIs explicit, and keep modeled scores separate from annotation.
