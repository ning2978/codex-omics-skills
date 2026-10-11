---
name: xenium-scimap-neighborhood-analysis
description: Build or review Xenium phenotype-composition RCN and target-cell NMF neighborhoods with separate center/contributor pools, ROI isolation, rank/seed/scale sensitivity and qualified result namespaces.
---

# Phenotype neighborhoods: RCN and target NMF

Use one local phenotype-composition representation for RCN and NMF. These are spatial clusters, not cell-expression clusters or new cell identities. Read [methods and portable tools](references/neighborhood-methods.md) before execution and [stage contracts](../../../shared/workflow-st-analysis/references/stage-contract.md) for handoff.

Require only unique cell key, explicit ROI boundary, coordinate unit, biological-unit nesting, trusted phenotype/evidence grade and selection masks. No expression matrix is needed for composition analysis. Resolve actual project label/RCN roles through the adapter; do not join by row order.

Declare the scientific question, centers and neighbor-contributor pool separately. Default centers are all eligible cells. Excluding epithelial centers keeps every eligible non-epithelial center and the full eligible contributor pool unless the user also excludes epithelial contributors. Epithelial membership comes from declared trusted labels, not a substring. Target NMF selects center rows after counting against the declared contributor pool.

Count within each ROI only, explicitly excluding self. Preserve raw counts, labeled-neighbor denominator and proportions. Missing neighbor phenotypes contribute no identity; report them. Zero labeled neighbors yield missing proportions and assignments, with explicit status. Edge/low-support centers stay flagged.

- RCN: fit the chosen clustering method to declared counts/proportions.
- NMF: fit nonnegative neighborhood phenotype features, compare rank and seed using reconstruction, component correspondence and support; lock rank before spatial-scale sensitivity. Optional clustering of NMF coordinates uses an explicitly selected method.
- Sensitivity: compare radius/KNN definitions, keeping eligible sets, labels and model settings fixed. Equal cluster numbers are never identities across runs.

The bundled tool provides an explicitly selected KMeans baseline and NMF grid:

```bash
python <skill-dir>/scripts/neighborhood.py --config <json> --out-dir <new-output>
```

It does not silently replace a requested scimap/Leiden method. For that method, inspect the installed API, use the project's validated runner, and apply the same contracts. New results use explicit assignment ID and namespace, remain exploratory, and do not become current automatically.

Verify overlapping-coordinate ROIs never share neighbors, centers/contributors can differ, zero-neighbor status persists and NMF columns are phenotypes. Export composition, assignments, loadings/weights, rank/seed/scale comparisons, support by biological unit and manifest. Cells and nested ROIs are not independent biological replicates; group inference routes to a separate sample-level design.
