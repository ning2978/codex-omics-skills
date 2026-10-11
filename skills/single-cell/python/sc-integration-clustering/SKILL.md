---
name: sc-integration-clustering
description: Normalize scRNA-seq counts, select matrix-appropriate HVGs, assess integration against a baseline and build explicit neighbors, Leiden candidates and UMAP. Use for initial clustering, retained-count reruns or unchanged-graph layout requests.
---

# Representation, integration and clustering

Read [../../workflow/workflow-sc-analysis/references/data-contract.md](../../workflow/workflow-sc-analysis/references/data-contract.md). Keep all-gene counts/log expression separate from scaled HVG work. Numerical helpers live in `workflow-sc-analysis`'s [sc_contracts.py](../../workflow/workflow-sc-analysis/scripts/sc_contracts.py).

## Normalization and HVG/PCA

- Normalize declared counts (default target sum 1e4), natural log1p into `layers["log1p"]`; preserve original counts. Corrected matrices have separate provenance; `.raw` has no implicit role.
- Default `seurat_v3` HVG uses **raw counts** with a justified batch key. On failure record exception/parameters and use `cell_ranger` on **normalized log1p**, with the same batch key. Degenerate bins/tiny batches may fail. Removing batch stratification requires explicit configured downgrade and rationale; never silently drop it. Report failure if no defensible method works.
- Record filter, HVG list/count, batch and matrix for every attempt. HVGs constrain representations, not markers. Scale HVG work only, avoid all-gene densification, and choose usable PCs/neighbors with variance/robustness and seeds recorded.
- Do not regress cell cycle, stress, interferon, donor or mitochondrial signals by default. Justify and compare requested regression.

## Integration decision

Preserve unintegrated PCA/graph/diagnostics. Distinguish technical batch from donor, sample and condition; examine contingency/confounding and shared biological types. `sample_id` alone does not justify correction. Full batch-condition confounding prevents separating technical/biological effects; correction cannot recover that inference.

Choose no integration or a justified method. Harmony changes a representation, BBKNN the graph, scVI requires its count/model runtime contract. Record fitted population, inputs/parameters and graph/`use_rep` scope. Compare baseline and candidates with within-type batch mixing (iLISI/kBET/neighbor composition where meaningful), type conservation (cLISI/silhouette), canonical markers/states, rare populations, donor support and held-out/reference evidence where available. UMAP mixing alone cannot select a method; distinct condition-specific populations need not mix. Record selected/rejected candidates and retain unintegrated expression for markers/DEG.

## Graph, clusters and embedding

- Neighbors explicitly selects representation/key/seed; graph integrators name their own graph. Leiden uses chosen adjacency, UMAP its explicit neighbor key. Save graph shape/hash and parameters.
- Output several coarse-to-fine resolution candidates (configurable `0.2/0.5/1.0` is an example), adapted to population size/hierarchy/continuous states. Compare marker coherence, competing lineages, sizes, sample dominance and seed stability (ARI/VI/coassignment). Select/justify final resolution; there is no universal 1.2 default. Stability alone is not identity evidence.
- UMAP distances/islands are visualization, not independent biological evidence. Record graph, seed/min_dist/spread and labels.
- Retained-cell reruns rebuild normalization through Leiden/UMAP from counts. `mode=layout_only` verifies unchanged keys/graph, adjusts only layout, preserves clusters/markers and records hashes. Representation/graph/cell changes invalidate this exception. Resolution-only sweeps may reuse unchanged graphs with evidence.

## Completion

Full counts/log H5AD with selected representation, graph, labels and finite two-dimensional UMAP; optional HVG work H5AD; HVG attempts, batch-condition table, integration assessment, resolution/seed table/rationale and run/validation manifest. Reopen to verify count/key integrity, graph alignment and sources. Matrix fingerprints/upstream evidence accompany recomputation flags. Hand off the **all-gene** object to `sc-marker-annotation`.

## Implementation handoff

Use an available validated project entrypoint through its repository adapter, or the shared numerical helpers. Require explicit counts origin, retained full-gene counts/log expression, matrix-appropriate HVG attempts, centered variance-scaled PCA, a baseline for correction, stage input/count/key/graph hashes, runtime and actual parameters. Resume may reuse only validated unchanged stages; changed cells invalidate the complete expression run. Historical presets are provenance, not defaults. The engine runs the declared method; scientific integration/resolution selection still needs the evidence above.

## Internal completion and user delivery

Read [the shared delivery contract](../../workflow/workflow-sc-analysis/references/delivery-contract.md)
when preparing a complete expression-analysis handoff. Completion artifacts and
validation above stay in the working directory outside `delivery/`; select final
UMAP PNGs, marker DotPlot PNGs, per-cluster annotation review CSV and README only.
Standalone plotting, marker export and review deliver only requested content.
Delivery selection does not authorize new analysis or reduce evidence requirements.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
