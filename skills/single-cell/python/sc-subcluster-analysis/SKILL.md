---
name: sc-subcluster-analysis
description: Use this skill when extracting a single-cell subset for reclustering and marker-driven subcluster annotation, especially for immune/myeloid/macrophage/monocyte subsets where raw counts must be preserved, HVG selection needs fallback logic, contamination must be judged from top markers, and labels should be proposed from current marker evidence rather than reused blindly. Also use it for requests involving new UMAP, re-embedding, rerun dimensionality reduction, recluster after filtering/removal/cleaning, filter then rerun, 重新降维, 重新聚类, or 新的 UMAP.
---

# Skill: sc-subcluster-analysis

## Purpose

Guide focused scRNA-seq subclustering from an existing AnnData object through:

- subset extraction
- counts-preserving preprocessing
- HVG fallback
- PCA/integration/neighbors/UMAP/Leiden
- top marker evidence
- contamination review
- candidate biological naming

## Core Workflow

**Hard contract for new embeddings:** When the user asks for a new UMAP, re-embedding, rerun dimensionality reduction, reclustering after filtering/removal/cleaning, or equivalent language, recompute the analysis from the retained cells' `layers["counts"]`: HVG, PCA, Harmony or the selected batch correction, neighbors, UMAP, and Leiden. Do not satisfy these requests by plotting, copying, or relabeling global UMAP coordinates or an existing `X_umap`.

Apply the same contract to gene-expression profiles from individually
segmented Xenium cells. They are cell-level expression observations, not
spot-level mixtures. Clear stale `X_umap`, neighbor graphs, and the target
cluster key before recomputing them.

**Old-coordinate handling:** Keep prior/global `X_umap` only as a historical field or as a separately named projection/overlay figure when the user explicitly requests projection or overlay. It must not become the primary `X_umap`, the main UMAP output, or evidence that the cleaned subset was re-embedded.

**Resolution Strategy:** For post-filtering subset reruns or high-resolution spatial/scRNA subsets, sweep at least `1.0/1.2/1.5` and use `1.2` as the default final Leiden resolution unless the user specifies otherwise. After contamination removal, always recompute the clean UMAP/Leiden from the retained cells' counts rather than reusing the pre-removal embedding or clusters.

**Exemption boundary:** A resolution-only sweep may reuse an existing graph
only when both the cells and graph are unchanged and the exception is recorded.
scimap RCN/CN, target-cell NMF neighborhood states, and CellCharter niche
clusters are not expression reclustering and do not require a new expression
UMAP.

1. **Confirm input contract**
   - Identify the source h5ad, subset key, subset labels, condition key, sample/batch key, and existing counts source.
   - Preserve counts in `layers["counts"]`. If no counts layer exists, document that `X` is being used and whether it appears raw/count-like.
   - Record if the input object has already been gene-subsetted; marker interpretation must state that marker ranking is limited to available genes.

2. **Build the subset**
   - Extract only the requested biological class unless the user explicitly wants boundary/contamination classes included.
   - Add basic QC from counts: `total_counts`, `n_genes_by_counts`.
   - Keep a full subset object separate from the HVG working copy.

3. **Recluster with fallback**
   - Start from counts.
   - Run `highly_variable_genes(..., flavor="seurat_v3", batch_key=batch_key)`.
   - If it fails, try `flavor="cell_ranger"` with the batch key.
   - If that fails, use `flavor="cell_ranger"` without batch key.
   - Normalize, log1p, subset to HVGs, scale, PCA.
   - Use Harmony when a suitable batch key exists; otherwise use PCA directly.
   - Run neighbors, UMAP, and Leiden with explicit resolution and seed.
   - For any post-filtering rerun, perform all steps in this item again from the retained cells' counts and write the recomputed embedding as the main `X_umap`.

4. **Compute marker evidence**
   - By default, run marker ranking on the subset-specific HVG working copy, where HVGs were recalculated from the extracted subset counts.
   - Do not reuse global/integration HVGs for marker ranking unless the user explicitly asks for that narrower scope.
   - Only run marker ranking on the full-gene subset object when the user explicitly requests full-gene markers; document the larger scope and runtime/memory tradeoff.
   - Use Wilcoxon with `pts=True`.
   - Export all markers with `gene`, `score`, `logfoldchanges`, `pvals`, `pvals_adj`, `pct_in`, and `pct_out`.
   - Export filtered top20 per cluster using:
     `pvals_adj < 0.05`, `pct_in > 0.3`, sorted by `logfoldchanges` descending.
   - Add a marker-scope note to marker tables or the run summary, such as `subset_recomputed_HVG`, `full_gene_subset`, or `global_integration_HVG`.

5. **Name cautiously**
   - Do not reuse old subtype names by default.
   - Generate candidate labels from the current cluster's strongest top markers plus a biological suffix supported by marker evidence.
   - Treat labels as candidates until reviewed. Prefer names like `GENE1_GENE2_function_like`.
   - If markers indicate contamination, label the cluster as requiring review rather than forcing a subtype name.

6. **Contamination review**
   - Check filtered top markers for non-target signals such as epithelial, T/NK, B/plasma, endothelial, stromal, RBC, neutrophil, and mast/basophil markers.
   - Export `potential_contamination_clusters.csv` with cluster, candidate label, evidence genes, and recommended action.
   - Do not delete clusters automatically unless the user explicitly asks.
   - When the user approves removal of contamination clusters, remove cells by the previous cluster key/value, then rerun the subclustering workflow from the remaining cells' counts.
   - Record removed cluster values, cells before removal, excluded cell count, and cells after removal in the run summary.
   - Record post-removal recomputation flags in the run summary: `umap_source=recomputed_clean_subset`, `recomputed_hvg=true`, `recomputed_neighbors=true`, and `recomputed_leiden=true`.
   - If reclustering reveals a new non-target cluster, repeat removal only when it is clearly marker-supported or the user confirms it.

## Required Outputs

- Annotated full subset h5ad with counts, UMAP, cluster key, and candidate label key.
- HVG working-copy h5ad.
- `all_markers` and filtered `top20` marker tables.
- Candidate annotation/evidence table.
- Potential contamination table.
- Counts by candidate label, dataset/condition/sample.
- UMAPs for cluster, candidate label, dataset, condition, sample, and marker module scores.
- Dotplot or violin marker panels.
- Run summary CSV documenting input, subset, genes available for markers, marker scope, HVG method, batch key, and marker filter.
- For any new UMAP or post-filtering rerun, the run summary must include `umap_source=recomputed_clean_subset`, `recomputed_hvg=true`, `recomputed_neighbors=true`, and `recomputed_leiden=true`.
- Also record `recomputed_pca=true` and `recomputed_umap=true`; retain all five
  flags in machine-readable H5AD/run provenance.
- If projection/overlay plots are requested, export them with explicit projection naming and document that they are not the main recomputed UMAP.

## Validation Checklist

- Re-open the output h5ad and confirm `X_umap`, cluster key, candidate label key, and `layers["counts"]`.
- Confirm the HVG working copy was recalculated from the extracted subset counts and that marker ranking used the intended scope.
- When the user requested a new UMAP or post-filtering rerun, confirm `X_umap` was recomputed from the retained subset and that the summary has `umap_source=recomputed_clean_subset`, `recomputed_hvg=true`, `recomputed_neighbors=true`, and `recomputed_leiden=true`.
- Require a two-dimensional current `X_umap`, neighbor graph, cluster key,
  `layers["counts"]`, and true
  `recomputed_hvg/pca/neighbors/umap/leiden` provenance before declaring the
  expression recluster complete.
- Confirm every cluster has filtered top20 marker evidence or explicitly report clusters without sufficient markers.
- Confirm contamination table exists even if empty; an empty table should still retain column headers.
- Confirm summary documents marker gene scope and HVG fallback result.
- If contamination was removed, confirm the summary records the excluded cluster values and cell counts.
- Treat validation as failed if the user requested a new UMAP and the summary contains `no_umap_recomputed=true`, or if `umap_source` points to global coordinates, an existing h5ad embedding, or any reused `X_umap`.

## Example Pattern

For a myeloid subset, extract monocyte/macrophage-like labels from a full-count object, recalculate HVGs from the subset counts, recluster with the project batch key, and generate candidate labels from subset-HVG top markers rather than reusing broad parent labels blindly.
