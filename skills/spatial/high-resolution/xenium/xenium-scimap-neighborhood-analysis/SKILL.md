---
name: xenium-scimap-neighborhood-analysis
description: Use for Xenium scimap neighborhood analysis based on local cell-type composition, including whole-population recurrent cellular neighborhoods, target-cell NMF neighborhood states, radius or KNN sensitivity analysis, ROI-level comparisons, and validation of spatial_count outputs.
---

# Xenium scimap Neighborhood Analysis

## Choose a Mode

| Question | Mode |
|---|---|
| Which recurrent neighborhood types occur across all cells? | RCN: `spatial_count -> spatial_cluster` |
| Which neighborhood states surround one target phenotype? | Target NMF: `spatial_count -> target subset -> NMF -> Leiden` |
| Does a conclusion survive another spatial scale? | Radius/KNN sensitivity |

These modes share one input contract and one local phenotype-count matrix. Do not maintain separate skills for them.

## Shared Input Contract

Require:

- an h5ad with `obsm["spatial"]`;
- an ROI/FOV/library key that prevents neighborhoods spanning images;
- a trusted phenotype column or a cell/cluster sidecar joined by stable identifiers;
- consistent coordinate units within every ROI;
- recorded radius or KNN definition, seed, exclusions, and label provenance.

Attach and audit labels before neighborhood counting. Never join sidecars by row position unless row identity is explicitly guaranteed.

## Shared Preflight

1. Confirm the requested phenotypes exist after label harmonization.
2. Exclude removed or untrusted labels explicitly.
3. Inspect cell density and coordinate units before selecting a radius.
4. Run `spatial_count` within each ROI/FOV only.
5. Verify the resulting feature columns are phenotype names, not genes.
6. Record whether features are raw neighbor counts or transformed proportions and whether the installed scimap version includes the focal cell.

Avoid full-cell downsampling before spatial counting except for a clearly labeled smoke test; it changes every local neighborhood.

When a wrapper or scimap version changes, run a synthetic boundary test with two ROIs that share coordinates but have disjoint phenotypes. Counts for each ROI must contain no phenotype contributed only by the other ROI.

## RCN Mode

1. Build the local phenotype-count matrix with `sm.tl.spatial_count`.
2. Cluster `spatial_count` with Leiden or KMeans.
3. Map raw clusters to stable labels such as `RCN1..RCNn`.
4. Export per-cell assignments, phenotype composition, ROI fractions, and optional group comparisons.

Use ROI or sample fractions as replicates for group-level inference, not individual cells.

## Target-Cell NMF Mode

1. Build neighborhoods using the full eligible cell universe.
2. Subset target cells only after `spatial_count` is complete.
3. Set `X` to the target cells' neighborhood phenotype counts.
4. Fit NMF, then build neighbors and Leiden clusters from NMF coordinates.
5. Rank phenotype features by neighborhood cluster.

Interpret ranked features as enriched neighbor phenotypes. NMF-Leiden labels are target-cell neighborhood states, not new cell types.

Choose NMF rank from a small grid using reconstruction error, repeated-seed stability, component support, and interpretability. Lock the selected rank before radius sensitivity. Record minimum target cells and minimum cells per reported state.

## Sensitivity Mode

- Compare several biologically plausible radii in the coordinate unit, or compare radius with KNN.
- Keep labels, eligible cells, seed, NMF components, and clustering parameters fixed when testing spatial scale.
- Compare composition, assignment correspondence, repeated-seed stability, group-effect direction, support, and NMF reconstruction error.
- State explicitly when a conclusion is radius-dependent.
- Flag tissue-edge cells or unusually low neighbor support instead of interpreting them as ordinary neighborhoods.

## Environment

Use the project's validated scimap environment. When creating one, use conda plus pip only and prefer a single compatible Python/package stack for all three modes. Record the versions of Python, scimap, anndata, scanpy, leidenalg, igraph, and scikit-learn. Verify API names against the installed scimap version; the current repository adapter is tested against scimap 2.3.6.

## Output Contract

RCN outputs:

- per-cell assignment and raw-to-stable label map;
- phenotype-by-RCN composition;
- RCN counts/fractions by ROI and optional group comparison;
- parameter record and spatial/composition plots.

Target-NMF outputs:

- target-cell neighborhood assignments;
- NMF loadings and ranked neighbor phenotypes;
- cluster counts/fractions by ROI and group;
- reconstruction error, parameter record, and spatial plots.

Sensitivity outputs:

- per-scale run summary;
- cross-scale composition and support tables;
- stability or correspondence table for assignments/components;
- explicit list of failed or under-supported scales.

## Validation

- Assignment row count matches the recorded analyzed or target-cell count.
- ROI/FOV boundaries were respected and no cross-ROI neighborhood was created.
- Fraction tables sum to approximately one within each ROI.
- NMF columns are neighbor phenotypes, not genes.
- Group comparisons use independent biological samples. If multiple ROIs are nested within a sample, aggregate at sample level or use a hierarchical model rather than treating ROIs as independent.
- Small target sets or single-ROI results remain exploratory.
- A biological claim based on one radius includes a sensitivity result or a stated limitation.
- A synthetic overlapping-coordinate test proves that the configured image/ROI key blocks cross-ROI neighbors.
