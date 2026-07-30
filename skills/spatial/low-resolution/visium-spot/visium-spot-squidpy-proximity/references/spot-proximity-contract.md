# Spot-Level Squidpy Proximity Contract

## Required AnnData Inputs

Minimum inputs for a reusable Visium spot-level proximity script:

- `adata.obs`: sample metadata and optional region/niche labels.
- `adata.obsm[abundance_key]`: deconvolution abundance matrix.
- `adata.uns[column_key]`: abundance column names, unless the abundance matrix is a DataFrame.
- `adata.obsp[connectivities_key]`: existing spatial graph.
- `adata.obsp[distances_key]`: existing spatial distances, if topology diagnostics use it.
- `adata.obsm[spatial_key]`: spatial coordinates for physical rings.

Prefer configurable CLI keys:

```text
--abundance-obsm-key
--celltype-columns-uns-key
--sample-key
--donor-key
--condition-key
--timepoint-key
--niche-key / --niche-value
--connectivities-key
--distances-key
--spatial-key
```

## Full-Mode Output Tables

Use explicit filenames that encode the distance definition:

```text
topology_<analysis>_1hop_celltype_proximity.csv
topology_<analysis>_2hop_celltype_proximity.csv
physical_<analysis>_ring1_<low>_<high>_celltype_proximity.csv
physical_<analysis>_ring2_<low>_<high>_celltype_proximity.csv
anchor_spot_summary_by_sample.csv
neighbor_spot_summary_by_sample.csv
README.md
```

Recommended row columns:

```text
analysis
distance_method
hop or ring_low/ring_high
neighbor_region
sample
donor
condition
timepoint
celltype
neighbor_mean_abundance
background_mean_abundance
fold_change
log2_fold_change
p_value
p_value_adj_bh
anchor_spot_count
neighbor_spot_count
background_spot_count
neighbor_niche_or_region_spot_count
neighbor_non_niche_or_region_spot_count
```

## Support-Aware Plotting

When summarizing by timepoint or condition:

1. Aggregate by `timepoint`, distance ring/hop, and `celltype`.
2. Weight mean abundance by `neighbor_spot_count`.
3. Keep `n_samples_with_neighbors` and `total_neighbor_spots`.
4. Filter before ranking plotted cell types:
   - minimum samples with neighbors
   - minimum total neighbor spots
   - optional minimum background abundance for log2FC
5. Write both:
   - filtered rank table for interpretation
   - filter summary table with failed rows and reasons

Use this interpretation split:

- abundance rank: "what is most present near anchors"
- log2FC rank: "what is specifically enriched relative to background"

Disagreement is expected. A cell type can have high abundance near anchors but low or negative log2FC if the broader region is also rich in that cell type.

## Common Failure Modes

- **False strong log2FC**: background mean is nearly zero. Add minimum background abundance and support filters.
- **Overstated D0/D30 or rare group patterns**: only one sample or very few neighbor spots. Mask low-support groups in plots.
- **Cross-sample graph edges**: stop and inspect graph construction before continuing.
- **Anchor dominates its own neighborhood**: exclude anchor spots from neighbor sets; optionally exclude the anchor cell type from "nearest cell type" ranking.
- **Wrong label space**: verify that region/niche labels come from the same workflow as the AnnData object.
- **Overinterpretation**: keep the result framed as Visium spot-level proximity.
