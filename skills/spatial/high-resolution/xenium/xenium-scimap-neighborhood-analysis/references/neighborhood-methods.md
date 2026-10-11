# Neighborhood methods and portable baseline

`spatial_count` implementations differ in self inclusion, row normalization and boundary handling. Record the installed scimap version/signatures and validate with two ROIs having overlapping coordinates and disjoint labels before trusting a wrapper. Composition counting can be equivalent across implementations, but KMeans and Leiden assignments are not interchangeable. Do not promise model parity from count parity.

Use radius with physical units for a fixed reach, KNN for a fixed neighbor count with density-dependent reach. The portable KNN baseline breaks equal-distance ties by stable cell key, excludes self, and uses all available contributors when fewer than k exist. Radius includes distances <= radius. Tissue edges, sparse cells and missing phenotypes affect denominators.

RCN centers and contributors are independently declared masks. A center may lie outside the contributor pool. NMF input columns are counts or proportions of neighbor phenotypes, never genes. Retain raw counts and proportions regardless of representation used for fitting. A zero-labeled-neighbor center has no fit/assignment; preserve its row and status. RCN/NMF remain exploratory spatial composition states.

Portable config:

```json
{
  "cells_csv":"/path/cells.csv", "cell_key":"obs_name", "roi_key":"roi",
  "biological_unit_key":"donor", "x":"x", "y":"y", "coordinate_unit":"um",
  "phenotype_key":"phenotype", "label_evidence_grade":"reviewed_operational",
  "center_mask_key":"center_eligible", "neighbor_mask_key":"neighbor_eligible",
  "mode":"nmf", "cluster_method":"kmeans", "representation":"proportions",
  "radii":[50,80,100], "ranks":[2,3], "seeds":[0,1], "n_states":3,
  "min_centers":20, "max_iter":1000,
  "assignment_id":"my_target_nmf_v1", "namespace":"my_target_nmf_v1",
  "upstream":[]
}
```

Use `mode=rcn` without ranks for KMeans RCN. Use `knn:[10,20]` instead of radii. Omit a mask key to include every row; supplied masks must be nonmissing bool/0/1. `center_phenotypes` optionally restricts centers after the independent mask. Coordinate units are declared, not inferred. ROI keys must be unique across the declared cohort and nested in exactly one biological unit. Do not downsample the contributor pool for an ordinary run.

Outputs include composition per scale, a long keyed assignment table (unique cell+namespace), fit summary, state composition, NMF loadings/weights, comparison table (ARI and matched component cosine) and biological-unit support. Namespace records scale/rank/seed; assignments are qualified, e.g. `my_run/radius80/rank2/seed0::NMFstate1`. Component indices also need their namespace. Compare profiles with matching/correspondence, never equate same numbers.

Choose ranks by reconstruction, repeated-seed correspondence, component/state support and interpretability. More components trivially improve reconstruction. Nonconverged/under-supported fits are explicitly marked and are not selected automatically. Lock rank/clustering settings before radius sensitivity. Cross-scale matching is descriptive; it is not a new formal annotation.

The tool uses dense center-by-phenotype counts and sklearn NMF/KMeans; memory scales with centers times phenotype count. Bound the pilot before large runs. Distribution includes `<skills-root>/_shared/spatial/*.py`. No project package, registry or count matrix is required. Geometry/composition validation does not prove biological cluster stability or group inference. Real-data rank/scale tuning remains task-specific.
