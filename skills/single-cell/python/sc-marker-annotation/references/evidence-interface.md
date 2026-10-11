# Portable evidence collector

Use `scripts/collect_evidence.py --input full.h5ad --config config.json --out new-directory`. It reuses the shared all-gene t-test/log-expression routines. It does not choose identities, thresholds for exclusion or biological stopping. It does not alter input H5ADs or write registries.

Configuration example (replace example feature IDs with sourced, assay-appropriate programs):

```json
{
  "cluster_key": "leiden", "sample_key": "sample_id", "unit_key": "donor_id",
  "species": "declared species", "tissue": "declared tissue", "platform": "xenium_segmented",
  "gene_scope": "all available panel genes; no HVG restriction",
  "counts_origin": {"kind":"raw_counts","source":"identified original artifact"},
  "min_detected_cells": 3, "target_sum": 10000, "method": "t-test",
  "qc_fields": ["cell_area", "nucleus_area", "nucleus_count", "segmentation_method"],
  "programs": [
    {"id":"target","role":"target","lineage":"target lineage","source":"reference/version","genes":["target_feature_1","target_feature_2"]},
    {"id":"alternative","role":"competing","lineage":"one foreign lineage","source":"reference/version","genes":["foreign_feature_1","foreign_feature_2"]}
  ]
}
```

Programs use explicit unique `var_names` feature IDs; ambiguous symbols must be mapped/audited first. Roles are `target`, `competing`, `state`. No tissue/species default panel. `unit_key: null` explicitly marks unavailable unit support. Samples must map to one declared biological unit; nested captures are not independent replicates. Choose sc/snRNA QC fields such as mitochondrial fraction and doublet scores when available. Supplied log expression must declare unintegrated natural-log1p provenance; if absent, a disposable in-memory log layer is prepared from the declared counts.

Outputs:

- `all_markers.csv.gz`, `top20.csv`: full available expression-filtered statistics and separate display; skips are explicit in manifest, never absence evidence. The default top20 display uses adjusted p < .05, pct_in > .3 and logFC order. Default t-test; alternate method requires an explicit justified config.
- `gene_availability.csv`: measured availability, detected cells and test eligibility are separate. Missing feature values are NA, never zero.
- `gene_evidence.csv`: cluster, sample and biological-unit detection and mean log expression for all panel genes, including low-detection genes outside marker filtering. `pct_out` always refers to the full retained cluster complement, including in unit-scoped rows.
- `program_support.csv`: any-gene and two-or-more distinct measured-gene detection by cluster/sample/unit. Two-marker support is NA if fewer than two genes are measured.
- `codetection.csv`: target >=1 and competing >=2 measured genes **in the same cell**, with evaluability; not a doublet classifier. Compare this with separate program support to distinguish potential mixtures; shared marker programs and sparse detection need expert interpretation.
- `unit_composition.csv`: counts and fractions within cluster and sample/unit; descriptive support only, no pseudoreplicated condition testing.
- `qc_evidence.csv`: median/quantiles/missing values or categorical counts, with unavailable fields explicit.
- `manifest.json`: config, hashes, matrix scope, detection source, skipped tests and technical/evidence status.

Review every requested group against complete markers, sourced program availability, specificity, target/competing co-detection, sample/unit recurrence and QC. An all-unavailable program is not negative evidence. The helper validates structural inputs; it cannot authenticate program sources or determine whether proposed coordinated markers are biologically independent. Record limitations and alternatives in the annotation/review table.
