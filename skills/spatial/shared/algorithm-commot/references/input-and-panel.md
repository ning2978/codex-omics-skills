# Explicit expression and measured-panel audit

Preparation config:

```json
{
  "input_h5ad":"/path/cells.h5ad", "matrix_source":"layers/counts", "matrix_kind":"counts",
  "output_transform":"normalize_log1p", "target_sum":10000, "counts_source":"layers/counts",
  "feature_key":"gene_symbol", "species":"human", "database_species":"human",
  "lr_database_csv":"/path/lr.csv", "database_name":"CellChat", "database_version":"declared_version",
  "complex_delimiter":"_", "min_cell_fraction":0.05,
  "roi_key":"roi", "biological_unit_key":"donor", "spatial_key":"spatial",
  "coordinate_unit":"um", "distance_threshold":100, "upstream":[]
}
```

Matrices are explicitly `X` or `layers/<name>`, with kind counts/normalized/log1p. Count input requires positive target_sum and normalization+log1p; normalized/log1p input uses `output_transform=preserve`. Negative/scaled/integrated matrices and fractional counts fail. The helper cannot infer matrix semantics from values alone; upstream provenance must establish them. It never changes the source H5AD.

`feature_key` selects `var_names` or a gene-symbol var column. Stable source var IDs are retained as `source_gene_id` while a new working object's var_names use audited canonical symbols for COMMOT lookup. Optional `aliases_csv` uses species-specific alias/canonical/species rows; ambiguous or many-to-one maps fail. Duplicate gene/cell keys fail rather than becoming suffixed names. No uppercase conversion or gene imputation is implicit.

LR CSV has ligand,receptor,pathway, one row per LR/pathway. Choose the complex delimiter from that database (COMMOT commonly uses underscore). Audit every subunit, not just one member; preserve unavailable and below-detection reasons. Global panel audit is descriptive, while eligible LR output is per ROI. Expression eligibility uses each subunit's measured positive fraction, not evidence that all subunits coexpress in the same cell. If the biology/method needs same-cell complex assembly, audit that separately. A filtered database must still preserve COMMOT's joint LR model; do not split LR pairs to accelerate it.

Coordinate unit and distance threshold are inseparable; convert pixels using documented image calibration before preparation. Independent ROI H5ADs prevent cross-boundary modeling. No model or formal identity is created. The bundled runtime probe in the execution reference remains the model-execution smoke check. A synthetic success does not validate a real tissue's LR biology or memory budget.

Method selection: COMMOT when collective spatial transport and sender/receiver geometry match the question; LIANA/CellChat/CellPhoneDB for expression/database evidence or biological-unit communication contrasts; Squidpy/stLearn for spatial association alternatives; NicheNet for a declared target-gene response hypothesis. Preserve species, measured panel/complex audit, matrices, ROI/unit design, DB/version, distance semantics and normalization in every handoff. Similar pathway names across methods do not make scores or null hypotheses equivalent. Use `sc-cell-communication` for nonspatial expression-based inference. Do not call association an inferred causal response without independent evidence.
