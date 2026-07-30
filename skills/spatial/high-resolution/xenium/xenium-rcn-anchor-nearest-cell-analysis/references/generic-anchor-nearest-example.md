# Generic Xenium RCN anchor example

Use this example only after resolving the current annotation and RCN tracks in
the target repository. Replace every placeholder with a registered input or a
new output path; do not infer current files from version suffixes.

## Input roles

- `registered_rcn_sidecar.csv.gz`: one row per globally unique observation,
  including the trusted RCN phenotype, ROI key, and neighbor-eligibility flag.
- `frozen_rcn_assignments.csv.gz`: frozen RCN/CN assignments keyed by the same
  observation identifier.
- `thin_spatial.h5ad`: unique `obs_names` and two-dimensional
  `obsm["spatial"]`.
- `RCN_A`: the selected frozen RCN label.
- `Endothelial`: an example anchor phenotype; choose the phenotype required by
  the analysis question.

## Reproducible command

Run from the repository root:

```bash
python /path/to/skill/scripts/rcn_anchor_nearest_cells.py \
  --sidecar inputs/registered_rcn_sidecar.csv.gz \
  --assignments inputs/frozen_rcn_assignments.csv.gz \
  --spatial-h5ad inputs/thin_spatial.h5ad \
  --anchor-rcn RCN_A \
  --anchor-phenotype Endothelial \
  --obs-key obs_name \
  --roi-column sample_roi \
  --phenotype-column rcn_phenotype \
  --eligible-column include_neighbor_pool \
  --rcn-column rcn_label \
  --track-id registered_track_identifier \
  --out-dir results/anchor_nearest_cell_analysis
```

Add the expected-count options only after computing the expected values from
an independent input audit:

```bash
  --expected-neighbor-count "${EXPECTED_NEIGHBORS}" \
  --expected-anchor-count "${EXPECTED_ANCHORS}" \
  --expected-phenotype-count "${EXPECTED_PHENOTYPES}"
```

## Acceptance checks

- The sidecar, assignments, and H5AD keys are unique as declared.
- Every join preserves the left-table row count.
- Anchors match both the selected RCN and anchor phenotype.
- Targets remain within the same ROI and the eligible-neighbor universe.
- The focal cell is never its own target.
- The output directory is new and no registry, annotation, H5AD, or frozen RCN
  assignment is modified.
