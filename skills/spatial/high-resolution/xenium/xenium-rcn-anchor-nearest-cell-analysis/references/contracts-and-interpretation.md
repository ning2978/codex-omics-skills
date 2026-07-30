# Contracts and interpretation

## Output tables

### `nearest_any_eligible_cell_per_anchor.csv.gz`

One row per anchor. Required fields:

- anchor observation key, ROI, coordinates;
- target observation key, phenotype, coordinates;
- Euclidean distance;
- query scope and tie-breaking rule.

The target may share the anchor phenotype, but must not be the focal cell.

### `nearest_heterotypic_cell_per_anchor.csv.gz`

Use the same schema, with every anchor-phenotype target removed from the
candidate pool before querying.

### `nearest_category_summary.csv`

For each query scope and phenotype, retain:

- cell-pooled count and fraction;
- conditional distance median and IQR among anchors for which the phenotype
  was the winning nearest category;
- ROI-macro mean/median fraction;
- number of ROI and positive ROI.

Zero-fill missing categories within each ROI before computing macro summaries.

### `anchor_to_phenotype_nearest_distance.csv.gz`

One row per anchor × eligible phenotype. Preserve rows with no target in the
anchor ROI and set target identity/distance to missing. This table supports
phenotype-specific distance rankings.

### Distance summaries

Retain both:

- cell-pooled median and IQR across anchors with an available target;
- median and IQR of ROI-specific medians;
- anchor coverage and ROI coverage.

Do not rank a phenotype with zero coverage as if its missing distance were
large or infinite.

## Recommended figures

1. Two-panel nearest-category bar plot:
   nearest-any and nearest-heterotypic; bars for cell-pooled fractions and
   points for ROI-macro means.
2. Phenotype-specific distance ranking:
   pooled median/IQR and median of ROI medians/IQR.
3. Optional per-group or per-pair plot only when the user asks for biological
   group comparison.

## Interpretation examples

Safe:

- “Another Endothelial cell was the nearest eligible cell for 64% of anchors.”
- “After excluding Endothelial targets, Pericyte/mural was the most frequent
  nearest phenotype.”
- “This supports a descriptive vascular/perivascular organization.”

Unsafe without a spatial null:

- “Pericytes are significantly attracted to endothelial cells.”
- “Macrophages are depleted around the anchor.”
- “The nearest-category fraction is a colocalization p value.”

## Formal enrichment extension

Nearest-cell identity is abundance-sensitive. If the user asks whether a
phenotype is preferentially colocalized:

1. Freeze anchors and the eligible target universe.
2. Preserve ROI boundaries.
3. Prefer permutations within exact ROI × spatial blocks, or use a
   prespecified point-process null that preserves density.
4. Keep phenotype counts fixed within the permutation unit.
5. Compare observed nearest-category frequency and/or radius-limited pair
   counts to the null.
6. Adjust over all tested phenotypes and metrics.
7. Summarize direction and recurrence at sample/pair level.
8. Label small or non-recurrent results exploratory.

Do not silently turn a descriptive nearest-distance request into this
inferential analysis.

## Common failures

- Using row order rather than observation keys.
- Allowing the focal cell to be its own nearest Endothelial target.
- Querying across ROI that reuse coordinate systems.
- Restricting the target pool to RCN-assigned center cells.
- Filling absent phenotype targets with a finite distance.
- Treating cells as biological replicates.
- Calling cell-center distances membrane-to-membrane contact.
- Interpreting a pooled fraction without showing ROI heterogeneity.
