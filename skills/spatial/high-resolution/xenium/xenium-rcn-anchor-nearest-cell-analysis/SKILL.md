---
name: xenium-rcn-anchor-nearest-cell-analysis
description: Analyze which Xenium cell phenotypes are spatially closest to a fixed anchor phenotype inside an existing RCN/CN, using same-ROI nearest-any, nearest-heterotypic, and phenotype-specific distance queries. Use for anchor-cell colocalization summaries, ROI or biological-unit recurrence checks, and auditing nearest-cell distance outputs without recomputing RCNs.
---

# Xenium RCN Anchor Nearest-Cell Analysis

## Scope

Use frozen RCN assignments and a trusted RCN phenotype sidecar to describe the
cells closest to a selected anchor population. Treat this as a downstream
spatial-distance analysis, not a new RCN run or a new cell annotation.

Use the bundled `scripts/rcn_anchor_nearest_cells.py` for deterministic
same-ROI queries. Read `references/contracts-and-interpretation.md` when
designing outputs, reviewing claims, or extending the workflow to formal
enrichment tests. Read `references/generic-anchor-nearest-example.md` for a
portable example with generic input names and expected-count gates.

## Resolve the registered track first

1. Read the repository annotation policy and current annotation registry.
2. Resolve exactly one RCN/CN use case and its single current track.
3. Join Xenium tables with the registered globally unique key, normally
   `obs_name`; never use row order or a non-unique `cell_id`.
4. Preserve candidate, exploratory, unresolved, and RCN-only label semantics.
5. Do not update the registry merely because this downstream analysis passes.

Stop if the human-readable guide and machine registry disagree.

## Define the question

Always produce these complementary views unless the user narrows the request:

- **Nearest any eligible cell:** exclude the focal cell but allow another cell
  of the anchor phenotype. This detects homotypic clustering.
- **Nearest heterotypic cell:** exclude all cells of the anchor phenotype.
  This answers which other phenotype most directly contacts the anchor.
- **Phenotype-specific nearest distance:** for every anchor, find the nearest
  target of each eligible phenotype. Rank pooled medians and ROI-macro medians,
  and report target coverage.

Do not confuse:

- the distance conditional on a phenotype winning the nearest-cell query; and
- the distance from every anchor to the nearest target of that phenotype.

## Input contract

Require:

- a phenotype sidecar with a unique observation key, ROI/FOV key, trusted
  phenotype, and eligible-neighbor flag;
- frozen RCN assignments keyed by the same observation key;
- an H5AD with unique `obs_names` and two-dimensional `obsm["spatial"]`;
- consistent coordinate units within every ROI;
- a selected RCN label and anchor phenotype.

Select anchors from the frozen center assignments first, then match the anchor
phenotype from the registered sidecar. Build the candidate target pool from
the complete eligible-neighbor universe. Never restrict targets to cells that
are themselves assigned to the selected RCN unless the user explicitly asks
that different question.

## Run the bundled workflow

Use a validated Python environment containing `anndata`, `pandas`, `numpy`,
`scipy`, and `matplotlib`.

```bash
python /path/to/skill/scripts/rcn_anchor_nearest_cells.py \
  --sidecar /path/to/registered_sidecar.csv.gz \
  --assignments /path/to/frozen_assignments.csv.gz \
  --spatial-h5ad /path/to/thin_spatial.h5ad \
  --anchor-rcn RCN_A \
  --anchor-phenotype Endothelial \
  --out-dir /path/to/new_result
```

Override column names only when the registered contract uses different names:

```bash
python /path/to/skill/scripts/rcn_anchor_nearest_cells.py \
  ... \
  --obs-key obs_name \
  --roi-column sample_roi \
  --phenotype-column rcn_phenotype \
  --eligible-column include_neighbor_pool \
  --rcn-column rcn_aligned_80
```

Pass optional expected counts as hard validation gates when they have been
established independently:

```bash
  --expected-neighbor-count "${EXPECTED_NEIGHBORS}" \
  --expected-anchor-count "${EXPECTED_ANCHORS}" \
  --expected-phenotype-count "${EXPECTED_PHENOTYPES}"
```

## Distance rules

- Query only within the same ROI/FOV/library.
- Exclude the focal `obs_name` from every nearest-cell query.
- Break exact distance ties by the lexicographically smallest target
  observation key.
- Use Euclidean cell-center distance and state the coordinate unit.
- Do not impose the RCN construction radius on exact nearest distances unless
  the user explicitly requests a radius-limited analysis.
- Retain missing phenotype-specific distances when a target phenotype is
  absent from an ROI; report both anchor coverage and ROI coverage.

## Summarize at two levels

Report cell-pooled results for descriptive composition and ROI-macro results
for recurrence:

- nearest-category count and fraction;
- mean or median fraction across ROI, with zero-filled absent categories;
- number of positive ROI;
- per-pair or per-sample fractions when biological-unit comparisons matter;
- pooled median nearest distance and IQR;
- median of ROI medians and IQR;
- target presence/coverage.

Do not use individual cells as independent replicates for group inference.
Aggregate nested ROI at sample/pair level or use a hierarchical model.

## Interpret conservatively

State that nearest-category frequency depends on cell abundance, tissue
density, segmentation geometry, and spatial clustering. A common nearest
phenotype is not automatically enriched or attracted.

Use terms such as “nearest”, “adjacent by cell-center distance”, or
“descriptive perivascular pattern”. Reserve “enriched”, “attracted”, or
“preferentially colocalized” for an explicit spatial null that preserves ROI
boundaries and local structure.

For formal enrichment, add a prespecified block-stratified permutation or
appropriate point-process null, adjust across tested phenotypes, and require
recurrence across biological units. Keep that extension separate from the
descriptive nearest-cell output.

## Validate before reporting

Require all of the following:

- sidecar, assignments, and H5AD observation keys are unique as declared;
- every join preserves the left-table row count;
- all anchors belong to the selected frozen RCN and anchor phenotype;
- every target is in the eligible pool;
- nearest output has exactly one row per anchor;
- focal cell is never its own target;
- heterotypic results contain no anchor phenotype;
- anchor and target ROI keys match;
- stored distances equal distances recomputed from saved coordinates;
- category fractions sum to one within each ROI and query scope;
- phenotype-specific output covers every anchor × phenotype combination;
- missing targets remain missing rather than being imputed;
- synthetic overlapping-coordinate ROI and tie-breaking tests pass;
- PNGs open and output hashes match the manifest;
- no annotation, H5AD, RCN assignment, or registry input changed.

Run the repository reclustering guard when present, even though this workflow
must not invoke expression or spatial reclustering.
