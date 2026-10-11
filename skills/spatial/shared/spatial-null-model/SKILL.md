---
name: spatial-null-model
description: Test fixed-anchor spatial proximity, cooccurrence and local phenotype enrichment with ROI/block/stratum label permutations, explicit null hypotheses and biological-unit support. Group-level permutation is a separate handoff.
---

# Spatial conditional label nulls

Use this for a declared spatial relationship on a fixed segmented-cell/point set. Read [null design and tool schema](references/null-contract.md) before execution. Geometry needs no expression matrix. Resolve label and any RCN-selection provenance through a project adapter.

Declare unique cell key, ROI boundary, units, biological-unit nesting, fixed anchor mask, complete candidate-pool mask and target-label field. Lock spatial block size/origin, optional strata, radius, statistic, direction, seed, permutation count and test families. Coordinate proximity does not imply direct contact or signaling.

Within ROI x spatial block x optional strata, permute known target labels while preserving each unit's phenotype counts. Keep coordinates, anchors and candidate membership fixed. Missing labels stay missing; missing strata cannot supply conditional evidence. Constant/too-small strata and invariant hypotheses are explicit `NOT_TESTABLE`, not negative findings.

```bash
python <skill-dir>/scripts/spatial_null.py --config <json> --out-dir <new-output>
```

Export observed/null summaries and distributions, observed-minus-null effects, `(b+1)/(B+1)` empirical p values, BH for the predeclared full family, margins/stratum diagnostics and biological-unit support. Repeat plausible block/radius scales with the same design; do not choose the smallest p value.

The tool implements conditional within-tissue label randomization, not a biological-group permutation test. It summarizes ROI means within each biological unit and gives units equal weight. Cells or nested ROIs do not become independent biological replicates. Group inference requires an estimable sample-level design and a separate method.

Same-data anchor/hypothesis selection is exploratory; record the selection lineage. This fixed-selection null does not correct data-dependent discovery. A selection-aware null would need to rerun the complete selection algorithm per permutation or independent validation. New spatial results do not promote upstream labels/RCN to current.

Resume accepts only matching complete outputs. Validate margin preservation, fixed missing positions, ROI isolation, reproducible seeds, family accounting and untestable status before using real data.
