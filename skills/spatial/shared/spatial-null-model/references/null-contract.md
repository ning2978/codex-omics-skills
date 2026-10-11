# Null contract and portable runner

The conditional null asks whether target-label arrangement differs from arrangements with the same label counts within declared ROI/block/strata, conditional on fixed anchors and geometry. Choose blocks larger than the tested local reach when scientifically justified and compare plausible scales. Blocks too small freeze labels; blocks too large may erase tissue heterogeneity. Stratification is justified by the null question, not by a desired p value. Preserve the entire declared candidate pool; selecting only observed target cells leaves nothing meaningful to permute.

Config:

```json
{
  "cells_csv":"/path/cells.csv", "cell_key":"obs_name", "roi_key":"roi",
  "biological_unit_key":"donor", "x":"x", "y":"y", "coordinate_unit":"um",
  "label_evidence_grade":"reviewed_operational", "target_label_key":"phenotype",
  "anchor_mask_key":"anchor", "candidate_mask_key":"candidate",
  "block_size":200, "block_origin":[0,0], "strata_keys":["zone"], "radius":50,
  "permutations":999, "seed":0,
  "families":{"primary":["near_T","enrich_T"]},
  "hypotheses":[
    {"id":"near_T","family":"primary","metric":"proximity","target":"T","alternative":"greater"},
    {"id":"enrich_T","family":"primary","metric":"local_enrichment","target":"T","alternative":"greater"}
  ],
  "upstream":[], "selection_provenance":{"status":"exploratory_same_data"}
}
```

Masks must be nonmissing bool/0/1; ROI IDs uniquely identify geometry and nest in a biological unit. Pool members with missing target labels retain missing positions and do not supply known-label denominators. Members with missing strata are recorded and excluded from conditional evidence. Candidate membership itself remains recorded/fixed. Anchors can be outside the contributor pool; self is excluded. If anchor phenotype must also remain fixed, declare anchor cells outside the permuted candidate pool.

Statistics per anchor:

- `proximity`: any known target within radius, 0/1.
- `cooccurrence`: number of known target neighbors within radius (an edge/neighbor-count summary, not direct contact).
- `nearest`: nearest known target distance capped at radius; absent target contributes radius. Less indicates stronger proximity.
- `local_enrichment`: known-target fraction among local known neighbors minus known-target fraction among same-ROI candidates (excluding self). No known local neighbors gives missing, not zero enrichment.

Mean within ROI, equal ROI means within biological unit, then equal unit means gives the reported observed/null statistic. `biological_unit_support.csv` gives unit observed values, evaluable anchors and nested ROI counts; it is descriptive support, not independent confirmation of selected hypotheses. Missing all support yields `NOT_TESTABLE`.

For greater/less tails, count null statistics at least as extreme including ties. For `two-sided`, require a finite predeclared `null_center` and count absolute distance from it; do not estimate a desired center after seeing results. Every empirical p is `(b+1)/(B+1)`. BH uses the entire declared family size, keeping untestable hypotheses missing (equivalent to conservative p=1 family accounting). Never silently shrink families to successful tests. Constant target margins or an invariant simulated statistic yield `NOT_TESTABLE` with reason.

The null distribution NPZ preserves hypothesis order/observed values; CSV summaries include effect, null mean/SD/quantiles, extreme count, p and q. Margins and permutation-unit diagnostics are exported. The runner validates margins and missing positions during every permutation. Complete runs replay by fingerprint; partial/mismatched/corrupt outputs fail. It holds the null array in memory (B x hypotheses) and has no incremental checkpoint; use a bounded pilot and a fresh directory after interruption.

Group inference handoff: genuine independent units, condition/exchangeability blocks, nesting, replication, confounding/paired design, sample-level statistic, predeclared hypotheses/families, unit-level permutations or suitable hierarchical model, support and sensitivity. No group-permutation tool is implemented here. Tissue label randomization p values do not test between-condition differences.
