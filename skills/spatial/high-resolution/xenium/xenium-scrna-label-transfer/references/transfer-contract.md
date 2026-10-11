# Gene mapping, predictions and review

Use a species/version-specific table with stable source gene ID, source symbol and canonical symbol. Alias maps have `alias,canonical,species`. Case folding is not universal (human/mouse conventions differ). Audit ambiguous alias-to-many and many-to-one source-gene collisions. Do not aggregate counts without an explicit biological rule; never collapse unrelated genes to repair a collision. Shared features are selected only after the one-to-one map passes. Keep stable IDs and unmodified full measured expression outside the transfer representation.

Declare reference and target matrix sources (X/layer), counts vs normalized/log1p, library-size transform, batch key and embedding provenance. Counts must be nonnegative integer-valued; integrated or scaled values are not measured marker evidence. Harmony/KNN and scVI are method choices, not automatic prerequisites. Honor project-specific scVI authorization/runtime gates. For an existing prediction audit, no new embedding or reclustering is needed.

Weighted KNN class sums divided by total positive neighbor weight yield `knn_vote_fraction`; top-two difference yields `knn_vote_margin`. Neighbors are not independent draws. Rare-class coverage, training label quality and reference distance affect votes. Calibrated confidence requires labels not used in fitting, held-out donors/batches, a specified calibration method and reliability/error by class. A score alone cannot establish formal annotation.

Review retains `predicted_label`, votes/margin, optional calibrated confidence/status, `review_label`, `review_status`, program sources/availability, target and competing support, biological units, applicable QC and alternatives. Broader labels require measured support. `unknown` is a valid conclusion; never fill it from the prediction or rename it as formal. Boundary-stain alone is not exclusion evidence.

Helper config for `python <skill-dir>/scripts/audit_transfer.py --config <json> --out-dir <new>`:

```json
{
  "species":"human", "stable_gene_key":"gene_id", "symbol_key":"symbol",
  "reference_genes_csv":"/path/reference_genes.csv", "target_genes_csv":"/path/panel.csv",
  "aliases_csv":"/path/aliases.csv", "cell_key":"obs_name",
  "votes_csv":"/path/votes.csv", "evidence_csv":"/path/evidence.csv",
  "min_vote_fraction":0.7, "min_vote_margin":0.2,
  "upstream":[], "label_evidence_grade":"prediction_review"
}
```

Aliases and supplied votes/evidence are optional paired modes. Gene tables contain stable ID and symbol columns; votes contain unique cell key, `reference_label,weight`. Evidence has every target cell exactly once, `program_audit_status` (supported/conflict/unavailable) and optional `supported_broad_label`. `supported` here records independently supplied program evidence, not biological validation by the helper. Source evidence columns are retained alongside derived review fields. The helper validates mapping collisions and unknown/broader handling without inferring new curated annotations.
