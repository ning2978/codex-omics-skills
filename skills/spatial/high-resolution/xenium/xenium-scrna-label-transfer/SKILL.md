---
name: xenium-scrna-label-transfer
description: Transfer or audit scRNA reference labels on segmented Xenium cells using explicitly mapped measured genes, separate KNN votes/calibration, and measured-program review with broader or unknown conclusions.
---

# Xenium label transfer and audit

Keep predicted labels, measured-expression review and curated annotation in separate fields. Marker-only cluster review routes to `sc-marker-annotation`. Read [the joint evidence review](../../../../single-cell/python/sc-marker-annotation/references/annotation-evidence-review.md) for biological decisions and [transfer mechanics](references/transfer-contract.md) for gene mapping, votes, calibration and handoff.

Require unique target/reference cell keys, species, explicit reference labels and evidence grade, explicit source matrices/transforms, stable gene identifiers and an audited alias table if used. Geometry is required only for spatial audit/maps. Project adapters resolve current native/reference/cross-modality roles separately; a paused role remains paused.

1. Audit gene IDs/symbols before intersection. Preserve stable IDs and full measured target counts; do not repair collisions with `var_names_make_unique`, capitalization or automatic alias collapse.
2. Choose an integration/transfer representation appropriate to measured panel coverage. Benchmark against a simple baseline and held-out biological units. Do not use a minimum shared-gene count as a universal quality gate.
3. Record weighted KNN predicted label, vote fraction, top-two margin, K and weighting. A vote fraction is not calibrated probability. Store calibrated confidence only with documented out-of-unit calibration and reliability/error estimates.
4. Combine predictions with full measured canonical/competing programs, biological-unit support and segmentation QC. Conflicts remain review flags; unsupported fine labels may receive a supported broader label or `unknown`. Missing panel genes are unavailable evidence.
5. Export keyed predictions, mapping audit, review decisions, coverage/conflict/unknown summaries and provenance. A new labeled H5AD is optional and preserves cell keys/full matrices. Curated fields are never overwritten.

The optional [audit helper](scripts/audit_transfer.py) audits gene mapping and supplied KNN votes/evidence; it does not fit embeddings or validate the biological evidence supplied to it. Its thresholds must be chosen on held-out units rather than on the target's desired composition.

Before downstream use, verify unchanged cell keys, no gene-mapping conflicts, finite vote fractions/margins in [0,1], explicit calibration status, and persistence of conflicts/unknown. Spatial summaries remain exploratory where replication or independent validation is missing.
