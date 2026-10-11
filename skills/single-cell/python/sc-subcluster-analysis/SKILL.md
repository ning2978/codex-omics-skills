---
name: sc-subcluster-analysis
description: Recluster scRNA/snRNA or segmented-cell subsets from retained counts, or coordinate iterative subclustering, marker review and program-supported automatic or requested manual cleanup rounds. Use for filter-and-rerun or an annotation-cleaning loop; preserve marker-only and unchanged-graph layout-only modes. Spots/bins use mixture interpretation.
---

# Subclustering and reviewed cleanup

Read [the data contract](../../workflow/workflow-sc-analysis/references/data-contract.md), [integration/clustering](../sc-integration-clustering/SKILL.md) and [marker annotation](../sc-marker-annotation/SKILL.md) completely. This skill owns stage transitions, exact-list decision records and stopping; those leaves own numerical methods and evidence. Repository paths, cohort codes and commands belong in a project adapter.

## Select the requested mode

- **Single expression run:** extract the requested subset, preserve counts and historical labels, rebuild the representation and review all available genes. One run does not authorize iterative exclusion.
- **Existing-cluster review / marker-only:** enter at evidence review, accurately record that clustering/UMAP were not recomputed. A plotting or marker-export request does not request new identities or deletion.
- **Layout-only:** use the integration leaf's unchanged-cell/unchanged-graph layout function, retaining clusters and markers. Finish with its validation; do not force a cleanup round.
- **Reviewed loop:** for a requested subclustering–annotation–cleanup loop, read [the round protocol and helper interfaces](references/review-loop.md). An explicit cleanup-loop request authorizes `cleanup_mode=auto_supported`: record the initial request, freeze exact lists/evidence/impact, then record an independent agent review with `auto-review` and apply every supported proposal. This is an automatic decision, never a fabricated user approval or simulated approval. Use `manual_per_round` when the user requests round-by-round confirmation; missing mode in older records remains manual. Single clustering, marker review and plotting do not activate deletion. Changed sources, lists or evidence invalidate decisions and require fresh review. Preserve unresolved/unselected cells, original artifacts, cumulative attrition and reasons; every removal requires a full retained-count Harmony rerun and fresh review. Stop at `no_supported_exclusions` only with no unresolved groups; otherwise retain doubts with `evidence_insufficient`; technical failure is separate.
- **Spot/bin data:** route to mixture/anatomy interpretation through spatial workflow. A small bin is not evidence of an individually segmented cell and cannot enter this cell-exclusion loop.

## Expression run contract

Identify count artifact/provenance, target lineage, species/tissue, platform, unique cell/feature keys, sample/capture, genuine biological unit, condition and selected batch. Validate counts and alignment; count-like X alone is not proof of raw counts. Preserve the source H5AD and state prior gene restrictions. Missing unit metadata permits descriptive review with an explicit limitation, not invented replication.

For every changed cell set or graph, start from retained `layers["counts"]`. Rebuild normalization, HVGs, PCA, mandatory Harmony, neighbors, UMAP and Leiden. Clear old `X_umap`, graphs and the target cluster key first. Old/global coordinates may appear only as explicitly requested historical projections, never the main result. Keep full available-gene counts/log expression separate from the scaled HVG work object.

Use the integration leaf's matrix-aware HVG fallback, baseline/integration comparison, explicit graph/representation and resolution/seed assessment. Harmony is required for this skill’s single expression runs and every cleanup rerun; retain unintegrated PCA/graph diagnostics, but use X_harmony for the main neighbors, UMAP and Leiden. Declare batch_key per analysis (never automatically use a biological condition). Require complete nonempty values and at least two retained batches; failure stops as technical_failure with no PCA/other-method fallback. Refit Harmony after every removal. No tissue-specific marker panel, fixed cluster IDs or universal resolution. A resolution-only sweep requires unchanged cells **and** graph; RCN/CN, neighborhood NMF and CellCharter are spatial composition/niche clustering, not expression reclustering.

Review complete available-gene markers (default t-test, unintegrated/unscaled log layer, explicit detection filter) and independent sourced positive/competing programs, including measured genes outside top20. Top20 is display only. Inspect panel availability, detection/specificity, same-cell co-detection, biological-unit recurrence and platform QC. Keep identity, subtype, state, display label and evidence grade separate; preserve candidate/review/unresolved outcomes and historical names.

## Cleanup decisions

Propose one action per observed cluster: retain, supported foreign-lineage exclusion, mixed whole-cluster exclusion, or unresolved review. Never decide removal from one foreign marker, shared cytotoxicity, cycle/stress, single-sample dominance or boundary-stain segmentation alone.

For segmented cells distinguish segmentation doublet/boundary error, broad-compartment misassignment, adjacent transcripts and genuine mixed biology. Two same-lineage foreign markers are screening evidence only. Removal requires a sourced, lineage-specific foreign program with multiple complementary markers, measured coverage and audited feature mapping; there is no universal gene-count or score cutoff. Review target programs, cell-level distributions and coordinated coexpression, biological-unit support/limitations and platform QC. Neither averages, module scores, UMAP nor seed instability alone justify deletion. Do not require every cell to express every program gene or treat panel/detection gaps as negative. Dual program positivity alone is not contamination. If genuine target subtype, biological mixed state or potential target-cell loss remains unresolved, retain as `review`. Whole-cluster removal needs supported technical mixture or an explicit outside-target-scope interpretation and quantified target-positive count/fraction and potential loss; partial-cell evidence cannot be expanded to a whole cluster. Retain or use the applicable at-most-one mixed-parent split and review a supported foreign child.

A supported mixed **whole-cluster** proposal can be applied automatically in an explicitly requested loop after its frozen round, namespace, cell list and target-loss assessment pass automatic review. It is a deliberate scope decision, not proof that every removed cell is foreign. Outside this explicit branch retain applicable project rules, including at most one mixed-parent subclustering and only a supported foreign child when required. No silent repeated splitting or cleaning.

Freeze evidence and action lists before the automatic review or requested manual approval; present cluster, sample and cumulative attrition. Apply all automatically supported IDs, or only manually approved IDs in manual mode, into a new counts derivative. Unselected/unresolved cells remain. Any source, cluster allocation, list or evidence change invalidates the old decision. After removal, recompute and review again before evaluating convergence.

## Completion and stopping

For expression runs reopen and verify immutable counts/keys, full log expression, actual HVG/PCA/X_harmony with parameters and representation fingerprints, named Harmony graph, finite current 2D UMAP and clusters. Retain `uns["expression_recluster_provenance"]` with `recomputed_hvg/pca/harmony/neighbors/umap/leiden=true`, supported by count/key/graph fingerprints, stage inputs, commands, parameters and runtime; flags alone are insufficient. Run project guards where supplied.

Each loop round retains complete markers/top view, panel/co-detection/unit/QC evidence, proposed actions, exact cell manifest, automatic decision or manual approval when applicable, all-gene annotated H5AD, HVG work H5AD, keyed sidecar, PNGs/palette and a validation manifest. An existing-cluster review may retain its source artifact by immutable reference and declare missing historical stages instead of claiming recomputation. A marker-only task needs only its requested evidence outputs.

The cumulative ledger covers every initial cell, historical round namespaces, excluded round/reason and final destination. Technical acceptance and scientific stopping are separate. Stop with `no_supported_exclusions` only after review finds no supported removal and no unresolved groups. Keep unresolved groups with `evidence_insufficient`; distinguish `user_stopped`, `awaiting_approval`, and `technical_failure`. Never turn a PASS or empty top20 into biological convergence. Default new labels to candidate/review; no registry promotion is implied.

## Internal completion and user delivery

Read [the shared delivery contract](../../workflow/workflow-sc-analysis/references/delivery-contract.md)
when preparing a complete expression-analysis handoff. Completion artifacts and
validation above stay in the working directory outside `delivery/`; select final
UMAP PNGs, marker DotPlot PNGs, per-cluster annotation review CSV and README only.
Standalone plotting, marker export and review deliver only requested content.
Delivery selection does not authorize new analysis or reduce evidence requirements.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
