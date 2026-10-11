---
name: sc-marker-annotation
description: Discover all-available-gene scRNA-seq cluster markers from explicit unintegrated log expression and build versioned identity, subtype and state annotations with competing-lineage and sample evidence. Use for marker review, mapping or annotation updates.
---

# Markers and annotation

Read [../../workflow/workflow-sc-analysis/references/data-contract.md](../../workflow/workflow-sc-analysis/references/data-contract.md). Cluster marker ranking is exploratory; group inference routes to `sc-sample-level-comparison`.

For biological annotation or review, read [references/annotation-evidence-review.md](references/annotation-evidence-review.md). It connects differential ranking, canonical/competing programs, biological-unit support and assay-appropriate QC. Marker export or plotting alone does not expand into a new annotation or exclusion task.

## Matrix and complete table

- Use **all available genes** passing a declared detection filter (e.g. >=3 retained cells, adjusted for tiny subsets). Audit canonical panel availability before filtering and preserve low-detection panel evidence separately. State prior gene restrictions; HVGs only construct representations.
- Explicit `rank_genes_groups(..., layer="log1p", use_raw=False, method="t-test", pts=True, n_genes=n_available)` or a disposable all-gene log object. Never use centered/scaled `X`, integrated expression or implicit `.raw` for markers/detection. Shared `rank_markers` enforces declared unscaled expression.
- t-test is a fast exploratory default; Wilcoxon when requested/justified. Neither makes cells biological replicates. Scanpy logFC is an approximation from mean log expression; flag undefined/low-detection statistics.
- Export every tested gene/cluster: stable ID/original symbol, score, logFC, p/FDR, pct_in/pct_out, group sizes, method, matrix, scope and filter. Detection uses `>0` in declared nonnegative unscaled expression.
- Optional legacy top20 (`pvals_adj<0.05`, `pct_in>0.3`, logFC order) is a configurable display. Review specificity/pct_out/effect and canonical low-expression markers beyond this view.
- Use FindAllMarkers-style cluster-versus-rest comparison with a complete tested-gene table and separate candidate filters. This does not change the existing Python default or require a second backend; native Seurat prefilters and method/version must be declared when requested.

## Evidence and labels

Review coherent positive and competing programs independently of the DE candidate cutoffs, recording their tissue/species/source and measured availability. Combine complete markers, detection/specificity, biological-unit recurrence and applicable QC under the linked review contract. Shared lineage markers need not be differential to support identity; unmeasured genes are not negative. Stress/cycle/interferon/activation are states, not new lineages by themselves.

Keep `cell_identity`, `cell_subtype`, `cell_state`, `annotation_status` and `annotation_version` separate. Permit `unknown/ambiguous`, candidate labels and competing explanations. Do not force coverage or reuse old labels blindly. Reference mapping records version/species/tissue/coverage, feature overlap, calibrated confidence/thresholds and unseen-type handling; compare with de novo evidence. Mapping score alone cannot establish identity.

Save per-cluster positive/support/competing genes and program sources, availability/specificity/detection, unit counts/support, QC, confidence/status, alternatives and decision. Missing checks remain explicit; single-unit dominance is a limitation, not automatic rejection. Annotation updates create a new version preserving history. Registry mutation must be within requested annotation scope and satisfy its own contract.

## Contamination and plots

Combine coherent foreign markers with target-positive, unit and platform-specific QC evidence. When needed, distinguish a cluster containing separate populations from same-cell co-detection. Segmented Xenium requires explicit segmentation/misassignment/adjacent-transcript/mixed-state alternatives; foreign expression alone does not establish ambient contamination. Export proposed actions/exact keys; evidence does not authorize iterative deletion. Authorized removal routes to `sc-subcluster-analysis` and recomputation.

Expression panels explicitly use log layer/`use_raw=False`; duplicate symbols require audited feature mapping. Publication DotPlots use `render-scrna-dotplot` (Xenium its renderer) with Seurat statistics/parity. Scaled heatmaps are labeled relative-signal visualization, separate from detection/statistics.

Retain internally complete markers/top view, panel audit, versioned evidence and sample/unit summaries. Include requested plots and a full-count annotated H5AD when annotation writing is in scope; review-only work can finish with evidence tables. Reopen/check written cell alignment, counts, labels/status and fraction provenance; inadequate clusters remain unresolved. Hand off identity/version/design for condition comparisons.

## Evidence helper and reviewed-loop handoff

For a portable structured review, read [the evidence interface](references/evidence-interface.md) and use [collect_evidence.py](scripts/collect_evidence.py). It reuses shared numerical tools to export full markers, panel availability, expression/detection, same-cell co-detection, sample/unit support and QC. It never labels or removes cells. Save display name separately from identity/subtype/state and evidence grade.

For an explicitly requested cleanup loop, `sc-subcluster-analysis` owns freezing, automatic program-supported review, filtering, recomputation and stopping. Use manual per-round confirmation when requested; old records lacking a mode stay manual. Ordinary marker review does not authorize deletion. Apply the complete foreign-program and target-loss criteria in the linked joint evidence review: two screening markers or dual positivity alone are insufficient. Supported mixed whole clusters may be removed within that loop; unresolved target biology or loss stays `review`. Evidence collection itself never deletes cells or equates technical PASS with scientific convergence.

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
