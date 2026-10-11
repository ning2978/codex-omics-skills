---
name: sc-sample-level-comparison
description: Compare scRNA-seq conditions using genuine biological units, raw-count edgeR pseudobulk, sample cell proportions and documented enrichment or scores. Use for group inference and descriptive results when replication or estimability is missing.
---

# Biological-unit comparisons

Read [../../workflow/workflow-sc-analysis/references/data-contract.md](../../workflow/workflow-sc-analysis/references/data-contract.md) and [references/statistical-handoff.md](references/statistical-handoff.md) for R models. Python prepares data/design; R runs edgeR/propeller. Cluster markers are not group DEG.

## Design gate

Declare experimental unit/contrast before tests: donor/specimen/pair/timepoint, technical captures and covariates. Collapse technical captures from the same unit/condition; repeated donor samples need blocking/correlation. Tabulate units, pairing, batches, missingness and type-specific eligibility; check batch-condition confounding, actual design rank and residual df.

Missing independent replication (including one donor with many cells), nonestimable/full-confounded contrast or no residual df -> `inference_status=descriptive_only` with reasons, summaries and no inferential p/FDR/group claims. Two independent units/group is a minimal technical gate, not sufficient power. Missing biological metadata permits exploration. Paired models include donor and require supported pairs per type; longitudinal/repeated designs may need other models. Never drop confounded covariates silently to make a test run.

## DEG default: edgeR pseudobulk

Sum **original raw counts** by cell type x biological unit x condition, collapsing technical captures only where the unit definition allows. Preserve stable IDs/symbols and cells per aggregate. Declare minimum cells/units and gene filters. A sample without cells of a type is not an all-zero transcriptome replicate. Check each type's design separately.

Use `DGEList`, `filterByExpr(design=...)`, TMM `calcNormFactors`, `estimateDisp`, robust `glmQLFit` and explicit `glmQLFTest` contrast. Check installed APIs/orientation. Export all tested genes, logFC direction, logCPM/statistic/p/BH FDR, aggregate counts, factors, design/contrast/df and exclusions. Include unit expression plots, model-supported effect uncertainty and donor consistency; inspect influential samples. Integrated/scaled/log expression is forbidden pseudobulk count input.

The bundled `scripts/edger_runtime_smoke.R` runs a real synthetic paired QL fit using the exact `Rscript --vanilla`; it checks APIs/library paths and contrast direction, not dataset power or validity.

## Proportions default: unit plots + propeller

Count each type by independently sampled unit with a declared denominator (all retained cells or defined parent). Include zeros for types absent from genuinely sampled units; distinguish no sampling/denominator. Plot individual units and paired lines.

Use speckle/propeller with recorded arcsin/logit transform and limma design, donor blocking for pairs. Convenience `propeller()` serves simple independent groups; use transformation helper + limma for paired/covariate designs. Export denominators/counts, eligibility, raw-proportion effects/uncertainty, transformed inference and p/FDR. No pooled-cell Fisher/chi-square replacement. Proportions are compositional and capture-dependent; evaluate denominator/rare-type sensitivity or remain descriptive.

## Enrichment and scores

- ORA: directional selected-gene rule, background = actually tested eligible genes.
- GSEA: full signed sample-model ranking, no significance prefilter; record statistic, ties, duplicate-ID mapping and direction.
- Scores: gene-set/species/version/availability, algorithm and control/background parameters; compare independent-unit summaries with suitable paired/repeated models.

Record ID conversion, database/version/date, set sizes, universe, contrast, multiple-testing family/BH FDR and complete results. Enrichment alone does not establish causality/pathway activation.

Deliver design status/reasons, aggregate/unit counts, complete tables/effects/plots, software/parameters and machine-readable checks. Descriptive runs still deliver unit plots/limitations. Formal manuscript claims additionally use `publication-grade-omics-analysis`.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
