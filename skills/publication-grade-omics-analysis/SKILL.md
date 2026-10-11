---
name: publication-grade-omics-analysis
description: Apply reviewer-defensible gates to omics planning, annotation, statistics, exclusions, validation, figures, and reporting. Use for formal biological conclusions or submission/publication deliverables. Ordinary exploratory cross-modality comparisons do not trigger the full gate.
---

# Publication-grade Omics Analysis

## Core contract

Prioritize scientific validity over positive-result yield or annotation coverage. Do not imply that publication is guaranteed. Make every formal claim traceable to data, prespecified rules, biological replication, reproducible code, and explicit limitations.

Read [references/reviewer-readiness-checklist.md](references/reviewer-readiness-checklist.md) before finalizing an analysis plan, launching a manuscript-intended run, or promoting a candidate result to a formal label or claim.

## Apply claim gates

Classify each result as `formal`, `candidate`, `exploratory`, or `unresolved`.

- Promote a result to `formal` only when its assay can measure the required evidence, its prespecified quantitative gate passes, and it recurs across the declared biological units.
- Use `candidate` for coherent but incompletely replicated, platform-limited, or single-case evidence.
- Use `exploratory` for data-derived hypotheses without confirmatory evidence.
- Preserve `unresolved` when evidence is insufficient or contradictory. Never relax a threshold solely to reduce this category.
- Record plausible competing explanations and evidence against them.

## Protect statistical validity

- Declare the biological unit before testing. Use patient, donor, specimen, pair, or ROI summaries as appropriate; do not treat cells, spots, or transcripts as independent replicates.
- Distinguish technical replicates from biological replicates.
- Prespecify primary outcomes, thresholds, covariates, contrasts, and exclusion rules.
- Report effect sizes, uncertainty, multiple-testing correction, and biological-unit recurrence. Do not rely on p-values alone.
- Separate discovery from confirmation. If no independent cohort exists, describe cross-sample recurrence as internal replication, not external validation.
- Avoid causal language for observational results unless the design supports it.

## Make annotation defensible

- Freeze a versioned marker panel before the formal run and record primary-literature support or an accepted atlas reference for each formal subtype.
- Record marker availability per modality; missing genes are unmeasured, not negative.
- Require coherent positive evidence. Do not let a single ambient, stress, state, or broad negative marker erase a lineage.
- Keep identity, differentiation subtype, tumor program, and transient state in separate fields.
- Treat UMAP separation and label transfer as supportive visualization/evidence, never as sufficient annotation evidence.
- Export full marker statistics, configured thresholds, positive/support/negative hits, competing programs, recurrence, and the final decision for every cluster.

## Control exclusions and threshold changes

- Never delete cells or clusters silently. Export exact identifiers, evidence, rule, before/after counts, and a retained-data sensitivity analysis.
- Require at least two coherent markers from the same foreign lineage before excluding a target-lineage cluster; subcluster mixed populations once before removing only a supported contaminant child.
- Keep a strict analysis whenever a relaxed operational threshold is used. Report strict-versus-relaxed confusion and whether conclusions change.
- Fail closed when a required field, marker, biological unit, input hash, or validation artifact is missing.

## Preserve reproducibility

Record input paths and hashes, counts provenance, software environment, package versions, random seeds, all parameters, commands, resumable stage status, output hashes, and validation results. Protect source data, prior versions, registries, and spatial coordinates from mutation.

Produce a manuscript handoff containing:

- a claim-evidence ledger;
- methods-ready parameters;
- biological-unit summaries;
- strict/relaxed sensitivity results;
- exclusions and attrition;
- limitations and alternative explanations;
- a reviewer-risk register;
- a pass/fail publication-readiness checklist.

Do not call the analysis publication-ready while any required gate fails.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
