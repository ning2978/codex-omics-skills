# Reviewer-readiness checklist

Use the sections relevant to the analysis. Record each item as `PASS`, `FAIL`, or `NOT_APPLICABLE`, with an artifact path and concise evidence.

## Study design and claims

- State the biological question, cohort, inclusion/exclusion criteria, biological unit, technical unit, and primary comparison.
- Mark each result `formal`, `candidate`, `exploratory`, or `unresolved`.
- Match wording to evidence: association, enrichment, recurrence, prediction, or validation.
- Avoid post hoc threshold selection driven only by a desired result.

## Data provenance and preprocessing

- Hash immutable inputs and record counts provenance.
- Preserve raw counts and identifiers.
- Document QC, filtering, normalization, feature selection, batch handling, dimensionality reduction, clustering, seeds, and software versions.
- Demonstrate that batch correction does not erase biology or create labels.
- Keep source objects, registries, coordinates, and prior result versions unchanged.

## Annotation

- Version and freeze marker panels before the formal run.
- Record literature/atlas support and marker availability by assay.
- Require coherent positive markers and examine competing lineages.
- Separate identity, subtype, tumor program, and state.
- Export all marker statistics and manual decisions.
- Treat transfer, UMAP appearance, and spatial proximity as supporting evidence only.

## Replication and statistics

- Perform inference at patient/donor/specimen/pair/ROI level, not cell level.
- Report the number of biological units supporting each claim.
- Provide effect size, uncertainty, multiple-testing correction, and denominator.
- Label internal recurrence and external validation accurately.
- Include sensitivity analyses for thresholds, integration, filtering, and rare groups.

## Exclusions and missingness

- Export exact excluded identifiers and coherent exclusion evidence.
- Report attrition and group-specific exclusion rates.
- Show that main conclusions are not created by the exclusion.
- Treat assay-missing markers as unmeasured.
- Retain negative and unresolved findings in the record.

## Deliverables

- Provide machine-readable sidecars, full marker tables, methods parameters, execution logs, manifests, validation summaries, and figure source data.
- Build a claim-evidence ledger linking every manuscript claim to its figures, tables, thresholds, biological units, and limitations.
- Maintain a reviewer-risk register with proposed wording or additional validation for each material risk.
- Fail publication readiness if any required evidence or provenance item is absent.
