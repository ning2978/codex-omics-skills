---
name: sc-cell-communication
description: Infer scRNA-seq ligand-receptor communication on demand using species databases, expression eligibility and sample-aware support. Use for CellChat, CellPhoneDB, LIANA or NicheNet; spatial geometry requires a spatial workflow.
---

# Question-driven communication

Read [../../workflow/workflow-sc-analysis/references/data-contract.md](../../workflow/workflow-sc-analysis/references/data-contract.md). This is an optional biological-question route, not mandatory QC.

Require reviewed identity/version, sample/donor/condition metadata, explicit unintegrated nonnegative expression, species and feature coverage. Targeted panels may omit ligand/receptor subunits; audit availability. Choose CellChat/CellPhoneDB/LIANA for ligand-receptor/consensus hypotheses, or NicheNet for receiver-target evidence with a declared response/background. R methods use the sparse bridge and exact API/environment checks.

Record database/version, ID mapping, complexes/cofactors and expression rules; resolve duplicate symbols explicitly. Missing subunits are not measured negatives. Analyze within samples with sufficient groups, summarize recurring edges and use an estimable sample model for condition differences. Record cell/expressed-fraction thresholds and missing type-sample pairs. Pooled-cell permutations are exploratory, not replicated group inference.

Record normalization/null/permutations/seeds, score definitions, tested-edge universe and FDR; examine expression/denominator/annotation sensitivity. Deliver all edges with source/target/complexes, eligibility, detection/cell counts, sample support, parameters and validation. Expression predictions do not prove contact, binding, causality or flux. Distinguish receptor predictions, receiver target evidence and orthogonal validation. Spatial proximity requires coordinates/platform methods; never infer it from UMAP.
