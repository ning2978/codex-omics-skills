---
name: xenium-scrna-label-transfer
description: Use when transferring scRNA-seq reference labels onto Xenium cells by shared genes, Harmony/scVI/KNN embeddings, or when auditing projected Xenium cell-type labels before downstream spatial neighborhood analysis.
---

# Xenium scRNA Label Transfer

## Purpose

Use this skill when Xenium cell labels should be derived or audited against a single-cell RNA-seq reference. It is based on the course `label_transfer.py` pattern, but with project-safe contracts and validation.

Do not use this for marker-only Leiden review; use `xenium-annotation-review` or a project annotation skill for that.

## Input Contract

Require:

- `reference.h5ad` with raw or normalized counts and `obs[label_key]`;
- `xenium.h5ad` with counts, `obsm["spatial"]`, and cell ids preserved;
- shared gene symbols after harmonizing capitalization and aliases;
- explicit output path for a labeled Xenium copy and confidence table.

Recommended minimum checks:

- at least several hundred shared genes, or an explicit warning if panel size is small;
- no duplicate `var_names` after `var_names_make_unique()`;
- reference labels have enough cells per class for KNN voting.

## Workflow

1. Intersect genes and keep counts/normalization choices explicit.
2. Preprocess reference and Xenium with the same transform where possible.
3. Integrate in a shared latent space:
   - `X_harmony` is acceptable for a light template;
   - use `scvi` if the project already has a stable scVI environment and model contract.
4. Transfer labels by weighted KNN from reference cells to Xenium cells.
5. Store:
   - `obs[predicted_label_key]`;
   - `obs[predicted_label_confidence]` as max class probability;
   - optional top-2 label margin.
6. Plot UMAP/embedding and spatial maps for transferred labels.

## Guardrails

- Never replace curated `obs["cell_type"]` without writing a separate key first.
- Do not impute/enhance Xenium gene expression from scRNA neighbors unless the user explicitly asks; keep label transfer separate from expression imputation.
- Treat rare labels cautiously; KNN voting often collapses rare cell states into abundant parents.
- If the Xenium panel is targeted, report missing reference marker genes for labels that look unstable.

## Output Contract

- labeled Xenium `.h5ad` copy;
- `label_transfer_summary.csv` with per-label counts, mean confidence, and low-confidence fraction;
- `low_confidence_cells.csv` for cells below the chosen threshold;
- UMAP and spatial PNGs colored by predicted label and confidence;
- parameters JSON recording reference path, Xenium path, gene intersection count, embedding method, K, and confidence cutoff.

## Validation

- Row count in output equals input Xenium `n_obs`.
- Predicted label categories are a subset of reference labels plus any explicit `unknown` bucket.
- Confidence values are in `[0, 1]`.
- Spatial plots show no obvious per-ROI/library artifact that dominates label assignment.
- Downstream neighborhood skills use the predicted label key only after this audit passes.
