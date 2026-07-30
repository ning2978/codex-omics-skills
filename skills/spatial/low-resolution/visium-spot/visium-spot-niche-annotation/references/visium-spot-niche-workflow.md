# Visium Spot Niche Annotation Workflow

## Minimal Procedure

1. Load the clustered h5ad or spot table.
2. Confirm `obs[label_key]`, `obsm["spatial"]`, `obs["sample"]`, and sample/timepoint metadata exist.
3. Export one-vs-rest markers per niche.
4. Rank positive markers by:
   `pvals_adj` ascending, `logfoldchanges` descending, `pct_ratio` descending, `scores` descending.
5. Export `top100` and `top200` marker tables.
6. Score a signature panel appropriate to the tissue.
7. Inspect the all-sample and per-sample spatial maps on HE.
8. Assign broad anatomical compartments before refined biological labels.
9. Write a manual-review table rather than overwriting h5ad labels unless the user explicitly approves.
10. Update the project context document after result generation.

## Mixed-Niche Comparison

For comparing fixed K values, compute:

- same-label spatial neighbor fraction from `obsp["spatial_connectivities"]`
- signature margin: best signature score minus second-best signature score
- signature entropy across tissue signatures
- dominant timepoint/condition fraction
- number of low-confidence or mixed/interface niches

Recommended decision rule:

- Prefer a K that separates major anatomical compartments such as epithelium, stroma, gland/appendage, vascular/perivascular, immune-rich, tumor/lesion, and artifact regions without creating many low-confidence fragments.
- Penalize a K if major anatomical regions are split into multiple nearly identical signatures without a clear spatial boundary.
- Penalize a K if one cluster contains strong competing epidermal, fibroblast, vascular, and blood signatures.

## Common Compartment Labels

- `normal epidermis`
- `activated epithelium`
- `dermis-interface`
- `fibroblast-rich dermis`
- `collagenous matrix dermis`
- `contractile/perivascular dermis`
- `vascular/inflammatory dermis`
- `appendage/sweat-gland epithelium`
- `artifact/low-confidence`
