# Final expression-analysis delivery contract

Read when preparing the final user handoff for a complete scRNA/snRNA or
segmented-cell Xenium reduction, clustering, annotation or subclustering flow.
This contract governs future delivery selection only. Preserve historical
results. It changes no annotation source, evidence gate, reclustering rule or
method-specific validation. Spatial workflows apply it only to segmented-cell
Xenium expression analysis, not spots/bins, geometry, RCN/NMF or niche analysis.

## User delivery

Create a separate `delivery/` within the new analysis working directory. By
default it contains only these four content types (filenames are examples):

| Content | Example | Required meaning |
|---|---|---|
| UMAP PNG | `umap_final.png` | Final cluster and annotation panels; a combined two-panel image is acceptable. Split lineages into separate PNGs when useful. |
| DotPlot PNG | `dotplot_markers.png` | Representative measured markers supporting the annotation; cluster grouping matches the review CSV. Use the assay's renderer and existing expression/statistics contract. |
| Cluster annotation review CSV | `annotation_review.csv` | One row per final cluster in its namespace, including exact cell count, identity/subtype, state, annotation grade, positive/competing markers, sample/unit support, doubts and final judgment. |
| Explanation | `README.md` | Data/cell scope, method and key parameters, annotation rationale, unresolved issues and paths to internal artifacts. |

Do not include spatial plots, proportions, heatmaps, extra statistical tables,
logs or duplicate reports by default. Add extra deliverables only on explicit
request. Standalone plotting, marker export or review tasks deliver only what
was requested; do not force the four-content bundle or start new analysis to
complete it. A failed or incomplete flow must be reported as such, never filled
with fabricated annotations or a claimed final bundle.

Recommended CSV columns: `cluster_namespace`, `cluster`, `n_cells`,
`cell_identity`, `cell_subtype`, `cell_state`, `annotation_status`,
`supporting_markers`, `competing_markers`, `sample_support`, `uncertainties`,
`final_judgment`. Equivalent explicit columns are acceptable. Namespace includes
the lineage/round/cluster key needed to disambiguate reused cluster numbers.
Identity, subtype, state and evidence grade remain separate. Preserve candidate,
review, exploratory, unresolved, unknown/ambiguous and pending-approval values;
do not rename or fill them into formal labels. Use explicit unavailable evidence
and support limitations, not invented negative evidence or biological replicates.

README records the exact final cell set/cluster key/namespace and grouping,
counts/log-expression source, normalization, HVG/PCA, integration/batch key,
neighbors, resolution/seed and UMAP parameters as applicable. Name annotation
sources/versions and supporting/competing evidence without bypassing registry
resolution. Explain detection/panel gaps, unresolved alternatives and technical
versus scientific stopping. For reviewed cleanup loops summarize retained and
excluded counts, cumulative attrition/reasons and stopping state; link internal
round evidence instead of delivering duplicate reports. List internal paths with
an explicit working-directory base or usable relative paths.

## Internal completion artifacts

Keep all-gene counts/annotated H5AD, HVG work objects, exact keyed per-cell
sidecars, complete markers, panel/evidence audits, sample/design tables,
exclusion and cumulative ledgers, baseline/integration/resolution comparisons,
palettes, logs, run/validation manifests and round decisions in the working
directory outside `delivery/`. Every existing method completion check still
applies; a concise delivery is not evidence that those checks passed.

Retain the explicit union label-to-hex palette beside the internal figure
collection and reference it from README. Delivered PNGs may be copies of those
validated figures; do not put palette or manifest JSON inside `delivery/`.
Packaging selects final validated artifacts and never triggers normalization,
marker discovery, clustering, annotation promotion or historical migration.
Do not erase intermediate artifacts to make a directory look concise.

## Delivery acceptance

- Inventory `delivery/`: only final UMAP PNGs, marker DotPlot PNGs, review CSV(s)
  and one README, unless extras were explicitly requested. PNGs must be nonempty
  and decode; obey the existing PNG-only and paired-PDF policy.
- Validate against the final keyed cell allocation, not filenames: cluster
  namespaces/labels match across UMAP, DotPlot and CSV; exact CSV cell counts
  sum to the final scope. For split lineage panels/tables reconcile the union;
  declare intentional subsets. DotPlot annotations must distinguish clusters
  even when several clusters share one identity. Record exact cell-set hashes,
  grouping and figure input provenance in the internal validation manifest.
- Check palette coverage and unique normalized colors over the union of labels;
  use the same mapping for plotted points and legend handles, with separate
  namespaces for cluster numbers. Inspect PNGs for readability/legend agreement.
- Check that unresolved/candidate/pending values and material limitations appear
  consistently in the review CSV and README. Structural acceptance does not
  validate biological identity or promote an annotation.
- Record the inventory, input/output hashes, checks and any explicit exceptions
  internally. Final replies link only these four delivery content types;
  internal artifacts are traceable through README.
