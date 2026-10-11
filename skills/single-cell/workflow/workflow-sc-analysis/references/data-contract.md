# Shared scRNA-seq data and stage contract

Read for the conventional workflow, its core leaves and subset reruns. Full-transcriptome scRNA-seq is the default; snRNA-seq uses a nuclear QC branch. Targeted spatial assays retain platform QC and gene-availability rules.

## Identity and experimental design

- Preserve species, reference assembly/annotation, chemistry, assay, count-generating software/version and filtered/unfiltered origin in `input_manifest.json`.
- `obs_names` are unique stable keys from **physical capture library + original barcode**, with a reversible map. Samples/donors can have multiple captures; never join naked barcodes across libraries.
- `var_names` are unique stable feature IDs. Keep `gene_id`, original `gene_symbol`, feature type and source row. Duplicate symbols remain distinct features; do not sum/select one silently. Record explicit aggregation.
- `sample_table.tsv` distinguishes capture/library, specimen/sample, independent biological unit, donor, technical batch, condition, pairing/time and covariates. Technical captures are not biological replicates. Declare the experimental unit for each contrast; donor/specimen differ in repeated designs.
- Missing biological-unit metadata allows exploration with `inference_status=descriptive_only`, but blocks between-condition inference. Do not invent donor IDs from cells.

## Matrix contract

| Location | Meaning and consumers |
|---|---|
| `layers["counts"]` | Original finite nonnegative raw counts from an identified artifact. Integer-like for UMI count models; numerical resemblance alone is not provenance. Never overwrite. |
| `layers["counts_corrected"]` | Optional ambient correction with method/version/parameters; may be fractional. Never round into purported raw counts. Choose a supported model or original counts explicitly. |
| `layers["log1p"]` | Unintegrated library-size normalized natural log1p, default target sum 1e4, from a declared count source. Explicit marker/expression plot input. Other normalization needs its own named contract. |
| HVG work object's `X` | Scaled/centered expression for PCA only, separate from the all-gene completion artifact. Never use for markers, detection fractions, counts or fold changes. |
| `obsm` / `obsp` | Representations and neighbor graphs; integrated representations are not expression matrices for DEG. |
| `.raw` | Optional snapshot; never rely on implicit `use_raw` selection. |

Store `uns["matrix_sources"]` and a serializable manifest. Each matrix declares kind, source artifact/hash and upstream matrix, cell/gene key hashes, normalization/log base/correction and gene filtering. Count recovery requires verified key alignment. Unknown raw-count provenance prevents declaring a counts-first rerun complete.

Marker scope is **all available genes passing a declared expression filter**, independent of HVGs. State prior gene restrictions. A top20 table is a display view, not annotation's gene universe. Detection uses expression **>0** in the declared nonnegative unscaled log layer (or explicitly declared raw counts with separate provenance). It equals original detection only when preprocessing preserves zeros; correction can change it.

## Stage evidence

Each stage writes a run/validation manifest with schema version, stage, `PASS/WARN/FAIL/SKIP/NOT_APPLICABLE`, actual commands, exact executables/package paths/versions, parameters, seeds, input/output paths and hashes, cell/gene key hashes, matrix sources, upstream stage dependencies, graph/representation and check evidence. Missing/ skipped checks never become PASS. A `recomputed=true` flag is insufficient without matrix identity and upstream evidence.

Reopen outputs and verify immutable counts, key order, dimensions, finite embeddings, graph shape/symmetry, labels and parameters. After cell deletion rebuild normalization, HVG, PCA, chosen integration, neighbors, Leiden and UMAP from retained counts; invalidate prior analysis or retain under historical names.

A **layout-only** UMAP change reuses an unchanged graph, preserves clusters and gene statistics, and records `mode=layout_only` plus graph/cell hashes. A resolution-only sweep likewise requires unchanged cells/graph. Neither is a counts-first rerun.

## Numerical helpers

[../scripts/sc_contracts.py](../scripts/sc_contracts.py) provides collision-safe keys, strict table alignment, log preparation, matrix-aware HVG fallback, full-gene markers, design gates, pseudobulk aggregation and unchanged-graph UMAP layout. Import by resolved skill path. It does not select biological labels, QC thresholds, integration or experimental design. Runners populate provenance and stage evidence. Behavioral regressions are bundled alongside it; cross-language tests additionally require an exact Seurat R environment.

Large H5AD preflight uses direct HDF5 chunks/tiles. `read_h5ad(backed="r")` is not a memory guarantee because versions may load layers eagerly. Sampled numerical checks never prove full-file integrity.

## Delivery boundary

Read [delivery-contract.md](delivery-contract.md) for the final user handoff.
Stage manifests and full analysis objects here are internal completion artifacts,
not default `delivery/` contents. Their validation requirements remain unchanged.
