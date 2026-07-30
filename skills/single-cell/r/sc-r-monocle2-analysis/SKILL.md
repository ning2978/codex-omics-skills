---
name: sc-r-monocle2-analysis
description: Use when running, adapting, scripting, or validating single-cell RNA-seq Monocle2 pseudotime or trajectory analysis in R from Seurat objects, especially when selecting target cell states, filtering samples/timepoints/conditions, preserving raw counts from Seurat v4/v5 assays, choosing ordering genes, setting root states, exporting plots/tables, or rerunning a subset without overwriting prior results.
---

# scRNA R Monocle2 Analysis

## Purpose

Build reusable, reviewable R workflows for Monocle2 trajectory analysis from annotated single-cell objects. Keep biological labels, sample/timepoint filters, output tags, and validation explicit so a rerun is reproducible and does not silently overwrite earlier results.

## Input Contract

Before running or writing code, identify:

- input object path and object type, usually Seurat RDS
- assay and counts source
- target label column, such as `cell_type`, `celltype_annotation`, `myeloid_fine`, or `subtype`
- target labels to include
- optional subset filters: sample, donor, condition, disease, treatment, timepoint
- root-state rule, usually "state enriched for the starting cell label"
- output tag and output directories

Do not infer labels from an older object. Inspect the actual metadata columns and label values in the current object.

## Script Pattern

Prefer a reusable CLI-style R script over notebook-only execution. A robust Monocle2 runner should support:

```text
--input-object=PATH
--label-col=COL
--target-labels=LABEL1,LABEL2
--start-label=LABEL1
--subset-col=COL
--subset-values=VALUE1,VALUE2
--output-tag=TAG
--qval-cutoff=0.01
--max-ordering-genes=1000
--min-ordering-genes=100
--cores=1
--seed=1234
--check-only
```

If the existing project script uses narrower names such as `--timepoints=`, keep them when they match the project, but preserve the generic concept: subset cells by explicit metadata values and write to a tagged output tree.

## Counts Handling

Use raw counts for `newCellDataSet()`.

For Seurat v4:

```r
counts <- Seurat::GetAssayData(object, assay = assay, slot = "counts")
```

For Seurat v5, counts may be split across layers such as `counts.1`, `counts.PWH26D1`, or `counts.sampleA`:

1. list count layers with `SeuratObject::Layers(assay_obj, search = "^counts")`
2. collect only requested cells from each layer
3. verify every requested cell was found
4. align features by union gene set when layers differ
5. `cbind` sparse matrices and restore requested cell order

After subsetting, remove all-zero genes and convert to `dgCMatrix`.

## Monocle2 Workflow

1. Build `phenoData` from selected metadata columns and `featureData` with `gene_short_name`.
2. Create the CDS with `VGAM::negbinomial.size()`.
3. Run `estimateSizeFactors()` and `detectGenes()`.
4. Select expressed genes with a minimum cell threshold, commonly `num_cells_expressed >= 10`.
5. Run `monocle::differentialGeneTest()` using the label column, unless the analysis requires a different model.
6. Sort by `qval`, keep genes passing the cutoff, cap at `max_ordering_genes`, and refuse under-supported runs below `min_ordering_genes`.
7. Run `setOrderingFilter()`, `estimateDispersions()`, `reduceDimension(method = "DDRTree")`, and `orderCells()`.
8. Choose the root state by a documented rule, usually highest fraction or count of the `start-label`.
9. Re-run `orderCells(cds, root_state = ...)`.

Keep warnings from Monocle2 or ggplot2 in the log, but distinguish non-fatal deprecation warnings from real failures.

## Output Contract

Write outputs into a dedicated tagged directory, for example:

```text
data/<analysis_name>_<tag>/cds.rds
tables/<analysis_name>_<tag>/cell_pseudotime.csv
tables/<analysis_name>_<tag>/ordering_gene_differential_test.csv
tables/<analysis_name>_<tag>/ordering_genes.csv
tables/<analysis_name>_<tag>/root_state_table.csv
tables/<analysis_name>_<tag>/pseudotime_by_<label_col>.csv
tables/<analysis_name>_<tag>/pseudotime_by_<subset_col>.csv
figures/<analysis_name>_<tag>/trajectory_by_pseudotime.png/.pdf
figures/<analysis_name>_<tag>/trajectory_by_<label_col>.png/.pdf
figures/<analysis_name>_<tag>/trajectory_by_state.png/.pdf
reports/<analysis_name>_<tag>_summary.md
reports/run_<analysis_name>_<tag>.log
```

Use no tag only for the all-cells or original baseline run. For filtered reruns, always use a clear tag such as `D1_D7`, `healthy_stone`, `treated_only`, or `labelA_labelB`.

## Check-Only Mode

Every reusable script should have a `--check-only` path that:

- loads packages and the input object
- verifies required metadata columns
- verifies requested labels and subset values exist
- prints target cell counts by label and subset
- exits before creating trajectory outputs

Run check-only before the full Monocle2 job, especially after changing labels or filters.

## Validation

Before saying the analysis finished, verify:

```bash
tail -n 80 reports/run_<analysis_name>_<tag>.log
find data/<analysis_name>_<tag> tables/<analysis_name>_<tag> figures/<analysis_name>_<tag> -maxdepth 1 -type f | sort
wc -l tables/<analysis_name>_<tag>/cell_pseudotime.csv
wc -l tables/<analysis_name>_<tag>/ordering_genes.csv
du -h data/<analysis_name>_<tag>/*.rds figures/<analysis_name>_<tag>/*.png figures/<analysis_name>_<tag>/*.pdf
```

Expected checks:

- cell pseudotime rows equal selected cells plus one header
- ordering genes rows equal selected ordering genes plus one header
- CDS RDS exists and is non-trivial in size
- trajectory PNG/PDF files are non-empty
- log ends with a finished message or report path
- summary records input object, filters, target counts, ordering gene count, and root state

For subset-filtered runs, also verify excluded groups are absent from `cell_pseudotime.csv` and from grouped pseudotime summaries.

## Interpretation Guardrails

- Do not claim lineage direction from Monocle2 alone; report the operational root-state rule.
- Treat pseudotime differences as trajectory summaries that need marker, composition, or perturbation support for biological interpretation.
- If a subgroup has very low cell count, flag it before comparing pseudotime distributions.
- If DDRTree geometry looks distorted, fix layout/readability first; do not over-prune states or draw misleading pseudotime-sorted lines.
- Preserve previous outputs unless the user explicitly asks for cleanup.

## Runner Pattern

A robust Monocle2 runner should parameterize labels, q-value thresholds, ordering gene caps, seed, check-only mode, group filters, output tags, Seurat v5 count layers, root-state selection, table/figure exports, and a Markdown summary. Reuse that structure for new project-specific runners, but do not hard-code one project's labels into a generic workflow.
