# Joint marker and annotation evidence review

Read when assigning or auditing biological identities, subtypes, states, or spatial compartments. Combine differential marker ranking, canonical and competing programs, biological-unit support, and assay-appropriate QC. A request only to export markers, plot supplied labels, or read a sidecar does not authorize a new annotation, clustering, or cell exclusion.

## Differential ranking and the candidate view

- The FindAllMarkers-style comparison is each declared cluster versus the remaining retained observations. Record the cluster key, cell scope, comparator, matrix, method, feature filter, group sizes, and candidate filters; a restricted comparator must be explicit.
- Use all available genes passing the declared expression filter, independently of HVGs. Preserve prior feature restrictions in the scope statement. An HVG-only object cannot provide all-gene annotation evidence; retrieve an aligned all-gene source when available, otherwise declare the limitation.
- Use explicit unintegrated, unscaled normalized log expression for ordinary marker ranking. Detection fractions come from `>0` in the declared nonnegative unscaled expression or separately declared counts; record whether correction changed zero patterns. PCA, batch-corrected embeddings, BANKSY features, and neighborhood-composition matrices cannot supply gene-expression effects or detection.
- Preserve complete tested-gene statistics and a separate filtered candidate view. Include effects, detection in/out, test statistics, raw/adjusted p-values and group sizes; record inadequate groups and undefined statistics rather than treating them as negative evidence. Candidate cutoffs and top-N size are configurable displays, not the annotation gene universe.
- Keep the existing Python t-test default. Use a requested or justified alternative explicitly; do not routinely run both methods, both backends, or HVG-only plus all-gene marker analyses. For a user-selected native Seurat workflow, declare the unintegrated assay/data layer, `test.use`, and version-specific prefilters. A default prefiltered FindAllMarkers result is a candidate view, not an exhaustive statistics table. Preserve every tested result before candidate filtering when producing a complete table. Scanpy and Seurat methods, fold changes, and adjusted p-values are not assumed numerically interchangeable.
- Cell-level cluster p-values are exploratory identity evidence. Between-condition inference requires genuine biological units and the sample-level comparison workflow.

## Canonical and competing programs

Choose programs appropriate to species, tissue, disease context, assay, and requested hierarchy. Record their source/version and marker-to-feature mapping. Include expected positive programs and plausible competing lineages or subtypes; distinguish identity programs from stress, cell cycle, interferon, and activation states. Do not turn one tissue's example panel into a universal list.

Audit panel availability before DE filtering. Keep availability, detection, and test eligibility separate: an unmeasured panel gene, a gene absent from a restricted object, an ambiguous feature mapping, a measured low-detection gene, and a gene skipped by a test have different meanings. Preserve measured low-detection panel summaries outside the candidate view. Do not silently resolve duplicate symbols or interpret an unavailable gene as a zero.

Review measured panel genes whether or not they are HVGs, significant, or in the top-N table. Use detection, expression/effect, specificity against plausible alternatives, and coherence across the program. A shared lineage marker can support identity without being differentially expressed between two related clusters. Conversely, a highly ranked state gene or one foreign-lineage gene cannot establish a new identity. An absent program can weaken an interpretation only when the assay measures informative genes with adequate detection power.

When two programs conflict, inspect cell-level detection if needed: separate subsets within a cluster differ from coordinated co-detection in the same cells. Cluster averages or a module score alone cannot distinguish them. Sparse detection limits both conclusions; avoid requiring every cell to express every marker. Module scores are optional supporting summaries, with their feature/control choices recorded.

## Biological-unit support and QC

For each reviewed group, summarize observation counts and fractions by sample and genuine biological unit, and the support for positive/competing programs within those units. Use per-unit detection/expression summaries where informative; a second DE run per unit is not mandatory. Record recurrence, sparse-unit limitations, and single-unit dominance. Multiple captures or ROIs from one specimen do not create independent replicates. Missing unit metadata is unavailable evidence, not invented replication.

A rare, tissue-specific, or condition-restricted group need not appear in every sample. Single-unit support limits generalizability but does not automatically invalidate a coherent identity. Condition enrichment cannot substitute for identity evidence or establish a condition effect from cell-level p-values.

Use existing platform-appropriate QC and the relevant QC skill:

| Observation model | QC context and interpretation |
|---|---|
| scRNA/snRNA cells | Count and detected-gene distributions, applicable mitochondrial/nuclear metrics, doublet evidence, and any independently supported ambient-RNA signal. State unavailable metrics and correction provenance. |
| Segmented Xenium cells | Counts/detected genes, cell and nucleus area, nucleus count, segmentation method, and geometry when relevant and available. Boundary stain is accepted; foreign expression is not automatically ambient contamination. Distinguish segmentation doublet/boundary error, broad-compartment mismatch, immediately adjacent transcripts, and biological mixed state. |
| Visium spots or HD bins | Counts/detected genes, tissue coverage, image alignment, edges/folds and composition where available. Coherent programs from multiple lineages can reflect a real mixed region. A small bin is not automatically a segmented single cell. |

QC supports an interpretation; one metric does not decide identity or justify deletion. Avoid importing universal thresholds across assays. Two same-lineage foreign markers are screening evidence only. Removal requires a sourced, lineage-specific foreign program with multiple complementary markers, measured coverage and audited feature mapping; there is no universal gene-count or score cutoff. Review target programs, cell-level distributions and coordinated coexpression, biological-unit support/limitations and platform QC. Neither averages, module scores, UMAP nor seed instability alone justify deletion. Do not require every cell to express every program gene or treat panel/detection gaps as negative. Dual program positivity alone is not contamination. If genuine target subtype, biological mixed state or potential target-cell loss remains unresolved, retain as `review`. Whole-cluster removal needs supported technical mixture or an explicit outside-target-scope interpretation and quantified target-positive count/fraction and potential loss; partial-cell evidence cannot be expanded to a whole cluster. Retain or use the applicable at-most-one mixed-parent split and review a supported foreign child. An explicitly requested cleanup loop may use automatic program-supported review, including supported mixed whole clusters, under its frozen exact-list contract; manual confirmation remains available when requested. Review evidence alone does not authorize removal or another exclusion cycle.

## Spatial and transfer evidence

For spot/bin compartments, retain anatomy-first review: compare markers and coherent programs with H&E location, sample support, QC, and optional deconvolution composition. Describe mixed/interface regions when supported. BANKSY domains, RCN and CellCharter niches describe spatial/compositional states and stay separate from cell-expression identities. Documented segmented-cell profiles may use cell-identity interpretation.

For label transfer, retain reference/version, feature overlap, predicted label, confidence and competing prediction separately. Audit measured target expression, canonical/competing programs, sample support and QC; transferred or imputed expression is not an independent target measurement. Ranking genes by predicted groups can be descriptive, but is not independent validation of those predictions. Use existing appropriate groups where available; do not create clusters merely to complete a transfer audit. A confident prediction conflicting with measured programs remains flagged, not promoted to a curated identity.

## Evidence record and decision

Extend existing review outputs rather than requiring a new file format. Keep enough information to recover:

- reviewed group/key and scope; candidate identity, subtype and state;
- complete marker table and candidate-view provenance;
- positive and competing programs, sources, informative genes, availability, detection and specificity;
- sample/biological-unit counts, recurrence and limitations;
- applicable QC, spatial/transfer evidence and missing checks;
- competing explanations, confidence rationale, status, decision and annotation version.

Assign the narrowest identity supported by coherent programs, specificity, unit support and QC together. Preserve unknown, ambiguous, candidate and unresolved outcomes when appropriate; a strong broad identity with weak subtype evidence remains broad. Neither differential significance nor complete label coverage is an annotation acceptance criterion. Keep stress/activation state separate from lineage and phenotype labels.

Before using project annotations, resolve the requested semantic role under the project's guide/registry contract and preserve its join, scope and grade. Native identities, transfer predictions, EPI subtypes and RCN fields cannot replace each other. Annotation writes require the requested scope and applicable validation; review does not itself update a registry or formalize a candidate. For review-only work, evidence tables suffice without creating an annotated H5AD.

Before accepting a review, check that each requested group has marker, panel, unit and QC evidence or an explicit limitation, and that the recorded decision addresses material conflicts. Preserve source objects, keys and existing evidence grades. Generate requested figures as PNG under the repository output policy.
