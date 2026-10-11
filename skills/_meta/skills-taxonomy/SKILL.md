---
name: skills-taxonomy
description: Route biomedical work to the smallest relevant skill. Use for scRNA-seq, spatial omics, Xenium, Visium/Visium HD, CODEX/IMC, deconvolution, annotation, clustering, or cross-modality analysis.
---

# Biomedical Skills Taxonomy

## On-Demand Leaf Loading

Treat every canonical path below as relative to the skills root, resolved as
`<skills-taxonomy-dir>/../..`.

After choosing the smallest route:

1. Read each selected leaf `SKILL.md` completely before task actions, even if
   that leaf is disabled from initial skill discovery.
2. Follow the selected leaf's resource-reading instructions exactly.
3. Load only the leaf skills required for the request.
4. If the user names a leaf skill, route to that leaf directly.
5. If a listed file is missing, stop and report the stale taxonomy entry.

Disabled discovery is a context-budget optimization, not retirement. The leaf
file remains authoritative for its method or platform.

## Routing Order

0. For cross-modality runtime, environment, CUDA, or GPU-training readiness, pair the selected modality route with `omics-analysis-runtime-preflight`.
1. Choose modality:
   - single-cell RNA-seq -> `single-cell/`
   - spatial transcriptomics or spatial proteomics -> `spatial/`
2. For spatial work, choose resolution:
   - spot-level -> `spatial/low-resolution/visium-spot/`
   - bin/cell-level -> `spatial/high-resolution/<platform>/`
   - cross-platform workflow or algorithm -> `spatial/shared/`
   - segmented-cell Xenium expression reclustering -> pair
     `spatial/shared/workflow-st-analysis` with
     `single-cell/python/sc-subcluster-analysis`; apply its counts-first
     HVG/PCA/integration/neighbors/UMAP/Leiden contract to the retained cells.
   - scRNA or segmented-cell Xenium marker/identity review ->
     `single-cell/python/sc-marker-annotation`, retaining platform QC and
     targeted-panel availability rules; existing-cluster review does not
     require new clustering.
   - Visium spot niche annotation -> `visium-spot-niche-annotation`; HD
     bin/domain interpretation -> `visiumhd-spatial-domain-neighborhood`.
     Preserve mixture/anatomy semantics unless inputs are documented
     segmented-cell profiles. Label-transfer audits use the matching
     transfer skill and measured-expression evidence.
   - requested subclustering–marker annotation–cleanup loop (scRNA/snRNA or
     segmented cells) -> `single-cell/python/sc-subcluster-analysis` owns
     rounds, program-supported automatic review (manual confirmation when requested) and stopping; pair with
     `single-cell/python/sc-marker-annotation` for evidence. Mixed whole-cluster
     proposals require sourced complete-program evidence, target-cell loss assessment and a frozen exact-list decision. Spot/bin
     mixtures do not enter this exclusion loop. Single-run, marker-only and
     layout-only modes remain distinct.
3. Distinguish existing segmentation inspection from requested new segmentation/transcript assignment; route the latter to `xenium-resegmentation`. Phenotype-composition RCN and target NMF share `xenium-scimap-neighborhood-analysis`; conditional spatial label nulls route to `spatial-null-model`. Communication panel/complex audit and transport modeling route to `algorithm-commot`. Geometry/composition routes require only their actual inputs, not count matrices.
4. Choose ownership:
   - unclear or multi-stage request -> `workflow-*`
   - named method -> `algorithm-*`
   - platform-specific data contract -> platform skill
5. When workflow and algorithm both matter, use both:
   - workflow skill owns stage, launch, and handoff contracts;
   - algorithm skill owns method parameters and result semantics.
6. Use a repository-local adapter only for concrete paths, cohort labels, or project scripts.
7. For formal conclusions or submission/publication deliverables, load
   `publication-grade-omics-analysis/SKILL.md` on demand. Ordinary exploration
   does not trigger its full reviewer-readiness gate.

Do not route by tissue or project name unless historical reproduction requires an archived adapter.

## Canonical Map

### Single-cell

- `render-scrna-dotplot`
- `single-cell/workflow/workflow-sc-analysis`
- `single-cell/python/sc-preprocessing-qc`
- `single-cell/python/sc-integration-clustering`
- `single-cell/python/sc-marker-annotation`
- `single-cell/python/sc-sample-level-comparison`
- `single-cell/python/sc-cell-communication`
- `single-cell/python/sc-subcluster-analysis`
- `single-cell/python/sc-multiomics-tcr-vdj-analysis`
- `single-cell/python/sc-pyscenic-regulon-analysis`
- `single-cell/python/sc-velocity-paga-cellrank-analysis`
- `single-cell/r/sc-r-monocle2-analysis`
- `single-cell/shared/h5ad-to-r-bridge`

### Spatial shared

- `spatial/shared/workflow-st-analysis`
- `spatial/shared/workflow-deconvolution-run`
- `spatial/shared/algorithm-cellcharter`
- `spatial/shared/algorithm-commot`
- `spatial/shared/spatial-null-model`
- `spatial/shared/algorithm-tangram`
- `spatial/shared/squidpy-spatial-image-features`

### Visium spot

- `spatial/low-resolution/visium-spot/st-data-preparation`
- `spatial/low-resolution/visium-spot/algorithm-cell2location`
- `spatial/low-resolution/visium-spot/algorithm-rctd`
- `spatial/low-resolution/visium-spot/cell2location-result-visualization`
- `spatial/low-resolution/visium-spot/stlearn-spatial-trajectory`
- `spatial/low-resolution/visium-spot/visium-spot-squidpy-proximity`
- `spatial/low-resolution/visium-spot/visium-spot-niche-annotation`

### High-resolution spatial

- Visium HD:
  - `spatial/high-resolution/visium-hd/visiumhd-qc-analysis`
  - `spatial/high-resolution/visium-hd/visiumhd-spatial-domain-neighborhood`
  - `spatial/high-resolution/visium-hd/visiumhd-cellcharter-visualization`
- Xenium:
  - `render-xenium-dotplot`
  - `spatial/high-resolution/xenium/xenium-segmentation-new`
  - `spatial/high-resolution/xenium/xenium-resegmentation`
  - `spatial/high-resolution/xenium/xenium-roi-h5ad-qc-segmentation`
  - `spatial/high-resolution/xenium/xenium-roi-dapi-overlays`
  - `spatial/high-resolution/xenium/xenium-scrna-label-transfer`
  - `spatial/high-resolution/xenium/xenium-scimap-neighborhood-analysis`
  - `spatial/high-resolution/xenium/xenium-rcn-anchor-nearest-cell-analysis`
- CODEX/IMC:
  - `spatial/high-resolution/codex-imc/codex-imc-neighborhood-analysis`

## Consolidation Rules

- Keep one stage router per modality.
- Keep method parameters in the matching algorithm skill.
- Merge branches that share the same inputs and intermediate representation.
- Keep visualization separate only when it has a distinct input/output contract or bundled implementation.
- Move retired aliases outside active discovery and record their replacements in `references/migration-map.md`.

Biological annotation leaves share [the joint marker evidence review](../../single-cell/python/sc-marker-annotation/references/annotation-evidence-review.md): complete differential ranking plus canonical/competing programs, biological-unit support and appropriate QC. Read it for annotation decisions, not automatically for every plot or marker export. Spatial domains/niches remain separate from cell identities; method owners retain their statistical defaults.

## Clustering Boundary

- Conventional scRNA/snRNA stages route to QC, integration/clustering,
  marker/annotation and sample-level comparison leaves respectively. They
  share `workflow-sc-analysis/references/data-contract.md`; load only needed
  leaves. Communication/velocity/SCENIC/TCR/Monocle are question-driven.
- Missing biological replication permits basic exploration, not group
  inference. Sample-level contrasts require design/replication/confounding
  gates, with edgeR pseudobulk and unit-level proportions as defaults.
- Pure UMAP layout requests with unchanged cells/graph route directly to
  `sc-integration-clustering`'s layout-only mode; preserve clusters/statistics.
- Route scRNA-seq and segmented-cell Xenium gene-expression reclustering to
  the counts-first expression contract, including a newly computed UMAP.
- Do not route scimap RCN/CN, target-cell NMF neighborhood states, or
  CellCharter niche clustering through the expression-reclustering UMAP gate;
  they cluster spatial composition or neighborhood features.
- General Xenium expression clustering may use scVI only on explicit user
  request, after exact counts/data/GPU preflight and a real scVI CUDA pilot
  pass, with CPU fallback disabled. Record batch key, model/training
  parameters, trainer device, peak CUDA allocation, training history, saved
  model, and complete expression-recluster provenance. CellCharter latent
  mode remains a compatible explicit-request alias, not the only allowed use.

## Maintenance

Run the bundled audit after adding, renaming, merging, or retiring biomedical skills:

```bash
python <skills-taxonomy-dir>/scripts/audit_bio_skills.py --project-root /path/to/project
```

Fix structural errors and stale active references before treating the catalog as canonical.

The audit validates structure and applicable project contracts, not analysis-method correctness. A project without `metadata/current_annotation_registry.json` receives `ANNOTATION_REGISTRY=NOT_APPLICABLE` and continues generic checks. An existing registry is validated and malformed/missing required assets fail; it is never silently skipped.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
