---
name: skills-taxonomy
description: Use as the global router for biomedical analysis skills, especially single-cell RNA-seq and spatial omics. Route first by modality, then spatial resolution, workflow stage, platform, and algorithm; avoid project-specific skills unless a repository adapter is explicitly required.
---

# Biomedical Skills Taxonomy

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
3. Choose ownership:
   - unclear or multi-stage request -> `workflow-*`
   - named method -> `algorithm-*`
   - platform-specific data contract -> platform skill
4. When workflow and algorithm both matter, use both:
   - workflow skill owns stage, launch, and handoff contracts;
   - algorithm skill owns method parameters and result semantics.
5. Use a repository-local adapter only for concrete paths, cohort labels, or project scripts.

Do not route by tissue or project name unless historical reproduction requires an archived adapter.

## Canonical Map

### Single-cell

- `single-cell/workflow/workflow-sc-analysis`
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
  - `spatial/high-resolution/xenium/xenium-segmentation-new`
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

## Clustering Boundary

- Route scRNA-seq and segmented-cell Xenium gene-expression reclustering to
  the counts-first expression contract, including a newly computed UMAP.
- Do not route scimap RCN/CN, target-cell NMF neighborhood states, or
  CellCharter niche clustering through the expression-reclustering UMAP gate;
  they cluster spatial composition or neighborhood features.
- Use scVI only when the user explicitly requests CellCharter expression
  latent mode. Do not select it as a general Xenium refinement default.

## Maintenance

Run the bundled audit after adding, renaming, merging, or retiring biomedical skills:

```bash
python <skills-taxonomy-dir>/scripts/audit_bio_skills.py --project-root /path/to/project
```

Fix structural errors and stale active references before treating the catalog as canonical.
