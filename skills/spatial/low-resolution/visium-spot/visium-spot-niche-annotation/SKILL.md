---
name: visium-spot-niche-annotation
description: Use when annotating or comparing 10x Visium spot-level spatial niches or CellCharter clusters from marker genes, spatial maps, histology anatomy, sample/group distribution, and deconvolution composition, especially when deciding epithelial, stromal, vascular, immune, glandular, artifact, or mixed/interface labels.
---

# Skill: Visium Spot Niche Annotation

## Purpose

Use this skill for 10x Visium spot-level niche annotation after a clustering workflow such as CellCharter, Squidpy neighborhood clustering, Banksy, Leiden on spatial features, or deconvolution-derived neighborhood clustering.

The core rule is anatomy first, marker support second. A spot-level niche is not a pure cell type. Do not name a niche from top markers alone when the spatial map shows a different anatomical compartment.

Read [the joint evidence review](../../../../single-cell/python/sc-marker-annotation/references/annotation-evidence-review.md) for full marker statistics, independent canonical/competing program checks, biological-unit support and spot QC. Apply its spot-mixture branch while retaining the anatomy-first order below.

## Required Inputs

- Clustered spot object or cluster table with a stable niche label column.
- Spatial cluster map on HE, ideally per sample plus an all-sample overview.
- Complete available-gene expression-filtered marker statistics plus a candidate view for each niche, ordinarily one-vs-rest, with explicit unintegrated expression/method/filter provenance. If only restricted historical tables are available, declare the review exploratory and record the missing gene scope.
- Sample, donor, condition, and timepoint distributions.
- Available spot QC and marker-panel sources/availability; missing evidence is reported explicitly.
- Optional but strongly recommended: deconvolution composition, HE overlay, or manual histology review.

## Annotation Order

1. Inspect the spatial map first.
   - Decide whether the cluster sits in surface epithelium, stromal tissue, gland/appendage regions, vessel/perivascular regions, immune-rich regions, fat/subcutis, tumor/lesion regions, or artifact/edge regions.
   - Record the image used for the decision.
2. Inspect complete markers and positive/competing programs.
   - Use top100 and top200 as candidate displays when available; measured canonical program genes are checked independently of HVG/DE/top-N selection.
   - Sort markers by `pvals_adj -> logfoldchanges -> pct_ratio -> score`, not score alone.
   - Treat HVG-limited marker tables as exploratory.
   - Review detection, effects, specificity and program coherence against anatomy and optional composition; one marker or a signature score alone does not establish a niche name.
3. Inspect sample/group/timepoint distribution.
   - Group enrichment can support a biological state, but should not override clear anatomy.
   - Timepoint or condition-specific enrichment should be stated as a hypothesis unless sample support is adequate.
   - Record genuine biological-unit support and applicable count/tissue/edge QC. Multiple sections or libraries from one specimen do not create independent replicates. Single-unit dominance limits generalizability without automatically rejecting a coherent region.
   - ECM-rich regions often need stromal, fibroblast, or remodeling labels rather than generic "matrix" labels.
4. Assign a broad anatomical compartment, then a refined label.
5. Flag mixed/interface labels instead of forcing a binary answer.

## Recommended Broad Compartments

- `epithelium`: surface epithelium, glandular epithelium, activated/stressed epithelium, tumor epithelium when applicable.
- `stroma`: fibroblast-rich, collagenous, smooth-muscle, matrix-rich, or reactive stromal regions.
- `gland_or_appendage`: gland, duct, hair follicle, sweat gland, or appendage-like epithelium where relevant to tissue.
- `vascular_or_perivascular`: endothelial, pericyte, smooth-muscle, or blood-rich regions.
- `immune_or_inflammatory`: immune-enriched or inflammatory interface regions.
- `mixed_or_interface`: epithelial-stromal, tumor-stroma, vascular-stromal, immune-stromal, or other boundary regions.
- `artifact_or_low_confidence`: tissue edge, folds, blood-only, low-quality, or inconsistent clusters.

## Marker Interpretation Guardrails

- The example genes below are tissue-dependent examples, not a universal panel; select and source programs for the actual species/tissue and retain measured availability. Unavailable genes are not negative evidence.
- Epithelial markers in a stromal-looking spatial region should become an interface or epithelial-stress hypothesis, not automatically a pure epithelial niche.
- Gland/appendage epithelial markers such as `DCD`, `PIP`, `MUCL1`, `SCGB2A2`, `AQP5`, `KRT7`, and `KRT8/18` should usually be gland/appendage-like when anatomy supports it.
- `POSTN`, `COL1A1`, `COL3A1`, `LUM`, `CTHRC1`, `FN1`, `MMP2`, and `COL6A3` support fibroblast/ECM dermis.
- `MYL9`, `TAGLN`, `ACTA2`, `MYH11`, `RGS5`, `PDGFRB`, `MCAM`, and `DES` support contractile or perivascular dermis.
- `PECAM1`, `VWF`, `PLVAP`, `RAMP2`, `CDH5`, and `CLDN5` support vascular/endothelial regions.
- `HBA1`, `HBA2`, and `HBB` support a blood/erythroid component; interpret it with vessels, tissue coverage and QC. Blood-rich anatomy or coherent mixed programs do not automatically imply contamination.

## Output Contract

Create or update a review table with at least:

- `niche`
- `broad_compartment`
- `subcompartment`
- `refined_label`
- `confidence`
- `spatial_observation`
- `marker_support`
- `program_sources_and_availability`
- `competing_programs`
- `biological_unit_support`
- `qc_evidence`
- `alternatives_and_status`
- `caution`

When comparing K values, also generate a comparison table with:

- `k`
- `niche`
- `n_spots`
- `same_niche_neighbor_fraction`
- `dominant_signature`
- `dominant_signature_score`
- `second_signature`
- `signature_margin`
- `signature_entropy`
- `dominant_timepoint`
- `dominant_timepoint_fraction`
- `mixing_flag`

Higher same-niche neighbor fraction, higher signature margin and lower signature entropy can describe more homogeneous niches. Interpret them alongside anatomy and unit support; a supported interface/mixed niche is not automatically worse, and these metrics are not calibrated identity confidence.

See `references/visium-spot-niche-workflow.md` for a compact reusable workflow.
