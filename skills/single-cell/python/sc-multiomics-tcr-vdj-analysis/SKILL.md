---
name: sc-multiomics-tcr-vdj-analysis
description: Use when integrating single-cell RNA-seq with scTCR, scBCR, bulk TCR, or VDJ repertoire data in Python, including 10x/Singleron/TRUST4 VDJ parsing, clonotype definition, clone expansion, diversity metrics, repertoire overlap, scirpy analysis, TCR distance clustering, and joining clonotypes back to AnnData.
---

# Skill: sc Multiomics TCR/VDJ Analysis

## Purpose

Use this skill for Python workflows that combine expression annotations with immune repertoire data. It covers scRNA+scTCR/scBCR, bulk TCR, TRUST4-derived TCR from scRNA, and VDJ clustering/summary analyses.

## Environment

Use a dedicated environment for repertoire packages:

```bash
conda create -n tcr_vdj_analysis python=3.10
conda activate tcr_vdj_analysis
pip install scanpy scirpy pandas scipy python-igraph tcrdist3 clustcr
```

Some packages have compiled dependencies; pin versions only after checking the project imports.

## Inputs

Typical inputs:

- AnnData `.h5ad` with cell states, sample IDs, and embeddings.
- 10x `filtered_contig_annotations.csv` or Singleron VDJ outputs.
- TRUST4 output tables for TCR reconstructed from scRNA.
- bulk TCR tables with sample, CDR3 nucleotide/amino acid, V/J genes, count, and frequency.

Normalize column names early:

- `sample`
- `barcode` for single-cell data
- `cdr3aa`, `cdr3nt`
- `v_gene`, `j_gene`
- `count`, `freq`
- chain labels such as TRA/TRB or IGH/IGK/IGL

## Workflow

1. **Read and standardize VDJ tables**
   - Keep sample labels explicit.
   - Strip or add barcode suffixes only with an overlap check against AnnData.
   - Remove nonproductive contigs unless the question requires them.

2. **Join to AnnData**

```python
adata.obs = adata.obs.join(vdj_metadata.set_index("barcode"), how="left")
```

Check the fraction of cells with VDJ annotations by sample and cell type.

3. **scirpy repertoire analysis**

```python
import scirpy as ir

ir.pp.merge_with_ir(adata, vdj_adata)
ir.tl.chain_qc(adata)
ir.pp.ir_dist(adata, metric="identity", sequence="aa")
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonal_expansion(adata)
```

4. **Diversity and overlap**
   - Report clonotype count, reads/cell count, mean frequency, CDR3 length, and convergence when available.
   - Use repertoire overlap matrices only after sample library sizes are clear.
   - For bulk TCR, compute diversity per sample before comparing groups.

5. **TCR distance clustering**
   - Use `tcrdist3` or `clustcr` for sequence-similarity clustering.
   - Standardize V/J gene names before distance calculation.
   - Record distance cutoff, chain, organism, and clustering method.

6. **Expression-linked interpretation**
   - Compare clone expansion across cell states.
   - Plot clonotype or clone-size categories on UMAP.
   - Test whether expanded clones are enriched in target states with sample-aware statistics.

## Output Checks

- barcode overlap report exists.
- VDJ coverage is summarized by sample and cell type.
- clonotype IDs are stable and not recomputed with hidden parameters.
- clone expansion metrics are stored in `adata.obs`.
- diversity/overlap outputs include sample-level denominators.
- sequence-distance clustering records organism, chain, V/J normalization, cutoff, and package.

## Boundaries

- Do not infer antigen specificity from clonotype expansion alone.
- Do not merge bulk and single-cell TCR frequencies without normalizing sample/library sizes.
- CNV or ordinary Scanpy integration utilities from course notebooks are supporting context, not part of this skill unless VDJ interpretation depends on them.
