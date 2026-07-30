---
name: sc-pyscenic-regulon-analysis
description: Use when running or reviewing pySCENIC/AUCell regulon analysis in Python, including loom input preparation, regulon activity matrices, AUCell thresholding, regulon binarization, RSS ranking, regulon annotation, and joining regulon activity back to Scanpy AnnData.
---

# Skill: sc pySCENIC Regulon Analysis

## Purpose

Use this skill for Python regulatory-network workflows based on pySCENIC outputs. It is intended for regulon activity interpretation and integration with Scanpy annotations, not for generic gene set scoring.

## Environment

Use a dedicated environment because pySCENIC and loom dependencies are version-sensitive:

```bash
conda create -n pyscenic_analysis python=3.10
conda activate pyscenic_analysis
pip install scanpy pyscenic loompy arboreto
```

If reproducing an older pySCENIC run, pin package versions from the project log.

## Inputs

Common inputs:

- expression matrix or `.loom` used for pySCENIC.
- pySCENIC output `.loom` containing regulons and AUCell activity.
- Scanpy `.h5ad` with cell annotations and embeddings.
- transcription factor database and motif ranking files used upstream.

Check:

- cell barcodes overlap between loom and AnnData.
- regulon names are unique and preserve activation/repression suffixes if present.
- species-specific TF/motif databases match the organism.

## Workflow

1. **Read pySCENIC loom**

```python
import scanpy as sc
import loompy
from pyscenic.export import add_scenic_metadata

adata = sc.read_h5ad("annotated.h5ad")
```

2. **Extract AUCell / regulon activity**
   - Prefer structured pySCENIC helpers when available.
   - Keep both continuous AUCell scores and binarized regulon calls if present.
   - Store matrices under explicit keys, e.g. `obsm['X_aucell']`.

3. **Align cells**

```python
shared = adata.obs_names.intersection(aucell.index)
if len(shared) == 0:
    raise ValueError("No shared cells between AnnData and pySCENIC AUCell matrix.")
adata = adata[shared].copy()
adata.obsm["X_aucell"] = aucell.loc[shared]
```

4. **Rank regulons by cell type**
   - Use regulon specificity score (RSS) or group-level AUCell summaries.
   - Rank within a trusted cell-type or cluster column.
   - Report both top regulons and their target gene counts.

5. **Visualize**
   - UMAP/spatial embedding overlays for selected regulons.
   - heatmaps of regulon activity by cell type.
   - RSS dot plots or ranked bars.

6. **Annotate interpretation**
   - Separate TF activity from TF expression.
   - Use regulon target genes to support annotation, but validate against known markers.

## Output Checks

- AUCell matrix shape matches AnnData cells after alignment.
- regulon metadata includes target counts and database/species provenance.
- binarization thresholds are recorded if binary regulon calls are used.
- top regulons are summarized by annotation group.
- selected regulon activity plots were generated on the intended embedding.

## Boundaries

- pySCENIC regulon activity is inferred regulatory activity, not direct chromatin accessibility.
- For true scATAC+scRNA integration, use a dedicated ATAC/multiome workflow; do not force pySCENIC to stand in for peak-gene linkage.
- Course notebooks with hard-coded loom paths should be treated as pattern sources, not canonical runnable scripts.
