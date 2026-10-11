---
name: sc-pyscenic-regulon-analysis
description: Run or review the pySCENIC GRN-to-motif-pruning-to-AUCell pipeline, including species database contracts, activity/RSS interpretation and exact cell alignment to AnnData. Use for regulon inference or validated existing outputs.
---

# Skill: sc pySCENIC Regulon Analysis

## Purpose

Use for complete regulatory-network inference or existing-output interpretation; record which stages were supplied versus executed. Generic gene-set scoring belongs to the sample comparison route.

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

0. **Execute or validate the upstream pipeline**
   - Specify cells x genes expression, normalization, stable feature-to-symbol mapping and species. No centered/scaled expression, no silent HVG-only restriction; record any gene filtering and TF coverage. Loom orientation/attributes must match CLI expectations (`Gene` and `CellID`), checked using `pyscenic --help` in the exact environment.
   - GRN: run `pyscenic grn expression.loom tf_list.txt --output adjacencies.tsv --method grnboost2 --seed <seed> --num_workers <n>` (or justified GENIE3). Save TF-target importance table and runtime/seed.
   - Motif pruning: run `pyscenic ctx adjacencies.tsv <ranking.feather> --annotations_fname <motifs.tbl> --expression_mtx_fname expression.loom --output regulons.csv --num_workers <n>` with compatible species/genome, database format, gene IDs and motif annotation version/hashes. Record database coverage/empty regulons and pruning thresholds; raw coexpression is not a motif-supported regulon.
   - AUCell: run `pyscenic aucell expression.loom regulons.csv --output scenic.loom --seed <seed> --num_workers <n>`, saving activity, regulon targets/weights and coverage. Check installed CLI flags before launching (0.12.x can differ); run a small real pipeline before a costly job. Do not claim upstream stages were validated merely from an activity loom.
   - Each stage has input/output hashes, command/version, parameters and pass/fail evidence. Supplied upstream results need compatible provenance; otherwise interpretation is exploratory with missing stages stated.

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
from sc_contracts import align_frame  # shared workflow's resolved scripts path
aligned = align_frame(aucell, adata.obs_names, exact=True)
if not aligned.columns.is_unique:
    raise ValueError("Duplicate regulon names")
adata.obsm["X_aucell"] = aligned
```

Default is exact cell set with explicit reordering. Missing/extra cells fail. Partial coverage must be an explicit configured decision: export missing/extra keys and attrition, then preserve all cells with clearly marked missing activity or create a named restricted derivative. Never silently intersect/delete cells. Use capture-aware keys for repeated barcodes, check finite continuous values and save regulon column order.

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
   - Compare condition activity on independent biological units, not cell-level p-values. RSS/thresholding is descriptive and data-dependent; record threshold method/seed and sensitivity. Motif-supported GRNs still do not establish causal TF regulation.

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

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
