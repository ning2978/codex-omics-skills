---
name: xenium-roi-h5ad-qc-segmentation
description: >-
  QC-filter per-ROI Xenium h5ad exports (n_genes threshold stats, anucleate removal,
  global cell_area percentiles, keep boundary/interior stain segmentation only).
  Use when working with split h5ad per sample/ROI, filtering cells before downstream
  analysis, or when the user mentions n_genes at least 50, morphology
  segmentation filter, or excluding nucleus-expansion fallback cells.
---

# Xenium ROI h5ad QC + Segmentation Filter

## Purpose

Filter **pre-split per-ROI AnnData** (one `{roi_id}.h5ad` per Qupath/Xenium Explorer selection) before downstream scRNA-style or spatial analysis. This is **not** re-segmentation (no Cellpose/Baysor/ProSeg); it filters cells already segmented by 10x multi-stain Xenium.

**Use when:**

- Input is a directory of `{ROI}.h5ad` from a sample split (e.g. `<sample>_output_h5ad/`)
- User wants `n_genes >= N` counts, simple QC, and **morphology-based** segmentation retention (boundary + interior stain only)

**Do NOT use when:**

- Input is raw `output-XETG*/cells.parquet` only → use `xenium-segmentation-new` for polygon/QC on parquet
- User wants to **re-run** segmentation from images or `transcripts.parquet`
- Task is full-slide single `cell_feature_matrix.h5` without ROI splits → use `xenium_qc_two_samples.py` pattern on `split_by_modality/xenium/{sample}/RAW_DATA/`

## Input contract

| Item | Requirement |
|------|-------------|
| Directory | `{sample}_output_h5ad/` containing `*.h5ad` (one file per ROI) |
| `adata.X` | Sparse counts, `n_vars` ~5051 for Xenium 5K panels |
| `adata.obs` | Must include `cell_area`, `nucleus_area`, `segmentation_method` |
| `adata.obsm` | `spatial` (2D µm, cell centroids) |
| ROI ID | File stem, e.g. `5191579.h5ad` → `roi_id=5191579` |

Optional sidecar (not required for filtering): `{sample}_coordinates.csv` (ROI polygon vertices), `{sample}_cells_stats/` (per-cell CSV mirrors).

## Filter pipeline (order matters)

Apply in sequence; log retention after **each** step:

1. **`n_genes_by_counts >= min_genes`** (default 50)
   - Compute via `sc.pp.calculate_qc_metrics(adata, expr_type="counts", log1p=False)` if column missing.
2. **Has nucleus**: `nucleus_area > 0` and finite (drops anucleate).
3. **`cell_area` global percentile band**: keep cells with `cell_area` in **[p_low, p_high]** computed across **all ROIs pooled** (default p1–p99), not per-ROI.
4. **Segmentation method** (morphology): keep only:
   - `Segmented by interior stain (18S)`
   - `Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)`
   Drop: `Segmented by nucleus expansion of 5.0µm`

Presets for CLI: `interior_stain`, `boundary_stain` (comma-separated).

## Deterministic implementation contract

A repository-local adapter may be named `xenium_roi_h5ad_qc_filter.py`.
Do not assume it is bundled with this skill; inspect the adapter before use.

Two-pass: (1) scan all ROIs for multi-threshold gene counts + global `cell_area` percentiles; (2) filter, write per-ROI h5ad, optional merged h5ad with `obs['roi_id']`.

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate scanpy

python /path/to/project/xenium_roi_h5ad_qc_filter.py \
  --input-dir "/path/to/{sample}_output_h5ad" \
  --sample <sample_id> \
  --min-genes 50 \
  --area-low-pct 1 --area-high-pct 99 \
  --keep-seg "interior_stain,boundary_stain" \
  --out-h5ad-dir "results/xenium_roi_h5ad_qc/<sample_id>/filtered_h5ad" \
  --out-report-dir "results/xenium_roi_h5ad_qc/<sample_id>/reports" \
  --merge
```

Use `--no-merge` to skip concatenated h5ad.

## Output contract

| Output | Location (project-relative) |
|--------|----------------------------|
| Per-ROI filtered h5ad | `{out-h5ad-dir}/{roi_id}.h5ad` |
| Merged filtered h5ad | `{out-h5ad-dir}/{sample}_filtered.h5ad` with `obs['roi_id']` |
| Gene threshold table | `{out-report-dir}/n_genes_threshold_summary.csv` (rows per ROI + `TOTAL`; columns `n_genes_ge_30/40/50/100`) |
| Step retention | `{out-report-dir}/filter_retention_per_roi.csv` |
| Segmentation breakdown | `{out-report-dir}/segmentation_method_breakdown.csv` (`before_filter` / `after_filter`) |
| Area thresholds | `{out-report-dir}/global_area_thresholds.json` |
| Run log | `{out-report-dir}/qc_filter_summary.log` |

## Validation checklist

After run, confirm:

- [ ] Per-ROI h5ad count matches input ROI count
- [ ] Merged `n_obs` equals sum of per-ROI `n_obs` (if `--merge`)
- [ ] All retained cells: `n_genes_by_counts >= min_genes`, `nucleus_area > 0`, area within global band
- [ ] **Zero** cells with `segmentation_method == nucleus expansion` in filtered objects
- [ ] Spot-read 1–2 ROI files + merged: `n_vars` unchanged, `obsm['spatial']` present

## Interpreting Retention

Low-retention ROIs often reflect few cells passing the selected `min_genes` threshold in that ROI, not necessarily a script bug. Inspect `filter_retention_per_roi.csv` and the `n_genes_threshold_summary.csv` columns before changing segmentation filters.

Write filtered outputs under the active analysis repository or results root, not back into raw download/storage directories.

## Relation to other skills

| Skill | Role |
|-------|------|
| `xenium-segmentation-new` | Parquet boundaries, polygon QC, vertex plots on **full slide** `cells.parquet` |
| **this skill** | **Per-ROI h5ad** gene QC + stain-based segmentation filter |
| project full-slide QC scripts | Full-slide `cell_feature_matrix.h5` + optional cells.parquet join |

Downstream after filtering: normalization/UMAP/Leiden on `{sample}_filtered.h5ad` with `roi_id` as batch covariate, or per-ROI analyses.

## Agent workflow

```
Task Progress:
- [ ] Confirm input dir has {sample}_output_h5ad/*.h5ad
- [ ] Run xenium_roi_h5ad_qc_filter.py (scanpy env)
- [ ] Read n_genes_threshold_summary.csv TOTAL row for user-facing stats
- [ ] Read filter_retention_per_roi.csv if user asks about specific ROI losses
- [ ] Verify nucleus_expansion count = 0 in filtered data
- [ ] Point user to merged h5ad for downstream
```

## Troubleshooting

| Issue | Action |
|-------|--------|
| `_CSRDataset` comparison error on backed X | Load full `read_h5ad` for pass1 (script already does) |
| Very low ROI retention | Check `n_genes_ge_50` column for that ROI in threshold summary |
| User wants stricter genes | Increase `--min-genes` (e.g. 100) |
| User wants boundary-only | `--keep-seg boundary_stain` |
| User wants to include nucleus expansion | `--keep-seg interior_stain,boundary_stain,nucleus_expansion` (not default) |
