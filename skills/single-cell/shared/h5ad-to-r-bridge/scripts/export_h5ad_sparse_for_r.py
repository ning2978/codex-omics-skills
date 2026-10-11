#!/usr/bin/env python
"""Export AnnData h5ad to a sparse Matrix Market bridge for R/Seurat/RCTD."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export an AnnData h5ad to sparse R-readable files without densifying expression matrices."
    )
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=["sc", "st"], required=True)
    parser.add_argument("--counts-layer", required=True, help="Raw count layer to export, e.g. counts or raw_counts.")
    parser.add_argument("--label-key", default=None, help="Optional cell-type label key expected in obs.")
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--export-x", action="store_true", help="Also export adata.X as x.mtx.gz.")
    parser.add_argument("--export-reductions", action="store_true", help="Export obsm reductions as CSV files.")
    parser.add_argument("--allow-dense-input", action="store_true", help="Allow dense matrices by converting them to sparse.")
    parser.add_argument("--no-round-integer-like", action="store_true", help="Do not round float counts that are integer-like.")
    parser.add_argument("--no-strict-spatial-check", action="store_true", help="Do not fail if obsm spatial disagrees with obs pixel columns.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def as_sparse(matrix, name: str, allow_dense: bool) -> sparse.spmatrix:
    if sparse.issparse(matrix):
        return matrix.tocsc()
    if not allow_dense:
        raise TypeError(f"{name} is dense. Refusing to export without --allow-dense-input.")
    return sparse.csc_matrix(matrix)


def matrix_stats(matrix: sparse.spmatrix) -> dict[str, object]:
    data = matrix.data
    if data.size == 0:
        return {
            "nnz": 0,
            "min_nonzero": None,
            "max_nonzero": None,
            "nonnegative": True,
            "all_finite": True,
            "integer_like": True,
            "max_integer_abs_diff": 0.0,
        }
    rounded = np.rint(data)
    max_diff = float(np.max(np.abs(data - rounded)))
    return {
        "nnz": int(matrix.nnz),
        "min_nonzero": float(np.min(data)),
        "max_nonzero": float(np.max(data)),
        "nonnegative": bool(np.min(data) >= 0),
        "all_finite": bool(np.isfinite(data).all()),
        "integer_like": bool(max_diff <= 1e-6),
        "max_integer_abs_diff": max_diff,
    }


def maybe_round_counts(matrix: sparse.spmatrix, stats: dict[str, object], no_round: bool) -> sparse.spmatrix:
    if no_round or not stats["integer_like"] or np.issubdtype(matrix.dtype, np.integer):
        return matrix
    rounded = matrix.copy()
    rounded.data = np.rint(rounded.data).astype(np.int64)
    return rounded


def write_mtx_gz(matrix: sparse.spmatrix, path: Path) -> None:
    tmp = path.with_suffix("")
    io.mmwrite(tmp, matrix)
    with open(tmp, "rb") as src, gzip.open(path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    tmp.unlink()


def write_lines_gz(values: list[str], path: Path) -> None:
    with gzip.open(path, "wt") as handle:
        for value in values:
            handle.write(f"{value}\n")


def write_features(adata: ad.AnnData, path: Path) -> None:
    gene_names = adata.var["gene_symbol"].astype(str).tolist() if "gene_symbol" in adata.var else adata.var_names.astype(str).tolist()
    gene_ids = None
    for key in ["gene_ids", "gene_id", "feature_id", "ensembl_id"]:
        if key in adata.var.columns:
            gene_ids = adata.var[key].astype(str).tolist()
            break
    if gene_ids is None:
        gene_ids = adata.var_names.astype(str).tolist()
    if len(set(gene_ids)) != len(gene_ids) or any(not g for g in gene_ids):
        raise ValueError("Export needs unique nonempty feature IDs; duplicate symbols cannot be used as keys.")
    feature_types = adata.var["feature_types"].astype(str).tolist() if "feature_types" in adata.var.columns else ["Gene Expression"] * len(gene_names)
    with gzip.open(path, "wt") as handle:
        for gene_id, gene_name, feature_type in zip(gene_ids, gene_names, feature_types, strict=True):
            handle.write(f"{gene_id}\t{gene_name}\t{feature_type}\n")


def metadata_frame(adata: ad.AnnData) -> pd.DataFrame:
    meta = adata.obs.copy()
    if "barcode" in meta.columns:
        existing = meta["barcode"].astype(str).to_numpy()
        obs_names = adata.obs_names.astype(str).to_numpy()
        if np.array_equal(existing, obs_names):
            meta = meta.drop(columns=["barcode"])
        else:
            if "source_barcode" in meta.columns:
                raise ValueError("source_barcode already exists; preserve the original barcode map explicitly")
            meta = meta.rename(columns={"barcode": "source_barcode"})
    meta.insert(0, "barcode", adata.obs_names.astype(str))
    for col in meta.columns:
        if isinstance(meta[col].dtype, pd.CategoricalDtype):
            meta[col] = meta[col].astype(str)
    return meta


def export_spatial(adata: ad.AnnData, output_dir: Path, strict: bool) -> dict[str, object]:
    if "spatial" not in adata.obsm:
        raise KeyError("ST mode requires adata.obsm['spatial'].")
    spatial = np.asarray(adata.obsm["spatial"])
    if spatial.shape != (adata.n_obs, 2):
        raise ValueError(f"Expected obsm['spatial'] shape {(adata.n_obs, 2)}, got {spatial.shape}.")

    pxl_col = spatial[:, 0].astype(float)
    pxl_row = spatial[:, 1].astype(float)
    validation = {
        "obsm_spatial_semantics": "x_y_equals_pxl_col_pxl_row",
        "obs_pixel_columns_checked": False,
        "obs_pixel_columns_match": None,
    }
    if {"pxl_col_in_fullres", "pxl_row_in_fullres"}.issubset(adata.obs.columns):
        obs_col = adata.obs["pxl_col_in_fullres"].astype(float).to_numpy()
        obs_row = adata.obs["pxl_row_in_fullres"].astype(float).to_numpy()
        col_match = bool(np.allclose(pxl_col, obs_col))
        row_match = bool(np.allclose(pxl_row, obs_row))
        validation.update(
            {
                "obs_pixel_columns_checked": True,
                "obs_pixel_columns_match": bool(col_match and row_match),
                "pxl_col_match": col_match,
                "pxl_row_match": row_match,
            }
        )
        if strict and not (col_match and row_match):
            raise ValueError("obsm['spatial'] does not match obs pixel columns.")

    spatial_dir = output_dir / "spatial"
    spatial_dir.mkdir(parents=True, exist_ok=True)
    positions = pd.DataFrame(
        {
            "barcode": adata.obs_names.astype(str),
            "in_tissue": adata.obs["in_tissue"].to_numpy() if "in_tissue" in adata.obs.columns else 1,
            "array_row": adata.obs["array_row"].to_numpy() if "array_row" in adata.obs.columns else np.nan,
            "array_col": adata.obs["array_col"].to_numpy() if "array_col" in adata.obs.columns else np.nan,
            "pxl_row_in_fullres": pxl_row,
            "pxl_col_in_fullres": pxl_col,
            "x": pxl_col,
            "y": pxl_row,
        }
    )
    positions.to_csv(spatial_dir / "tissue_positions.csv", index=False)
    return validation


def export_reductions(adata: ad.AnnData, output_dir: Path) -> list[str]:
    reductions_dir = output_dir / "reductions"
    written = []
    for key in sorted(adata.obsm.keys()):
        if key == "spatial":
            continue
        values = adata.obsm[key]
        values = values.toarray() if sparse.issparse(values) else np.asarray(values)
        if values.ndim != 2 or values.shape[0] != adata.n_obs:
            raise ValueError(f"Cannot preserve obsm[{key!r}] as a reduction; unsupported shape")
        if values.dtype.kind not in "fiu" or not np.isfinite(values).all():
            raise ValueError(f"Reduction {key!r} must be a finite numeric matrix")
        if not key or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in key):
            raise ValueError(f"Unsafe reduction filename key: {key!r}")
        reductions_dir.mkdir(parents=True, exist_ok=True)
        cols = [f"{key}_{i + 1}" for i in range(values.shape[1])]
        df = pd.DataFrame(values, index=adata.obs_names.astype(str), columns=cols)
        df.index.name = "barcode"
        path = reductions_dir / f"{key}.csv"
        df.to_csv(path)
        written.append(str(path.relative_to(output_dir)))
    return written


def main() -> None:
    args = parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{args.output_dir} already exists and is not empty. Pass --overwrite.")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Unique cell and feature IDs required; do not silently make duplicates unique.")
    if args.counts_layer not in adata.layers:
        raise KeyError(f"Missing layer {args.counts_layer!r}; available layers: {list(adata.layers.keys())}")
    if args.label_key and args.label_key not in adata.obs.columns:
        raise KeyError(f"Missing label key {args.label_key!r} in obs.")

    counts = as_sparse(adata.layers[args.counts_layer], f"layers[{args.counts_layer!r}]", args.allow_dense_input)
    stats_before = matrix_stats(counts)
    if not stats_before.get("all_finite", True):
        raise ValueError("Counts contain nonfinite values")
    if not stats_before["nonnegative"]:
        raise ValueError("Counts contain negative values.")
    if not stats_before["integer_like"]:
        raise ValueError(f"Counts are not integer-like; max abs diff {stats_before['max_integer_abs_diff']}.")
    counts = maybe_round_counts(counts, stats_before, args.no_round_integer_like)
    counts_gene_by_obs = counts.T.tocoo()

    write_mtx_gz(counts_gene_by_obs, args.output_dir / "counts.mtx.gz")
    write_lines_gz(adata.obs_names.astype(str).tolist(), args.output_dir / "barcodes.tsv.gz")
    write_features(adata, args.output_dir / "features.tsv.gz")
    feature_meta = adata.var.copy()
    if "source_var_name" in feature_meta:
        raise ValueError("source_var_name is reserved for the bridge feature mapping")
    feature_meta.insert(0, "source_var_name", adata.var_names.astype(str))
    feature_meta.to_csv(args.output_dir / "feature_metadata.csv", index=False)
    metadata_frame(adata).to_csv(args.output_dir / "metadata.csv", index=False)

    x_exported = False
    if args.export_x:
        x = as_sparse(adata.X, "adata.X", args.allow_dense_input)
        write_mtx_gz(x.T.tocoo(), args.output_dir / "x.mtx.gz")
        x_exported = True

    spatial_validation = None
    if args.mode == "st":
        spatial_validation = export_spatial(adata, args.output_dir, strict=not args.no_strict_spatial_check)

    reductions = export_reductions(adata, args.output_dir) if args.export_reductions else []

    manifest = {
        "schema_version": 2,
        "input_h5ad": str(args.input_h5ad),
        "mode": args.mode,
        "counts_layer": args.counts_layer,
        "shape_obs_by_vars": [int(adata.n_obs), int(adata.n_vars)],
        "matrix_orientation": "genes_by_observations",
        "barcodes_n": int(adata.n_obs),
        "features_n": int(adata.n_vars),
        "counts_dtype_in_h5ad": str(adata.layers[args.counts_layer].dtype),
        "counts_stats": stats_before,
        "counts_rounded_if_integer_like": bool(not args.no_round_integer_like and stats_before["integer_like"]),
        "label_key": args.label_key,
        "label_n": int(adata.obs[args.label_key].nunique()) if args.label_key else None,
        "sample_key": args.sample_key,
        "sample_values": sorted(map(str, adata.obs[args.sample_key].unique())) if args.sample_key in adata.obs.columns else [],
        "x_exported": x_exported,
        "reductions_exported": reductions,
        "spatial_validation": spatial_validation,
        "files": {
            "counts": "counts.mtx.gz",
            "barcodes": "barcodes.tsv.gz",
            "features": "features.tsv.gz",
            "metadata": "metadata.csv",
            "feature_metadata": "feature_metadata.csv",
            "spatial_positions": "spatial/tissue_positions.csv" if args.mode == "st" else None,
        },
    }
    with open(args.output_dir / "source_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    validation = pd.DataFrame(
        [
            {"metric": "n_obs", "value": adata.n_obs},
            {"metric": "n_vars", "value": adata.n_vars},
            {"metric": "counts_nnz", "value": stats_before["nnz"]},
            {"metric": "counts_nonnegative", "value": stats_before["nonnegative"]},
            {"metric": "counts_integer_like", "value": stats_before["integer_like"]},
            {"metric": "label_n", "value": manifest["label_n"] if manifest["label_n"] is not None else ""},
            {"metric": "spatial_obs_pixel_match", "value": spatial_validation.get("obs_pixel_columns_match", "") if spatial_validation else ""},
        ]
    )
    validation.to_csv(args.output_dir / "validation_summary.csv", index=False)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
