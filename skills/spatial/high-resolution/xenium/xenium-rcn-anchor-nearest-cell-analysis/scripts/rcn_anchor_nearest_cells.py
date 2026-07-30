#!/usr/bin/env python3
"""Run same-ROI nearest-cell analysis for a fixed phenotype in a frozen RCN.

The script is annotation-agnostic: callers must first resolve a trusted RCN
track and supply its sidecar and assignments. All joins use the requested
observation key. No input is modified.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version as package_version
import json
from pathlib import Path
import platform
import sys
from typing import Any, Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import scipy
from scipy.spatial import cKDTree


ANALYSIS_ID = "rcn_anchor_nearest_cells"


class ContractError(RuntimeError):
    """Raised when an input, spatial, or output contract fails."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sidecar", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--spatial-h5ad", type=Path, required=True)
    parser.add_argument("--anchor-rcn", required=True)
    parser.add_argument("--anchor-phenotype", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--obs-key", default="obs_name")
    parser.add_argument("--roi-column", default="sample_roi")
    parser.add_argument("--pair-column", default="pair_id")
    parser.add_argument("--group-column", default="group_merged")
    parser.add_argument("--phenotype-column", default="rcn_phenotype")
    parser.add_argument("--eligible-column", default="include_neighbor_pool")
    parser.add_argument("--source-class-column", default="source_cell_class")
    parser.add_argument("--rcn-column", default="rcn_aligned_80")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--coordinate-unit", default="um")
    parser.add_argument("--use-case", default="")
    parser.add_argument("--track-id", default="")
    parser.add_argument("--expected-neighbor-count", type=int)
    parser.add_argument("--expected-anchor-count", type=int)
    parser.add_argument("--expected-phenotype-count", type=int)
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format="%.17g")


def write_json(value: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def normalize_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    values = series.astype(str).str.strip().str.lower()
    observed = set(values.unique())
    if not observed.issubset({"true", "false"}):
        raise ContractError(
            f"Invalid values in eligible flag: {sorted(observed - {'true', 'false'})}"
        )
    return values.eq("true")


def require_columns(path: Path, columns: Iterable[str]) -> list[str]:
    if not path.is_file():
        raise ContractError(f"Missing input: {path}")
    header = pd.read_csv(path, nrows=0).columns.astype(str).tolist()
    missing = sorted(set(columns) - set(header))
    if missing:
        raise ContractError(f"{path} lacks columns: {missing}")
    return header


def load_coordinates(
    path: Path,
    obs_key: str,
    spatial_key: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not path.is_file():
        raise ContractError(f"Missing spatial H5AD: {path}")
    adata = ad.read_h5ad(path, backed="r")
    try:
        if not adata.obs_names.is_unique:
            raise ContractError("Spatial H5AD obs_names are not unique")
        if spatial_key not in adata.obsm:
            raise ContractError(f'H5AD lacks obsm["{spatial_key}"]')
        coordinates = np.asarray(adata.obsm[spatial_key], dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 2:
            raise ContractError(f"Unexpected spatial shape: {coordinates.shape}")
        if not np.isfinite(coordinates).all():
            raise ContractError("Spatial coordinates contain non-finite values")
        frame = pd.DataFrame(
            {
                obs_key: adata.obs_names.astype(str),
                "x": coordinates[:, 0],
                "y": coordinates[:, 1],
            }
        )
        metadata = {
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "spatial_key": spatial_key,
            "spatial_shape": list(map(int, coordinates.shape)),
        }
    finally:
        adata.file.close()
    return frame, metadata


def load_inputs(args: argparse.Namespace) -> dict[str, Any]:
    sidecar_required = [
        args.obs_key,
        args.roi_column,
        args.pair_column,
        args.group_column,
        args.phenotype_column,
        args.eligible_column,
    ]
    sidecar_header = require_columns(args.sidecar, sidecar_required)
    source_class_available = args.source_class_column in sidecar_header
    sidecar_usecols = list(sidecar_required)
    if source_class_available:
        sidecar_usecols.append(args.source_class_column)
    sidecar = pd.read_csv(args.sidecar, usecols=sidecar_usecols, low_memory=False)
    sidecar[args.obs_key] = sidecar[args.obs_key].astype(str)
    if not sidecar[args.obs_key].is_unique:
        raise ContractError("Sidecar observation key is not unique")
    sidecar[args.eligible_column] = normalize_bool(sidecar[args.eligible_column])
    if source_class_available:
        sidecar["source_class"] = sidecar[args.source_class_column].astype(str)
    else:
        sidecar["source_class"] = sidecar[args.phenotype_column].astype(str)
    sidecar = sidecar.rename(
        columns={
            args.obs_key: "obs_key",
            args.roi_column: "roi",
            args.pair_column: "pair",
            args.group_column: "group",
            args.phenotype_column: "phenotype",
            args.eligible_column: "eligible",
        }
    )
    for column in ["obs_key", "roi", "pair", "group", "phenotype"]:
        sidecar[column] = sidecar[column].astype(str)

    require_columns(args.assignments, [args.obs_key, args.rcn_column])
    assignments = pd.read_csv(
        args.assignments,
        usecols=[args.obs_key, args.rcn_column],
        low_memory=False,
    ).rename(columns={args.obs_key: "obs_key", args.rcn_column: "rcn"})
    assignments["obs_key"] = assignments["obs_key"].astype(str)
    assignments["rcn"] = assignments["rcn"].astype(str)
    if not assignments["obs_key"].is_unique:
        raise ContractError("Assignment observation key is not unique")

    coordinates, spatial_metadata = load_coordinates(
        args.spatial_h5ad, args.obs_key, args.spatial_key
    )
    coordinates = coordinates.rename(columns={args.obs_key: "obs_key"})
    coordinates["obs_key"] = coordinates["obs_key"].astype(str)

    pool = sidecar.loc[
        sidecar["eligible"],
        [
            "obs_key",
            "roi",
            "pair",
            "group",
            "phenotype",
            "source_class",
        ],
    ].copy()
    pool_before = len(pool)
    pool = pool.merge(
        coordinates,
        on="obs_key",
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if (
        len(pool) != pool_before
        or not pool["_merge"].eq("both").all()
        or pool[["x", "y"]].isna().any().any()
    ):
        raise ContractError("Eligible-pool obs-key join to spatial H5AD failed")
    pool = pool.drop(columns="_merge")

    assigned = assignments.loc[
        assignments["rcn"].eq(str(args.anchor_rcn)), ["obs_key"]
    ].copy()
    center_candidates = assigned.merge(
        sidecar,
        on="obs_key",
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if not center_candidates["_merge"].eq("both").all():
        raise ContractError("Frozen-center obs-key join to sidecar failed")
    center_candidates = center_candidates.drop(columns="_merge")
    centers = center_candidates.loc[
        center_candidates["phenotype"].eq(str(args.anchor_phenotype))
    ].copy()
    if centers.empty:
        raise ContractError("No anchor matched the selected RCN and phenotype")
    if not centers["eligible"].all():
        raise ContractError("One or more selected anchors are outside the eligible pool")
    centers_before = len(centers)
    centers = centers.merge(
        coordinates,
        on="obs_key",
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if (
        len(centers) != centers_before
        or not centers["_merge"].eq("both").all()
        or centers[["x", "y"]].isna().any().any()
    ):
        raise ContractError("Anchor obs-key join to spatial H5AD failed")
    centers = centers.drop(columns="_merge")

    phenotypes = sorted(pool["phenotype"].unique().tolist())
    if str(args.anchor_phenotype) not in phenotypes:
        raise ContractError("Anchor phenotype is absent from the eligible target pool")
    if (
        args.expected_neighbor_count is not None
        and len(pool) != args.expected_neighbor_count
    ):
        raise ContractError(
            f"Eligible pool={len(pool)}, expected={args.expected_neighbor_count}"
        )
    if (
        args.expected_anchor_count is not None
        and len(centers) != args.expected_anchor_count
    ):
        raise ContractError(
            f"Anchors={len(centers)}, expected={args.expected_anchor_count}"
        )
    if (
        args.expected_phenotype_count is not None
        and len(phenotypes) != args.expected_phenotype_count
    ):
        raise ContractError(
            f"Phenotypes={len(phenotypes)}, expected={args.expected_phenotype_count}"
        )
    return {
        "sidecar": sidecar,
        "assignments": assignments,
        "pool": pool,
        "centers": centers,
        "phenotypes": phenotypes,
        "spatial_metadata": spatial_metadata,
    }


def stable_nearest(
    anchors: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    query_scope: str,
    exclude_same_obs: bool,
    allow_missing_roi: bool,
) -> pd.DataFrame:
    anchors = anchors.sort_values("obs_key", kind="mergesort").reset_index(drop=True)
    targets = targets.sort_values("obs_key", kind="mergesort").reset_index(drop=True)
    n_anchor = len(anchors)
    target_name = np.full(n_anchor, None, dtype=object)
    target_phenotype = np.full(n_anchor, None, dtype=object)
    target_class = np.full(n_anchor, None, dtype=object)
    target_x = np.full(n_anchor, np.nan)
    target_y = np.full(n_anchor, np.nan)
    distance = np.full(n_anchor, np.nan)
    tie_count = np.zeros(n_anchor, dtype=np.int64)
    target_groups = {
        str(roi): group.sort_values("obs_key", kind="mergesort").reset_index(
            drop=True
        )
        for roi, group in targets.groupby("roi", observed=True, sort=False)
    }
    for roi, anchor_group in anchors.groupby("roi", observed=True, sort=True):
        roi_name = str(roi)
        target_group = target_groups.get(roi_name)
        anchor_indices = anchor_group.index.to_numpy(dtype=np.int64)
        if target_group is None or target_group.empty:
            if allow_missing_roi:
                continue
            raise ContractError(f"No target in anchor ROI {roi_name}")
        anchor_xy = anchor_group[["x", "y"]].to_numpy(dtype=float)
        anchor_names = anchor_group["obs_key"].astype(str).to_numpy()
        target_xy = target_group[["x", "y"]].to_numpy(dtype=float)
        target_names = target_group["obs_key"].astype(str).to_numpy()
        tree = cKDTree(target_xy)
        k = min(len(target_group), 2 if exclude_same_obs else 1)
        approximate_distance, approximate_index = tree.query(
            anchor_xy, k=k, workers=1
        )
        approximate_distance = np.asarray(approximate_distance, dtype=float)
        approximate_index = np.asarray(approximate_index, dtype=np.int64)
        if k == 1:
            approximate_distance = approximate_distance[:, None]
            approximate_index = approximate_index[:, None]
        for local_anchor, global_anchor in enumerate(anchor_indices):
            initial = [
                (float(value), int(index))
                for value, index in zip(
                    approximate_distance[local_anchor],
                    approximate_index[local_anchor],
                )
                if not (
                    exclude_same_obs
                    and target_names[int(index)] == anchor_names[local_anchor]
                )
            ]
            if not initial:
                if allow_missing_roi:
                    continue
                raise ContractError(
                    f"No non-self target for {anchor_names[local_anchor]}"
                )
            minimum = min(value[0] for value in initial)
            tolerance = max(1e-9, abs(minimum) * 1e-12)
            exact: list[tuple[str, int, float]] = []
            for index in tree.query_ball_point(
                anchor_xy[local_anchor], minimum + tolerance
            ):
                index = int(index)
                if exclude_same_obs and target_names[index] == anchor_names[local_anchor]:
                    continue
                exact_distance = float(
                    np.linalg.norm(target_xy[index] - anchor_xy[local_anchor])
                )
                if exact_distance <= minimum + tolerance:
                    exact.append((target_names[index], index, exact_distance))
            if not exact:
                _, index = min(initial, key=lambda value: value[0])
                exact = [
                    (
                        target_names[index],
                        index,
                        float(
                            np.linalg.norm(
                                target_xy[index] - anchor_xy[local_anchor]
                            )
                        ),
                    )
                ]
            chosen_name, chosen_index, chosen_distance = min(
                exact, key=lambda value: value[0]
            )
            target_name[global_anchor] = chosen_name
            target_phenotype[global_anchor] = str(
                target_group["phenotype"].iloc[chosen_index]
            )
            target_class[global_anchor] = str(
                target_group["source_class"].iloc[chosen_index]
            )
            target_x[global_anchor] = target_xy[chosen_index, 0]
            target_y[global_anchor] = target_xy[chosen_index, 1]
            distance[global_anchor] = chosen_distance
            tie_count[global_anchor] = len(exact)

    result = anchors[
        ["obs_key", "roi", "pair", "group", "x", "y"]
    ].rename(
        columns={
            "obs_key": "anchor_obs_key",
            "x": "anchor_x",
            "y": "anchor_y",
        }
    )
    result["target_obs_key"] = pd.array(target_name, dtype="string")
    result["target_phenotype"] = pd.array(target_phenotype, dtype="string")
    result["target_source_class"] = pd.array(target_class, dtype="string")
    result["target_x"] = target_x
    result["target_y"] = target_y
    result["distance"] = distance
    result["equal_distance_target_count"] = tie_count
    result["target_present_in_roi"] = result["target_obs_key"].notna()
    result["query_scope"] = query_scope
    result["tie_break_rule"] = (
        "lexicographically smallest target observation key at equal distance"
    )
    return result


def nearest_category_outputs(
    nearest_any: pd.DataFrame,
    nearest_heterotypic: pd.DataFrame,
    phenotypes: list[str],
    anchor_phenotype: str,
) -> dict[str, pd.DataFrame]:
    nearest = pd.concat([nearest_any, nearest_heterotypic], ignore_index=True)
    pooled = (
        nearest.groupby(["query_scope", "target_phenotype"], observed=True)
        .agg(
            n_anchors=("anchor_obs_key", "size"),
            median_conditional_distance=("distance", "median"),
            q25_conditional_distance=("distance", lambda x: float(x.quantile(0.25))),
            q75_conditional_distance=("distance", lambda x: float(x.quantile(0.75))),
        )
        .reset_index()
    )
    pooled["n_anchors_total"] = pooled.groupby("query_scope")[
        "n_anchors"
    ].transform("sum")
    pooled["cell_pooled_fraction"] = (
        pooled["n_anchors"] / pooled["n_anchors_total"]
    )
    roi_totals = (
        nearest.groupby(["query_scope", "roi", "pair", "group"], observed=True)
        .size()
        .rename("n_anchors_roi")
        .reset_index()
    )
    roi_frames = []
    for scope, scope_data in nearest.groupby("query_scope", observed=True):
        allowed = list(phenotypes)
        if scope == "nearest_heterotypic":
            allowed = [value for value in allowed if value != anchor_phenotype]
        counts = (
            scope_data.groupby(
                ["roi", "pair", "group", "target_phenotype"], observed=True
            )
            .size()
            .rename("n_nearest")
            .reset_index()
        )
        base = roi_totals.loc[roi_totals["query_scope"].eq(scope)].drop(
            columns="query_scope"
        )
        grid = base.assign(_key=1).merge(
            pd.DataFrame({"target_phenotype": allowed, "_key": 1}),
            on="_key",
        ).drop(columns="_key")
        grid = grid.merge(
            counts,
            on=["roi", "pair", "group", "target_phenotype"],
            how="left",
            validate="one_to_one",
        )
        grid["n_nearest"] = grid["n_nearest"].fillna(0).astype(int)
        grid["fraction_within_roi"] = grid["n_nearest"] / grid["n_anchors_roi"]
        grid["query_scope"] = scope
        roi_frames.append(grid)
    by_roi = pd.concat(roi_frames, ignore_index=True)
    macro = (
        by_roi.groupby(["query_scope", "target_phenotype"], observed=True)
        .agg(
            n_rois=("roi", "nunique"),
            n_positive_rois=("n_nearest", lambda x: int((x > 0).sum())),
            mean_roi_fraction=("fraction_within_roi", "mean"),
            median_roi_fraction=("fraction_within_roi", "median"),
            q25_roi_fraction=("fraction_within_roi", lambda x: float(x.quantile(0.25))),
            q75_roi_fraction=("fraction_within_roi", lambda x: float(x.quantile(0.75))),
        )
        .reset_index()
    )
    summary = pooled.merge(
        macro,
        on=["query_scope", "target_phenotype"],
        validate="one_to_one",
    ).sort_values(
        ["query_scope", "cell_pooled_fraction"],
        ascending=[True, False],
        kind="mergesort",
    )
    by_pair = (
        by_roi.groupby(
            ["query_scope", "pair", "group", "target_phenotype"], observed=True
        )
        .agg(
            n_nearest=("n_nearest", "sum"),
            n_anchors_pair_group=("n_anchors_roi", "sum"),
            n_rois=("roi", "nunique"),
        )
        .reset_index()
    )
    by_pair["fraction_within_pair_group"] = (
        by_pair["n_nearest"] / by_pair["n_anchors_pair_group"]
    )
    return {"summary": summary, "by_roi": by_roi, "by_pair_group": by_pair}


def phenotype_distance_outputs(
    centers: pd.DataFrame,
    pool: pd.DataFrame,
    phenotypes: list[str],
    anchor_phenotype: str,
) -> dict[str, pd.DataFrame]:
    frames = []
    for phenotype in phenotypes:
        frame = stable_nearest(
            centers,
            pool.loc[pool["phenotype"].eq(phenotype)],
            query_scope=f"nearest_{phenotype}",
            exclude_same_obs=phenotype == anchor_phenotype,
            allow_missing_roi=True,
        )
        frame["requested_phenotype"] = phenotype
        frames.append(frame)
    long = pd.concat(frames, ignore_index=True)
    pooled = (
        long.groupby("requested_phenotype", observed=True)
        .agg(
            n_anchors_total=("anchor_obs_key", "size"),
            n_anchors_with_target=("target_present_in_roi", "sum"),
            median_distance=("distance", "median"),
            q25_distance=("distance", lambda x: float(x.quantile(0.25))),
            q75_distance=("distance", lambda x: float(x.quantile(0.75))),
        )
        .reset_index()
    )
    pooled["anchor_coverage"] = (
        pooled["n_anchors_with_target"] / pooled["n_anchors_total"]
    )
    by_roi = (
        long.groupby(
            ["requested_phenotype", "roi", "pair", "group"], observed=True
        )
        .agg(
            n_anchors=("anchor_obs_key", "size"),
            n_anchors_with_target=("target_present_in_roi", "sum"),
            median_distance=("distance", "median"),
        )
        .reset_index()
    )
    by_roi["target_present_in_roi"] = by_roi["n_anchors_with_target"].gt(0)
    roi_macro = (
        by_roi.groupby("requested_phenotype", observed=True)
        .agg(
            n_rois=("roi", "nunique"),
            n_rois_with_target=("target_present_in_roi", "sum"),
            median_of_roi_medians=("median_distance", "median"),
            q25_of_roi_medians=("median_distance", lambda x: float(x.quantile(0.25))),
            q75_of_roi_medians=("median_distance", lambda x: float(x.quantile(0.75))),
        )
        .reset_index()
    )
    roi_macro["roi_coverage"] = (
        roi_macro["n_rois_with_target"] / roi_macro["n_rois"]
    )
    summary = pooled.merge(
        roi_macro,
        on="requested_phenotype",
        validate="one_to_one",
    ).sort_values(
        ["median_distance", "requested_phenotype"], kind="mergesort"
    )
    return {
        "long": long,
        "summary": summary,
        "by_roi": by_roi,
    }


def plot_categories(summary: pd.DataFrame, path: Path, dpi: int) -> None:
    scopes = [
        ("nearest_any_eligible", "Nearest eligible cell"),
        ("nearest_heterotypic", "Nearest heterotypic cell"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    for ax, (scope, title) in zip(axes, scopes):
        data = (
            summary.loc[summary["query_scope"].eq(scope)]
            .nlargest(12, "cell_pooled_fraction")
            .sort_values("cell_pooled_fraction")
        )
        ax.barh(
            data["target_phenotype"],
            100 * data["cell_pooled_fraction"],
            color="#5B8FA8",
            label="cell-pooled",
        )
        ax.scatter(
            100 * data["mean_roi_fraction"],
            data["target_phenotype"],
            color="black",
            s=25,
            label="mean across ROI",
            zorder=3,
        )
        ax.set_xlabel("anchors (%)")
        ax.set_title(title)
        ax.grid(axis="x", color="#dddddd", linewidth=0.7)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Identity of nearest cells to fixed RCN anchors", fontweight="bold")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_distances(summary: pd.DataFrame, path: Path, dpi: int, unit: str) -> None:
    data = summary.nsmallest(15, "median_distance").sort_values(
        "median_distance", ascending=False
    )
    median = data["median_distance"].to_numpy()
    low = median - data["q25_distance"].to_numpy()
    high = data["q75_distance"].to_numpy() - median
    roi_median = data["median_of_roi_medians"].to_numpy()
    roi_low = roi_median - data["q25_of_roi_medians"].to_numpy()
    roi_high = data["q75_of_roi_medians"].to_numpy() - roi_median
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(10.5, 8), constrained_layout=True)
    ax.errorbar(
        median,
        y + 0.13,
        xerr=np.vstack([low, high]),
        fmt="o",
        capsize=2,
        label="cell-pooled median [IQR]",
    )
    ax.errorbar(
        roi_median,
        y - 0.13,
        xerr=np.vstack([roi_low, roi_high]),
        fmt="s",
        capsize=2,
        label="median of ROI medians [IQR]",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(data["requested_phenotype"])
    ax.set_xlabel(f"nearest phenotype-specific distance ({unit})")
    ax.set_title("Closest phenotype-specific targets")
    ax.grid(axis="x", color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def synthetic_test() -> bool:
    anchors = pd.DataFrame(
        {
            "obs_key": ["a1", "b1"],
            "roi": ["A", "B"],
            "pair": ["P1", "P2"],
            "group": ["X", "Y"],
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
        }
    )
    targets = pd.DataFrame(
        {
            "obs_key": ["a2", "a0", "b2"],
            "roi": ["A", "A", "B"],
            "pair": ["P1", "P1", "P2"],
            "group": ["X", "X", "Y"],
            "phenotype": ["T1", "T2", "T3"],
            "source_class": ["T1", "T2", "T3"],
            "x": [1.0, -1.0, 0.25],
            "y": [0.0, 0.0, 0.0],
        }
    )
    result = stable_nearest(
        anchors,
        targets,
        query_scope="synthetic",
        exclude_same_obs=False,
        allow_missing_roi=False,
    )
    return result["target_obs_key"].astype(str).tolist() == ["a0", "b2"]


def validate(
    args: argparse.Namespace,
    pool: pd.DataFrame,
    centers: pd.DataFrame,
    phenotypes: list[str],
    nearest_any: pd.DataFrame,
    nearest_heterotypic: pd.DataFrame,
    category: dict[str, pd.DataFrame],
    distance: dict[str, pd.DataFrame],
) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []

    def add(check: str, passed: bool, observed: Any, expected: Any) -> None:
        checks.append(
            {
                "check": check,
                "status": "PASS" if passed else "FAIL",
                "observed": observed,
                "expected": expected,
            }
        )

    pool_index = pool.set_index("obs_key")
    for scope, frame in [
        ("nearest_any", nearest_any),
        ("nearest_heterotypic", nearest_heterotypic),
    ]:
        target_roi = frame["target_obs_key"].map(pool_index["roi"].astype(str))
        target_label = frame["target_obs_key"].map(
            pool_index["phenotype"].astype(str)
        )
        recomputed = np.hypot(
            frame["target_x"] - frame["anchor_x"],
            frame["target_y"] - frame["anchor_y"],
        )
        passed = (
            len(frame) == len(centers)
            and frame["anchor_obs_key"].is_unique
            and frame["target_obs_key"].notna().all()
            and ~frame["anchor_obs_key"].astype(str).eq(
                frame["target_obs_key"].astype(str)
            ).any()
            and frame["roi"].astype(str).eq(target_roi).all()
            and frame["target_phenotype"].astype(str).eq(target_label).all()
            and np.allclose(frame["distance"], recomputed, atol=1e-9, rtol=1e-12)
        )
        if scope == "nearest_heterotypic":
            passed = passed and ~frame["target_phenotype"].eq(
                args.anchor_phenotype
            ).any()
        add(scope, passed, len(frame), len(centers))

    pooled_sum = category["summary"].groupby("query_scope")[
        "cell_pooled_fraction"
    ].sum()
    roi_sum = category["by_roi"].groupby(["query_scope", "roi"])[
        "fraction_within_roi"
    ].sum()
    add(
        "category_fractions_sum_to_one",
        np.allclose(pooled_sum, 1.0, atol=1e-12)
        and np.allclose(roi_sum, 1.0, atol=1e-12),
        {"pooled": pooled_sum.to_dict(), "roi_min": float(roi_sum.min()), "roi_max": float(roi_sum.max())},
        1.0,
    )
    long = distance["long"]
    finite = long.loc[long["target_present_in_roi"]]
    target_roi = finite["target_obs_key"].map(pool_index["roi"].astype(str))
    target_label = finite["target_obs_key"].map(
        pool_index["phenotype"].astype(str)
    )
    recomputed = np.hypot(
        finite["target_x"] - finite["anchor_x"],
        finite["target_y"] - finite["anchor_y"],
    )
    expected_long = len(centers) * len(phenotypes)
    add(
        "phenotype_distance_grid",
        len(long) == expected_long
        and long["requested_phenotype"].nunique() == len(phenotypes)
        and finite["roi"].astype(str).eq(target_roi).all()
        and finite["requested_phenotype"].astype(str).eq(target_label).all()
        and np.allclose(finite["distance"], recomputed, atol=1e-9, rtol=1e-12),
        len(long),
        expected_long,
    )
    add("synthetic_roi_tie_test", synthetic_test(), "passed", "passed")
    failed = [row["check"] for row in checks if row["status"] != "PASS"]
    if failed:
        raise ContractError("Validation failed: " + ", ".join(failed))
    return checks


def build_manifest(out_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(value for value in out_dir.rglob("*") if value.is_file()):
        if path.name == "output_manifest.csv":
            continue
        rows.append(
            {
                "relative_path": path.relative_to(out_dir).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    if args.dpi < 72:
        raise ContractError("--dpi must be at least 72")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs(args)
    pool = inputs["pool"]
    centers = inputs["centers"]
    phenotypes = inputs["phenotypes"]
    nearest_any = stable_nearest(
        centers,
        pool,
        query_scope="nearest_any_eligible",
        exclude_same_obs=True,
        allow_missing_roi=False,
    )
    nearest_heterotypic = stable_nearest(
        centers,
        pool.loc[~pool["phenotype"].eq(args.anchor_phenotype)],
        query_scope="nearest_heterotypic",
        exclude_same_obs=False,
        allow_missing_roi=False,
    )
    category = nearest_category_outputs(
        nearest_any, nearest_heterotypic, phenotypes, args.anchor_phenotype
    )
    distance = phenotype_distance_outputs(
        centers, pool, phenotypes, args.anchor_phenotype
    )

    write_csv(centers, args.out_dir / "tables/anchor_centers.csv.gz")
    write_csv(
        nearest_any,
        args.out_dir / "tables/nearest_any_eligible_cell_per_anchor.csv.gz",
    )
    write_csv(
        nearest_heterotypic,
        args.out_dir / "tables/nearest_heterotypic_cell_per_anchor.csv.gz",
    )
    write_csv(category["summary"], args.out_dir / "tables/nearest_category_summary.csv")
    write_csv(category["by_roi"], args.out_dir / "tables/nearest_category_by_roi.csv")
    write_csv(
        category["by_pair_group"],
        args.out_dir / "tables/nearest_category_by_pair_group.csv",
    )
    write_csv(
        distance["long"],
        args.out_dir / "tables/anchor_to_phenotype_nearest_distance.csv.gz",
    )
    write_csv(
        distance["summary"],
        args.out_dir / "tables/phenotype_distance_summary.csv",
    )
    write_csv(
        distance["by_roi"],
        args.out_dir / "tables/phenotype_distance_by_roi.csv",
    )
    plot_categories(
        category["summary"],
        args.out_dir / "figures/nearest_celltype_summary.png",
        args.dpi,
    )
    plot_distances(
        distance["summary"],
        args.out_dir / "figures/phenotype_distance_ranking.png",
        args.dpi,
        args.coordinate_unit,
    )
    checks = validate(
        args,
        pool,
        centers,
        phenotypes,
        nearest_any,
        nearest_heterotypic,
        category,
        distance,
    )
    for path in sorted((args.out_dir / "figures").glob("*.png")):
        with Image.open(path) as image:
            image.verify()
        checks.append(
            {
                "check": f"png_open_{path.stem}",
                "status": "PASS",
                "observed": path.stat().st_size,
                "expected": "readable non-empty PNG",
            }
        )
    write_csv(pd.DataFrame(checks), args.out_dir / "validation/validation.csv")
    input_manifest = pd.DataFrame(
        [
            {
                "role": role,
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for role, path in [
                ("sidecar", args.sidecar),
                ("assignments", args.assignments),
                ("spatial_h5ad", args.spatial_h5ad),
            ]
        ]
    )
    write_csv(input_manifest, args.out_dir / "input_manifest.csv")
    write_json(
        {
            "analysis_id": ANALYSIS_ID,
            "completed_at_utc": utc_now(),
            "command": " ".join(map(str, sys.argv)),
            "anchor_rcn": args.anchor_rcn,
            "anchor_phenotype": args.anchor_phenotype,
            "n_anchors": len(centers),
            "n_neighbor_pool": len(pool),
            "n_phenotypes": len(phenotypes),
            "n_rois": int(centers["roi"].nunique()),
            "n_pair_group_units": int(
                centers[["pair", "group"]].drop_duplicates().shape[0]
            ),
            "coordinate_unit": args.coordinate_unit,
            "use_case": args.use_case,
            "track_id": args.track_id,
            "same_roi_only": True,
            "focal_cell_excluded": True,
            "rcn_recomputed": False,
            "annotation_modified": False,
            "python": platform.python_version(),
            "anndata": package_version("anndata"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "validation_passed": len(checks),
            "validation_total": len(checks),
            "status": "PASS",
        },
        args.out_dir / "run_manifest.json",
    )
    manifest = build_manifest(args.out_dir)
    write_csv(manifest, args.out_dir / "output_manifest.csv")
    print(
        json.dumps(
            {
                "status": "PASS",
                "out_dir": str(args.out_dir),
                "anchors": len(centers),
                "neighbor_pool": len(pool),
                "phenotypes": len(phenotypes),
                "validation": f"{len(checks)}/{len(checks)}",
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
