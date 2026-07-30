#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map 10x Visium HD CellCharter niche labels from h5ad onto embedded H&E images."
    )
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--cluster-key", default="cluster_cellcharter")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--image-key", choices=("hires", "lowres"), default="hires")
    parser.add_argument("--output-prefix", default="visiumhd_cellcharter")
    parser.add_argument("--spot-size", type=float, default=0.9)
    parser.add_argument("--single-niche-spot-size", type=float, default=0.55)
    parser.add_argument("--alpha", type=float, default=0.82)
    parser.add_argument("--single-niche-alpha", type=float, default=0.9)
    parser.add_argument("--dpi", type=int, default=360)
    parser.add_argument("--write-per-sample", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--write-per-niche", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def sorted_categories(values: pd.Series) -> list[str]:
    if isinstance(values.dtype, pd.CategoricalDtype):
        return [str(value) for value in values.cat.categories]
    as_str = values.astype(str)
    try:
        return [str(value) for value in sorted({int(value) for value in as_str})]
    except ValueError:
        return sorted(as_str.unique())


def palette_for_niches(niches: list[str]) -> dict[str, Any]:
    cmap = plt.get_cmap("tab20" if len(niches) <= 20 else "hsv")
    return {niche: cmap(i / max(1, len(niches) - 1)) for i, niche in enumerate(niches)}


def legend_handles(palette: dict[str, Any]) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=6,
            label=f"niche {niche}",
        )
        for niche, color in palette.items()
    ]


def require_inputs(adata: Any, args: argparse.Namespace) -> None:
    if args.sample_key not in adata.obs:
        raise ValueError(f"Missing obs[{args.sample_key!r}].")
    if args.cluster_key not in adata.obs:
        raise ValueError(f"Missing obs[{args.cluster_key!r}].")
    if args.spatial_key not in adata.obsm:
        raise ValueError(f"Missing obsm[{args.spatial_key!r}].")
    if "spatial" not in adata.uns or not isinstance(adata.uns["spatial"], dict):
        raise ValueError("Expected sample-keyed H&E images under uns['spatial'].")


def get_image_and_scale(adata: Any, sample_id: str, image_key: str) -> tuple[np.ndarray, float, dict[str, Any]]:
    entry = adata.uns["spatial"].get(sample_id)
    if not isinstance(entry, dict):
        raise ValueError(f"Missing uns['spatial'][{sample_id!r}].")
    images = entry.get("images", {})
    scalefactors = entry.get("scalefactors", {})
    if image_key not in images:
        available = ", ".join(sorted(images)) or "none"
        raise ValueError(f"{sample_id}: image {image_key!r} not found; available: {available}.")
    scale = float(scalefactors.get(f"tissue_{image_key}_scalef", 1.0))
    return np.asarray(images[image_key]), scale, dict(scalefactors)


def build_mapping_table(adata: Any, args: argparse.Namespace) -> pd.DataFrame:
    coords = np.asarray(adata.obsm[args.spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(f"obsm[{args.spatial_key!r}] must contain at least two columns.")

    frame = pd.DataFrame(
        {
            "obs_name": adata.obs_names,
            "sample_id": adata.obs[args.sample_key].astype(str).to_numpy(),
            "niche": adata.obs[args.cluster_key].astype(str).to_numpy(),
            "fullres_col": coords[:, 0],
            "fullres_row": coords[:, 1],
        }
    )

    rows = []
    for sample_id, sample_frame in frame.groupby("sample_id", sort=False):
        _, scale, scalefactors = get_image_and_scale(adata, sample_id, args.image_key)
        sample_out = sample_frame.copy()
        sample_out["image_key"] = args.image_key
        sample_out["tissue_image_scale"] = scale
        sample_out["he_col"] = sample_out["fullres_col"] * scale
        sample_out["he_row"] = sample_out["fullres_row"] * scale
        sample_out["bin_size_um"] = scalefactors.get("bin_size_um", np.nan)
        sample_out["spot_diameter_fullres"] = scalefactors.get("spot_diameter_fullres", np.nan)
        sample_out["microns_per_pixel"] = scalefactors.get("microns_per_pixel", np.nan)
        rows.append(sample_out)
    return pd.concat(rows, axis=0, ignore_index=True)


def draw_all_niches(
    ax: Any,
    image: np.ndarray,
    sample_frame: pd.DataFrame,
    palette: dict[str, Any],
    spot_size: float,
    alpha: float,
) -> None:
    ax.imshow(image)
    colors = sample_frame["niche"].map(palette).to_numpy()
    ax.scatter(
        sample_frame["he_col"],
        sample_frame["he_row"],
        c=list(colors),
        s=spot_size,
        linewidths=0,
        alpha=alpha,
        rasterized=True,
    )
    ax.set_axis_off()


def save_overview(
    adata: Any,
    mapping: pd.DataFrame,
    samples: list[str],
    palette: dict[str, Any],
    args: argparse.Namespace,
    outdir: Path,
) -> Path:
    ncols = min(2, max(1, len(samples)))
    nrows = int(np.ceil(len(samples) / ncols))
    output = outdir / f"{args.output_prefix}_all_samples_niches_on_he_overview.png"
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 6 * nrows), constrained_layout=True)
    axes_arr = np.asarray(axes).reshape(-1)
    for ax, sample_id in zip(axes_arr, samples):
        image, _, _ = get_image_and_scale(adata, sample_id, args.image_key)
        sample_frame = mapping.loc[mapping["sample_id"] == sample_id]
        draw_all_niches(ax, image, sample_frame, palette, args.spot_size, args.alpha)
        ax.set_title(sample_id, fontsize=10)
    for ax in axes_arr[len(samples) :]:
        ax.set_axis_off()
    fig.legend(
        handles=legend_handles(palette),
        title="CellCharter",
        loc="center right",
        bbox_to_anchor=(1.08, 0.5),
        fontsize=7,
    )
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def save_sample_overlay(
    image: np.ndarray,
    sample_frame: pd.DataFrame,
    sample_id: str,
    palette: dict[str, Any],
    args: argparse.Namespace,
    outdir: Path,
) -> Path:
    output = outdir / "per_sample_all_niches" / f"{sample_id}_{args.output_prefix}_all_niches_on_he.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    draw_all_niches(ax, image, sample_frame, palette, args.spot_size, args.alpha)
    ax.set_title(f"{sample_id}: CellCharter niches on H&E", fontsize=10)
    fig.legend(handles=legend_handles(palette), title="CellCharter", loc="center right", bbox_to_anchor=(1.24, 0.5), fontsize=7)
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def save_single_niche_overlay(
    image: np.ndarray,
    sample_frame: pd.DataFrame,
    sample_id: str,
    niche: str,
    color: Any,
    args: argparse.Namespace,
    outdir: Path,
) -> Path:
    output = outdir / "per_sample_per_niche" / sample_id / f"{sample_id}_{args.output_prefix}_niche_{niche}_on_he.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    selected = sample_frame.loc[sample_frame["niche"] == niche]
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    ax.imshow(image)
    ax.scatter(
        sample_frame["he_col"],
        sample_frame["he_row"],
        c="white",
        s=max(args.single_niche_spot_size * 0.45, 0.12),
        linewidths=0,
        alpha=0.16,
        rasterized=True,
    )
    ax.scatter(
        selected["he_col"],
        selected["he_row"],
        c=[color],
        s=args.single_niche_spot_size,
        linewidths=0,
        alpha=args.single_niche_alpha,
        rasterized=True,
    )
    ax.set_title(f"{sample_id}: niche {niche} on H&E (n={selected.shape[0]})", fontsize=10)
    ax.set_axis_off()
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def main() -> None:
    args = parse_args()
    if not args.input_h5ad.exists():
        raise FileNotFoundError(args.input_h5ad)
    outdir = args.outdir or args.input_h5ad.parent / "plots" / "he_niche_mapping"
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"Output directory is not empty; pass --overwrite: {outdir}")
    if outdir.exists() and args.overwrite:
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input_h5ad)
    require_inputs(adata, args)
    samples = sorted_categories(adata.obs[args.sample_key])
    niches = sorted_categories(adata.obs[args.cluster_key])
    palette = palette_for_niches(niches)
    mapping = build_mapping_table(adata, args)

    written: list[Path] = []
    mapping_path = outdir / f"{args.output_prefix}_bin_to_he_coordinates.csv"
    mapping.to_csv(mapping_path, index=False)
    written.append(mapping_path)
    written.append(save_overview(adata, mapping, samples, palette, args, outdir))

    if args.write_per_sample or args.write_per_niche:
        for sample_id in samples:
            image, _, _ = get_image_and_scale(adata, sample_id, args.image_key)
            sample_frame = mapping.loc[mapping["sample_id"] == sample_id]
            if args.write_per_sample:
                written.append(save_sample_overlay(image, sample_frame, sample_id, palette, args, outdir))
            if args.write_per_niche:
                for niche in niches:
                    written.append(save_single_niche_overlay(image, sample_frame, sample_id, niche, palette[niche], args, outdir))

    counts = pd.crosstab(mapping["sample_id"], mapping["niche"]).reindex(index=samples, columns=niches, fill_value=0)
    counts_path = outdir / f"{args.output_prefix}_niche_bin_counts_by_sample.csv"
    counts.to_csv(counts_path)
    written.append(counts_path)

    params_path = outdir / f"{args.output_prefix}_he_mapping_parameters.json"
    output_paths = [*written, params_path]
    write_json(
        params_path,
        {
            "input_h5ad": str(args.input_h5ad.resolve()),
            "outdir": str(outdir.resolve()),
            "sample_key": args.sample_key,
            "cluster_key": args.cluster_key,
            "spatial_key": args.spatial_key,
            "image_key": args.image_key,
            "output_prefix": args.output_prefix,
            "n_samples": len(samples),
            "samples": samples,
            "n_niches": len(niches),
            "niches": niches,
            "n_bins": int(mapping.shape[0]),
            "dpi": args.dpi,
            "spot_size": args.spot_size,
            "write_per_sample": args.write_per_sample,
            "write_per_niche": args.write_per_niche,
            "outputs": [str(path.resolve()) for path in output_paths],
        },
    )
    written.append(params_path)

    print(f"Wrote {len(written)} files under {outdir}")
    print(f"Coordinate table: {mapping_path}")


if __name__ == "__main__":
    main()
