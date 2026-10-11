"""QC and PNG rendering of existing segmentation; no segmentation models are run."""
from __future__ import annotations

from pathlib import Path
import argparse
import sys
_SHARED = next(p / "_shared/spatial" for p in Path(__file__).resolve().parents if (p / "_shared/spatial/run_contract.py").is_file())
sys.path.insert(0, str(_SHARED))
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from run_contract import assert_new_outputs, file_sha256, write_json
from rings import load_rings, ring_vertex_counts


def read_cells(path):
    cells = pd.read_parquet(path)
    if "cell_id" not in cells or cells.cell_id.isna().any() or cells.cell_id.duplicated().any():
        raise ValueError("Unique nonmissing cell_id required within this single segmentation sample")
    return cells


def flatten_rings(rings, cells):
    vertices, ids, labels = [], [], []
    for cid, subrings in rings.items():
        if cid not in cells.index:
            raise ValueError(f"Boundary cell {cid!r} missing from cells.parquet")
        for label, ring in subrings.items():
            if len(ring) < 3:
                raise ValueError(f"Truncated ring {cid!r}/{label!r}: {len(ring)} vertices")
            vertices.append(ring)
            ids.append(cid)
            labels.append(label)
    return vertices, ids, labels


def ring_collection(rings, cells, color_by=None, log_color=False, nuclei=False):
    vertices, ids, labels = flatten_rings(rings, cells)
    kwargs = {"edgecolors":"black" if nuclei else "#555555", "linewidths":0.4,
              "facecolors":"none"}
    legend, scalar = [], None
    if color_by:
        if color_by not in cells:
            raise ValueError(f"Unknown color column: {color_by}")
        values = cells.loc[ids, color_by]
        if pd.api.types.is_numeric_dtype(values):
            arr = values.to_numpy(float)
            finite = arr[np.isfinite(arr) & ((arr > 0) if log_color else True)]
            if not finite.size:
                raise ValueError("No finite valid values for polygon color")
            low, high = float(finite.min()), float(finite.max())
            norm = LogNorm(low, max(high, low*1.001)) if log_color else Normalize(low, max(high,low+1e-6))
            kwargs.pop("facecolors")
            kwargs.update(array=np.ma.masked_invalid(arr),cmap="viridis",norm=norm)
            scalar = plt.cm.ScalarMappable(norm=norm,cmap="viridis")
        else:
            categories = sorted(values.astype(str).unique())
            colors = {v:plt.cm.tab20(i/max(len(categories)-1,1)) for i,v in enumerate(categories)}
            kwargs["facecolors"] = [colors[v] for v in values.astype(str)]
            legend = [Patch(facecolor=colors[v],label=v) for v in categories]
    collection = PolyCollection(vertices, **kwargs)
    return collection, scalar, legend, list(zip(ids, labels))


def render_current(input_dir, output_png, *, bbox=None, color_by=None, log_color=False,
                   overlay_nucleus=True, dpi=160):
    if Path(output_png).suffix.lower() != ".png":
        raise ValueError("segmentation-plot-current outputs PNG only")
    assert_new_outputs([output_png, Path(output_png).with_suffix(".json")])
    cells = read_cells(input_dir/"cells.parquet").set_index("cell_id")
    # Read complete rings. The axes perform clipping, including rings crossing a crop
    # whose centroids/all vertices lie outside it. No vertex-level crop distortion.
    cell_rings = load_rings(input_dir/"cell_boundaries.parquet")
    nucleus_rings = load_rings(input_dir/"nucleus_boundaries.parquet") if overlay_nucleus else {}
    pc, scalar, legend, cell_keys = ring_collection(cell_rings,cells,color_by,log_color)
    npc, _, _, nucleus_keys = ring_collection(nucleus_rings,cells,nuclei=True)
    fig, ax = plt.subplots(figsize=(8,8))
    ax.add_collection(pc)
    if nucleus_keys:
        ax.add_collection(npc)
    if bbox is None:
        ax.autoscale_view()
        ax.invert_yaxis()
    else:
        xmin,xmax,ymin,ymax = bbox
        if not np.isfinite(bbox).all() or not xmin < xmax or not ymin < ymax:
            plt.close(fig)
            raise ValueError("bbox must be finite with xmin<xmax and ymin<ymax")
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymax,ymin)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    if scalar:
        fig.colorbar(scalar,ax=ax,label=color_by)
    if legend:
        ax.legend(handles=legend,fontsize=7)
    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output_png,dpi=dpi,bbox_inches="tight")
    plt.close(fig)
    from PIL import Image
    with Image.open(output_png) as img:
        img.verify()
    manifest = {"status":"PASS", "entrypoint":"segmentation-plot-current",
        "cell_rings":len(cell_keys), "nucleus_rings":len(nucleus_keys),
        "loaded_ring_scope":"all cells/rings; complete polygons clipped by axes",
        "bbox":list(bbox) if bbox else None,"overlay_nucleus":overlay_nucleus,
        "output_png":str(output_png.resolve()),"png_sha256":file_sha256(output_png),
        "inputs":{name:file_sha256(input_dir/name) for name in
                  ["cells.parquet","cell_boundaries.parquet"] + (["nucleus_boundaries.parquet"] if overlay_nucleus else [])}}
    write_json(output_png.with_suffix(".json"),manifest)
    return manifest


def qc_current(input_dir, out_dir):
    assert_new_outputs([out_dir/"segmentation_summary.json"])
    cells = read_cells(input_dir/"cells.parquet")
    tables = {}
    for kind in ("cell", "nucleus"):
        table = ring_vertex_counts(input_dir/f"{kind}_boundaries.parquet")
        if not table.cell_id.isin(cells.cell_id).all():
            raise ValueError("Boundary identifiers missing from cells.parquet")
        table["truncated_ring"] = table.n_vertices < 3
        tables[kind] = table
    summary = {"status":"PASS", "entrypoint":"segmentation-qc-current", "n_cells":len(cells),
        "cell_rings":len(tables["cell"]),"nucleus_rings":len(tables["nucleus"]),
        "truncated_cell_rings":int(tables["cell"].truncated_ring.sum()),
        "truncated_nucleus_rings":int(tables["nucleus"].truncated_ring.sum()),
        "cells_without_cell_rings":int((~cells.cell_id.isin(tables["cell"].cell_id)).sum()),
        "cells_without_nucleus_rings":int((~cells.cell_id.isin(tables["nucleus"].cell_id)).sum()),
        "segmentation_rerun":False,
        "inputs":{name:file_sha256(input_dir/name) for name in
                  ["cells.parquet","cell_boundaries.parquet","nucleus_boundaries.parquet"]}}
    if summary["truncated_cell_rings"] or summary["truncated_nucleus_rings"] or summary["cells_without_cell_rings"]:
        summary["status"] = "WARN"
    out_dir.mkdir(parents=True,exist_ok=True)
    for kind,table in tables.items():
        table.to_csv(out_dir/f"{kind}_ring_vertex_counts.csv",index=False)
    if "segmentation_method" in cells:
        cells.segmentation_method.value_counts(dropna=False).rename("n_cells").to_csv(out_dir/"segmentation_method_breakdown.csv")
    columns = [c for c in ["cell_area","nucleus_area","nucleus_count","total_counts"] if c in cells]
    fig,axes = plt.subplots(1,max(len(columns),1),figsize=(4*max(len(columns),1),4),squeeze=False)
    for col,ax in zip(columns,axes.ravel()):
        values = cells[col].to_numpy(float)
        ax.hist(values[np.isfinite(values)],bins=30)
        ax.set_title(col)
    fig.savefig(out_dir/"segmentation_qc.png",bbox_inches="tight")
    plt.close(fig)
    write_json(out_dir/"segmentation_summary.json",summary)
    return summary


def qc_main(argv=None):
    parser = argparse.ArgumentParser(description="Inspect existing segmented cells and all boundary rings")
    parser.add_argument("--input-dir",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    args = parser.parse_args(argv)
    result = qc_current(args.input_dir,args.out_dir)
    print(f"{result['status']} segmentation-qc-current: {result['cell_rings']} cell rings; {result['nucleus_rings']} nucleus rings")
    return result


def plot_main(argv=None):
    parser = argparse.ArgumentParser(description="Render all cell/nucleus rings from existing segmentation")
    parser.add_argument("--input-dir",type=Path,required=True)
    parser.add_argument("--output-png",type=Path,required=True)
    scope = parser.add_mutually_exclusive_group(required=True)
    scope.add_argument("--bbox",nargs=4,type=float)
    scope.add_argument("--full",action="store_true")
    parser.add_argument("--color-by")
    parser.add_argument("--log-color",action="store_true")
    parser.add_argument("--no-nucleus",action="store_true")
    parser.add_argument("--dpi",type=int,default=160)
    args = parser.parse_args(argv)
    result = render_current(args.input_dir,args.output_png,bbox=args.bbox,color_by=args.color_by,
        log_color=args.log_color,overlay_nucleus=not args.no_nucleus,dpi=args.dpi)
    print(f"PASS segmentation-plot-current: {result['cell_rings']} cell rings; {result['nucleus_rings']} nucleus rings")
    return result


def main():
    import json
    from run_contract import Run
    parser = argparse.ArgumentParser(description="Portable complete-ring QC and PNG render")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bbox", nargs=4, type=float)
    parser.add_argument("--color-by")
    parser.add_argument("--no-nucleus", action="store_true")
    parser.add_argument("--upstream-json", type=Path)
    args = parser.parse_args()
    inputs = [args.input_dir/n for n in ["cells.parquet", "cell_boundaries.parquet", "nucleus_boundaries.parquet"]]
    if args.upstream_json:
        inputs.append(args.upstream_json)
    config = {"coordinate_unit":"um", "orientation":"top_left_y_down", "bbox":args.bbox,
              "color_by":args.color_by,"overlay_nucleus":not args.no_nucleus,
              "upstream":json.loads(args.upstream_json.read_text()) if args.upstream_json else []}
    run = Run(args.out_dir, config, inputs, __file__)
    if run.resumed:
        print("PASS complete-artifact resume")
        return
    # Validate every coordinate and ring before writing QC outputs.
    cells = read_cells(args.input_dir/"cells.parquet").set_index("cell_id")
    for name in ("cell", "nucleus"):
        flatten_rings(load_rings(args.input_dir/f"{name}_boundaries.parquet"), cells)
    summary = qc_current(args.input_dir, args.out_dir)
    render = render_current(args.input_dir, args.out_dir/"segmentation.png", bbox=args.bbox,
                            color_by=args.color_by, overlay_nucleus=not args.no_nucleus)
    render["orientation"] = "top_left_y_down"
    render["coordinate_unit"] = "um"
    write_json(args.out_dir/"segmentation.json", render)
    run.finish({"status":summary["status"], "scope":"existing_ring_integrity_and_png_decode",
                "biological_segmentation_accuracy":"NOT_VALIDATED"})
    print(summary["status"] + " portable segmentation inspection")

if __name__ == "__main__":
    main()
