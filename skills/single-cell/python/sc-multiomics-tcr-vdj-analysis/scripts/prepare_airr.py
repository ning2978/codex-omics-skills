"""Scirpy >=0.13 AIRR layout, validated with 0.22.4; retain all RNA cells."""
from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd


def key(library, barcode):
    return json.dumps([str(library), str(barcode)], ensure_ascii=False, separators=(",", ":"))


def attach_airr(rna, library_paths, *, library_col="capture_library", barcode_col="barcode"):
    import scirpy as ir
    from mudata import MuData
    required = [library_col, barcode_col, "sample_id"]
    if any(col not in rna.obs or rna.obs[col].isna().any() for col in required):
        raise ValueError("Complete capture_library, original barcode and sample_id required")
    rna = rna.copy()
    original_keys = rna.obs_names.copy()
    join_keys = pd.Index([key(a, b) for a, b in zip(rna.obs[library_col], rna.obs[barcode_col], strict=True)])
    if not original_keys.is_unique or not join_keys.is_unique:
        raise ValueError("Duplicate RNA cell or capture-barcode key")
    frames, report = [], []
    for library, path in library_paths.items():
        if Path(path).suffix == ".csv":
            contigs = pd.read_csv(path, dtype={"barcode": str, "contig_id": str})
            if contigs[["barcode", "contig_id"]].isna().any().any() or contigs.duplicated(["barcode", "contig_id"]).any():
                raise ValueError("Missing/duplicated contig key; distinct chains share cells, not contig IDs")
        # read_10x_vdj stores every contig for one cell in an AIRR list, not extra obs rows.
        airr = ir.io.read_10x_vdj(Path(path))
        if not airr.obs_names.is_unique:
            raise ValueError("Reader produced duplicate cell IDs; inspect contig records")
        airr.obs_names = pd.Index([key(library, b) for b in airr.obs_names])
        wanted = join_keys[rna.obs[library_col].astype(str).to_numpy() == str(library)]
        shared = wanted.intersection(airr.obs_names, sort=False)
        report.append({"capture_library": str(library), "rna_cells": len(wanted),
                       "vdj_cells": airr.n_obs, "matched": len(shared),
                       "rna_without_vdj": len(wanted.difference(airr.obs_names)),
                       "vdj_without_rna": len(airr.obs_names.difference(wanted))})
        if len(shared):
            frames.append(airr[shared].copy())
    if not frames:
        raise ValueError("No capture-aware RNA/VDJ match")
    airr = ad.concat(frames, axis=0, index_unique=None, merge="same")
    if not airr.obs_names.is_unique:
        raise ValueError("VDJ library mapping produced duplicate keys")
    positions = join_keys.get_indexer(airr.obs_names)
    if (positions < 0).any():
        raise ValueError("Unexpected VDJ-only cell after filtering")
    airr.obs_names = original_keys[positions]
    for col in required:
        airr.obs[col] = rna.obs.loc[airr.obs_names, col].to_numpy()
    mdata = MuData({"rna": rna, "airr": airr})
    if mdata.mod["rna"].n_obs != rna.n_obs or set(mdata.obs_names) != set(original_keys):
        raise ValueError("RNA cell set changed during modality assembly")
    return mdata, pd.DataFrame(report)


def export_airr_metadata(mdata, *, columns=("clone_id", "clonal_expansion", "chain_pairing")):
    rna, airr = mdata.mod["rna"], mdata.mod["airr"]
    if not airr.obs_names.is_unique:
        raise ValueError("Duplicate AIRR cell keys")
    for col in columns:
        if col not in airr.obs:
            raise KeyError(col)
        # Reindex only, never duplicate rows or drop RNA cells with missing receptors.
        rna.obs[col] = airr.obs[col].reindex(rna.obs_names)
    return rna
