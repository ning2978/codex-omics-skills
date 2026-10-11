"""Ring-preserving Xenium parquet readers, separate from historical polygon helpers."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

COLUMNS = ["cell_id", "label_id", "vertex_x", "vertex_y"]


def _parquet(path):
    pf = pq.ParquetFile(path)
    missing = set(COLUMNS) - set(pf.schema_arrow.names)
    if missing:
        raise ValueError(f"Boundary parquet missing columns {sorted(missing)}")
    return pf


def _validate(frame):
    if frame[["cell_id", "label_id"]].isna().any().any():
        raise ValueError("Missing cell_id/label_id ring key")
    if not np.isfinite(frame[["vertex_x", "vertex_y"]].to_numpy(float)).all():
        raise ValueError("Non-finite ring vertices")


def load_rings(path: Path, cell_ids=None):
    """cell_id -> label_id -> complete vertices; supports interleaved/unsorted ring rows.

    Vertex order within each ring is the source row order. Vertices are never
    cropped or sorted by coordinate. Crop renderers clip complete polygons.
    """
    pf = _parquet(path)
    selected = set(cell_ids) if cell_ids is not None else None
    parts = []
    for i in range(pf.num_row_groups):
        frame = pf.read_row_group(i, columns=COLUMNS).to_pandas()
        _validate(frame)
        if selected is not None:
            frame = frame[frame.cell_id.isin(selected)]
        if not frame.empty:
            parts.append(frame)
    if not parts:
        return {}
    frames = pd.concat(parts, ignore_index=True)
    rings = {}
    for (cid, label), group in frames.groupby(["cell_id", "label_id"], sort=False):
        rings.setdefault(cid, {})[label] = group[["vertex_x", "vertex_y"]].to_numpy(float)
    return rings


def iter_rings(path: Path, chunk_row_groups=1):
    """Yield (cell_id, label_id, vertices); fail when a completed ring reappears.

    Rings may span any number of row groups. A cell may reappear with a
    different label; the continuity contract applies to the complete ring key.
    """
    if chunk_row_groups < 1:
        raise ValueError("chunk_row_groups must be positive")
    pf = _parquet(path)
    finished, current, buffers = set(), None, []
    for start in range(0, pf.num_row_groups, chunk_row_groups):
        frames = [pf.read_row_group(i, columns=COLUMNS).to_pandas()
                  for i in range(start, min(start + chunk_row_groups, pf.num_row_groups))]
        frame = pd.concat(frames, ignore_index=True)
        _validate(frame)
        keys = list(zip(frame.cell_id, frame.label_id))
        begin = 0
        while begin < len(frame):
            key = keys[begin]
            end = begin + 1
            while end < len(frame) and keys[end] == key:
                end += 1
            if key != current:
                if key in finished:
                    raise ValueError(f"Non-contiguous ring {key!r} reappeared after completion; use load_rings for unsorted data")
                if current is not None:
                    finished.add(current)
                    yield (*current, np.concatenate(buffers))
                current, buffers = key, []
            buffers.append(frame.iloc[begin:end][["vertex_x", "vertex_y"]].to_numpy(float))
            begin = end
    if current is not None:
        yield (*current, np.concatenate(buffers))


def ring_vertex_counts(path: Path):
    """Sum ring sizes across all row groups, even for unsorted rows, without loading vertices."""
    pf = _parquet(path)
    counts = defaultdict(int)
    for i in range(pf.num_row_groups):
        frame = pf.read_row_group(i, columns=["cell_id", "label_id"]).to_pandas()
        if frame.isna().any().any():
            raise ValueError("Missing cell_id/label_id ring key")
        for key, n in frame.groupby(["cell_id", "label_id"], sort=False).size().items():
            counts[key] += int(n)
    return pd.DataFrame([{"cell_id":k[0], "label_id":k[1], "n_vertices":v}
                         for k,v in counts.items()], columns=["cell_id","label_id","n_vertices"])
