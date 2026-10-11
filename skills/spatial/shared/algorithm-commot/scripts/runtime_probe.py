#!/usr/bin/env python3
"""Benchmark process-level parallelism for independent COMMOT tasks.

This probe uses deterministic synthetic AnnData objects and does not touch
project result files.  It compares serial execution with a spawn-based process
pool, verifies exact sparse-matrix hashes, and records worker process IDs.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=500)
    parser.add_argument("--lr-pairs", type=int, default=4)
    parser.add_argument("--tasks", type=int, default=4)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads-per-worker", type=int, default=1)
    parser.add_argument("--python-hash-seed", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path)
    return parser


def sparse_digest(matrix: sparse.spmatrix) -> str:
    csr = sparse.csr_matrix(matrix)
    digest = hashlib.sha256()
    digest.update(np.asarray(csr.shape, dtype=np.int64).tobytes())
    digest.update(csr.indptr.astype(np.int64, copy=False).tobytes())
    digest.update(csr.indices.astype(np.int64, copy=False).tobytes())
    digest.update(csr.data.astype(np.float64, copy=False).tobytes())
    return digest.hexdigest()


def run_task(payload: tuple[int, int, int, int]) -> dict[str, Any]:
    task_index, cells, lr_pairs, seed = payload
    from threadpoolctl import threadpool_info, threadpool_limits

    thread_limit = int(os.environ.get("COMMOT_PROBE_THREADS", "1"))
    with threadpool_limits(limits=thread_limit):
        import commot as ct

        rng = np.random.default_rng(seed + task_index)
        genes = [f"L{i}" for i in range(lr_pairs)] + [f"R{i}" for i in range(lr_pairs)]
        expression = rng.gamma(shape=1.5, scale=1.0, size=(cells, len(genes)))
        expression[rng.random(expression.shape) < 0.65] = 0.0
        # Guarantee nonzero mass for every ligand and receptor.
        expression[:lr_pairs, :lr_pairs] += np.eye(lr_pairs) * 2.0
        expression[lr_pairs : 2 * lr_pairs, lr_pairs:] += np.eye(lr_pairs) * 2.0
        side = int(np.ceil(np.sqrt(cells)))
        grid = np.column_stack(np.unravel_index(np.arange(cells), (side, side))).astype(float)
        grid += rng.normal(0.0, 0.02, size=grid.shape)
        adata = ad.AnnData(
            X=sparse.csr_matrix(expression),
            obs=pd.DataFrame(index=[f"task{task_index}_cell{i}" for i in range(cells)]),
            var=pd.DataFrame(index=genes),
        )
        adata.obsm["spatial"] = grid
        database = pd.DataFrame(
            {
                "ligand": [f"L{i}" for i in range(lr_pairs)],
                "receptor": [f"R{i}" for i in range(lr_pairs)],
                "pathway": [f"P{i % 2}" for i in range(lr_pairs)],
            }
        )
        started = time.perf_counter()
        ct.tl.spatial_communication(
            adata,
            database_name="probe",
            df_ligrec=database,
            dis_thr=3.0,
            heteromeric=True,
            pathway_sum=True,
            smooth=False,
        )
        elapsed = time.perf_counter() - started
        keys = sorted(key for key in adata.obsp if key.startswith("commot-probe-"))
        digests = {key: sparse_digest(adata.obsp[key]) for key in keys}
        return {
            "task": task_index,
            "pid": os.getpid(),
            "elapsed_seconds": elapsed,
            "threads": (
                sum(1 for _ in Path("/proc/self/task").glob("*"))
                if Path("/proc/self/task").is_dir()
                else None
            ),
            "threadpools": [
                {
                    "internal_api": item.get("internal_api"),
                    "num_threads": item.get("num_threads"),
                }
                for item in threadpool_info()
            ],
            "keys": keys,
            "digests": digests,
        }


def pool_initializer(threads_per_worker: int) -> None:
    value = str(int(threads_per_worker))
    os.environ["COMMOT_PROBE_THREADS"] = value
    os.environ["OMP_NUM_THREADS"] = value
    os.environ["OPENBLAS_NUM_THREADS"] = value
    os.environ["MKL_NUM_THREADS"] = value


def run_mode(
    payloads: list[tuple[int, int, int, int]], workers: int, threads_per_worker: int
) -> tuple[float, list[dict[str, Any]]]:
    started = time.perf_counter()
    # Use a spawned pool for both modes so interpreter hash initialization,
    # imports, and thread controls are comparable.
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
        initializer=pool_initializer,
        initargs=(threads_per_worker,),
    ) as executor:
        results = list(executor.map(run_task, payloads))
    return time.perf_counter() - started, results


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if min(args.cells, args.lr_pairs, args.tasks, args.workers, args.threads_per_worker) < 1:
        raise ValueError("All integer parameters must be positive")
    payloads = [
        (task, int(args.cells), int(args.lr_pairs), int(args.seed))
        for task in range(int(args.tasks))
    ]
    # COMMOT 0.0.3 constructs ligand/receptor orders from Python sets.  This
    # must be set before spawn creates either the serial or parallel workers.
    os.environ["PYTHONHASHSEED"] = str(int(args.python_hash_seed))
    thread_value = str(int(args.threads_per_worker))
    # Spawn imports NumPy/Scipy before the worker initializer runs, so these
    # controls also have to be present in the parent environment pre-spawn.
    os.environ["OMP_NUM_THREADS"] = thread_value
    os.environ["OPENBLAS_NUM_THREADS"] = thread_value
    os.environ["MKL_NUM_THREADS"] = thread_value
    serial_seconds, serial = run_mode(payloads, 1, int(args.threads_per_worker))
    parallel_seconds, parallel = run_mode(
        payloads, min(int(args.workers), len(payloads)), int(args.threads_per_worker)
    )
    serial_by_task = {int(row["task"]): row["digests"] for row in serial}
    parallel_by_task = {int(row["task"]): row["digests"] for row in parallel}
    exact_match = serial_by_task == parallel_by_task
    result = {
        "status": "PASS" if exact_match else "FAIL",
        "scope": "synthetic_independent_commot_tasks",
        "python": sys.executable,
        "platform": platform.platform(),
        "parameters": {
            "cells": int(args.cells),
            "lr_pairs": int(args.lr_pairs),
            "tasks": int(args.tasks),
            "workers": min(int(args.workers), len(payloads)),
            "threads_per_worker": int(args.threads_per_worker),
            "seed": int(args.seed),
            "python_hash_seed": int(args.python_hash_seed),
        },
        "serial_wall_seconds": serial_seconds,
        "parallel_wall_seconds": parallel_seconds,
        "speedup": serial_seconds / parallel_seconds if parallel_seconds else None,
        "exact_sparse_output_match": exact_match,
        "serial_worker_pids": sorted({int(row["pid"]) for row in serial}),
        "parallel_worker_pids": sorted({int(row["pid"]) for row in parallel}),
        "serial_tasks": serial,
        "parallel_tasks": parallel,
        "limitations": [
            "Synthetic timing does not predict full-dataset wall time or peak memory.",
            "The probe validates process-level task parallelism, not GPU execution.",
        ],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(rendered + "\n", encoding="utf-8")
        temporary.replace(args.output)
    print(rendered)
    return 0 if exact_match else 2


if __name__ == "__main__":
    raise SystemExit(main())
