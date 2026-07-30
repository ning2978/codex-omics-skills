#!/usr/bin/env python3
"""Run a real, minimal scVI optimization loop on a requested accelerator."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def history_summary(history: Any) -> dict[str, Any]:
    import numpy as np

    result: dict[str, Any] = {}
    finite_across_history = True
    numeric_series = 0
    for key, value in dict(history).items():
        raw_values = value.to_numpy().reshape(-1) if hasattr(value, "to_numpy") else np.asarray(value, dtype=object).reshape(-1)
        converted: list[float] = []
        for item in raw_values:
            try:
                if hasattr(item, "detach"):
                    item = item.detach().cpu()
                if hasattr(item, "item"):
                    item = item.item()
                converted.append(float(item))
            except (TypeError, ValueError, RuntimeError):
                continue
        values = np.asarray(converted, dtype=float)
        if values.size == 0:
            continue
        numeric_series += 1
        finite = bool(np.isfinite(values).all())
        finite_across_history = finite_across_history and finite
        result[str(key)] = {
            "values": int(values.size),
            "first": float(values[0]),
            "last": float(values[-1]),
            "all_finite": finite,
        }
    return {
        "series": result,
        "numeric_series": numeric_series,
        "all_finite": finite_across_history and numeric_series > 0,
    }


def write_report(report: dict[str, Any], output: Path | None) -> None:
    rendered = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accelerator", choices=("gpu", "cpu"), default="gpu")
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--cells", type=int, default=256)
    parser.add_argument("--genes", type=int, default=96)
    parser.add_argument("--batches", type=int, default=2)
    parser.add_argument("--latent", type=int, default=5)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    started = time.monotonic()
    report: dict[str, Any] = {
        "schema_version": 1,
        "check": "scvi_training_smoke",
        "status": "FAIL",
        "generated_at_utc": utc_now(),
        "requested_accelerator": args.accelerator,
        "parameters": {
            "epochs": args.epochs,
            "cells": args.cells,
            "genes": args.genes,
            "batches": args.batches,
            "latent": args.latent,
            "seed": args.seed,
        },
        "python": {"executable": str(Path(sys.executable).resolve()), "version": platform.python_version()},
        "packages": {
            name: version(name)
            for name in ("scvi-tools", "torch", "lightning", "pytorch-lightning", "anndata", "numpy", "scipy")
        },
        "scope": "Synthetic integer-count plumbing test; not a full-data memory, convergence, or biology test.",
    }

    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scvi
        import torch
        from scipy import sparse

        if args.epochs < 1 or args.cells < 32 or args.genes < 16 or args.batches < 1:
            raise ValueError("epochs>=1, cells>=32, genes>=16, and batches>=1 are required")

        cuda_available = bool(torch.cuda.is_available())
        device_count = int(torch.cuda.device_count())
        report["framework_cuda"] = {
            "torch_version": torch.__version__,
            "torch_cuda_build": torch.version.cuda,
            "cuda_available": cuda_available,
            "device_count": device_count,
        }
        if args.accelerator == "gpu" and (not cuda_available or device_count < 1):
            report["failure_stage"] = "cuda_precondition"
            raise RuntimeError("GPU requested but torch reports no CUDA device")

        rng = np.random.default_rng(args.seed)
        batch_codes = np.arange(args.cells) % args.batches
        base_rate = rng.gamma(shape=2.0, scale=0.8, size=(args.cells, args.genes))
        batch_scale = 1.0 + 0.25 * batch_codes[:, None]
        counts = rng.poisson(base_rate * batch_scale).astype(np.int32)
        counts[:, : min(4, args.genes)] += (batch_codes[:, None] + 1).astype(np.int32)
        obs = pd.DataFrame(
            {"batch": pd.Categorical([f"batch_{value}" for value in batch_codes])},
            index=[f"cell_{index}" for index in range(args.cells)],
        )
        var = pd.DataFrame(index=[f"gene_{index}" for index in range(args.genes)])
        adata = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()
        input_sum = int(adata.layers["counts"].sum())

        scvi.settings.seed = args.seed
        scvi.model.SCVI.setup_anndata(adata, batch_key="batch", layer="counts")
        model = scvi.model.SCVI(
            adata,
            n_layers=1,
            n_latent=args.latent,
            gene_likelihood="nb",
        )
        if args.accelerator == "gpu":
            torch.cuda.reset_peak_memory_stats(0)

        training_started = time.monotonic()
        model.train(
            max_epochs=args.epochs,
            accelerator=args.accelerator,
            devices=1,
            train_size=0.9,
            batch_size=min(64, args.cells),
            early_stopping=False,
            enable_progress_bar=False,
            enable_model_summary=False,
            enable_checkpointing=False,
        )
        training_seconds = time.monotonic() - training_started
        trainer = model.trainer
        root_device = str(trainer.strategy.root_device)
        latent = np.asarray(model.get_latent_representation())
        history = history_summary(model.history)
        output_sum = int(adata.layers["counts"].sum())
        latent_ok = latent.shape == (args.cells, args.latent) and bool(np.isfinite(latent).all())
        count_preserved = input_sum == output_sum
        expected_device = root_device.startswith("cuda") if args.accelerator == "gpu" else root_device.startswith("cpu")

        peak_bytes = 0
        if args.accelerator == "gpu":
            torch.cuda.synchronize(0)
            peak_bytes = int(torch.cuda.max_memory_allocated(0))
        gpu_allocation_ok = peak_bytes > 0 if args.accelerator == "gpu" else True

        report["training"] = {
            "completed": True,
            "elapsed_seconds": round(training_seconds, 4),
            "trainer_root_device": root_device,
            "trainer_accelerator_class": type(trainer.accelerator).__name__,
            "peak_cuda_allocated_bytes": peak_bytes,
            "history": history,
            "latent_shape": list(map(int, latent.shape)),
            "latent_all_finite": bool(np.isfinite(latent).all()),
            "input_count_sum": input_sum,
            "output_count_sum": output_sum,
            "counts_preserved": count_preserved,
        }
        report["gates"] = {
            "training_completed": True,
            "observed_requested_device": expected_device,
            "finite_training_history": history["all_finite"],
            "finite_expected_latent_shape": latent_ok,
            "raw_counts_preserved": count_preserved,
            "nonzero_gpu_allocation_if_requested": gpu_allocation_ok,
        }
        report["status"] = "PASS" if all(report["gates"].values()) else "FAIL"
        if report["status"] == "FAIL":
            report["failure_stage"] = "post_training_validation"
    except Exception as exc:
        report.setdefault("failure_stage", "import_setup_or_training")
        report["exception"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }

    report["total_elapsed_seconds"] = round(time.monotonic() - started, 4)
    write_report(report, args.output)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
