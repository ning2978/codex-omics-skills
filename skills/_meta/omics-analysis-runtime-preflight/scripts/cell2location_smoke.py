#!/usr/bin/env python3
"""Train both cell2location model stages on synthetic counts and one requested device."""

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


def trainer_summary(model: Any) -> dict[str, Any]:
    trainer = model.trainer
    return {
        "root_device": str(trainer.strategy.root_device),
        "accelerator_class": type(trainer.accelerator).__name__,
        "history": history_summary(model.history),
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
    parser.add_argument("--reference-cells", type=int, default=180)
    parser.add_argument("--spots", type=int, default=48)
    parser.add_argument("--genes", type=int, default=48)
    parser.add_argument("--cell-types", type=int, default=3)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    started = time.monotonic()
    report: dict[str, Any] = {
        "schema_version": 1,
        "check": "cell2location_training_smoke",
        "status": "FAIL",
        "generated_at_utc": utc_now(),
        "requested_accelerator": args.accelerator,
        "parameters": {
            "epochs_per_stage": args.epochs,
            "reference_cells": args.reference_cells,
            "spots": args.spots,
            "genes": args.genes,
            "cell_types": args.cell_types,
            "N_cells_per_location": 5,
            "detection_alpha": 20,
            "seed": args.seed,
        },
        "python": {"executable": str(Path(sys.executable).resolve()), "version": platform.python_version()},
        "packages": {
            name: version(name)
            for name in (
                "cell2location",
                "scvi-tools",
                "torch",
                "lightning",
                "pytorch-lightning",
                "pyro-ppl",
                "anndata",
                "numpy",
                "scipy",
            )
        },
        "scope": (
            "Synthetic integer-count test of RegressionModel and spatial Cell2location training. "
            "It does not validate final parameter choices, signature export, posterior export, scale, or biology."
        ),
    }

    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scvi
        import torch
        from cell2location.models import Cell2location, RegressionModel
        from scipy import sparse

        if (
            args.epochs < 1
            or args.reference_cells < 60
            or args.spots < 16
            or args.genes < 24
            or args.cell_types < 2
        ):
            raise ValueError("epochs>=1, reference-cells>=60, spots>=16, genes>=24, and cell-types>=2 are required")

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
        torch.manual_seed(args.seed)
        scvi.settings.seed = args.seed
        genes = [f"gene_{index}" for index in range(args.genes)]
        factor_names = [f"cell_type_{index}" for index in range(args.cell_types)]

        signatures = rng.gamma(shape=2.0, scale=0.7, size=(args.genes, args.cell_types)) + 0.15
        marker_width = max(2, args.genes // (args.cell_types * 4))
        for factor in range(args.cell_types):
            start = factor * marker_width
            signatures[start : start + marker_width, factor] += 3.0
        signature_frame = pd.DataFrame(signatures, index=genes, columns=factor_names)

        labels = np.arange(args.reference_cells) % args.cell_types
        batches = np.arange(args.reference_cells) % 2
        ref_rate = signatures[:, labels].T * (1.0 + 0.15 * batches[:, None])
        ref_counts = rng.poisson(ref_rate).astype(np.int32)
        ref = ad.AnnData(
            X=sparse.csr_matrix(ref_counts),
            obs=pd.DataFrame(
                {
                    "cell_type": pd.Categorical([factor_names[value] for value in labels]),
                    "batch": pd.Categorical([f"batch_{value}" for value in batches]),
                },
                index=[f"ref_cell_{index}" for index in range(args.reference_cells)],
            ),
            var=pd.DataFrame(index=genes),
        )
        ref.layers["counts"] = ref.X.copy()

        proportions = rng.dirichlet(np.ones(args.cell_types), size=args.spots)
        expected_cells = rng.integers(3, 8, size=args.spots)
        spatial_rate = (proportions * expected_cells[:, None]) @ signatures.T
        spatial_rate *= rng.uniform(0.7, 1.3, size=(args.spots, 1))
        spatial_counts = rng.poisson(spatial_rate).astype(np.int32)
        spatial = ad.AnnData(
            X=sparse.csr_matrix(spatial_counts),
            obs=pd.DataFrame(
                {"sample": pd.Categorical(["sample_0"] * args.spots)},
                index=[f"spot_{index}" for index in range(args.spots)],
            ),
            var=pd.DataFrame(index=genes),
        )
        spatial.layers["counts"] = spatial.X.copy()
        spatial.obsm["spatial"] = rng.uniform(0, 1000, size=(args.spots, 2)).astype(np.float32)

        if args.accelerator == "gpu":
            torch.cuda.reset_peak_memory_stats(0)

        RegressionModel.setup_anndata(ref, layer="counts", batch_key="batch", labels_key="cell_type")
        reference_model = RegressionModel(ref)
        pyro_device = "0" if args.accelerator == "gpu" else "auto"
        reference_started = time.monotonic()
        reference_model.train(
            max_epochs=args.epochs,
            batch_size=min(64, args.reference_cells),
            train_size=1,
            accelerator=args.accelerator,
            device=pyro_device,
            enable_progress_bar=False,
            enable_model_summary=False,
            enable_checkpointing=False,
        )
        reference = trainer_summary(reference_model)
        reference["elapsed_seconds"] = round(time.monotonic() - reference_started, 4)

        Cell2location.setup_anndata(spatial, layer="counts", batch_key="sample")
        spatial_model = Cell2location(
            spatial,
            cell_state_df=signature_frame,
            N_cells_per_location=5,
            detection_alpha=20,
        )
        spatial_started = time.monotonic()
        spatial_model.train(
            max_epochs=args.epochs,
            batch_size=None,
            train_size=1,
            accelerator=args.accelerator,
            device=pyro_device,
            enable_progress_bar=False,
            enable_model_summary=False,
            enable_checkpointing=False,
        )
        spatial_training = trainer_summary(spatial_model)
        spatial_training["elapsed_seconds"] = round(time.monotonic() - spatial_started, 4)

        peak_bytes = 0
        if args.accelerator == "gpu":
            torch.cuda.synchronize(0)
            peak_bytes = int(torch.cuda.max_memory_allocated(0))
        gpu_allocation_ok = peak_bytes > 0 if args.accelerator == "gpu" else True
        reference_device_ok = (
            reference["root_device"].startswith("cuda")
            if args.accelerator == "gpu"
            else reference["root_device"].startswith("cpu")
        )
        spatial_device_ok = (
            spatial_training["root_device"].startswith("cuda")
            if args.accelerator == "gpu"
            else spatial_training["root_device"].startswith("cpu")
        )

        report["training"] = {
            "reference_model": reference,
            "spatial_model": spatial_training,
            "peak_cuda_allocated_bytes": peak_bytes,
            "reference_count_sum": int(ref.layers["counts"].sum()),
            "spatial_count_sum": int(spatial.layers["counts"].sum()),
            "signature_shape": list(map(int, signature_frame.shape)),
        }
        report["gates"] = {
            "reference_training_completed": True,
            "spatial_training_completed": True,
            "reference_observed_requested_device": reference_device_ok,
            "spatial_observed_requested_device": spatial_device_ok,
            "reference_history_finite": reference["history"]["all_finite"],
            "spatial_history_finite": spatial_training["history"]["all_finite"],
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
