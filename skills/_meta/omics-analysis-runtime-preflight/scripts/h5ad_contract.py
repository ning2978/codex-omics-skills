#!/usr/bin/env python3
"""Validate a representative AnnData contract using a backed, bounded read."""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def add_check(
    checks: list[dict[str, Any]],
    name: str,
    passed: bool,
    evidence: Any,
    *,
    required: bool = True,
) -> None:
    checks.append(
        {
            "name": name,
            "status": "PASS" if passed else ("FAIL" if required else "WARN"),
            "required": required,
            "evidence": evidence,
        }
    )


def materialize_sample(matrix: Any, n_rows: int, n_cols: int) -> tuple[Any, str]:
    from scipy import sparse

    sampled = matrix[:n_rows, :n_cols]
    if hasattr(sampled, "to_memory"):
        sampled = sampled.to_memory()
    if sparse.issparse(sampled):
        return sampled.data, f"{type(sampled).__name__}.data"
    import numpy as np

    return np.asarray(sampled).reshape(-1), type(sampled).__name__


def resolve_count_matrix(adata: Any, source: str) -> tuple[Any, str]:
    if source == "auto":
        if "counts" in adata.layers:
            return adata.layers["counts"], "layer:counts"
        if adata.raw is not None:
            return adata.raw.X, "raw"
        return adata.X, "X"
    if source == "X":
        return adata.X, "X"
    if source == "raw":
        if adata.raw is None:
            raise KeyError("adata.raw is absent")
        return adata.raw.X, "raw"
    if source.startswith("layer:"):
        key = source.split(":", 1)[1]
        if key not in adata.layers:
            raise KeyError(f"layers[{key!r}] is absent")
        return adata.layers[key], f"layer:{key}"
    raise ValueError(f"Unsupported count source: {source}")


def numeric_summary(values: Any) -> dict[str, Any]:
    import numpy as np

    values = np.asarray(values)
    finite = np.isfinite(values)
    finite_values = values[finite]
    if finite_values.size:
        minimum = float(finite_values.min())
        maximum = float(finite_values.max())
        mean = float(finite_values.mean())
        integer_like = bool(np.allclose(finite_values, np.rint(finite_values), atol=1e-6, rtol=0.0))
        nonnegative = minimum >= 0.0
    else:
        minimum = maximum = mean = None
        integer_like = False
        nonnegative = False
    return {
        "sampled_values": int(values.size),
        "finite_values": int(finite.sum()),
        "all_finite": bool(finite.all()),
        "nonnegative": nonnegative,
        "integer_like": integer_like,
        "minimum": minimum,
        "maximum": maximum,
        "mean": mean,
    }


def write_report(report: dict[str, Any], output: Path | None) -> None:
    rendered = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--kind", choices=("scrna", "xenium", "visium", "generic"), required=True)
    parser.add_argument(
        "--count-source",
        default="auto",
        help="One of auto, X, raw, or layer:<name>.",
    )
    parser.add_argument("--obs-key", action="append", default=[], help="Required obs key; repeatable.")
    parser.add_argument("--obsm-key", action="append", default=[], help="Additional required obsm key; repeatable.")
    parser.add_argument("--sample-rows", type=int, default=256)
    parser.add_argument("--sample-cols", type=int, default=512)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    checks: list[dict[str, Any]] = []
    report: dict[str, Any] = {
        "schema_version": 1,
        "check": "h5ad_contract",
        "generated_at_utc": utc_now(),
        "input": str(args.input.resolve()),
        "kind": args.kind,
        "requested_count_source": args.count_source,
        "python": {"executable": str(Path(sys.executable).resolve()), "version": platform.python_version()},
        "sampling_limit": {"rows": args.sample_rows, "columns": args.sample_cols},
        "sampling_disclaimer": "Matrix integrity checks cover a bounded sample, not every stored value.",
        "checks": checks,
    }

    adata = None
    try:
        import anndata as ad
        import numpy as np

        add_check(checks, "input_exists", args.input.is_file(), {"path": str(args.input), "size_bytes": args.input.stat().st_size if args.input.is_file() else None})
        if not args.input.is_file():
            raise FileNotFoundError(args.input)

        adata = ad.read_h5ad(args.input, backed="r")
        add_check(checks, "backed_reopen", True, {"isbacked": bool(adata.isbacked)})
        report["shape"] = [int(adata.n_obs), int(adata.n_vars)]
        report["keys"] = {
            "obs": list(map(str, adata.obs.columns)),
            "var": list(map(str, adata.var.columns)),
            "layers": list(map(str, adata.layers.keys())),
            "obsm": list(map(str, adata.obsm.keys())),
            "uns": list(map(str, adata.uns.keys())),
            "raw_present": adata.raw is not None,
        }
        add_check(checks, "nonzero_shape", adata.n_obs > 0 and adata.n_vars > 0, report["shape"])
        add_check(checks, "unique_obs_names", bool(adata.obs_names.is_unique), {"n_obs": int(adata.n_obs)})
        add_check(checks, "unique_var_names", bool(adata.var_names.is_unique), {"n_vars": int(adata.n_vars)})

        for key in args.obs_key:
            present = key in adata.obs
            evidence: dict[str, Any] = {"key": key, "present": present}
            if present:
                values = adata.obs[key]
                evidence.update(
                    {
                        "missing": int(values.isna().sum()),
                        "unique": int(values.nunique(dropna=True)),
                        "examples": list(map(str, values.dropna().astype(str).unique()[:8])),
                    }
                )
            add_check(checks, f"obs_key:{key}", present and evidence.get("missing", 1) == 0, evidence)

        required_obsm = list(dict.fromkeys((*args.obsm_key, *(("spatial",) if args.kind in {"xenium", "visium"} else ()))))
        for key in required_obsm:
            present = key in adata.obsm
            evidence = {"key": key, "present": present}
            valid = False
            if present:
                matrix = adata.obsm[key]
                rows = min(args.sample_rows, adata.n_obs)
                sample = matrix[:rows]
                if hasattr(sample, "to_memory"):
                    sample = sample.to_memory()
                sample = np.asarray(sample)
                evidence.update(
                    {
                        "shape": list(map(int, matrix.shape)),
                        "sample_rows": int(rows),
                        "sample_all_finite": bool(np.isfinite(sample).all()),
                    }
                )
                valid = matrix.shape[0] == adata.n_obs and matrix.shape[1] >= 2 and evidence["sample_all_finite"]
            add_check(checks, f"obsm_key:{key}", valid, evidence)

        try:
            count_matrix, resolved_source = resolve_count_matrix(adata, args.count_source)
            report["resolved_count_source"] = resolved_source
            n_rows = min(max(1, args.sample_rows), adata.n_obs)
            n_cols = min(max(1, args.sample_cols), count_matrix.shape[1])
            values, representation = materialize_sample(count_matrix, n_rows, n_cols)
            summary = numeric_summary(values)
            summary.update(
                {
                    "resolved_source": resolved_source,
                    "matrix_shape": list(map(int, count_matrix.shape)),
                    "sample_shape": [int(n_rows), int(n_cols)],
                    "sample_representation": representation,
                }
            )
            report["count_sample"] = summary
            count_ok = summary["all_finite"] and summary["nonnegative"] and summary["integer_like"]
            add_check(checks, "sampled_raw_count_contract", count_ok, summary)
            add_check(
                checks,
                "sample_contains_nonzero_values",
                bool(summary["sampled_values"] and summary["maximum"] is not None and summary["maximum"] > 0),
                {"maximum": summary["maximum"], "sampled_values": summary["sampled_values"]},
                required=False,
            )
        except Exception as exc:
            add_check(
                checks,
                "sampled_raw_count_contract",
                False,
                {"error": f"{type(exc).__name__}: {exc}"},
            )

    except Exception as exc:
        report["exception"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        if not any(check["name"] == "backed_reopen" for check in checks):
            add_check(checks, "backed_reopen", False, report["exception"])
    finally:
        if adata is not None and getattr(adata, "file", None) is not None:
            try:
                adata.file.close()
            except Exception:
                pass

    required_failures = [check["name"] for check in checks if check["required"] and check["status"] == "FAIL"]
    warnings = [check["name"] for check in checks if check["status"] == "WARN"]
    report["required_failures"] = required_failures
    report["warnings"] = warnings
    report["status"] = "FAIL" if required_failures else ("WARN" if warnings else "PASS")
    write_report(report, args.output)
    return 1 if required_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
