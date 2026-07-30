#!/usr/bin/env python3
"""Probe one exact omics Python environment and the host GPU stack."""

from __future__ import annotations

import argparse
import glob
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


DEFAULT_PACKAGES = (
    "anndata",
    "scanpy",
    "scvi-tools",
    "cell2location",
    "cellcharter",
    "scimap",
    "squidpy",
    "torch",
    "lightning",
    "pytorch-lightning",
    "numpy",
    "scipy",
    "pandas",
    "pyarrow",
    "harmonypy",
    "bbknn",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def decode_stream(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def run_command(command: Sequence[str], timeout: float) -> dict[str, Any]:
    started = time.monotonic()
    result: dict[str, Any] = {
        "command": list(command),
        "timeout_seconds": timeout,
    }
    try:
        completed = subprocess.run(
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
            check=False,
        )
        result.update(
            {
                "returncode": completed.returncode,
                "stdout": completed.stdout.strip(),
                "stderr": completed.stderr.strip(),
                "timed_out": False,
            }
        )
    except FileNotFoundError as exc:
        result.update(
            {
                "returncode": None,
                "stdout": "",
                "stderr": str(exc),
                "timed_out": False,
                "not_found": True,
            }
        )
    except subprocess.TimeoutExpired as exc:
        result.update(
            {
                "returncode": None,
                "stdout": decode_stream(exc.stdout).strip(),
                "stderr": decode_stream(exc.stderr).strip(),
                "timed_out": True,
            }
        )
    result["elapsed_seconds"] = round(time.monotonic() - started, 4)
    return result


def read_text(path: Path) -> dict[str, Any]:
    try:
        return {"exists": True, "text": path.read_text(encoding="utf-8", errors="replace").strip()}
    except OSError as exc:
        return {"exists": path.exists(), "error": f"{type(exc).__name__}: {exc}"}


def package_versions(extra: Sequence[str]) -> dict[str, str | None]:
    names = list(dict.fromkeys((*DEFAULT_PACKAGES, *extra)))
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def parse_child_json(probe: dict[str, Any]) -> dict[str, Any] | None:
    if probe.get("returncode") != 0 or not probe.get("stdout"):
        return None
    lines = str(probe["stdout"]).splitlines()
    for line in reversed(lines):
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    return None


def import_probe(module: str, timeout: float) -> dict[str, Any]:
    code = f"""
import importlib
import json
import time
started = time.monotonic()
try:
    mod = importlib.import_module({module!r})
    out = {{
        "module": {module!r},
        "ok": True,
        "version": getattr(mod, "__version__", None),
        "elapsed_seconds": round(time.monotonic() - started, 4),
    }}
except Exception as exc:
    out = {{
        "module": {module!r},
        "ok": False,
        "error": type(exc).__name__ + ": " + str(exc),
        "elapsed_seconds": round(time.monotonic() - started, 4),
    }}
print(json.dumps(out, ensure_ascii=False), flush=True)
raise SystemExit(0 if out["ok"] else 1)
"""
    probe = run_command([sys.executable, "-u", "-c", code], timeout)
    probe["result"] = parse_child_json(probe)
    return probe


def torch_probe(timeout: float) -> dict[str, Any]:
    code = r'''
import json
import time
out = {"ok": False}
started = time.monotonic()
try:
    import torch
    out.update({
        "ok": True,
        "torch_version": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "cudnn_version": torch.backends.cudnn.version(),
    })
    available = bool(torch.cuda.is_available())
    count = int(torch.cuda.device_count())
    out.update({"cuda_available": available, "device_count": count})
    devices = []
    for index in range(count):
        props = torch.cuda.get_device_properties(index)
        devices.append({
            "index": index,
            "name": torch.cuda.get_device_name(index),
            "capability": list(torch.cuda.get_device_capability(index)),
            "total_memory_bytes": int(props.total_memory),
        })
    out["devices"] = devices
    out["cuda_tensor_probe"] = False
    if available and count > 0:
        torch.cuda.reset_peak_memory_stats(0)
        x = torch.arange(16384, dtype=torch.float32, device="cuda:0").reshape(128, 128)
        y = x @ x.T
        torch.cuda.synchronize(0)
        out["cuda_tensor_probe"] = bool(torch.isfinite(y).all().item())
        out["peak_allocated_bytes"] = int(torch.cuda.max_memory_allocated(0))
except Exception as exc:
    out.update({"ok": False, "error": type(exc).__name__ + ": " + str(exc)})
out["elapsed_seconds"] = round(time.monotonic() - started, 4)
print(json.dumps(out, ensure_ascii=False), flush=True)
'''
    probe = run_command([sys.executable, "-u", "-c", code], timeout)
    probe["result"] = parse_child_json(probe)
    return probe


def memory_info() -> dict[str, int]:
    values: dict[str, int] = {}
    path = Path("/proc/meminfo")
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^(MemTotal|MemAvailable|SwapTotal|SwapFree):\s+(\d+)\s+kB$", line)
        if match:
            values[match.group(1) + "_bytes"] = int(match.group(2)) * 1024
    return values


def write_report(report: dict[str, Any], output: Path | None) -> None:
    rendered = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", action="append", default=[], help="Extra distribution to inventory.")
    parser.add_argument("--import-module", action="append", default=[], help="Module to import in an isolated child process.")
    parser.add_argument("--probe-timeout", type=float, default=30.0, help="Seconds per child or system probe.")
    parser.add_argument("--require-gpu", action="store_true", help="Return FAIL unless the full framework GPU gate passes.")
    parser.add_argument("--output", type=Path, help="Optional JSON report path.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    lspci_command = [shutil.which("lspci") or "lspci", "-nn"]
    lspci = run_command(lspci_command, args.probe_timeout)
    nvidia_lines = [
        line for line in str(lspci.get("stdout", "")).splitlines() if "NVIDIA" in line.upper()
    ]
    physical_gpu = bool(nvidia_lines)

    smi = run_command(
        [
            shutil.which("nvidia-smi") or "nvidia-smi",
            "--query-gpu=index,name,uuid,driver_version,memory.total,memory.used,compute_cap,pstate",
            "--format=csv,noheader,nounits",
        ],
        args.probe_timeout,
    )
    device_nodes = sorted(glob.glob("/dev/nvidia*"))
    compute_nodes = [path for path in device_nodes if re.fullmatch(r"/dev/nvidia\d+", path)]
    proc_gpu_dirs = sorted(glob.glob("/proc/driver/nvidia/gpus/*"))

    torch = torch_probe(args.probe_timeout)
    torch_result = torch.get("result") or {}
    gates = {
        "physical_nvidia_device": physical_gpu,
        "nvidia_smi_ok": smi.get("returncode") == 0 and bool(smi.get("stdout")),
        "compute_device_node_present": bool(compute_nodes),
        "torch_cuda_available": torch_result.get("cuda_available") is True
        and int(torch_result.get("device_count", 0)) > 0,
        "cuda_tensor_operation": torch_result.get("cuda_tensor_probe") is True,
    }
    gpu_ready = all(gates.values())

    imports = {module: import_probe(module, args.probe_timeout) for module in args.import_module}
    import_failures = [
        module
        for module, probe in imports.items()
        if probe.get("timed_out") or not (probe.get("result") or {}).get("ok", False)
    ]

    if import_failures or (args.require_gpu and not gpu_ready):
        overall_status = "FAIL"
    elif not gpu_ready:
        overall_status = "WARN"
    else:
        overall_status = "PASS"

    report: dict[str, Any] = {
        "schema_version": 1,
        "check": "omics_runtime_probe",
        "status": overall_status,
        "generated_at_utc": utc_now(),
        "require_gpu": bool(args.require_gpu),
        "host": {
            "hostname": platform.node(),
            "platform": platform.platform(),
            "kernel": platform.release(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "memory": memory_info(),
        },
        "python": {
            "executable": str(Path(sys.executable).resolve()),
            "version": platform.python_version(),
            "prefix": sys.prefix,
        },
        "packages": package_versions(args.package),
        "imports": imports,
        "import_failures": import_failures,
        "gpu": {
            "status": "PASS" if gpu_ready else "FAIL",
            "gates": gates,
            "physical_nvidia_lines": nvidia_lines,
            "lspci_probe": lspci,
            "nvidia_smi_probe": smi,
            "device_nodes": device_nodes,
            "compute_device_nodes": compute_nodes,
            "proc_gpu_directories": proc_gpu_dirs,
            "driver_version_file": read_text(Path("/proc/driver/nvidia/version")),
            "torch_probe": torch,
        },
        "limitations": [
            "This probe does not prove that a named model trains on CUDA; run its model smoke test.",
            "A passing tensor probe does not prove full-dataset memory fit or convergence.",
        ],
    }
    write_report(report, args.output)
    return 1 if overall_status == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
