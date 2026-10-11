#!/usr/bin/env python3
"""Read H5AD contracts directly with bounded HDF5 matrix chunks (metadata is O(keys))."""
from __future__ import annotations
import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


def shape_of(node):
    import h5py
    if isinstance(node, h5py.Dataset):
        return tuple(node.shape)
    if "shape" in node.attrs:
        return tuple(map(int, node.attrs["shape"]))
    if node.attrs.get("encoding-type") == "dataframe":
        return (len(node[node.attrs.get("_index", "_index")]), len(node.attrs["column-order"]))
    raise ValueError(f"Unsupported matrix encoding at {node.name}")


def matrix_chunks(node, rows, cols, chunk_values):
    import h5py
    import numpy as np
    if isinstance(node, h5py.Dataset):
        for row in range(rows):
            for col in range(0, cols, chunk_values):
                yield np.asarray(node[row:row+1, col:min(cols, col+chunk_values)]).ravel()
        return
    encoding = node.attrs.get("encoding-type")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding in {"csr_matrix", "csc_matrix"}:
        major, minor = (rows, cols) if encoding == "csr_matrix" else (cols, rows)
        start, end = int(node["indptr"][0]), int(node["indptr"][major])
        for offset in range(start, end, chunk_values):
            limit = min(end, offset+chunk_values)
            indices = node["indices"][offset:limit]
            values = node["data"][offset:limit]
            yield values[indices < minor]
        return
    if encoding == "dataframe":
        for key in node.attrs["column-order"][:cols]:
            col = node[str(key)]
            if not isinstance(col, h5py.Dataset):
                raise ValueError("Non-numeric dataframe embedding")
            for offset in range(0, rows, chunk_values):
                yield col[offset:min(rows, offset+chunk_values)]
        return
    raise ValueError(f"Unsupported matrix encoding: {encoding}")


def numeric_summary(chunks):
    import numpy as np
    count = 0
    finite = nonnegative = integer_like = True
    minimum = maximum = None
    for values in chunks:
        values = np.asarray(values)
        count += values.size
        valid = np.isfinite(values)
        finite &= bool(valid.all())
        good = values[valid]
        if good.size:
            low, high = float(good.min()), float(good.max())
            minimum = low if minimum is None else min(minimum, low)
            maximum = high if maximum is None else max(maximum, high)
            nonnegative &= bool(low >= 0)
            integer_like &= bool(np.allclose(good, np.rint(good), atol=1e-6, rtol=0))
    # Empty stored sparse data represents valid implicit zeros.
    return {"sampled_values": int(count), "all_finite": finite,
            "nonnegative": nonnegative, "integer_like": integer_like,
            "minimum": minimum, "maximum": maximum}


def resolve_count_matrix(handle, source):
    if source == "auto":
        for candidate in ("layer:counts", "raw", "X"):
            try:
                return resolve_count_matrix(handle, candidate)
            except KeyError:
                pass
        raise KeyError("No count candidate")
    path = "raw/X" if source == "raw" else "layers/"+source.split(":", 1)[1] if source.startswith("layer:") else "X" if source == "X" else None
    if path is None:
        raise ValueError(f"Unsupported count source: {source}")
    return handle[path], source


def read_index(group):
    import numpy as np
    key = group.attrs.get("_index", "_index")
    return np.asarray(group[key].asstr()[:])


def main(argv=None):
    import h5py
    import numpy as np
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--kind", choices=("scrna", "xenium", "visium", "generic"), required=True)
    parser.add_argument("--count-source", default="auto")
    parser.add_argument("--obs-key", action="append", default=[])
    parser.add_argument("--obsm-key", action="append", default=[])
    parser.add_argument("--sample-rows", type=int, default=256)
    parser.add_argument("--sample-cols", type=int, default=512)
    parser.add_argument("--chunk-values", type=int, default=65536)
    parser.add_argument("--full-count-scan", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if min(args.sample_rows, args.sample_cols, args.chunk_values) < 1:
        parser.error("Sampling/chunk limits must be positive")
    checks = []
    def check(name, ok, evidence, required=True):
        checks.append({"name": name, "status": "PASS" if ok else "FAIL" if required else "WARN",
                       "required": required, "evidence": evidence})
    report = {"schema_version": 2, "check": "h5ad_contract",
              "generated_at_utc": datetime.now(timezone.utc).isoformat(),
              "input": str(args.input.resolve()), "kind": args.kind,
              "python": {"executable": sys.executable, "version": platform.python_version()},
              "requested_count_source": args.count_source,
              "reader": "h5py_direct_chunked", "chunk_values": args.chunk_values,
              "sampling_limit": {"rows": args.sample_rows, "columns": args.sample_cols},
              "sampling_disclaimer": "Only requested count values are scanned; metadata/key uniqueness uses O(number of keys) memory. Numerics do not prove raw-count provenance.",
              "checks": checks}
    try:
        with h5py.File(args.input, "r") as handle:
            check("hdf5_reopen", True, {"reader": "h5py", "no_anndata_or_layer_materialization": True})
            obs, var = read_index(handle["obs"]), read_index(handle["var"])
            n_obs, n_vars = len(obs), len(var)
            report["shape"] = [n_obs, n_vars]
            report["keys"] = {k: list(handle[k].keys()) if k in handle else [] for k in ("obs", "var", "layers", "obsm", "uns")}
            check("nonzero_shape", n_obs > 0 and n_vars > 0, report["shape"])
            check("unique_obs_names", len(np.unique(obs)) == n_obs, {"n_obs": n_obs})
            check("unique_var_names", len(np.unique(var)) == n_vars, {"n_vars": n_vars})
            for key in args.obs_key:
                present = key in handle["obs"]
                evidence = {"key": key, "present": present}
                missing = 0
                if present:
                    node = handle["obs"][key]
                    if isinstance(node, h5py.Group) and "codes" in node:
                        missing = sum(int((node["codes"][i:i+args.chunk_values] < 0).sum()) for i in range(0, n_obs, args.chunk_values))
                    elif isinstance(node, h5py.Dataset):
                        for i in range(0, n_obs, args.chunk_values):
                            vals = node[i:i+args.chunk_values]
                            missing += int((~np.isfinite(vals)).sum()) if vals.dtype.kind in "fci" else int(np.isin(vals, [b"", ""]).sum())
                    else:
                        raise ValueError(f"Unsupported obs encoding: {key}")
                    evidence["missing"] = missing
                check("obs_key:"+key, present and missing == 0, evidence)
            required_obsm = list(dict.fromkeys(args.obsm_key + (["spatial"] if args.kind in {"xenium", "visium"} else [])))
            for key in required_obsm:
                path = "obsm/"+key
                if path not in handle:
                    check("obsm_key:"+key, False, {"present": False})
                    continue
                node = handle[path]
                shape = shape_of(node)
                if len(shape) != 2:
                    check("obsm_key:"+key, False, {"shape": list(shape)})
                    continue
                summary = numeric_summary(matrix_chunks(node, min(args.sample_rows, shape[0]), shape[1], args.chunk_values))
                check("obsm_key:"+key, shape[0] == n_obs and shape[1] >= 2 and summary["all_finite"], {"shape": list(shape), "sample": summary})
            node, resolved = resolve_count_matrix(handle, args.count_source)
            report["resolved_count_source"] = resolved
            check("explicit_count_source", args.count_source != "auto", {"resolved": resolved, "provenance_verified": False}, required=False)
            shape = shape_of(node)
            expected_vars = len(read_index(handle["raw/var"])) if resolved == "raw" else n_vars
            check("count_shape", shape == (n_obs, expected_vars), {"shape": list(shape)})
            rows, cols = shape if args.full_count_scan else (min(args.sample_rows, shape[0]), min(args.sample_cols, shape[1]))
            summary = numeric_summary(matrix_chunks(node, rows, cols, args.chunk_values))
            summary.update(matrix_shape=list(shape), sample_shape=[rows, cols], scope="full_stored_count_values" if args.full_count_scan else "sampled_count_values")
            report["count_sample"] = summary
            check("full_raw_count_numerics" if args.full_count_scan else "sampled_raw_count_contract",
                  summary["all_finite"] and summary["nonnegative"] and summary["integer_like"], summary)
            check("sample_contains_nonzero_values", summary["maximum"] is not None and summary["maximum"] > 0, summary, required=False)
    except Exception as exc:
        report["exception"] = {"type": type(exc).__name__, "message": str(exc)}
        check("read_contract", False, report["exception"])
    report["required_failures"] = [c["name"] for c in checks if c["required"] and c["status"] == "FAIL"]
    report["warnings"] = [c["name"] for c in checks if c["status"] == "WARN"]
    report["status"] = "FAIL" if report["required_failures"] else "WARN" if report["warnings"] else "PASS"
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered+"\n", encoding="utf-8")
    print(rendered)
    return int(bool(report["required_failures"]))


if __name__ == "__main__":
    raise SystemExit(main())
