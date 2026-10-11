#!/usr/bin/env python3
"""Render derivative-only Xenium DotPlots with native Seurat and strict parity."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from scipy.io import mmwrite
import yaml
from annotation_contract import resolve as resolve_annotation_contract


SKILL_DIR = Path(__file__).resolve().parents[1]
DEFAULT_ORACLE = SKILL_DIR.parent / "_shared" / "dotplot" / "seurat_dotplot_renderer.R"
MODALITY = globals().get("MODALITY", "xenium")
RENDERER_NAME = "render-" + MODALITY + "-dotplot"
MANIFEST_REL = Path("manifest/sha256_manifest.csv")
MANIFEST_HASH_REL = Path("manifest/sha256_manifest.sha256")
STYLE_DEFAULTS = json.loads((Path(__file__).resolve().parents[2] / "_shared/dotplot/style_defaults.json").read_text())
PROGRAM_PALETTE = tuple(STYLE_DEFAULTS["program_palette"])
ALLOWED_USE_CASES = {
    "scrna_top_marker",
    "xenium_top_marker",
    "cross_modality",
    "xenium_rcn",
    "epi_subtype",
}

DOT_DEFAULTS: dict[str, Any] = {
    "normalization_scale_factor": 10_000.0,
    "scaling_ddof": 1,
    "col_min": -2.5,
    "col_mid": 0.0,
    "col_max": 2.5,
    "dot_min": 0.0,
    "scale_by": "radius",
    "dot_scale": 6.0,
    "color_low": "#2166AC",
    "color_mid": "#FFFFFF",
    "color_high": "#B2182B",
}

X_AXIS_TEXT_ANGLE_DEGREES = 60.0
X_AXIS_TEXT_HJUST = 1.0
X_AXIS_TEXT_VJUST = 1.0
X_AXIS_TEXT_FONT_SIZE_PT = 9.0
FIGURE_BASE_WIDTH_INCHES = 6.0
MINIMUM_MARKER_SLOT_INCHES = 0.42
MARKER_SLOT_BASE_INCHES = 0.18
MARKER_SLOT_GENE_CHARACTER_INCHES = 0.025


def now_local() -> str:
    return datetime.now().astimezone().isoformat()


def load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return value


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def resolve_path(root: Path, value: str | Path) -> Path:
    candidate = Path(value).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def display_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def file_record(path: Path, root: Path) -> dict[str, Any]:
    return {
        "path": display_path(path, root),
        "absolute_path": str(path.resolve()),
        "bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def render_value(value: Any) -> str:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (dict, list, tuple, set)):
        if isinstance(value, set):
            value = sorted(value)
        return json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)
    return str(value)


def resolve_figure_width(
    requested_width_inches: float,
    genes: Sequence[str],
    auto_expand_width: bool,
) -> dict[str, Any]:
    marker_count = len(genes)
    if marker_count < 1:
        raise ValueError("figure width cannot be resolved for an empty marker panel")
    longest_gene_length = max(len(str(gene)) for gene in genes)
    marker_slot_width_inches = max(
        MINIMUM_MARKER_SLOT_INCHES,
        MARKER_SLOT_BASE_INCHES
        + MARKER_SLOT_GENE_CHARACTER_INCHES * longest_gene_length,
    )
    computed_minimum_width_inches = (
        FIGURE_BASE_WIDTH_INCHES + marker_count * marker_slot_width_inches
    )
    if not auto_expand_width and requested_width_inches < computed_minimum_width_inches:
        raise ValueError(
            "output.width_inches is below the marker-label layout minimum while "
            "output.auto_expand_width is false: "
            f"requested={requested_width_inches:.6g}, "
            f"minimum={computed_minimum_width_inches:.6g}"
        )
    final_width_inches = (
        max(requested_width_inches, computed_minimum_width_inches)
        if auto_expand_width
        else requested_width_inches
    )
    return {
        "auto_expand_width": auto_expand_width,
        "requested_width_inches": float(requested_width_inches),
        "base_width_inches": FIGURE_BASE_WIDTH_INCHES,
        "marker_count": marker_count,
        "longest_gene_length": longest_gene_length,
        "minimum_marker_slot_inches": MINIMUM_MARKER_SLOT_INCHES,
        "marker_slot_formula": "max(0.42, 0.18 + 0.025 * longest_gene_length)",
        "marker_slot_width_inches": float(marker_slot_width_inches),
        "computed_minimum_width_inches": float(computed_minimum_width_inches),
        "final_width_inches": float(final_width_inches),
    }


def add_check(
    checks: list[dict[str, str]],
    check_id: str,
    observed: Any,
    expected: Any,
    passed: bool,
    details: str = "",
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "observed": render_value(observed),
            "expected": render_value(expected),
            "status": "PASS" if passed else "FAIL",
            "details": details,
        }
    )


def prepare_config(config_path: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = load_yaml(config_path)
    if not str(config.get("analysis_id", "")).strip():
        raise ValueError("analysis_id is required")
    config["status"] = str(config.get("status", "exploratory_not_current"))
    project_root = resolve_path(config_path.parent, config.get("project_root", "."))
    if not project_root.is_dir():
        raise FileNotFoundError(f"project_root does not exist: {project_root}")
    config["_config_path"] = config_path
    config["_project_root"] = project_root

    inputs = config.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("inputs must be a mapping")
    for key in ("counts_h5ad", "sidecar", "marker_yaml"):
        if not str(inputs.get(key, "")).strip():
            raise ValueError(f"inputs.{key} is required")
        inputs[f"_{key}"] = resolve_path(project_root, inputs[key])
        if not inputs[f"_{key}"].is_file():
            raise FileNotFoundError(inputs[f"_{key}"])
    if inputs.get("counts_layer", "counts") != "counts":
        raise ValueError('canonical rendering requires inputs.counts_layer: counts')

    join = config.get("join")
    if not isinstance(join, dict):
        raise ValueError("join must be an explicit mapping")
    for key in ("sidecar_key", "target_key"):
        if not str(join.get(key, "")).strip():
            raise ValueError(f"join.{key} is required")
        join[key] = str(join[key])
    required_true = ("sidecar_key_unique", "target_key_unique")
    for key in required_true:
        if join.get(key) is not True:
            raise ValueError(f"join.{key} must be true")
    if join.get("row_order_join_forbidden") is not True:
        raise ValueError("join.row_order_join_forbidden must be true")
    if join["target_key"] != "obs_name":
        raise ValueError("the target key must be the AnnData obs_name index")
    if MODALITY == "xenium" and join["sidecar_key"] != "obs_name":
        if join["sidecar_key"] == "cell_id":
            raise ValueError("Xenium cell_id joins are forbidden; use obs_name -> obs_name")
        raise ValueError("Xenium sidecar joins must use obs_name -> obs_name")

    groupings = config.get("groupings")
    if not isinstance(groupings, dict) or not groupings:
        raise ValueError("groupings must be a non-empty mapping")
    for name, block in groupings.items():
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(name)):
            raise ValueError(
                f"grouping name {name!r} must contain only letters, digits, ., _, or -"
            )
        if not isinstance(block, dict) or not str(block.get("field", "")).strip():
            raise ValueError(f"groupings.{name}.field is required")
        order = [str(item) for item in block.get("order", [])]
        if len(order) < 1 or len(order) != len(set(order)):
            raise ValueError(f"groupings.{name}.order must contain >=1 unique labels")
        block["order"] = order
        label_map = block.get("label_map", {})
        if not isinstance(label_map, dict):
            raise ValueError(f"groupings.{name}.label_map must be a mapping")
        label_map = {str(key): str(value) for key, value in label_map.items()}
        allow_many_to_one = bool(block.get("allow_many_to_one_label_map", False))
        if label_map and set(label_map.values()) != set(order):
            raise ValueError(
                f"groupings.{name}.label_map values must cover the declared order exactly"
            )
        if (
            label_map
            and len(label_map.values()) != len(set(label_map.values()))
            and not allow_many_to_one
        ):
            raise ValueError(
                f"groupings.{name}.label_map is many-to-one; set "
                "allow_many_to_one_label_map: true to authorize the explicit merge"
            )
        block["label_map"] = label_map
        block["allow_many_to_one_label_map"] = allow_many_to_one
        block["title"] = str(block.get("title", name.replace("_", " ")))

    selection = config.get("marker_selection")
    if not isinstance(selection, dict):
        raise ValueError("marker_selection must be a mapping")
    for key in (
        "include_categories",
        "exclude_categories",
        "include_programs",
        "exclude_programs",
    ):
        selection[key] = [str(item) for item in selection.get(key, [])]
    threshold = float(selection.get("minimum_max_group_pct_exp", 0.0))
    if not 0.0 <= threshold <= 100.0:
        raise ValueError("minimum_max_group_pct_exp must be between 0 and 100")
    selection["minimum_max_group_pct_exp"] = threshold
    filter_groupings = [str(item) for item in selection.get("filter_groupings", groupings)]
    unknown_filter_groupings = sorted(set(filter_groupings) - set(groupings))
    if unknown_filter_groupings:
        raise ValueError(f"unknown filter_groupings: {unknown_filter_groupings}")
    selection["filter_groupings"] = filter_groupings
    deduplication = selection.get("deduplication", {"strategy": "first_occurrence"})
    if not isinstance(deduplication, dict):
        raise ValueError("marker_selection.deduplication must be a mapping")
    if deduplication.get("strategy", "first_occurrence") not in {
        "error",
        "first_occurrence",
        "move_selected_to_program",
    }:
        raise ValueError("unsupported marker deduplication strategy")
    selection["deduplication"] = deduplication

    dotplot = dict(DOT_DEFAULTS)
    supplied_dotplot = config.get("dotplot", {})
    if not isinstance(supplied_dotplot, dict):
        raise ValueError("dotplot must be a mapping")
    supplied_mid = str(supplied_dotplot.get("color_mid", "")).upper()
    if supplied_mid == "#F7F7F7":
        raise ValueError(
            "dotplot.color_mid #F7F7F7 is a retired pre-v7 value; "
            "migrate the run config to #FFFFFF"
        )
    dotplot.update(supplied_dotplot)
    for key, expected in DOT_DEFAULTS.items():
        if dotplot[key] != expected:
            raise ValueError(
                f"dotplot.{key} must equal {expected!r} for bundled Seurat parity"
            )
    config["dotplot"] = dotplot

    oracle = config.get("seurat_native", config.get("oracle"))
    if not isinstance(oracle, dict):
        raise ValueError("seurat_native must be a mapping")
    for key in (
        "rscript",
        "expected_r_version",
        "expected_seurat_version",
        "expected_seuratobject_version",
    ):
        if not str(oracle.get(key, "")).strip():
            raise ValueError(f"seurat_native.{key} is required")
    oracle_script = oracle.get("script", DEFAULT_ORACLE)
    oracle["_script"] = resolve_path(project_root, oracle_script)
    if not oracle["_script"].is_file():
        raise FileNotFoundError(oracle["_script"])
    for key, default in {
        "avg_exp_abs_tolerance": 1e-9,
        "avg_exp_log1p_abs_tolerance": 1e-9,
        "zscore_abs_tolerance": 1e-9,
        "scaled_abs_tolerance": 1e-9,
        "pct_exp_abs_tolerance": 1e-10,
    }.items():
        oracle[key] = float(oracle.get(key, default))
        if oracle[key] < 0:
            raise ValueError(f"seurat_native.{key} must be non-negative")
    config["seurat_native"] = oracle

    output = config.get("output")
    if not isinstance(output, dict) or not str(output.get("root", "")).strip():
        raise ValueError("output.root is required")
    output["_root"] = resolve_path(project_root, output["root"])
    output["dpi"] = int(output.get("dpi", 240))
    if output["dpi"] <= 0:
        raise ValueError("output.dpi must be positive")
    output["width_inches"] = float(output.get("width_inches", 20.0))
    output["height_inches"] = None if output.get("height_inches") is None else float(output["height_inches"])
    style = dict(STYLE_DEFAULTS)
    overrides = config.get("style", {})
    if not isinstance(overrides, dict) or set(overrides) - (set(style) - {"program_palette"}):
        raise ValueError("style contains unsupported keys; set program colors in marker YAML")
    style.update(overrides)
    for key in set(style) - {"program_palette"}:
        style[key] = float(style[key])
        if not np.isfinite(style[key]) or style[key] <= 0:
            raise ValueError(f"style.{key} must be finite and positive")
    config["style"] = style
    for block in groupings.values():
        block["height_inches"] = output["height_inches"] or max(3.2, 2.2 + 0.30 * len(block["order"]))
        block["style"] = style
        block["subtitle"] = str(block.get("subtitle", ""))
        if block["height_inches"] <= style["strip_height_inches"]:
            raise ValueError("figure height must exceed strip height")
    if not np.isfinite(output["width_inches"]) or output["width_inches"] <= 0 or (output["height_inches"] is not None and (not np.isfinite(output["height_inches"]) or output["height_inches"] <= 0)):
        raise ValueError("output width_inches and height_inches must be positive")
    auto_expand_width = output.get("auto_expand_width", True)
    if not isinstance(auto_expand_width, bool):
        raise ValueError("output.auto_expand_width must be true or false")
    output["auto_expand_width"] = auto_expand_width
    project = config.setdefault("project", {})
    if not isinstance(project, dict):
        raise ValueError("project must be a mapping")
    use_case = str(project.get("annotation_use_case", "")).strip()
    if use_case not in ALLOWED_USE_CASES:
        raise ValueError(
            "project.annotation_use_case must resolve exactly one of: "
            + ", ".join(sorted(ALLOWED_USE_CASES))
        )
    for key in ("annotation_guide", "annotation_registry"):
        if not str(project.get(key, "")).strip():
            raise ValueError(f"project.{key} is required")
        path = resolve_path(project_root, project[key])
        if not path.is_file():
            raise FileNotFoundError(path)
        project[f"_{key}"] = path
    track = resolve_annotation_contract(project["_annotation_guide"], project["_annotation_registry"], use_case)
    expected = project.get("expected_annotation_contract") or {}
    if not isinstance(expected, dict):
        raise ValueError("project.expected_annotation_contract must be a mapping")
    for key, value in expected.items():
        if track.get(key) != value:
            raise ValueError(f"annotation contract mismatch for {use_case}.{key}")
    project["_resolved_annotation_track"] = track
    bindings = {}
    roles = {use_case: track}
    for name, block in groupings.items():
        role = block.get("annotation_use_case", use_case)
        if role not in ALLOWED_USE_CASES:
            raise ValueError(f"unknown annotation role: {role}")
        if role not in roles:
            roles[role] = resolve_annotation_contract(project["_annotation_guide"], project["_annotation_registry"], role)
        current = roles[role]
        if str(current.get("status", "")).startswith(("paused", "pending", "historical")):
            raise ValueError(f"{role}: current source unavailable ({current.get('status')})")
        field = block["field"]
        assignment_id = block.get("assignment_id")
        if assignment_id:
            candidates = [s for s in current.get("assignments", []) if assignment_id in (s.get("assignment_id"), s.get("namespace")) and s.get("label_column") == field]
        else:
            candidates = [s for s in current.get("sources", []) if field in s.get("label_columns", {}).values()]
        requested_source = block.get("annotation_source_id")
        if requested_source:
            candidates = [s for s in candidates if s.get("source_id") == requested_source]
        if len(candidates) != 1:
            raise ValueError(f"{role}:{field}: resolve exactly one registered source or qualified assignment; found {len(candidates)}")
        source = candidates[0]
        source_join = source.get("join", {})
        if any(source_join.get(k) != join[k] for k in ("sidecar_key", "target_key", "sidecar_key_unique")) or source_join.get("row_order_join_forbidden") is not True:
            raise ValueError(f"{role}:{field}: join differs from selected source {source.get('source_id', assignment_id)}")
        source_modality = source.get("modality", current.get("modality"))
        if source_modality in ("scrna", "xenium") and source_modality != MODALITY:
            raise ValueError(f"{role}:{field}: source modality is {source_modality}, renderer is {MODALITY}")
        source_path = resolve_path(project_root, source["path"])
        if not source_path.is_file():
            raise FileNotFoundError(source_path)
        bindings[name] = {"use_case": role, "track_id": current["track_id"], "source": source, "path": str(source_path)}
    project["_annotation_bindings"] = bindings
    project["_resolved_annotation_roles"] = roles
    project["_registry_join_source_ids"] = [b["source"].get("source_id", b["source"].get("assignment_id")) for b in bindings.values()]
    return config



def validate_annotation_sources(config, aligned):
    """Verify derivative labels against actual registered rows, never row order."""
    key = config["join"]["sidecar_key"]
    cache = {}
    records = []
    def load(path, source_key):
        cache_key = (str(path), source_key)
        if cache_key not in cache:
            frame = read_sidecar(Path(path))
            if source_key not in frame or frame[source_key].isna().any() or not frame[source_key].astype(str).is_unique:
                raise ValueError(f"registered source key is missing or nonunique: {path}:{source_key}")
            cache[cache_key] = frame.set_index(frame[source_key].astype(str), drop=False)
        return cache[cache_key]
    for name, binding in config["project"]["_annotation_bindings"].items():
        source = binding["source"]
        field = config["groupings"][name]["field"]
        actual_hash = sha256_file(Path(binding["path"]))
        if source.get("sha256") and source["sha256"] != actual_hash:
            raise ValueError(f"{name}: registered source sha256 mismatch")
        full = load(binding["path"], key)
        missing = aligned.index.difference(full.index)
        if len(missing):
            raise ValueError(f"{name}: {len(missing)} obs_names outside registered source universe")
        rows = full.loc[aligned.index]
        if field not in rows:
            raise ValueError(f"{name}: source field absent: {field}")
        expected = rows[field].astype('string').fillna('<MISSING>')
        actual = aligned[field].astype('string').fillna('<MISSING>')
        if not expected.equals(actual):
            raise ValueError(f"{name}: derivative labels differ from registered source by obs_name")
        # Required operational gates are derived from the registered contract,
        # including source-specific prose until the registry gains typed predicates.
        gates = set()
        text = source.get('fallback', '')
        gates.update(re.findall(r'only when ([A-Za-z_][A-Za-z0-9_]*) is true', text))
        for exclusion in source.get('exclusion_gates', []):
            match = re.fullmatch(r'([A-Za-z_][A-Za-z0-9_]*) == false', exclusion)
            if match:
                gates.add(match[1])
        gates.update(re.findall(r'([A-Za-z_][A-Za-z0-9_]*) to omit remaining unresolved rows', text))
        if field == source.get('label_columns', {}).get('primary_phenotype') and 'primary_neighbor_eligible' in source.get('status_columns', []):
            gates.add('primary_neighbor_eligible')
        if field == source.get('label_columns', {}).get('strict_phenotype') and 'strict_neighbor_eligible' in source.get('status_columns', []):
            gates.add('strict_neighbor_eligible')
        flags = source.get('analysis_flags', {})
        for gate_role in ('operational_include', 'resolved_only_include'):
            if isinstance(flags.get(gate_role), str):
                gates.add(flags[gate_role])
        mode = config['groupings'][name].get('annotation_gate_mode', 'display')
        if mode not in ('display', 'fine_analysis'):
            raise ValueError(f'{name}: annotation_gate_mode must be display or fine_analysis')
        if mode == 'fine_analysis':
            gate = flags.get('fine_eligible') or flags.get('fine_analysis')
            if not gate:
                matches = re.findall(r'([A-Za-z_][A-Za-z0-9_]*) for fine-analysis eligibility', text)
                gate = matches[0] if len(matches) == 1 else None
            if not gate:
                raise ValueError(f'{name}: source has no typed fine-analysis gate; cannot infer eligibility')
            gates.add(gate)
        for gate in sorted(gates):
            if gate not in rows or not rows[gate].astype('string').str.lower().isin(['true', '1']).all():
                raise ValueError(f'{name}: registered applicability gate failed: {gate}')
        scope = source.get('cohort_scope', '')
        cohort = source.get('cohort_gate', {})
        if scope and not cohort:
            cohort = config['project'].get('cohort_scope_gates', {}).get(scope)
            if not cohort:
                raise ValueError(f'{name}: cohort scope requires explicit typed predicates: {scope}')
        predicates = dict(cohort.get('allowed_values', {}))
        if cohort.get('sample_ids'):
            predicates['sample_id'] = cohort['sample_ids']
        if not set(cohort).issubset({'sample_ids', 'allowed_values'}):
            raise ValueError(f'{name}: unsupported cohort gate keys')
        if (scope or cohort) and not predicates:
            raise ValueError(f'{name}: empty cohort gate predicates')
        for column, allowed in predicates.items():
            if not isinstance(allowed, list) or not allowed or column not in rows or not rows[column].astype(str).isin([str(v) for v in allowed]).all():
                raise ValueError(f'{name}: source cohort gate failed: {column}')
        # A lower-priority source cannot stand in for a registered higher-priority
        # overlap, even when the primary fields have different names.
        track = config['project']['_resolved_annotation_roles'][binding['use_case']]
        for higher in track.get('sources', []):
            if higher.get('priority', 0) <= source.get('priority', 0):
                continue
            if higher.get('modality', MODALITY) != MODALITY:
                continue
            other = load(resolve_path(config['_project_root'], higher['path']), higher['join']['sidecar_key'])
            if len(aligned.index.intersection(other.index)):
                raise ValueError(f"{name}: higher-priority source overlaps selected cells: {higher['source_id']}")
        records.append({**binding, 'field': field, 'gate_mode': mode,
                        'verified_gates': sorted(gates), 'rows_verified': len(rows),
                        'source_sha256': actual_hash})
    return records


def decode_strings(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind == "S":
        return np.char.decode(array, "utf-8")
    if array.dtype.kind == "U":
        return array
    return np.asarray(
        [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in array],
        dtype=object,
    )


def dataframe_index_key(group: h5py.Group) -> str:
    key = group.attrs.get("_index", "_index")
    return key.decode("utf-8") if isinstance(key, bytes) else str(key)


def read_h5ad_contract(path: Path, layer: str) -> tuple[pd.Index, list[str], tuple[int, int], str]:
    with h5py.File(path, "r") as handle:
        if "layers" not in handle or layer not in handle["layers"]:
            raise ValueError(f'H5AD is missing layers["{layer}"]')
        counts = handle["layers"][layer]
        encoding = counts.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8")
        if str(encoding) != "csr_matrix":
            raise ValueError('layers["counts"] must use AnnData CSR encoding')
        shape = tuple(int(value) for value in counts.attrs["shape"])
        obs_values = decode_strings(handle["obs"][dataframe_index_key(handle["obs"])][:])
        var_values = decode_strings(handle["var"][dataframe_index_key(handle["var"])][:])
    obs_names = pd.Index(obs_values.astype(str), name="obs_name")
    var_names = [str(value) for value in var_values]
    if len(obs_names) != shape[0] or len(var_names) != shape[1]:
        raise ValueError("H5AD axes do not match counts shape")
    if not obs_names.is_unique:
        raise ValueError("H5AD obs_names are not unique")
    if len(var_names) != len(set(var_names)):
        raise ValueError("H5AD var_names are not unique")
    return obs_names, var_names, shape, str(encoding)


def read_sidecar(path: Path, usecols: Sequence[str] | None = None) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith(".parquet"):
        return pd.read_parquet(path, columns=None if usecols is None else list(usecols))
    if name.endswith((".tsv", ".tsv.gz")):
        return pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    if name.endswith((".csv", ".csv.gz")):
        return pd.read_csv(path, usecols=usecols, low_memory=False)
    raise ValueError("sidecar must be CSV, CSV.GZ, TSV, TSV.GZ, or Parquet")


def align_sidecar_to_obs_names(
    obs_names: pd.Index,
    sidecar: pd.DataFrame,
    label_columns: Sequence[str],
    sidecar_key: str,
) -> pd.DataFrame:
    required = {sidecar_key, *label_columns}
    missing_columns = sorted(required - set(sidecar.columns))
    if missing_columns:
        raise ValueError(f"sidecar is missing columns: {missing_columns}")
    if sidecar[sidecar_key].isna().any():
        raise ValueError(f"sidecar {sidecar_key} contains missing values")
    source_names = pd.Index(sidecar[sidecar_key].astype(str), name="obs_name")
    if not source_names.is_unique:
        duplicates = source_names[source_names.duplicated()].unique()[:5].tolist()
        raise ValueError(f"sidecar {sidecar_key} is not unique: {duplicates}")
    target_names = pd.Index(obs_names.astype(str), name="obs_name")
    missing = target_names.difference(source_names)
    extra = source_names.difference(target_names)
    if len(missing) or len(extra):
        raise ValueError(
            f"{sidecar_key} -> obs_name exact-set mismatch: "
            f"missing_from_sidecar={len(missing)}, extra_in_sidecar={len(extra)}"
        )
    aligned = (
        sidecar.assign(__target_obs_name=source_names)
        .set_index("__target_obs_name")
        .loc[target_names]
    )
    aligned.index.name = "obs_name"
    if not aligned.index.equals(target_names):
        raise ValueError(f"{sidecar_key} -> obs_name alignment failed")
    return aligned.loc[:, list(label_columns)].copy()


def build_group_contract(
    aligned: pd.DataFrame,
    groupings: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, list[str]], pd.DataFrame]:
    codes_by_level: dict[str, np.ndarray] = {}
    counts_by_level: dict[str, np.ndarray] = {}
    groups_by_level: dict[str, list[str]] = {}
    rows: list[dict[str, Any]] = []
    for level, block in groupings.items():
        field = str(block["field"])
        order = [str(item) for item in block["order"]]
        series = aligned[field]
        if series.isna().any():
            raise ValueError(f"{level}:{field} contains missing labels")
        source_labels = series.astype(str)
        label_map = {str(key): str(value) for key, value in block.get("label_map", {}).items()}
        if label_map:
            observed_source = set(source_labels.unique())
            if observed_source != set(label_map):
                raise ValueError(
                    f"{level} source labels differ from label_map: "
                    f"missing={sorted(set(label_map) - observed_source)}, "
                    f"extra={sorted(observed_source - set(label_map))}"
                )
            labels = source_labels.map(label_map)
        else:
            labels = source_labels
        observed = set(labels.unique())
        if observed != set(order):
            raise ValueError(
                f"{level} labels differ from declared order: "
                f"missing={sorted(set(order) - observed)}, extra={sorted(observed - set(order))}"
            )
        categorical = pd.Categorical(labels, categories=order, ordered=True)
        codes = categorical.codes.astype(np.int32)
        if np.any(codes < 0):
            raise ValueError(f"{level} produced invalid group codes")
        counts = np.bincount(codes, minlength=len(order)).astype(np.int64)
        observed_counts = {group: int(counts[index]) for index, group in enumerate(order)}
        expected_counts = block.get("expected_counts")
        if expected_counts is not None:
            expected = {str(key): int(value) for key, value in expected_counts.items()}
            if observed_counts != expected:
                raise ValueError(
                    f"{level} group counts changed: observed={observed_counts}, expected={expected}"
                )
        codes_by_level[level] = codes
        counts_by_level[level] = counts
        groups_by_level[level] = order
        for group_order, group in enumerate(order):
            rows.append(
                {
                    "grouping_level": level,
                    "group_field": field,
                    "group_order": group_order,
                    "group_label": group,
                    "source_labels": json.dumps(
                        sorted(
                            source
                            for source, display in label_map.items()
                            if display == group
                        )
                        if label_map
                        else [group],
                        ensure_ascii=False,
                    ),
                    "n_cells": int(counts[group_order]),
                }
            )
    return codes_by_level, counts_by_level, groups_by_level, pd.DataFrame(rows)


def expand_marker_slots(marker_yaml: Path) -> pd.DataFrame:
    source = load_yaml(marker_yaml)
    programs = source.get("programs")
    rows: list[dict[str, Any]] = []
    source_slot = 0

    source_program_global_order = 0

    def add_program(category: str, program: Mapping[str, Any], program_order: int) -> None:
        nonlocal source_slot
        nonlocal source_program_global_order
        name = str(program.get("name", "")).strip()
        if not name:
            raise ValueError("every marker program requires a name")
        short_label = str(program.get("short_label", name)).strip()
        if not short_label:
            raise ValueError(f"{category}:{name} short_label cannot be empty")
        color_value = program.get("color")
        color = str(color_value).upper() if color_value is not None else ""
        if color and not re.fullmatch(r"#[0-9A-F]{6}", color):
            raise ValueError(f"{category}:{name} color must be a six-digit hex color")
        role_blocks: list[tuple[str, Sequence[Any]]] = []
        if "genes" in program:
            role_blocks.append(("gene", program.get("genes", [])))
        else:
            for role in ("core", "support"):
                role_blocks.append((role, program.get(role, [])))
        marker_order = 0
        for role, genes in role_blocks:
            if genes is None:
                continue
            if not isinstance(genes, list):
                raise ValueError(f"{category}:{name}:{role} markers must be a list")
            for role_order, gene_value in enumerate(genes):
                if isinstance(gene_value, dict):
                    gene = str(gene_value.get("gene", "")).strip()
                    marker_role = str(gene_value.get("role", role)).strip()
                else:
                    gene = str(gene_value).strip()
                    marker_role = role
                if not gene:
                    raise ValueError(f"{category}:{name} contains an empty gene")
                if not marker_role:
                    raise ValueError(f"{category}:{name}:{gene} role cannot be empty")
                rows.append(
                    {
                        "source_slot": source_slot,
                        "source_category": category,
                        "source_program": name,
                        "source_program_order": program_order,
                        "source_marker_role": marker_role,
                        "source_role_order": role_order,
                        "source_marker_order": marker_order,
                        "source_program_global_order": source_program_global_order,
                        "source_program_short_label": short_label,
                        "source_program_color": color,
                        "gene": gene,
                    }
                )
                source_slot += 1
                marker_order += 1
        source_program_global_order += 1

    if isinstance(programs, dict):
        for category, category_programs in programs.items():
            if not isinstance(category_programs, list):
                raise ValueError(f"programs.{category} must be a list")
            for program_order, program in enumerate(category_programs):
                if not isinstance(program, dict):
                    raise ValueError(f"programs.{category}[{program_order}] must be a mapping")
                add_program(str(category), program, program_order)
    elif isinstance(programs, list):
        category_orders: dict[str, int] = {}
        for program in programs:
            if not isinstance(program, dict):
                raise ValueError("flat marker programs must be mappings")
            category = str(program.get("category", "target"))
            program_order = category_orders.get(category, 0)
            add_program(category, program, program_order)
            category_orders[category] = program_order + 1
    else:
        raise ValueError("marker YAML programs must be a mapping or list")
    if not rows:
        raise ValueError("marker YAML contains no marker slots")
    result = pd.DataFrame(rows)
    if result[["source_category", "source_program"]].drop_duplicates().duplicated(
        "source_program"
    ).any():
        raise ValueError("marker program names must be unique across categories")
    return result


def _first_occurrence_assignments(
    eligible: pd.DataFrame,
    duplicate_policy: str,
) -> list[dict[str, Any]]:
    duplicated = eligible.loc[eligible["gene"].duplicated(keep=False), "gene"].unique().tolist()
    if duplicated and duplicate_policy == "error":
        raise ValueError(f"duplicate selected marker genes: {duplicated[:20]}")
    assignments: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in eligible.itertuples(index=False):
        if row.gene in seen:
            continue
        seen.add(row.gene)
        assignments.append(
            {
                "gene": row.gene,
                "category": row.source_category,
                "program": row.source_program,
                "program_short_label": row.source_program_short_label,
                "program_color": row.source_program_color,
                "marker_role": row.source_marker_role,
                "chosen_source_slot": int(row.source_slot),
            }
        )
    return assignments


def select_and_deduplicate_markers(
    slots: pd.DataFrame,
    selection: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    include_categories = set(selection["include_categories"])
    exclude_categories = set(selection["exclude_categories"])
    include_programs = set(selection["include_programs"])
    exclude_programs = set(selection["exclude_programs"])
    category_ok = (~slots["source_category"].isin(exclude_categories)) & (
        slots["source_category"].isin(include_categories) if include_categories else True
    )
    program_ok = (~slots["source_program"].isin(exclude_programs)) & (
        slots["source_program"].isin(include_programs) if include_programs else True
    )
    eligible_mask = category_ok & program_ok
    eligible = slots.loc[eligible_mask].sort_values("source_slot").copy()
    if eligible.empty:
        raise ValueError("marker selection removed every marker slot")

    dedupe = selection["deduplication"]
    strategy = str(dedupe.get("strategy", "first_occurrence"))
    if strategy in {"error", "first_occurrence"}:
        assignments = _first_occurrence_assignments(eligible, strategy)
        shared_block: list[dict[str, Any]] = []
        shared_position: str | int = "last"
    else:
        shared = dedupe.get("shared_program")
        if not isinstance(shared, dict):
            raise ValueError("move_selected_to_program requires shared_program")
        shared_name = str(shared.get("name", "")).strip()
        shared_genes = [str(gene) for gene in shared.get("genes", [])]
        if not shared_name or not shared_genes or len(shared_genes) != len(set(shared_genes)):
            raise ValueError("shared_program requires a name and unique genes")
        missing_shared = [gene for gene in shared_genes if gene not in set(eligible["gene"])]
        if missing_shared:
            raise ValueError(f"shared genes are absent after selection: {missing_shared}")
        remaining_policy = str(dedupe.get("remaining_duplicates", "first_occurrence"))
        if remaining_policy not in {"error", "first_occurrence"}:
            raise ValueError("remaining_duplicates must be error or first_occurrence")
        remaining = eligible.loc[~eligible["gene"].isin(shared_genes)]
        assignments = _first_occurrence_assignments(remaining, remaining_policy)
        shared_category = str(shared.get("category", eligible.iloc[0]["source_category"]))
        shared_short_label = str(shared.get("short_label", shared_name)).strip()
        shared_color_value = shared.get("color")
        shared_color = (
            str(shared_color_value).upper() if shared_color_value is not None else ""
        )
        if shared_color and not re.fullmatch(r"#[0-9A-F]{6}", shared_color):
            raise ValueError("shared_program.color must be a six-digit hex color")
        shared_role = str(shared.get("role", "shared_marker")).strip()
        shared_block = []
        for gene in shared_genes:
            chosen = eligible.loc[eligible["gene"].eq(gene)].sort_values("source_slot").iloc[0]
            shared_block.append(
                {
                    "gene": gene,
                    "category": shared_category,
                    "program": shared_name,
                    "program_short_label": shared_short_label,
                    "program_color": shared_color,
                    "marker_role": shared_role or str(chosen["source_marker_role"]),
                    "chosen_source_slot": int(chosen["source_slot"]),
                }
            )
        shared_position = shared.get("position", "first")

    program_blocks: list[tuple[tuple[str, str], list[dict[str, Any]]]] = []
    program_lookup: dict[tuple[str, str], int] = {}
    for row in assignments:
        key = (str(row["category"]), str(row["program"]))
        if key not in program_lookup:
            program_lookup[key] = len(program_blocks)
            program_blocks.append((key, []))
        program_blocks[program_lookup[key]][1].append(row)
    if shared_block:
        shared_key = (str(shared_block[0]["category"]), str(shared_block[0]["program"]))
        if shared_key in program_lookup:
            raise ValueError("shared program name collides with a selected program")
        if shared_position == "first":
            insert_at = 0
        elif shared_position == "last":
            insert_at = len(program_blocks)
        elif isinstance(shared_position, int) and 0 <= shared_position <= len(program_blocks):
            insert_at = shared_position
        else:
            raise ValueError("shared_program.position must be first, last, or a valid integer")
        program_blocks.insert(insert_at, (shared_key, shared_block))

    panel_rows: list[dict[str, Any]] = []
    provisional_order = 0
    for program_order, ((category, program), block) in enumerate(program_blocks):
        explicit_colors = {str(row.get("program_color", "")).upper() for row in block}
        explicit_colors.discard("")
        if len(explicit_colors) > 1:
            raise ValueError(f"program {program} declares inconsistent colors")
        program_color = (
            next(iter(explicit_colors))
            if explicit_colors
            else PROGRAM_PALETTE[program_order % len(PROGRAM_PALETTE)]
        )
        short_labels = {str(row.get("program_short_label", program)) for row in block}
        if len(short_labels) != 1:
            raise ValueError(f"program {program} declares inconsistent short labels")
        program_short_label = next(iter(short_labels))
        for program_gene_order, row in enumerate(block):
            panel_rows.append(
                {
                    "provisional_order": provisional_order,
                    "program_order": program_order,
                    "program": program,
                    "program_short_label": program_short_label,
                    "program_color": program_color,
                    "program_gene_order": program_gene_order,
                    "gene": row["gene"],
                    "marker_role": row["marker_role"],
                    "category": category,
                    "chosen_source_slot": int(row["chosen_source_slot"]),
                }
            )
            provisional_order += 1
    panel = pd.DataFrame(panel_rows)
    if panel.empty or panel["gene"].duplicated().any():
        raise ValueError("deduplicated marker panel is empty or non-unique")

    assignment_lookup = panel.set_index("gene").to_dict("index")
    slot_audit = slots.copy()
    decisions: list[str] = []
    assigned_programs: list[Any] = []
    assigned_orders: list[Any] = []
    for row, category_selected, program_selected in zip(
        slots.itertuples(index=False), category_ok.tolist(), program_ok.tolist()
    ):
        if not category_selected:
            decision = "excluded_category"
        elif not program_selected:
            decision = "excluded_program"
        else:
            target = assignment_lookup[row.gene]
            if int(target["chosen_source_slot"]) == int(row.source_slot):
                decision = "selected_primary_slot"
            elif target["program"] != row.source_program:
                decision = "deduplicated_to_shared_or_first_program"
            else:
                decision = "deduplicated_repeated_slot"
        decisions.append(decision)
        target = assignment_lookup.get(row.gene)
        assigned_programs.append(target["program"] if target is not None else pd.NA)
        assigned_orders.append(target["provisional_order"] if target is not None else pd.NA)
    slot_audit["slot_decision"] = decisions
    slot_audit["assigned_program"] = assigned_programs
    slot_audit["provisional_order"] = pd.array(assigned_orders, dtype="Int64")
    return panel, slot_audit


def scan_full_counts(
    h5ad_path: Path,
    layer: str,
    marker_var_indices: np.ndarray,
    group_codes: Mapping[str, np.ndarray],
    group_counts: Mapping[str, np.ndarray],
    chunk_rows: int,
    scale_factor: float,
) -> tuple[dict[str, dict[str, np.ndarray]], sparse.csr_matrix, sparse.csr_matrix, dict[str, Any]]:
    levels = list(group_codes)
    n_markers = len(marker_var_indices)
    sums = {
        level: np.zeros((len(group_counts[level]), n_markers), dtype=np.float64)
        for level in levels
    }
    detected = {
        level: np.zeros((len(group_counts[level]), n_markers), dtype=np.int64)
        for level in levels
    }
    count_chunks: list[sparse.csr_matrix] = []
    data_chunks: list[sparse.csr_matrix] = []
    audit: dict[str, Any] = {
        "negative_values": 0,
        "noninteger_values": 0,
        "nonfinite_values": 0,
        "stored_zero_values": 0,
        "empty_library_cells": 0,
        "nonzero_entries_scanned": 0,
        "max_raw_count": 0.0,
        "min_library_size": math.inf,
        "max_library_size": 0.0,
        "sum_library_size": 0.0,
    }
    with h5py.File(h5ad_path, "r") as handle:
        encoded = handle["layers"][layer]
        indptr = np.asarray(encoded["indptr"][:], dtype=np.int64)
        n_obs, n_vars = (int(value) for value in encoded.attrs["shape"])
        for level in levels:
            if len(group_codes[level]) != n_obs:
                raise ValueError(f"{level} group-code length differs from counts")
        print(
            f"[{RENDERER_NAME}] scanning {n_obs:,} cells × {n_vars:,} genes; "
            f"retaining {n_markers} measured candidate markers",
            flush=True,
        )
        for start in range(0, n_obs, chunk_rows):
            end = min(start + chunk_rows, n_obs)
            data_start, data_end = int(indptr[start]), int(indptr[end])
            data = np.asarray(encoded["data"][data_start:data_end], dtype=np.float64)
            indices = np.asarray(encoded["indices"][data_start:data_end], dtype=np.int64)
            local_indptr = indptr[start : end + 1] - data_start
            audit["negative_values"] += int(np.count_nonzero(data < 0))
            audit["noninteger_values"] += int(
                np.count_nonzero(np.abs(data - np.rint(data)) > 1e-8)
            )
            audit["nonfinite_values"] += int(np.count_nonzero(~np.isfinite(data)))
            audit["stored_zero_values"] += int(np.count_nonzero(data == 0))
            audit["nonzero_entries_scanned"] += int(len(data))
            if len(data):
                audit["max_raw_count"] = max(audit["max_raw_count"], float(data.max()))
            matrix = sparse.csr_matrix(
                (data, indices, local_indptr), shape=(end - start, n_vars)
            )
            library_size = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
            audit["empty_library_cells"] += int(np.count_nonzero(library_size <= 0))
            if len(library_size):
                audit["min_library_size"] = min(
                    audit["min_library_size"], float(library_size.min())
                )
                audit["max_library_size"] = max(
                    audit["max_library_size"], float(library_size.max())
                )
                audit["sum_library_size"] += float(library_size.sum(dtype=np.float64))
            if not np.isfinite(library_size).all() or np.any(library_size <= 0):
                raise ValueError("full-library sizes must be finite and positive")
            selected = matrix[:, marker_var_indices].astype(np.float64).tocsr()
            selected.eliminate_zeros()
            cp10k = selected.multiply((scale_factor / library_size)[:, None]).tocsr()
            lognormalized = cp10k.copy()
            lognormalized.data = np.log1p(lognormalized.data)
            count_chunks.append(selected)
            data_chunks.append(lognormalized)
            for level in levels:
                codes = group_codes[level][start:end]
                for code in np.unique(codes):
                    take = codes == code
                    sums[level][int(code)] += np.asarray(cp10k[take].sum(axis=0)).ravel()
                    detected[level][int(code)] += np.asarray(
                        selected[take].getnnz(axis=0)
                    ).ravel()
            if end == n_obs or end % 200_000 < chunk_rows:
                print(f"[{RENDERER_NAME}] processed {end:,}/{n_obs:,}", flush=True)
    aggregates: dict[str, dict[str, np.ndarray]] = {}
    for level in levels:
        n_cells = np.asarray(group_counts[level], dtype=np.int64)
        aggregates[level] = {
            "sum_cp10k": sums[level],
            "detected_cells": detected[level],
            "avg_exp_cp10k": sums[level] / n_cells[:, None],
            "pct_exp": 100.0 * detected[level] / n_cells[:, None],
        }
    marker_counts = sparse.vstack(count_chunks, format="csr")
    marker_data = sparse.vstack(data_chunks, format="csr")
    if marker_counts.shape != (len(next(iter(group_codes.values()))), n_markers):
        raise ValueError("retained marker matrix shape changed")
    if marker_data.shape != marker_counts.shape:
        raise ValueError("lognormalized marker matrix shape changed")
    return aggregates, marker_counts, marker_data, audit


def filter_marker_panel(
    provisional: pd.DataFrame,
    slot_audit: pd.DataFrame,
    var_lookup: Mapping[str, int],
    aggregates: Mapping[str, Mapping[str, np.ndarray]],
    measured_panel: pd.DataFrame,
    selection: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    threshold = float(selection["minimum_max_group_pct_exp"])
    filter_groupings = list(selection["filter_groupings"])
    measured_position = {
        gene: index for index, gene in enumerate(measured_panel["gene"].astype(str))
    }
    rows: list[dict[str, Any]] = []
    for marker in provisional.itertuples(index=False):
        available = marker.gene in var_lookup
        maxima: dict[str, float] = {}
        if available:
            index = measured_position[marker.gene]
            for level in aggregates:
                maxima[level] = float(aggregates[level]["pct_exp"][:, index].max())
            overall = max(maxima[level] for level in filter_groupings)
            if overall < threshold:
                decision = "excluded_below_max_pct_threshold"
                reason = (
                    f"maximum pct.exp across filter groupings is < {threshold:g}%"
                )
            else:
                decision = "retained"
                reason = (
                    f"maximum pct.exp across filter groupings is >= {threshold:g}%"
                )
        else:
            overall = np.nan
            decision = "excluded_unavailable_not_measured"
            reason = "absent from H5AD var_names; never filled with zero"
        row: dict[str, Any] = {
            "provisional_order": int(marker.provisional_order),
            "program_order": int(marker.program_order),
            "program": marker.program,
            "program_gene_order": int(marker.program_gene_order),
            "gene": marker.gene,
            "category": marker.category,
            "gene_available": available,
            "var_index": int(var_lookup[marker.gene]) if available else pd.NA,
            "max_pct_exp_across_filter_groupings": overall,
            "minimum_required_max_pct_exp": threshold,
            "gene_decision": decision,
            "decision_reason": reason,
        }
        for level in aggregates:
            row[f"{level}_max_pct_exp"] = maxima.get(level, np.nan)
        rows.append(row)
    gene_audit = pd.DataFrame(rows)
    retained = gene_audit.loc[gene_audit["gene_decision"].eq("retained"), "gene"].tolist()
    if not retained:
        raise ValueError("marker availability/detection filtering removed every gene")
    final_panel = provisional.loc[provisional["gene"].isin(retained)].copy()
    final_panel["final_order"] = range(len(final_panel))
    retained_programs = final_panel[["category", "program"]].drop_duplicates().reset_index(drop=True)
    program_order_lookup = {
        (row.category, row.program): index
        for index, row in enumerate(retained_programs.itertuples(index=False))
    }
    final_panel["program_order"] = [
        program_order_lookup[(row.category, row.program)]
        for row in final_panel.itertuples(index=False)
    ]
    final_panel["program_gene_order"] = final_panel.groupby(
        ["category", "program"], sort=False
    ).cumcount()
    program_sizes = final_panel.groupby("program_order", sort=False)["gene"].transform("size")
    program_starts = final_panel.groupby("program_order", sort=False)["final_order"].transform("min")
    program_ends = final_panel.groupby("program_order", sort=False)["final_order"].transform("max")
    final_panel["program_size"] = program_sizes.astype(int)
    final_panel["program_start_order"] = program_starts.astype(int)
    final_panel["program_end_order"] = program_ends.astype(int)
    final_panel["strip_xmin"] = final_panel["program_start_order"] + 0.5
    final_panel["strip_xmax"] = final_panel["program_end_order"] + 1.5
    final_panel["strip_x"] = (final_panel["strip_xmin"] + final_panel["strip_xmax"]) / 2.0
    final_panel = final_panel[
        [
            "final_order",
            "program_order",
            "program",
            "program_short_label",
            "program_color",
            "program_size",
            "program_start_order",
            "program_end_order",
            "strip_xmin",
            "strip_xmax",
            "strip_x",
            "program_gene_order",
            "gene",
            "marker_role",
            "category",
            "provisional_order",
            "chosen_source_slot",
        ]
    ]
    gene_lookup = gene_audit.set_index("gene")
    final_lookup = final_panel.set_index("gene")["final_order"].to_dict()
    slot_audit = slot_audit.copy()
    slot_audit["gene_decision"] = slot_audit["gene"].map(gene_lookup["gene_decision"])
    slot_audit["decision_reason"] = slot_audit["gene"].map(gene_lookup["decision_reason"])
    slot_audit["final_order"] = pd.array(
        slot_audit["gene"].map(final_lookup), dtype="Int64"
    )
    selection_excluded = slot_audit["slot_decision"].isin(
        ["excluded_category", "excluded_program"]
    )
    slot_audit.loc[selection_excluded, "gene_decision"] = "excluded_by_selection"
    slot_audit.loc[selection_excluded, "decision_reason"] = slot_audit.loc[
        selection_excluded, "slot_decision"
    ]
    return gene_audit, slot_audit, final_panel


def build_python_values(
    groups_by_level: Mapping[str, Sequence[str]],
    counts_by_level: Mapping[str, np.ndarray],
    final_panel: pd.DataFrame,
    measured_panel: pd.DataFrame,
    aggregates: Mapping[str, Mapping[str, np.ndarray]],
    dotplot: Mapping[str, Any],
    status: str,
    library_genes: int,
) -> pd.DataFrame:
    measured_position = {
        gene: index for index, gene in enumerate(measured_panel["gene"].astype(str))
    }
    rows: list[dict[str, Any]] = []
    for level, groups in groups_by_level.items():
        n_cells = np.asarray(counts_by_level[level], dtype=np.int64)
        for group_order, group in enumerate(groups):
            for marker in final_panel.itertuples(index=False):
                index = measured_position[marker.gene]
                avg = float(aggregates[level]["avg_exp_cp10k"][group_order, index])
                detected = int(aggregates[level]["detected_cells"][group_order, index])
                pct = float(aggregates[level]["pct_exp"][group_order, index])
                rows.append(
                    {
                        "result_status": status,
                        "category": marker.category,
                        "grouping_level": level,
                        "group_order": group_order,
                        "group_label": str(group),
                        "n_cells": int(n_cells[group_order]),
                        "final_order": int(marker.final_order),
                        "program_order": int(marker.program_order),
                        "program": marker.program,
                        "program_gene_order": int(marker.program_gene_order),
                        "gene": marker.gene,
                        "detected_cells": detected,
                        "avg_exp_cp10k": avg,
                        "avg_exp_log1p": float(np.log1p(avg)),
                        "avg_exp_zscore_unclipped": np.nan,
                        "avg_exp_scaled": np.nan,
                        "pct_exp": pct,
                        "count_source": "layers[counts]",
                        "library_size_source": "complete_layers_counts_library",
                        "library_genes": int(library_genes),
                        "normalization_scale_factor": float(
                            dotplot["normalization_scale_factor"]
                        ),
                        "scaling_ddof": int(dotplot["scaling_ddof"]),
                        "scaling_clip_min": float(dotplot["col_min"]),
                        "scaling_clip_max": float(dotplot["col_max"]),
                        "dot_scale_by": dotplot["scale_by"],
                        "dot_scale": float(dotplot["dot_scale"]),
                    }
                )
    values = pd.DataFrame(rows)
    for (level, gene), index in values.groupby(
        ["grouping_level", "gene"], sort=False
    ).groups.items():
        ordered_index = values.loc[index].sort_values("group_order").index
        expression = values.loc[ordered_index, "avg_exp_log1p"].to_numpy(np.float64)
        if len(expression) == 1:
            values.loc[ordered_index, "avg_exp_scaled"] = expression
            continue
        standard_deviation = float(expression.std(ddof=int(dotplot["scaling_ddof"])))
        if not np.isfinite(standard_deviation) or standard_deviation <= 0:
            continue
        zscore = (expression - expression.mean()) / standard_deviation
        values.loc[ordered_index, "avg_exp_zscore_unclipped"] = zscore
        values.loc[ordered_index, "avg_exp_scaled"] = np.clip(
            zscore, float(dotplot["col_min"]), float(dotplot["col_max"])
        )
    values["expression_scaling"] = values["grouping_level"].map({k: "native_log1p_mean_unscaled" if len(v) == 1 else "sample_zscore_ddof1" for k, v in groups_by_level.items()})
    return values


def resolve_rscript(value: str) -> Path:
    candidate = Path(value).expanduser()
    if candidate.is_absolute() or len(candidate.parts) > 1:
        resolved = candidate.resolve()
    else:
        found = shutil.which(value)
        if found is None:
            raise FileNotFoundError(f"configured Rscript is not on PATH: {value}")
        resolved = Path(found).resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise FileNotFoundError(f"configured Rscript is not executable: {resolved}")
    return resolved


def run_seurat_renderer(
    marker_counts: sparse.csr_matrix,
    marker_data: sparse.csr_matrix,
    obs_names: Sequence[str],
    group_codes: Mapping[str, np.ndarray],
    groups_by_level: Mapping[str, Sequence[str]],
    groupings: Mapping[str, Mapping[str, Any]],
    final_panel: pd.DataFrame,
    oracle_config: Mapping[str, Any],
    log_path: Path,
    temporary_root: Path | None,
    native_dir: Path,
    figure_dir: Path,
    width_inches: float,
    height_inches: float,
    dpi: int,
) -> tuple[pd.DataFrame, dict[str, str]]:
    rscript = resolve_rscript(str(oracle_config["rscript"]))
    oracle_script = Path(oracle_config["_script"]).resolve()
    genes = final_panel["gene"].astype(str).tolist()
    if marker_counts.shape != marker_data.shape or marker_counts.shape[1] != len(genes):
        raise ValueError("Seurat renderer matrix dimensions are inconsistent")
    if marker_counts.shape[0] != len(obs_names):
        raise ValueError("Seurat renderer obs_name length is inconsistent")
    temporary_parent = None if temporary_root is None else str(temporary_root.resolve())
    with tempfile.TemporaryDirectory(
        prefix=f"{RENDERER_NAME}_native_seurat.", dir=temporary_parent
    ) as temporary:
        temp = Path(temporary)
        counts_path = temp / "marker_counts_gene_by_cell.mtx"
        data_path = temp / "marker_data_log1p_cp10k_gene_by_cell.mtx"
        marker_contract_path = temp / "marker_contract.csv"
        cells_path = temp / "cells.csv"
        levels_path = temp / "group_levels.csv"
        print(f"[{RENDERER_NAME}] writing temporary Seurat matrices", flush=True)
        mmwrite(counts_path, marker_counts.T.tocoo(), symmetry="general")
        mmwrite(data_path, marker_data.T.tocoo(), symmetry="general")
        r_panel = final_panel.rename(
            columns={
                "program": "marker_program",
                "program_short_label": "marker_program_short_label",
                "program_color": "marker_program_color",
                "program_size": "marker_program_size",
                "marker_role": "marker_role",
                "final_order": "marker_order",
            }
        )
        r_panel.to_csv(marker_contract_path, index=False)
        cell_table: dict[str, Any] = {"obs_name": [str(value) for value in obs_names]}
        level_rows: list[dict[str, Any]] = []
        for level, groups in groups_by_level.items():
            codes = np.asarray(group_codes[level], dtype=np.int32)
            labels = np.asarray(list(groups), dtype=object)
            cell_table[level] = labels[codes]
            for order, group in enumerate(groups):
                level_rows.append(
                    {
                        "grouping_level": level,
                        "group_order": order,
                        "group_label": group,
                        "group_title": str(groupings[level]["title"]),
                        "group_subtitle": str(groupings[level].get("subtitle", "")),
                        "height_inches": groupings[level].get("height_inches", height_inches),
                        **{k: v for k, v in groupings[level].get("style", STYLE_DEFAULTS).items() if k != "program_palette"},
                    }
                )
        pd.DataFrame(cell_table).to_csv(cells_path, index=False)
        pd.DataFrame(level_rows).to_csv(levels_path, index=False)
        command = [
            str(rscript),
            "--vanilla",
            str(oracle_script),
            str(counts_path),
            str(data_path),
            str(marker_contract_path),
            str(cells_path),
            str(levels_path),
            str(native_dir),
            str(figure_dir),
            f"{width_inches:.17g}",
            f"{height_inches:.17g}",
            str(int(dpi)),
        ]
        environment = os.environ.copy()
        environment.update(
            {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        )
        print(f"[{RENDERER_NAME}] invoking native Seurat DotPlot renderer", flush=True)
        completed = subprocess.run(
            command,
            cwd=oracle_script.parent,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            "command: "
            + " ".join(command[:3])
            + " <temporary real-cell inputs>\n"
            + f"returncode: {completed.returncode}\n"
            + f"stdout:\n{completed.stdout}\n"
            + f"stderr:\n{completed.stderr}\n",
            encoding="utf-8",
        )
        if completed.returncode != 0:
            raise RuntimeError("native Seurat renderer failed; see logs/seurat_renderer.log")

    environment_path = native_dir / "seurat_render_environment.csv"
    if not environment_path.is_file():
        raise ValueError("native Seurat renderer did not write its environment record")
    environment = pd.read_csv(environment_path, dtype=str, keep_default_na=False)
    if len(environment) != 1:
        raise ValueError("unexpected Seurat render environment record")
    versions = {
        "r_version": str(environment.loc[0, "r_version"]),
        "seurat_version": str(environment.loc[0, "seurat_version"]),
        "seuratobject_version": str(environment.loc[0, "seuratobject_version"]),
        "patchwork_version": str(environment.loc[0, "patchwork_version"]),
    }

    oracle_parts: list[pd.DataFrame] = []
    required_native = {"id", "features.plot", "avg.exp", "avg.exp.scaled", "pct.exp"}
    for level, groups in groups_by_level.items():
        path = native_dir / f"{level}_seurat_native_pdata.csv"
        figure_path = figure_dir / f"{level}_seurat_dotplot.png"
        if not path.is_file() or not figure_path.is_file():
            raise ValueError(f"native Seurat output is missing for grouping {level}")
        native = pd.read_csv(path, low_memory=False)
        if not required_native.issubset(native.columns):
            raise ValueError(f"native Seurat p$data schema changed for {level}")
        if len(native) != len(groups) * len(genes):
            raise ValueError(f"native Seurat p$data row count changed for {level}")
        oracle = pd.DataFrame(
            {
                "grouping_level": level,
                "group_label": native["id"].astype(str),
                "gene": native["features.plot"].astype(str),
                "seurat_avg_exp_cp10k": native["avg.exp"].astype(float),
                "seurat_avg_exp_log1p": np.log1p(native["avg.exp"].astype(float)),
                "seurat_avg_exp_zscore_unclipped": np.nan,
                "seurat_avg_exp_scaled": native["avg.exp.scaled"].astype(float),
                "seurat_pct_exp": native["pct.exp"].astype(float),
                "seurat_version": versions["seurat_version"],
                "seuratobject_version": versions["seuratobject_version"],
                "r_version": versions["r_version"],
            }
        )
        for gene, index in oracle.groupby("gene", sort=False).groups.items():
            values = oracle.loc[index, "seurat_avg_exp_log1p"].to_numpy(np.float64)
            if len(values) == 1:
                continue
            standard_deviation = float(values.std(ddof=1))
            if not np.isfinite(standard_deviation) or standard_deviation <= 0:
                continue
            oracle.loc[index, "seurat_avg_exp_zscore_unclipped"] = (
                values - values.mean()
            ) / standard_deviation
        oracle_parts.append(oracle)
    oracle = pd.concat(oracle_parts, ignore_index=True)
    expected_rows = sum(len(groups) for groups in groups_by_level.values()) * len(genes)
    if len(oracle) != expected_rows or oracle.duplicated(
        ["grouping_level", "group_label", "gene"]
    ).any():
        raise ValueError("native Seurat p$data keys are incomplete or non-unique")
    return oracle, versions


def build_parity(
    python_values: pd.DataFrame,
    oracle: pd.DataFrame,
    oracle_config: Mapping[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    parity = python_values.merge(
        oracle,
        on=["grouping_level", "group_label", "gene"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    if len(parity) != len(python_values) or not parity["_merge"].eq("both").all():
        raise ValueError("Python and Seurat keys do not match exactly")
    parity = parity.drop(columns="_merge")
    pairs = {
        "avg_exp_cp10k": ("seurat_avg_exp_cp10k", "avg_exp_abs_tolerance"),
        "avg_exp_log1p": ("seurat_avg_exp_log1p", "avg_exp_log1p_abs_tolerance"),
        "avg_exp_zscore_unclipped": (
            "seurat_avg_exp_zscore_unclipped",
            "zscore_abs_tolerance",
        ),
        "avg_exp_scaled": ("seurat_avg_exp_scaled", "scaled_abs_tolerance"),
        "pct_exp": ("seurat_pct_exp", "pct_exp_abs_tolerance"),
    }
    pass_columns: list[str] = []
    summary_rows: list[dict[str, Any]] = []
    for python_column, (seurat_column, tolerance_key) in pairs.items():
        difference_column = f"abs_diff_{python_column}"
        pass_column = f"pass_{python_column}"
        parity[difference_column] = np.abs(
            parity[python_column].astype(float) - parity[seurat_column].astype(float)
        )
        tolerance = float(oracle_config[tolerance_key])
        parity[pass_column] = parity[difference_column].le(tolerance) | (parity[python_column].isna() & parity[seurat_column].isna())
        pass_columns.append(pass_column)
        maximum = float(parity[difference_column].max()) if parity[difference_column].notna().any() else 0.0
        summary_rows.append(
            {
                "metric": python_column,
                "max_abs_difference": maximum,
                "absolute_tolerance": tolerance,
                "status": "PASS" if parity[pass_column].all() else "FAIL",
            }
        )
    parity["row_pass"] = parity[pass_columns].all(axis=1)
    parity["avg.exp"] = parity["seurat_avg_exp_cp10k"]
    parity["avg.exp.log1p"] = parity["seurat_avg_exp_log1p"]
    parity["avg.exp.scaled"] = parity["seurat_avg_exp_scaled"]
    parity["pct.exp"] = parity["seurat_pct_exp"]
    parity = parity.sort_values(
        ["grouping_level", "group_order", "final_order"]
    ).reset_index(drop=True)
    versions = {
        "r_version": str(oracle["r_version"].iloc[0]),
        "seurat_version": str(oracle["seurat_version"].iloc[0]),
        "seuratobject_version": str(oracle["seuratobject_version"].iloc[0]),
    }
    expected_versions = {
        "r_version": str(oracle_config["expected_r_version"]),
        "seurat_version": str(oracle_config["expected_seurat_version"]),
        "seuratobject_version": str(oracle_config["expected_seuratobject_version"]),
    }
    for key, observed in versions.items():
        source_column = key
        if oracle[source_column].astype(str).nunique() != 1 or observed != expected_versions[key]:
            raise ValueError(
                f"oracle version mismatch for {key}: observed={observed}, "
                f"expected={expected_versions[key]}"
            )
    summary = pd.DataFrame(summary_rows)
    if not summary["status"].eq("PASS").all() or not parity["row_pass"].all():
        raise ValueError("Python–Seurat numeric parity failed")
    return parity, summary, versions


def image_info(path: Path) -> dict[str, Any]:
    with Image.open(path) as image:
        image.load()
        rgba = np.asarray(image.convert("RGBA"))
        return {
            "format": str(image.format),
            "width": int(image.width),
            "height": int(image.height),
            "bytes": int(path.stat().st_size),
            "corner_rgba": [int(value) for value in rgba[0, 0]],
        }


def program_strip_audit(panel: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "program_order",
        "program",
        "program_short_label",
        "program_color",
        "program_size",
        "program_start_order",
        "program_end_order",
        "strip_xmin",
        "strip_xmax",
        "strip_x",
    ]
    audit = panel[columns].drop_duplicates("program_order").sort_values("program_order")
    expected = np.arange(len(audit), dtype=np.int64)
    if not np.array_equal(audit["program_order"].to_numpy(np.int64), expected):
        raise ValueError("filtered marker-program order is not contiguous")
    if int(audit["program_size"].sum()) != len(panel):
        raise ValueError("filtered marker-program strip lengths do not cover the panel")
    if not np.allclose(
        audit["strip_xmax"].to_numpy(float) - audit["strip_xmin"].to_numpy(float),
        audit["program_size"].to_numpy(float),
    ):
        raise ValueError("filtered marker-program strip geometry is inconsistent")
    return audit.reset_index(drop=True)


def native_layout_contract(
    level: str,
    values: pd.DataFrame,
    panel: pd.DataFrame,
    strip_audit: pd.DataFrame,
    width_contract: Mapping[str, Any],
) -> dict[str, Any]:
    groups = (
        values[["group_order", "group_label", "n_cells"]]
        .drop_duplicates()
        .sort_values("group_order")
    )
    return {
        "grouping_level": level,
        "x_axis_index": panel["final_order"].astype(int).tolist(),
        "x_axis_genes": panel["gene"].astype(str).tolist(),
        "x_axis_text": {
            "angle_degrees": X_AXIS_TEXT_ANGLE_DEGREES,
            "hjust": X_AXIS_TEXT_HJUST,
            "vjust": X_AXIS_TEXT_VJUST,
            "font_size_pt": X_AXIS_TEXT_FONT_SIZE_PT,
        },
        "width_contract": dict(width_contract),
        "y_axis_index": groups["group_order"].astype(int).tolist(),
        "y_axis_labels": groups["group_label"].astype(str).tolist(),
        "program_spans": strip_audit.to_dict("records"),
        "size_legend_title": "Percent Expressed",
        "size_semantics": "native Seurat DotPlot scale.by=radius and dot.scale=6 defaults",
        "color_legend_title": "Average Expression",
        "render_palette": ["#2166AC", "#FFFFFF", "#B2182B"],
        "render_midpoint": 0.0,
        "render_limits": [-2.5, 2.5],
        "annotation_strip_legend": False,
    }


def write_artifact_reopen(stage: Path, path: Path, artifact_type: str) -> dict[str, Any]:
    relative = path.relative_to(stage).as_posix()
    if artifact_type == "png":
        info = image_info(path)
        return {"path": relative, "artifact_type": artifact_type, "status": "PASS", **info}
    if artifact_type == "csv":
        rows = len(pd.read_csv(path, low_memory=False))
        return {"path": relative, "artifact_type": artifact_type, "status": "PASS", "rows": rows}
    if artifact_type == "csv_gz":
        rows = len(pd.read_csv(path, compression="gzip", low_memory=False))
        return {"path": relative, "artifact_type": artifact_type, "status": "PASS", "rows": rows}
    if artifact_type == "json":
        json.loads(path.read_text(encoding="utf-8"))
        return {"path": relative, "artifact_type": artifact_type, "status": "PASS"}
    raise ValueError(f"unsupported reopen artifact type: {artifact_type}")


def protected_input_paths(config: Mapping[str, Any]) -> list[Path]:
    root = Path(config["_project_root"])
    inputs = config["inputs"]
    oracle = config["seurat_native"]
    paths = [
        Path(__file__).resolve(),
        Path(__file__).resolve().parents[2] / "_shared/dotplot/engine.py",
        Path(__file__).resolve().parents[2] / "_shared/dotplot/annotation_contract.py",
        Path(config["_config_path"]),
        Path(inputs["_counts_h5ad"]),
        Path(inputs["_sidecar"]),
        Path(inputs["_marker_yaml"]),
        Path(oracle["_script"]),
    ]
    project = config.get("project", {})
    for key in ("annotation_guide", "annotation_registry"):
        if project.get(key):
            paths.append(resolve_path(root, project[key]))
    for value in project.get("protected_files", []):
        paths.append(resolve_path(root, value))
    for track in project.get("_resolved_annotation_roles", {}).values():
        for source in track.get("sources", []):
            paths.append(resolve_path(root, source["path"]))
    for binding in project.get("_annotation_bindings", {}).values():
        paths.append(Path(binding["path"]))
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            if not resolved.is_file():
                raise FileNotFoundError(f"protected input is missing: {resolved}")
            seen.add(resolved)
            unique.append(resolved)
    return unique


def input_hash_table(paths: Sequence[Path], root: Path) -> pd.DataFrame:
    return pd.DataFrame([file_record(path, root) for path in paths])


def check_declared_input_hashes(config: Mapping[str, Any], before: pd.DataFrame) -> None:
    inputs = config["inputs"]
    lookup = before.set_index("absolute_path")["sha256"].to_dict()
    for key, path_key in (
        ("counts_h5ad_sha256", "_counts_h5ad"),
        ("sidecar_sha256", "_sidecar"),
        ("marker_yaml_sha256", "_marker_yaml"),
    ):
        if inputs.get(key):
            observed = lookup[str(Path(inputs[path_key]).resolve())]
            expected = str(inputs[key])
            if observed != expected:
                raise ValueError(f"declared SHA256 mismatch for inputs.{key}")


def run_recluster_guard(config: Mapping[str, Any], log_path: Path) -> dict[str, Any]:
    root = Path(config["_project_root"])
    guard_contract = root / "metadata" / "recluster_guard.json"
    project = config.get("project", {})
    enabled = bool(project.get("run_recluster_guard", guard_contract.is_file()))
    if guard_contract.is_file() and not enabled:
        reason = str(project.get("recluster_guard_disabled_reason", "")).strip()
        if not reason:
            raise ValueError("disabling a present recluster guard requires a non-empty reason")
        result = {"status": "NOT_RUN_WITH_REASON", "reason": reason, "returncode": None}
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result
    if not enabled:
        result = {"status": "NOT_APPLICABLE", "reason": "no project guard requested", "returncode": None}
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result
    command = [sys.executable, "-m", "xenium_downstream.recluster_guard", "--root", str(root)]
    completed = subprocess.run(
        command,
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        f"command: {' '.join(command)}\nreturncode: {completed.returncode}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}\n",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError("repository recluster guard failed")
    return {"status": "PASS", "reason": "project guard completed", "returncode": 0}


def environment_record(versions: Mapping[str, str]) -> dict[str, Any]:
    packages: dict[str, str] = {}
    for name in ("h5py", "numpy", "pandas", "Pillow", "PyYAML", "scipy"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "not_installed"
    return {
        "generated_at": now_local(),
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": packages,
        "oracle": dict(versions),
    }


def write_output_manifest(output: Path) -> None:
    excluded = {MANIFEST_REL.as_posix(), MANIFEST_HASH_REL.as_posix()}
    rows: list[dict[str, Any]] = []
    for path in sorted(item for item in output.rglob("*") if item.is_file()):
        relative = path.relative_to(output).as_posix()
        if relative in excluded:
            continue
        rows.append(
            {"path": relative, "bytes": int(path.stat().st_size), "sha256": sha256_file(path)}
        )
    manifest_path = output / MANIFEST_REL
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(manifest_path, index=False)
    digest = sha256_file(manifest_path)
    (output / MANIFEST_HASH_REL).write_text(
        f"{digest}  {MANIFEST_REL.name}\n", encoding="utf-8"
    )


def verify_output_package(output: Path, verify_inputs: bool) -> dict[str, Any]:
    output = output.resolve()
    failures: list[str] = []
    manifest_path = output / MANIFEST_REL
    manifest_hash_path = output / MANIFEST_HASH_REL
    if not manifest_path.is_file() or not manifest_hash_path.is_file():
        return {"status": "FAIL", "failures": ["manifest_missing"], "output": str(output)}
    expected_manifest_hash = manifest_hash_path.read_text(encoding="utf-8").split()[0]
    if sha256_file(manifest_path) != expected_manifest_hash:
        failures.append("manifest_sha256_mismatch")
    manifest = pd.read_csv(manifest_path)
    expected_paths = set(manifest["path"].astype(str))
    actual_paths = {
        path.relative_to(output).as_posix()
        for path in output.rglob("*")
        if path.is_file()
        and path.relative_to(output).as_posix()
        not in {MANIFEST_REL.as_posix(), MANIFEST_HASH_REL.as_posix()}
    }
    if expected_paths != actual_paths:
        failures.append("manifest_file_set_mismatch")
    for row in manifest.itertuples(index=False):
        path = output / str(row.path)
        if not path.is_file():
            failures.append(f"missing:{row.path}")
        elif int(path.stat().st_size) != int(row.bytes) or sha256_file(path) != str(row.sha256):
            failures.append(f"hash_or_size:{row.path}")
    pdfs = list(output.rglob("*.pdf"))
    if pdfs:
        failures.append(f"pdf_count:{len(pdfs)}")
    pngs = list(output.rglob("*.png"))
    if not pngs:
        failures.append("png_missing")
    for path in pngs:
        try:
            info = image_info(path)
            if info["format"] != "PNG" or info["bytes"] <= 0:
                failures.append(f"invalid_png:{path.relative_to(output)}")
            if info["corner_rgba"] != [255, 255, 255, 255]:
                failures.append(f"nonwhite_png_background:{path.relative_to(output)}")
        except Exception:
            failures.append(f"invalid_png:{path.relative_to(output)}")
    forbidden = [
        path
        for path in output.rglob("*")
        if path.is_file()
        and (
            path.suffix.lower() in {".h5ad", ".h5"}
            or "sidecar" in path.name.casefold()
            or "umap" in path.name.casefold()
        )
    ]
    if forbidden:
        failures.append(f"forbidden_derivative_artifact_count:{len(forbidden)}")
    completion_path = output / "completion.json"
    if not completion_path.is_file() or json.loads(completion_path.read_text())["status"] != "PASS":
        failures.append("completion_not_pass")
    validation_path = output / "03_validation/validation_summary.csv"
    if not validation_path.is_file() or not pd.read_csv(validation_path)["status"].eq("PASS").all():
        failures.append("validation_not_pass")
    if verify_inputs:
        contract_path = output / "00_contract/run_contract.json"
        if not contract_path.is_file():
            failures.append("run_contract_missing")
        else:
            contract = json.loads(contract_path.read_text(encoding="utf-8"))
            for record in contract.get("protected_inputs", []):
                path = Path(record["absolute_path"])
                if not path.is_file():
                    failures.append(f"live_input_missing:{path}")
                elif int(path.stat().st_size) != int(record["bytes"]) or sha256_file(path) != record["sha256"]:
                    failures.append(f"live_input_changed:{path}")
    return {
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "output": str(output),
        "manifest_rows": int(len(manifest)),
        "png_count": len(pngs),
        "pdf_count": len(pdfs),
        "forbidden_derivative_artifact_count": len(forbidden),
        "live_inputs_verified": bool(verify_inputs),
    }


def run_build(args: argparse.Namespace) -> Path:
    config = prepare_config(args.config)
    root = Path(config["_project_root"])
    output = resolve_path(root, args.output) if args.output else Path(config["output"]["_root"])
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"output already exists; choose a new derivative directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    checks: list[dict[str, str]] = []
    try:
        protected = protected_input_paths(config)
        before = input_hash_table(protected, root)
        check_declared_input_hashes(config, before)
        (stage / "00_contract").mkdir(parents=True, exist_ok=True)
        before.to_csv(stage / "00_contract/input_sha256_before.csv", index=False)
        write_json(
            stage / "00_contract/annotation_resolution.json",
            {
                "use_case": config["project"]["annotation_use_case"],
                "track_id": config["project"]["_resolved_annotation_track"].get(
                    "track_id"
                ),
                "status": config["project"]["_resolved_annotation_track"].get("status"),
                "registry": str(config["project"]["_annotation_registry"]),
                "guide": str(config["project"]["_annotation_guide"]),
                "selection_rule": "one_use_case_per_semantic_role",
                "bindings": config["project"]["_annotation_bindings"],
                "resolved_roles": config["project"]["_resolved_annotation_roles"],
                "annotation_changed": False,
            },
        )

        inputs = config["inputs"]
        obs_names, var_names, shape, encoding = read_h5ad_contract(
            Path(inputs["_counts_h5ad"]), str(inputs.get("counts_layer", "counts"))
        )
        expected_cells = inputs.get("expected_cells")
        expected_genes = inputs.get("expected_genes")
        if expected_cells is not None and shape[0] != int(expected_cells):
            raise ValueError(f"H5AD cell count changed: {shape[0]} != {expected_cells}")
        if expected_genes is not None and shape[1] != int(expected_genes):
            raise ValueError(f"H5AD gene count changed: {shape[1]} != {expected_genes}")
        add_check(checks, "h5ad_counts_encoding", encoding, "csr_matrix", encoding == "csr_matrix")
        add_check(checks, "h5ad_obs_name_unique", obs_names.is_unique, True, obs_names.is_unique)
        add_check(checks, "h5ad_var_name_unique", len(var_names) == len(set(var_names)), True, len(var_names) == len(set(var_names)))

        grouping_fields = list(
            dict.fromkeys(
                str(block["field"]) for block in config["groupings"].values()
            )
        )
        sidecar_key = str(config["join"]["sidecar_key"])
        sidecar = read_sidecar(
            Path(inputs["_sidecar"]), [sidecar_key, *grouping_fields]
        )
        aligned = align_sidecar_to_obs_names(
            obs_names, sidecar, grouping_fields, sidecar_key
        )
        source_checks = validate_annotation_sources(config, aligned)
        write_json(stage / "00_contract/annotation_source_validation.json", {"sources": source_checks})
        shuffled = align_sidecar_to_obs_names(
            obs_names,
            sidecar.sample(frac=1.0, random_state=17),
            grouping_fields,
            sidecar_key,
        )
        shuffle_invariant = aligned.equals(shuffled)
        write_json(
            stage / "00_contract/join_contract.json",
            {
                "modality": MODALITY,
                "sidecar_key": sidecar_key,
                "target_key": "obs_name",
                "sidecar_key_unique": True,
                "target_key_unique": True,
                "exact_set_equality": True,
                "row_order_join_forbidden": True,
                "cell_id_join_forbidden": MODALITY == "xenium",
                "registry_join_source_ids": config["project"][
                    "_registry_join_source_ids"
                ],
                "shuffle_invariant_verified": bool(shuffle_invariant),
                "rows": int(len(aligned)),
            },
        )
        add_check(
            checks,
            f"{sidecar_key}_to_obs_name_exact_set_join",
            len(aligned),
            shape[0],
            len(aligned) == shape[0],
        )
        add_check(
            checks,
            f"{sidecar_key}_to_obs_name_join_shuffle_invariant",
            shuffle_invariant,
            True,
            shuffle_invariant,
        )
        codes, group_counts, groups, group_table = build_group_contract(
            aligned, config["groupings"]
        )
        for level in groups:
            add_check(
                checks,
                f"{level}_total_cells",
                int(group_counts[level].sum()),
                shape[0],
                int(group_counts[level].sum()) == shape[0],
            )

        slots = expand_marker_slots(Path(inputs["_marker_yaml"]))
        provisional, slot_audit = select_and_deduplicate_markers(
            slots, config["marker_selection"]
        )
        var_lookup = {gene: index for index, gene in enumerate(var_names)}
        measured_panel = provisional.loc[provisional["gene"].isin(var_lookup)].copy()
        if measured_panel.empty:
            raise ValueError("none of the selected markers are measured")
        marker_indices = np.asarray(
            [var_lookup[gene] for gene in measured_panel["gene"]], dtype=np.int64
        )
        aggregates, marker_counts, marker_data, count_audit = scan_full_counts(
            Path(inputs["_counts_h5ad"]),
            str(inputs.get("counts_layer", "counts")),
            marker_indices,
            codes,
            group_counts,
            args.chunk_rows,
            float(config["dotplot"]["normalization_scale_factor"]),
        )
        for key in (
            "negative_values",
            "noninteger_values",
            "nonfinite_values",
            "stored_zero_values",
            "empty_library_cells",
        ):
            add_check(checks, f"counts_{key}", count_audit[key], 0, count_audit[key] == 0)
        if any(count_audit[key] != 0 for key in (
            "negative_values", "noninteger_values", "nonfinite_values", "stored_zero_values", "empty_library_cells"
        )):
            raise ValueError("raw-count integrity audit failed")

        gene_audit, slot_audit, final_panel = filter_marker_panel(
            provisional,
            slot_audit,
            var_lookup,
            aggregates,
            measured_panel,
            config["marker_selection"],
        )
        add_check(checks, "final_marker_unique", final_panel["gene"].is_unique, True, final_panel["gene"].is_unique)
        unavailable_count = int(gene_audit["gene_decision"].eq("excluded_unavailable_not_measured").sum())
        add_check(
            checks,
            "unavailable_markers_not_filled_with_zero",
            unavailable_count,
            "audited",
            True,
            "Every unavailable marker remains in marker_filter_gene_audit.csv",
        )
        measured_position = {
            gene: index for index, gene in enumerate(measured_panel["gene"].astype(str))
        }
        final_positions = np.asarray(
            [measured_position[gene] for gene in final_panel["gene"]], dtype=np.int64
        )
        final_counts = marker_counts[:, final_positions].tocsr()
        final_data = marker_data[:, final_positions].tocsr()
        python_values = build_python_values(
            groups,
            group_counts,
            final_panel,
            measured_panel,
            aggregates,
            config["dotplot"],
            config["status"],
            shape[1],
        )
        numeric_dir = stage / "02_numeric_sources"
        native_dir = numeric_dir / "native_seurat"
        validation_dir = stage / "03_validation"
        figure_dir = stage / "01_figures"
        numeric_dir.mkdir(parents=True, exist_ok=True)
        native_dir.mkdir(parents=True, exist_ok=True)
        validation_dir.mkdir(parents=True, exist_ok=True)
        figure_dir.mkdir(parents=True, exist_ok=True)
        strip_audit = program_strip_audit(final_panel)
        render_dpi = int(args.dpi or config["output"]["dpi"])
        width_contract = resolve_figure_width(
            float(config["output"]["width_inches"]),
            final_panel["gene"].astype(str).tolist(),
            bool(config["output"]["auto_expand_width"]),
        )
        effective_width_inches = float(width_contract["final_width_inches"])
        add_check(
            checks,
            "x_axis_label_rotation_degrees",
            X_AXIS_TEXT_ANGLE_DEGREES,
            60.0,
            X_AXIS_TEXT_ANGLE_DEGREES == 60.0,
        )
        add_check(
            checks,
            "x_axis_label_justification",
            [X_AXIS_TEXT_HJUST, X_AXIS_TEXT_VJUST],
            [1.0, 1.0],
            X_AXIS_TEXT_HJUST == 1.0 and X_AXIS_TEXT_VJUST == 1.0,
        )
        add_check(
            checks,
            "x_axis_label_font_size_pt",
            config["style"]["gene_font_pt"],
            "configured positive font size",
            config["style"]["gene_font_pt"] > 0,
        )
        add_check(
            checks,
            "figure_width_layout_contract",
            width_contract,
            "final width satisfies the computed marker-label minimum",
            effective_width_inches
            >= float(width_contract["computed_minimum_width_inches"]),
        )
        oracle, renderer_versions = run_seurat_renderer(
            final_counts,
            final_data,
            obs_names.astype(str).tolist(),
            codes,
            groups,
            config["groupings"],
            final_panel,
            config["seurat_native"],
            stage / "logs/seurat_renderer.log",
            args.oracle_temp_dir,
            native_dir,
            figure_dir,
            effective_width_inches,
            float(config["output"]["height_inches"] or 3.2),
            render_dpi,
        )
        parity, parity_summary, versions = build_parity(
            python_values, oracle, config["seurat_native"]
        )
        versions["patchwork_version"] = renderer_versions["patchwork_version"]
        expected_versions = {
            "r_version": str(config["seurat_native"]["expected_r_version"]),
            "seurat_version": str(config["seurat_native"]["expected_seurat_version"]),
            "seuratobject_version": str(
                config["seurat_native"]["expected_seuratobject_version"]
            ),
        }
        observed_versions = {key: versions[key] for key in expected_versions}
        add_check(
            checks,
            "seurat_environment_exact",
            observed_versions,
            expected_versions,
            observed_versions == expected_versions,
        )
        add_check(checks, "python_seurat_all_rows_pass", int(parity["row_pass"].sum()), len(parity), bool(parity["row_pass"].all()))

        final_panel.to_csv(numeric_dir / "final_marker_panel.csv", index=False)
        strip_audit.to_csv(numeric_dir / "marker_program_strip_audit.csv", index=False)
        gene_audit.to_csv(numeric_dir / "marker_filter_gene_audit.csv", index=False, float_format="%.12g")
        slot_audit.to_csv(numeric_dir / "marker_filter_slot_audit.csv", index=False, float_format="%.12g")
        group_table.to_csv(numeric_dir / "group_counts.csv", index=False)
        parity.to_csv(
            validation_dir / "python_seurat_dotplot_parity.csv.gz",
            index=False,
            compression="gzip",
            float_format="%.17g",
        )
        parity_summary.to_csv(
            validation_dir / "python_seurat_dotplot_parity_summary.csv", index=False
        )

        figure_rows: list[dict[str, Any]] = []
        layouts: dict[str, Any] = {}
        # R's png device truncates non-integer pixel dimensions supplied by
        # ggsave (for example, 23.22 in * 240 dpi -> 5572 px).
        expected_width = int(effective_width_inches * render_dpi)
        reopen_rows: list[dict[str, Any]] = []
        for level, block in config["groupings"].items():
            expected_height = int(block["height_inches"] * render_dpi)
            values = parity.loc[parity["grouping_level"].eq(level)].copy()
            numeric_rel = Path(f"02_numeric_sources/{level}_seurat_dotplot_values.csv.gz")
            figure_rel = Path(f"01_figures/{level}_seurat_dotplot.png")
            native_rel = Path(
                f"02_numeric_sources/native_seurat/{level}_seurat_native_pdata.csv"
            )
            values.to_csv(
                stage / numeric_rel,
                index=False,
                compression="gzip",
                float_format="%.17g",
            )
            layouts[level] = native_layout_contract(
                level, values, final_panel, strip_audit, width_contract
            )
            layouts[level]["style"] = config["style"]
            layouts[level]["height_inches"] = block["height_inches"]
            layouts[level]["title"] = block["title"]
            layouts[level]["subtitle"] = block["subtitle"]
            layouts[level]["x_axis_text"]["font_size_pt"] = config["style"]["gene_font_pt"]
            layouts[level]["color_legend_title"] = "Expression (unscaled)" if len(block["order"]) == 1 else "Average Expression"
            layouts[level]["expression_scaling_note"] = "Single group: native log1p mean expression; no across-group scaling." if len(block["order"]) == 1 else "Across-group sample z-score (ddof=1)."
            layouts[level]["expression_scaling"] = "native_log1p_mean_unscaled" if len(block["order"]) == 1 else "sample_zscore_ddof1"
            info = image_info(stage / figure_rel)
            add_check(checks, f"{level}_png_decodes", info["format"], "PNG", info["format"] == "PNG" and info["bytes"] > 0)
            add_check(
                checks,
                f"{level}_png_geometry",
                [info["width"], info["height"]],
                [expected_width, expected_height],
                info["width"] == expected_width and info["height"] == expected_height,
            )
            add_check(
                checks,
                f"{level}_png_white_background",
                info["corner_rgba"],
                [255, 255, 255, 255],
                info["corner_rgba"] == [255, 255, 255, 255],
            )
            reopen_rows.extend(
                [
                    write_artifact_reopen(stage, stage / figure_rel, "png"),
                    write_artifact_reopen(stage, stage / numeric_rel, "csv_gz"),
                    write_artifact_reopen(stage, stage / native_rel, "csv"),
                ]
            )
            figure_rows.append(
                {
                    "figure_id": level,
                    "figure_path": figure_rel.as_posix(),
                    "numeric_source_path": numeric_rel.as_posix(),
                    "native_seurat_pdata_path": native_rel.as_posix(),
                    "grouping_level": level,
                    "group_field": str(block["field"]),
                    "cells": shape[0],
                    "markers": len(final_panel),
                    "result_status": config["status"],
                }
            )
        write_json(stage / "00_contract/figure_layout_contract.json", layouts)
        pd.DataFrame(reopen_rows).to_csv(
            validation_dir / "artifact_reopen.csv", index=False
        )
        (stage / "manifest").mkdir(parents=True, exist_ok=True)
        pd.DataFrame(figure_rows).to_csv(stage / "manifest/figure_source_index.csv", index=False)
        pdf_count = len(list(stage.rglob("*.pdf")))
        add_check(checks, "output_pdf_count", pdf_count, 0, pdf_count == 0)

        guard = run_recluster_guard(config, validation_dir / "recluster_guard.log")
        add_check(
            checks,
            "repository_recluster_guard",
            guard["status"],
            "PASS or NOT_APPLICABLE",
            guard["status"] in {"PASS", "NOT_APPLICABLE", "NOT_RUN_WITH_REASON"},
            str(guard.get("reason", "")),
        )
        after = input_hash_table(protected, root)
        integrity = before.rename(
            columns={"bytes": "before_bytes", "sha256": "before_sha256"}
        ).merge(
            after.rename(columns={"bytes": "after_bytes", "sha256": "after_sha256"}),
            on=["path", "absolute_path"],
            how="outer",
            validate="one_to_one",
        )
        integrity["unchanged"] = integrity["before_bytes"].eq(
            integrity["after_bytes"]
        ) & integrity["before_sha256"].eq(integrity["after_sha256"])
        integrity.to_csv(validation_dir / "input_sha256_before_after.csv", index=False)
        add_check(checks, "all_protected_inputs_sha256_unchanged", bool(integrity["unchanged"].all()), True, bool(integrity["unchanged"].all()))

        write_json(stage / "00_contract/environment.json", environment_record(versions))
        run_contract = {
            "analysis_id": str(config["analysis_id"]),
            "renderer": RENDERER_NAME,
            "modality": MODALITY,
            "generated_at": now_local(),
            "status": config["status"],
            "canonical_status": "PASS",
            "annotation_use_case": config["project"]["annotation_use_case"],
            "resolved_annotation_track": {
                "track_id": config["project"]["_resolved_annotation_track"].get(
                    "track_id"
                ),
                "status": config["project"]["_resolved_annotation_track"].get("status"),
                "source_ids": [
                    source.get("source_id")
                    for source in config["project"]["_resolved_annotation_track"].get(
                        "sources", []
                    )
                ],
            },
            "counts_source": 'layers["counts"]',
            "cells": shape[0],
            "library_genes": shape[1],
            "join": {
                "sidecar_key": sidecar_key,
                "target_key": "obs_name",
                "sidecar_key_unique": True,
                "target_key_unique": True,
                "exact_set_equality": True,
                "row_order_join_forbidden": True,
                "cell_id_join_forbidden": MODALITY == "xenium",
                "registry_join_source_ids": config["project"][
                    "_registry_join_source_ids"
                ],
            },
            "groupings": {
                level: {
                    "field": str(config["groupings"][level]["field"]),
                    "label_map": dict(config["groupings"][level].get("label_map", {})),
                    "allow_many_to_one_label_map": bool(
                        config["groupings"][level].get(
                            "allow_many_to_one_label_map", False
                        )
                    ),
                    "order": list(groups[level]),
                    "counts": {
                        group: int(group_counts[level][index])
                        for index, group in enumerate(groups[level])
                    },
                }
                for level in groups
            },
            "marker_source_slots": len(slots),
            "provisional_unique_markers": len(provisional),
            "final_marker_count": len(final_panel),
            "final_markers": final_panel["gene"].astype(str).tolist(),
            "marker_selection_contract": dict(config["marker_selection"]),
            "minimum_max_group_pct_exp": float(config["marker_selection"]["minimum_max_group_pct_exp"]),
            "filter_groupings": list(config["marker_selection"]["filter_groupings"]),
            "dotplot_semantics": dict(config["dotplot"]),
            "figure_geometry": {
                "width_inches": effective_width_inches,
                "requested_width_inches": float(
                    width_contract["requested_width_inches"]
                ),
                "computed_minimum_width_inches": float(
                    width_contract["computed_minimum_width_inches"]
                ),
                "final_width_inches": effective_width_inches,
                "auto_expand_width": bool(width_contract["auto_expand_width"]),
                "base_width_inches": float(width_contract["base_width_inches"]),
                "marker_slot_width_inches": float(
                    width_contract["marker_slot_width_inches"]
                ),
                "minimum_marker_slot_inches": float(
                    width_contract["minimum_marker_slot_inches"]
                ),
                "longest_gene_length": int(width_contract["longest_gene_length"]),
                "height_inches": config["output"]["height_inches"],
                "height_inches_by_grouping": {k: v["height_inches"] for k, v in config["groupings"].items()},
                "style": config["style"],
                "dpi": render_dpi,
                "width_pixels": expected_width,
                "height_pixels_by_grouping": {k: int(v["height_inches"] * render_dpi) for k, v in config["groupings"].items()},
                "background": "white",
                "x_axis_text_angle_degrees": X_AXIS_TEXT_ANGLE_DEGREES,
                "x_axis_text_hjust": X_AXIS_TEXT_HJUST,
                "x_axis_text_vjust": X_AXIS_TEXT_VJUST,
                "x_axis_text_font_size_pt": config["style"]["gene_font_pt"],
            },
            "configured_seurat_native": {
                key: str(config["seurat_native"][key])
                for key in (
                    "rscript",
                    "expected_r_version",
                    "expected_seurat_version",
                    "expected_seuratobject_version",
                    "avg_exp_abs_tolerance",
                    "avg_exp_log1p_abs_tolerance",
                    "zscore_abs_tolerance",
                    "scaled_abs_tolerance",
                    "pct_exp_abs_tolerance",
                )
            },
            "seurat_native": versions,
            "native_dotplot_call": 'DotPlot(object, features=markers, assay="RNA")',
            "native_pdata_exported_before_render_override": True,
            "render_palette": ["#2166AC", "#FFFFFF", "#B2182B"],
            "render_midpoint": 0.0,
            "render_limits": [-2.5, 2.5],
            "marker_program_strip_legend": False,
            "python_seurat_parity_status": "PASS",
            "count_audit": count_audit,
            "protected_inputs": before.to_dict("records"),
            "annotation_changed": False,
            "h5ad_changed": False,
            "sidecar_changed": False,
            "registry_changed": False,
            "guide_changed": False,
            "formal_annotation_created": False,
            "recluster_performed": False,
            "figure_format": "PNG_only",
            "recluster_guard": guard,
        }
        write_json(stage / "00_contract/run_contract.json", run_contract)
        failed = [row for row in checks if row["status"] != "PASS"]
        pd.DataFrame(checks).to_csv(validation_dir / "validation_summary.csv", index=False)
        write_json(
            stage / "completion.json",
            {
                "analysis_id": str(config["analysis_id"]),
                "completed_at": now_local(),
                "status": "PASS" if not failed else "FAIL",
                "result_status": config["status"],
                "validation_checks": len(checks),
                "validation_passed": len(checks) - len(failed),
                "validation_failed": len(failed),
                "figures": [row["figure_path"] for row in figure_rows],
                "markers": len(final_panel),
                "pdf_count": pdf_count,
                "canonical_requires_real_native_seurat": True,
                "native_seurat_dotplot_rendered": True,
                "figure_geometry_pixels": {k: [expected_width, int(v["height_inches"] * render_dpi)] for k, v in config["groupings"].items()},
            },
        )
        if failed:
            raise ValueError(
                "validation failed: " + ", ".join(row["check_id"] for row in failed)
            )
        write_output_manifest(stage)
        verification = verify_output_package(stage, verify_inputs=True)
        if verification["status"] != "PASS":
            raise ValueError(f"staged package verification failed: {verification['failures']}")
        stage.rename(output)
        print(f"[{RENDERER_NAME}] PASS: {output}", flush=True)
        return output
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def run_self_test() -> None:
    obs = pd.Index(["c1", "c2", "c3"], name="obs_name")
    sidecar_key = "obs_name" if MODALITY == "xenium" else "cell_id"
    sidecar = pd.DataFrame(
        {sidecar_key: ["c3", "c1", "c2"], "group": ["B", "A", "B"]}
    )
    aligned = align_sidecar_to_obs_names(obs, sidecar, ["group"], sidecar_key)
    shuffled = align_sidecar_to_obs_names(
        obs, sidecar.sample(frac=1.0, random_state=3), ["group"], sidecar_key
    )
    if not aligned.equals(shuffled):
        raise AssertionError("obs_name join is not shuffle invariant")
    mapped = pd.DataFrame(
        {"sample_role": ["B", "P", "P2"]}, index=pd.Index(obs, name="obs_name")
    )
    _, mapped_counts, mapped_groups, _ = build_group_contract(
        mapped,
        {
            "conditions": {
                "field": "sample_role",
                "order": ["Control", "Treatment"],
                "label_map": {"B": "Control", "P": "Treatment", "P2": "Treatment"},
                "allow_many_to_one_label_map": True,
            }
        },
    )
    if mapped_groups["conditions"] != ["Control", "Treatment"] or mapped_counts["conditions"].tolist() != [1, 2]:
        raise AssertionError("explicit many-to-one grouping map changed")
    duplicate = pd.concat([sidecar, sidecar.iloc[[0]]], ignore_index=True)
    try:
        align_sidecar_to_obs_names(obs, duplicate, ["group"], sidecar_key)
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate obs_name did not fail")
    mismatched = sidecar.loc[sidecar[sidecar_key].ne("c3")].copy()
    try:
        align_sidecar_to_obs_names(obs, mismatched, ["group"], sidecar_key)
    except ValueError:
        pass
    else:
        raise AssertionError("obs_name set mismatch did not fail")
    try:
        resolve_rscript(f"/definitely/missing/{RENDERER_NAME}/Rscript")
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("missing configured Rscript did not fail closed")
    r_source = DEFAULT_ORACLE.read_text(encoding="utf-8")
    if r_source.count('plot <- DotPlot(object, features=markers, assay="RNA")') != 1:
        raise AssertionError("native DotPlot call contract changed")
    if "plot$data <-" in r_source:
        raise AssertionError("native plot$data must never be mutated")
    if "#2166AC->#FFFFFF->#B2182B" not in r_source:
        raise AssertionError("v7 palette contract changed")
    for r_contract in (
        "axis.text.x = ggplot2::element_text(",
        "angle = 60",
        "hjust = 1",
        "vjust = 1",
        "size = style$gene_font_pt",
    ):
        if r_contract not in r_source:
            raise AssertionError(f"R x-axis gene-label contract changed: {r_contract}")
    if DOT_DEFAULTS["color_mid"] != "#FFFFFF":
        raise AssertionError("v7 midpoint color changed")
    width = resolve_figure_width(20.0, ["G" for _ in range(41)], True)
    if not math.isclose(width["computed_minimum_width_inches"], 23.22):
        raise AssertionError("automatic width formula changed")
    if not math.isclose(width["final_width_inches"], 23.22):
        raise AssertionError("automatic width expansion changed")
    legacy_width = resolve_figure_width(30.0, ["G1", "LONGGENE"], True)
    if legacy_width["final_width_inches"] != 30.0:
        raise AssertionError("legacy numeric width compatibility changed")
    long_gene_width = resolve_figure_width(1.0, ["G" * 20], True)
    if not math.isclose(long_gene_width["marker_slot_width_inches"], 0.68):
        raise AssertionError("long-gene marker slot formula changed")
    try:
        resolve_figure_width(20.0, ["G" for _ in range(41)], False)
    except ValueError:
        pass
    else:
        raise AssertionError("undersized fixed width did not fail closed")
    if (
        X_AXIS_TEXT_ANGLE_DEGREES,
        X_AXIS_TEXT_HJUST,
        X_AXIS_TEXT_VJUST,
        X_AXIS_TEXT_FONT_SIZE_PT,
    ) != (60.0, 1.0, 1.0, 9.0):
        raise AssertionError("x-axis gene-label contract changed")
    with tempfile.TemporaryDirectory(prefix=f"{RENDERER_NAME}.markers.") as temporary:
        marker_path = Path(temporary) / "markers.yaml"
        marker_path.write_text(
            yaml.safe_dump(
                {
                    "programs": [
                        {
                            "name": "P1",
                            "short_label": "One",
                            "genes": [
                                {"gene": "G1", "role": "core"},
                                {"gene": "G2", "role": "support"},
                            ],
                        },
                        {
                            "name": "P2",
                            "short_label": "Two",
                            "color": "#E45756",
                            "genes": ["G2", "G3"],
                        },
                    ]
                },
                sort_keys=False,
            ),
            encoding="utf-8",
        )
        slots = expand_marker_slots(marker_path)
        panel, audit = select_and_deduplicate_markers(
            slots,
            {
                "include_categories": [],
                "exclude_categories": [],
                "include_programs": [],
                "exclude_programs": [],
                "deduplication": {"strategy": "first_occurrence"},
            },
        )
        if panel["gene"].tolist() != ["G1", "G2", "G3"]:
            raise AssertionError("deterministic marker deduplication changed")
        colors = panel.drop_duplicates("program_order")["program_color"].tolist()
        if colors != [PROGRAM_PALETTE[0], "#E45756"]:
            raise AssertionError("marker-program color fallback changed")
        repeated = audit.loc[
            audit["gene"].eq("G2") & audit["source_program"].eq("P2"),
            "slot_decision",
        ].tolist()
        if repeated != ["deduplicated_to_shared_or_first_program"]:
            raise AssertionError("duplicate marker audit changed")
    print(f"PASS {RENDERER_NAME} self-test")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    render_parser = subparsers.add_parser("render", help="render a new canonical package")
    render_parser.add_argument("--config", type=Path, required=True)
    render_parser.add_argument("--output", type=Path, default=None)
    render_parser.add_argument("--chunk-rows", type=int, default=20_000)
    render_parser.add_argument("--dpi", type=int, default=None)
    render_parser.add_argument("--oracle-temp-dir", type=Path, default=None)
    validate_parser = subparsers.add_parser("validate", help="verify an existing package")
    validate_parser.add_argument("--output", type=Path, required=True)
    validate_parser.add_argument("--verify-inputs", action="store_true")
    subparsers.add_parser("self-test", help="run lightweight deterministic checks")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "self-test":
        run_self_test()
        return 0
    if args.command == "validate":
        result = verify_output_package(args.output, args.verify_inputs)
        print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
        return 0 if result["status"] == "PASS" else 1
    if args.chunk_rows <= 0:
        raise ValueError("--chunk-rows must be positive")
    if args.dpi is not None and args.dpi <= 0:
        raise ValueError("--dpi must be positive")
    if args.oracle_temp_dir is not None:
        args.oracle_temp_dir = args.oracle_temp_dir.resolve()
        if not args.oracle_temp_dir.is_dir():
            raise FileNotFoundError(args.oracle_temp_dir)
    output = run_build(args)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
