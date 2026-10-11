"""Reusable numeric contracts; no biological labels or automatic data exclusions."""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd
from scipy import sparse


def cell_keys(libraries, barcodes) -> pd.Index:
    """Reversible tuple encoding avoids both cross-capture and delimiter collisions."""
    libraries, barcodes = list(libraries), list(barcodes)
    if len(libraries) != len(barcodes):
        raise ValueError("Library and barcode lengths differ")
    if any(pd.isna(v) or str(v) == "" for v in libraries + barcodes):
        raise ValueError("Missing library/barcode")
    keys = pd.Index([json.dumps([str(a), str(b)], ensure_ascii=False, separators=(",", ":"))
                     for a, b in zip(libraries, barcodes, strict=True)])
    if not keys.is_unique:
        raise ValueError("Duplicate library-barcode pair")
    return keys


def align_frame(frame: pd.DataFrame, keys, *, exact=True) -> pd.DataFrame:
    keys = pd.Index(keys)
    if not keys.is_unique or not frame.index.is_unique:
        raise ValueError("Duplicate alignment keys")
    missing, extra = keys.difference(frame.index), frame.index.difference(keys)
    if len(missing) or (exact and len(extra)):
        raise ValueError(f"Key mismatch: missing={len(missing)}, extra={len(extra)}")
    return frame.loc[keys].copy()


def validate_matrix(matrix, *, integer=False):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Expression must be finite and nonnegative, never centered/scaled")
    if integer and not np.allclose(values, np.rint(values), atol=1e-6, rtol=0):
        raise ValueError("Raw count model requires integer-like counts; do not round corrections")


def matrix_fingerprint(matrix):
    """Canonical CSR fingerprint; use artifact hashes for disk-scale runners."""
    mat = sparse.csr_matrix(matrix, copy=True)
    mat.sum_duplicates()
    mat.sort_indices()
    mat.eliminate_zeros()
    h = hashlib.sha256()
    h.update(np.asarray(mat.shape, dtype="<i8").tobytes())
    for values, dtype in ((mat.indptr, "<i8"), (mat.indices, "<i8"), (mat.data, "<f8")):
        h.update(np.asarray(values, dtype=dtype).tobytes())
    return h.hexdigest()


def key_fingerprint(keys):
    return hashlib.sha256(json.dumps(list(map(str, keys)), ensure_ascii=False,
                                     separators=(",", ":")).encode()).hexdigest()


def update_umap_layout(adata, *, neighbors_key, expected_cells, expected_graph,
                       random_state=0, min_dist=.5, spread=1.):
    """Only re-layout a validated unchanged graph; never run clustering/markers."""
    import scanpy as sc
    if key_fingerprint(adata.obs_names) != expected_cells:
        raise ValueError("Layout-only cell set/order changed; retained-count rerun required")
    if neighbors_key not in adata.uns:
        raise ValueError("Explicit existing neighbor graph required")
    graph_key = adata.uns[neighbors_key]["connectivities_key"]
    if matrix_fingerprint(adata.obsp[graph_key]) != expected_graph:
        raise ValueError("Layout-only graph changed; reclassify requested analysis")
    sc.tl.umap(adata, neighbors_key=neighbors_key, random_state=random_state,
               min_dist=min_dist, spread=spread)
    adata.uns["umap_layout_provenance"] = {
        "mode": "layout_only", "cells_sha256": expected_cells,
        "graph_sha256": expected_graph, "neighbors_key": neighbors_key,
        "random_state": random_state, "min_dist": min_dist, "spread": spread,
    }


def prepare_log(adata, *, source="counts", target_sum=1e4):
    import scanpy as sc
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Unique stable cell/feature IDs required")
    if source not in adata.layers:
        raise KeyError(source)
    if not adata.uns.get("matrix_sources", {}).get(source):
        raise ValueError("Declare count-source provenance in uns['matrix_sources'] first")
    validate_matrix(adata.layers[source], integer=(source == "counts"))
    if not np.isfinite(target_sum) or target_sum <= 0:
        raise ValueError("Positive finite normalization target required")
    # Disposable object prevents mutation of X, raw, counts or old normalizations.
    import anndata as ad
    work = ad.AnnData(X=adata.layers[source].astype(np.float64).copy())
    sc.pp.normalize_total(work, target_sum=target_sum)
    sc.pp.log1p(work)
    adata.layers["log1p"] = work.X.copy()
    adata.uns.setdefault("matrix_sources", {})["log1p"] = {
        "kind": "unintegrated_log1p", "source_layer": source,
        "target_sum": float(target_sum), "log_base": "natural",
        "source_fingerprint": matrix_fingerprint(adata.layers[source]),
    }


def validate_log(adata, layer="log1p"):
    meta = adata.uns.get("matrix_sources", {}).get(layer, {})
    if meta.get("kind") != "unintegrated_log1p" or meta.get("log_base") != "natural":
        raise ValueError("Explicit unintegrated natural-log1p matrix provenance required")
    validate_matrix(adata.layers[layer])


def select_hvg(adata, *, n_top_genes=2000, batch_key=None,
               allow_unbatched=False, downgrade_reason=None):
    import scanpy as sc
    if batch_key and (batch_key not in adata.obs or adata.obs[batch_key].isna().any()):
        raise ValueError("Missing/incomplete batch key")
    validate_matrix(adata.layers["counts"], integer=True)
    validate_log(adata)
    attempts = []
    choices = [("seurat_v3", "counts", batch_key), ("cell_ranger", "log1p", batch_key)]
    if allow_unbatched:
        if not batch_key or not downgrade_reason:
            raise ValueError("Unbatched downgrade requires batch_key and explicit rationale")
        choices.append(("cell_ranger", "log1p", None))
    for flavor, layer, batch in choices:
        record = {"flavor": flavor, "layer": layer, "batch_key": batch,
                  "n_top_genes": min(n_top_genes, adata.n_vars)}
        try:
            result = sc.pp.highly_variable_genes(
                adata, layer=layer, flavor=flavor, batch_key=batch,
                n_top_genes=record["n_top_genes"], inplace=False)
            if not result.highly_variable.any():
                raise ValueError("No usable HVGs")
        except Exception as exc:
            record.update(status="FAIL", error=f"{type(exc).__name__}: {exc}")
            attempts.append(record)
            continue
        record.update(status="PASS")
        if batch_key and batch is None:
            record["downgrade_reason"] = downgrade_reason
        attempts.append(record)
        for col in result:
            adata.var[col] = result[col]
        return attempts
    raise RuntimeError("All HVG attempts failed: " + json.dumps(attempts))


def rank_markers(adata, groupby, *, layer="log1p", min_cells=3, method="t-test"):
    """Full available gene scope, regardless of raw/HVG/scaled work matrices."""
    import anndata as ad
    import scanpy as sc
    validate_log(adata, layer)
    if not adata.var_names.is_unique or not adata.obs_names.is_unique:
        raise ValueError("Unique feature/cell keys required")
    if groupby not in adata.obs or adata.obs[groupby].isna().any():
        raise ValueError("Missing cluster assignments")
    mat = adata.layers[layer]
    detected = np.asarray((mat > 0).sum(axis=0)).ravel()
    keep = detected >= min_cells
    if min_cells < 1 or not keep.any():
        raise ValueError("No eligible marker genes or invalid detection filter")
    work = ad.AnnData(X=mat[:, keep].copy(), obs=adata.obs[[groupby]].copy(),
                     var=adata.var.loc[keep].copy())
    work.obs[groupby] = work.obs[groupby].astype("category").cat.remove_unused_categories()
    sizes = work.obs[groupby].value_counts()
    if len(sizes) < 2 or (sizes < 2).any():
        raise ValueError("Exploratory marker test needs >=2 clusters and >=2 cells/cluster")
    work.layers[layer] = work.X.copy()
    sc.tl.rank_genes_groups(work, groupby, layer=layer, use_raw=False, method=method,
                           pts=True, n_genes=work.n_vars)
    result = sc.get.rank_genes_groups_df(work, group=None)
    result = result.rename(columns={"names": "gene", "scores": "score",
                                     "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"})
    result["gene_symbol"] = result.gene.map(work.var.get("gene_symbol", pd.Series(work.var_names, index=work.var_names)))
    result["n_in"] = result.group.map(sizes).astype(int)
    result["n_out"] = work.n_obs - result.n_in
    result["matrix_source"] = f"layer:{layer}"
    result["marker_scope"] = "all_available_expression_filtered_genes"
    result["method"] = method
    result["min_detected_cells"] = min_cells
    return result


def assess_design(samples, *, unit="biological_unit", condition="condition", donor="donor",
                  batch="technical_batch", paired=False):
    """Conservative gate for two-condition independent/paired designs, not power proof."""
    reasons = []
    required = [unit, condition] + ([donor] if paired else [])
    if any(k not in samples or samples[k].isna().any() or (samples[k].astype(str) == "").any()
           for k in required):
        return {"inference_status": "descriptive_only", "reasons": ["Missing experimental-unit/condition/pair metadata"]}
    # A sample/unit may have technical rows, but covariates cannot conflict within it.
    for col in [donor, batch]:
        if col in samples and samples.groupby([unit, condition], observed=True)[col].nunique(dropna=False).gt(1).any():
            reasons.append(f"Conflicting {col} within unit-condition; specify a supported design")
    data = samples.drop_duplicates([unit, condition]).copy()
    groups = data[condition].unique()
    if len(groups) != 2:
        reasons.append("Helper requires exactly two contrast conditions")
    if data.groupby(condition, observed=True)[unit].nunique().min() < 2:
        reasons.append("Fewer than two independent units per condition")
    if not paired and data.groupby(unit, observed=True)[condition].nunique().gt(1).any():
        reasons.append("Repeated units need a paired/repeated-measures model")
    if not paired and donor in data and data[donor].duplicated().any():
        reasons.append("Repeated donor samples are not independent units")
    if paired:
        pairs = pd.crosstab(data[donor], data[condition])
        if len(groups) != 2 or len(pairs) < 2 or not (pairs == 1).all().all():
            reasons.append("Need at least two complete one-sample-per-condition donor pairs")
    columns = [np.ones((len(data), 1)), pd.get_dummies(data[condition], drop_first=True, dtype=float).to_numpy()]
    if paired:
        columns.append(pd.get_dummies(data[donor], drop_first=True, dtype=float).to_numpy())
    if batch and batch in data:
        if data[batch].isna().any():
            reasons.append("Missing technical batch metadata")
        columns.append(pd.get_dummies(data[batch], drop_first=True, dtype=float).to_numpy())
    design = np.column_stack(columns)
    rank = int(np.linalg.matrix_rank(design)) if len(data) else 0
    if rank < design.shape[1]:
        reasons.append("Nonestimable/confounded design; do not silently drop terms")
    if len(data) <= rank:
        reasons.append("No residual degrees of freedom")
    return {"inference_status": "descriptive_only" if reasons else "eligible_for_model",
            "reasons": reasons, "design_rank": rank, "residual_df": len(data)-rank,
            "units_per_condition": data.groupby(condition, observed=True)[unit].nunique().to_dict()}


def pseudobulk(adata, *, keys=("cell_identity", "biological_unit", "condition"), min_cells=10):
    """Sum technical captures in declared units; no zero columns for unobserved types."""
    validate_matrix(adata.layers["counts"], integer=True)
    if min_cells < 1 or any(k not in adata.obs or adata.obs[k].isna().any() for k in keys):
        raise ValueError("Complete aggregation keys and positive min_cells required")
    tuples = pd.MultiIndex.from_frame(adata.obs[list(keys)])
    codes, labels = pd.factorize(tuples, sort=True)
    labels = labels.set_names(list(keys))
    membership = sparse.csr_matrix((np.ones(adata.n_obs), (codes, np.arange(adata.n_obs))),
                                   shape=(len(labels), adata.n_obs))
    sums = (membership @ sparse.csr_matrix(adata.layers["counts"])).tocsr()
    meta = labels.to_frame(index=False)
    meta["n_cells"] = np.bincount(codes, minlength=len(labels))
    meta["eligible"] = meta.n_cells >= min_cells
    # Preserve ineligible observed aggregates in the ledger, not in model counts.
    meta.index = pd.Index([json.dumps(list(v), ensure_ascii=False) for v in labels], name="aggregate_id")
    return sums[meta.eligible.to_numpy()], meta.loc[meta.eligible].copy(), meta


def merge_velocity_layers(adata, velocity, *, allow_gene_subset=False):
    """Exact analyzed-cell conservation; optional gene restriction has an attrition report."""
    if any(not idx.is_unique for idx in (adata.obs_names, adata.var_names,
                                         velocity.obs_names, velocity.var_names)):
        raise ValueError("Unique stable cell/gene keys required, not duplicate symbols")
    missing = adata.obs_names.difference(velocity.obs_names)
    if len(missing):
        raise ValueError(f"Velocity missing {len(missing)} analyzed cells; choose an explicit subset first")
    genes = adata.var_names[adata.var_names.isin(velocity.var_names)]
    if not len(genes) or (len(genes) != adata.n_vars and not allow_gene_subset):
        raise ValueError("Velocity gene mismatch; explicit restricted derivative required")
    out = adata[:, genes].copy()
    aligned = velocity[out.obs_names, out.var_names]
    for layer in ("spliced", "unspliced", "ambiguous"):
        if layer not in velocity.layers:
            if layer == "ambiguous":
                continue
            raise ValueError(f"Missing velocity layer {layer}")
        validate_matrix(aligned.layers[layer])
        out.layers[layer] = aligned.layers[layer].copy()
    out.uns["velocity_alignment"] = {
        "analyzed_cells": out.n_obs, "extra_velocity_cells": len(velocity.obs_names.difference(out.obs_names)),
        "genes_before": adata.n_vars, "genes_after": out.n_vars,
        "omitted_genes": adata.var_names.difference(genes).tolist(),
    }
    return out


def validate_velocity(adata, mode):
    if mode not in {"stochastic", "dynamical"}:
        raise ValueError("Choose stochastic or dynamical")
    for layer in ("spliced", "unspliced"):
        if layer not in adata.layers:
            raise ValueError(f"Missing {layer}")
        validate_matrix(adata.layers[layer])
        mat = adata.layers[layer]
        positive = (mat > 0).nnz if sparse.issparse(mat) else np.count_nonzero(mat > 0)
        if mat.shape != adata.shape or not positive:
            raise ValueError(f"Empty/misaligned {layer}")
    if mode == "dynamical":
        alpha = adata.var.get("fit_alpha")
        if "recover_dynamics" not in adata.uns or alpha is None or not np.any(np.isfinite(alpha) & (alpha > 0)):
            raise ValueError("Dynamical mode/latent_time require successful recover_dynamics first")
