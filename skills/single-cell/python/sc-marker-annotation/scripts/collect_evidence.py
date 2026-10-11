"""Portable measured-expression review; emits evidence, never labels or exclusions."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

HELPER = Path(__file__).resolve().parents[3] / 'workflow/workflow-sc-analysis/scripts/sc_contracts.py'
spec = importlib.util.spec_from_file_location('sc_contracts_evidence', HELPER)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def collect(source, config, out):
    source, out = Path(source).resolve(), Path(out).resolve()
    if out.exists():
        raise FileExistsError(out)
    a = ad.read_h5ad(source)
    cluster, sample, unit = (config[k] for k in ('cluster_key', 'sample_key', 'unit_key'))
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError('Unique cell and feature keys required; map symbols explicitly')
    for key in [cluster, sample] + ([unit] if unit else []):
        if key not in a.obs or a.obs[key].isna().any() or a.obs[key].astype(str).eq('').any():
            raise ValueError(f'Missing/incomplete {key}')
    if unit and a.obs.groupby(sample, observed=True)[unit].nunique().gt(1).any():
        raise ValueError('A sample maps to multiple biological units')
    for key in ('species', 'tissue', 'platform', 'gene_scope', 'counts_origin'):
        if not config.get(key):
            raise ValueError(f'Declare {key}')
    c.validate_matrix(a.layers['counts'], integer=True)
    if 'log1p' not in a.layers:
        a.uns.setdefault('matrix_sources', {})['counts'] = config['counts_origin']
        c.prepare_log(a, target_sum=config.get('target_sum', 1e4))
    c.validate_log(a)
    programs = config['programs']
    if not programs or len({p['id'] for p in programs}) != len(programs):
        raise ValueError('Nonempty unique program IDs required')
    for p in programs:
        if p.get('role') not in {'target', 'competing', 'state'} or not all(p.get(k) for k in ('id', 'source', 'lineage', 'genes')):
            raise ValueError('Programs require id, role, lineage, source and explicit feature IDs')
        if len(set(p['genes'])) != len(p['genes']):
            raise ValueError('Duplicate feature in program')
    minimum = int(config.get('min_detected_cells', 3))
    if minimum < 1:
        raise ValueError('Positive detection filter required')
    counts = sparse.csr_matrix(a.layers['counts'])
    expr = sparse.csr_matrix(a.layers['log1p'])
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    mask = detected >= minimum
    availability = []
    program_hits = {}
    for p in programs:
        found = []
        for gene in p['genes']:
            measured = gene in a.var_names
            idx = a.var_names.get_loc(gene) if measured else None
            availability.append(dict(program=p['id'], role=p['role'], lineage=p['lineage'], source=p['source'],
                gene=gene, availability='measured' if measured else 'unavailable_in_input',
                detected_cells=int(detected[idx]) if measured else None,
                marker_test_eligible=bool(mask[idx]) if measured else None))
            if measured:
                found.append(idx)
        # A missing panel never becomes a zero-expression program.
        n = np.asarray((counts[:, found] > 0).sum(axis=1)).ravel() if found else None
        program_hits[p['id']] = (found, n)
    labels = a.obs[cluster].astype(str).to_numpy()
    units = a.obs[unit].astype(str).to_numpy() if unit else None
    samples = a.obs[sample].astype(str).to_numpy()
    genes = list(dict.fromkeys(g for p in programs for g in p['genes']))
    gene_rows, program_rows, pair_rows, qc_rows, composition = [], [], [], [], []
    qc_fields = list(dict.fromkeys(['total_counts', 'n_genes_by_counts'] + config.get('qc_fields', [])))
    a.obs['total_counts'] = np.asarray(counts.sum(axis=1)).ravel()
    a.obs['n_genes_by_counts'] = np.asarray((counts > 0).sum(axis=1)).ravel()
    for label in sorted(set(labels)):
        inside = labels == label
        scopes = [('cluster', 'all', inside)]
        scopes += [('sample', v, inside & (samples == v)) for v in sorted(set(samples[inside]))]
        if unit:
            scopes += [('biological_unit', v, inside & (units == v)) for v in sorted(set(units[inside]))]
        for level, name, rows in scopes:
            base = dict(cluster=label, scope=level, unit=name, n_cells=int(rows.sum()))
            if level != 'cluster':
                denominator = int((samples == name).sum()) if level == 'sample' else int((units == name).sum())
                composition.append(dict(**base, fraction_of_cluster=float(rows.sum()/inside.sum()), fraction_of_unit=float(rows.sum()/denominator)))
            for gene in genes:
                if gene not in a.var_names:
                    gene_rows.append(dict(**base, gene=gene, measured=False, pct_detected=None, mean_log1p=None, pct_out=None))
                    continue
                idx = a.var_names.get_loc(gene)
                gene_rows.append(dict(**base, gene=gene, measured=True,
                    pct_detected=float((counts[rows, idx] > 0).mean()), mean_log1p=float(expr[rows, idx].mean()),
                    pct_out=float((counts[~inside, idx] > 0).mean()) if (~inside).any() else None))
            for p in programs:
                found, n = program_hits[p['id']]
                program_rows.append(dict(**base, program=p['id'], role=p['role'], lineage=p['lineage'],
                    measured_genes=len(found), declared_genes=len(p['genes']),
                    pct_any=float((n[rows] >= 1).mean()) if found else None,
                    pct_two_or_more=float((n[rows] >= 2).mean()) if len(found) >= 2 else None))
            for target in [p for p in programs if p['role'] == 'target']:
                for foreign in [p for p in programs if p['role'] == 'competing']:
                    tf, tn = program_hits[target['id']]; ff, fn = program_hits[foreign['id']]
                    evaluable = bool(tf) and len(ff) >= 2
                    pair_rows.append(dict(**base, target=target['id'], competing=foreign['id'],
                        threshold='target>=1 AND competing>=2 distinct measured features', evaluable=evaluable,
                        pct_codetected=float(((tn[rows] >= 1) & (fn[rows] >= 2)).mean()) if evaluable else None))
            for field in qc_fields:
                if field not in a.obs:
                    qc_rows.append(dict(**base, field=field, status='unavailable', value=None, statistic='unavailable'))
                elif pd.api.types.is_numeric_dtype(a.obs[field]):
                    v = a.obs.loc[rows, field]
                    for stat, val in [('median', v.median()), ('q10', v.quantile(.1)), ('q90', v.quantile(.9)), ('missing', int(v.isna().sum()))]:
                        qc_rows.append(dict(**base, field=field, status='observed', value=val, statistic=stat))
                else:
                    for val, num in a.obs.loc[rows, field].astype('string').fillna('unavailable').value_counts().items():
                        qc_rows.append(dict(**base, field=field, status='observed', value=int(num), statistic='count:'+val))
    # Reuse the shared full-gene test, explicitly recording insufficient groups.
    sizes = a.obs[cluster].value_counts()
    skipped = []
    if len(sizes) < 2 or (sizes < 2).any() or not mask.any():
        skipped.append({'status':'SKIP', 'reason':'Insufficient groups/cells or no genes pass detection filter', 'sizes':sizes.to_dict()})
        markers = pd.DataFrame(columns=['group','gene','score','logfoldchanges','pvals','pvals_adj','pct_in','pct_out'])
    else:
        markers = c.rank_markers(a, cluster, min_cells=minimum, method=config.get('method', 't-test'))
    out.mkdir(parents=True)
    tables = {'gene_availability':availability, 'gene_evidence':gene_rows, 'program_support':program_rows,
              'codetection':pair_rows, 'unit_composition':composition, 'qc_evidence':qc_rows}
    for name, rows in tables.items():
        columns = {
            'codetection':['cluster','scope','unit','n_cells','target','competing','threshold','evaluable','pct_codetected'],
            'unit_composition':['cluster','scope','unit','n_cells','fraction_of_cluster','fraction_of_unit']
        }.get(name)
        pd.DataFrame(rows, columns=columns).to_csv(out/(name+'.csv'), index=False)
    markers.to_csv(out/'all_markers.csv.gz', index=False)
    top = markers[(markers.pvals_adj < .05) & (markers.pct_in > .3)].sort_values('logfoldchanges', ascending=False).groupby('group', observed=True).head(20)
    top.to_csv(out/'top20.csv', index=False)
    manifest = {'schema':'annotation_evidence_v1', 'technical_status':'PASS', 'biological_decision':'not_automated',
        'source':str(source), 'source_sha256':sha(source), 'config':config, 'actual_command':sys.argv,
        'runtime':{'executable':sys.executable, 'versions':{pkg:importlib.metadata.version(pkg) for pkg in ('anndata','scanpy','numpy','scipy','pandas')}},
        'implementation':{str(p):sha(p) for p in (Path(__file__).resolve(), HELPER)},
        'counts_sha256':c.matrix_fingerprint(counts), 'cells_sha256':c.key_fingerprint(a.obs_names),
        'genes_sha256':c.key_fingerprint(a.var_names), 'detection_matrix':'layers/counts > 0',
        'expression_matrix':'layers/log1p', 'matrix_sources':a.uns.get('matrix_sources', {}),
        'marker_scope':'all_available_expression_filtered_genes', 'marker_skipped':skipped,
        'unit_support_status':'available' if unit else 'unavailable_descriptive_only',
        'outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--config', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    args = p.parse_args()
    result = collect(args.input, json.loads(args.config.read_text()), args.out)
    print(json.dumps({'technical_status':result['technical_status'], 'out':str(args.out)}))

if __name__ == '__main__':
    main()
