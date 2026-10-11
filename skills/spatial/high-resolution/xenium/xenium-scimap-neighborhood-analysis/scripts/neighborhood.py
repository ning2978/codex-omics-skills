"""Portable phenotype composition RCN/KMeans and target NMF; no expression clustering."""
from pathlib import Path
import argparse
import itertools
import json
import sys

shared = next(p/'_shared/spatial' for p in Path(__file__).resolve().parents if (p/'_shared/spatial/run_contract.py').is_file())
sys.path.insert(0, str(shared))
from run_contract import Run, cell_table, mask, neighborhoods, write_json

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF
from sklearn.metrics import adjusted_rand_score


def count_features(table, cfg, centers, pool, radius=None, k=None):
    labels = table[cfg['phenotype_key']]
    features = sorted(labels.iloc[pool].dropna().astype(str).unique())
    if not features:
        raise ValueError('No measured/trusted phenotype features in contributor pool')
    lookup = {v:i for i,v in enumerate(features)}
    graph = neighborhoods(table, cfg, centers, pool, radius=radius, k=k)
    counts = np.zeros((len(centers), len(features)), dtype=np.int64)
    for row, center in enumerate(centers):
        for label in labels.iloc[graph[int(center)]].dropna().astype(str):
            counts[row, lookup[label]] += 1
    return counts, features, graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if cfg['mode'] not in {'rcn', 'nmf'} or not cfg.get('namespace') or not cfg.get('assignment_id'):
        raise ValueError('Explicit mode, assignment_id and namespace required')
    if cfg.get('cluster_method') != 'kmeans':
        raise ValueError('This portable baseline implements explicit KMeans only; scimap/Leiden require their own runner')
    if ('radii' in cfg) == ('knn' in cfg):
        raise ValueError('Choose radii OR knn grid')
    if cfg.get('representation') not in {'counts', 'proportions'}:
        raise ValueError('Choose counts or proportions representation explicitly')
    if not cfg.get('seeds') or len(set(cfg['seeds'])) != len(cfg['seeds']):
        raise ValueError('Distinct explicit seeds required')
    run = Run(args.out_dir, cfg, [args.config, cfg['cells_csv']], __file__)
    if run.resumed:
        print('PASS complete-artifact resume')
        return
    table = cell_table(cfg['cells_csv'], cfg, labels=True)
    centers = np.flatnonzero(mask(table, cfg.get('center_mask_key')))
    pool = np.flatnonzero(mask(table, cfg.get('neighbor_mask_key')))
    if cfg.get('center_phenotypes') is not None:
        centers = centers[table.iloc[centers][cfg['phenotype_key']].isin(cfg['center_phenotypes'])]
    if not len(centers):
        raise ValueError('No eligible centers')
    scales = [('radius', r) for r in cfg.get('radii', [])] + [('knn', k) for k in cfg.get('knn', [])]
    if not scales or len(set(scales)) != len(scales):
        raise ValueError('Distinct nonempty spatial-scale grid required')
    assignments, summaries, compositions, loadings, weights, stability = [], [], [], [], [], []
    results = []
    for method, value in scales:
        counts, features, graph = count_features(table, cfg, centers, pool,
            radius=value if method == 'radius' else None, k=value if method == 'knn' else None)
        scale = f'{method}{value}'
        sums = counts.sum(axis=1)
        valid = sums > 0
        proportions = np.divide(counts, sums[:, None], out=np.full(counts.shape, np.nan), where=sums[:, None] != 0)
        raw = table.iloc[centers][[cfg['cell_key'], cfg['roi_key'], cfg['biological_unit_key']]].reset_index(drop=True)
        raw['n_neighbors'] = [len(graph[int(c)]) for c in centers]
        raw['n_labeled_neighbors'] = sums
        raw['status'] = np.where(valid, 'eligible', 'zero_labeled_neighbors')
        raw.join(pd.DataFrame(counts, columns=features).add_prefix('count::')).join(
            pd.DataFrame(proportions, columns=features).add_prefix('proportion::')).to_csv(args.out_dir/f'composition_{scale}.csv', index=False)
        ranks = cfg['ranks'] if cfg['mode'] == 'nmf' else [None]
        for rank, seed in itertools.product(ranks, cfg['seeds']):
            namespace = f"{cfg['namespace']}/{scale}" + (f'/rank{rank}' if rank else '') + f'/seed{seed}'
            data = counts[valid] if cfg['representation'] == 'counts' else proportions[valid]
            n_states = cfg['n_states']
            support = len(data)
            fit_ok = support >= max(cfg.get('min_centers', 2), n_states) and (rank is None or 1 <= rank <= min(data.shape))
            frame = raw.copy()
            frame['namespace'] = namespace
            frame['assignment_id'] = cfg['assignment_id']
            frame['assignment'] = pd.Series([None]*len(frame), dtype=object)
            status = 'PASS' if fit_ok else 'UNDER_SUPPORTED'
            error, iterations, components = None, None, None
            n_states_found, smallest_state = None, None
            if fit_ok:
                if rank is not None:
                    model = NMF(n_components=rank, init='random', random_state=seed, max_iter=cfg.get('max_iter', 1000))
                    representation = model.fit_transform(data)
                    error, iterations = float(model.reconstruction_err_), int(model.n_iter_)
                    components = model.components_
                    if iterations >= model.max_iter:
                        status = 'NONCONVERGED'
                    for component, values in enumerate(components):
                        for phenotype, loading in zip(features, values):
                            loadings.append({'namespace':namespace,'component':component+1,'neighbor_phenotype':phenotype,'loading':float(loading)})
                    for key, row in zip(raw.loc[valid, cfg['cell_key']], representation):
                        for component, weight in enumerate(row):
                            weights.append({cfg['cell_key']:key,'namespace':namespace,'component':component+1,'weight':float(weight)})
                else:
                    representation = data
                model_cluster = KMeans(n_clusters=n_states, random_state=seed, n_init=10)
                labels = model_cluster.fit_predict(representation)
                state_sizes = np.unique(labels, return_counts=True)[1]
                n_states_found, smallest_state = len(state_sizes), int(state_sizes.min())
                if n_states_found < n_states or smallest_state < cfg.get('min_state_centers', 1):
                    status = 'UNDER_SUPPORTED_STATES'
                suffix = 'RCN' if cfg['mode'] == 'rcn' else 'NMFstate'
                frame.loc[valid, 'assignment'] = [f'{namespace}::{suffix}{label+1}' for label in labels]
                frame.loc[valid, 'status'] = 'exploratory' if status == 'PASS' else status
                results.append((namespace, scale, rank, frame['assignment'].to_numpy(), components))
                for label in np.unique(labels):
                    for phenotype, fraction in zip(features, proportions[valid][labels == label].mean(axis=0)):
                        compositions.append({'namespace':namespace,'assignment':f'{namespace}::{suffix}{label+1}',
                            'neighbor_phenotype':phenotype,'mean_proportion':float(fraction),'n_centers':int((labels == label).sum())})
            else:
                frame.loc[valid, 'status'] = status
            assignments.append(frame)
            summaries.append({'namespace':namespace,'scale':scale,'rank':rank,'seed':seed,'status':status,
                             'n_centers':len(centers),'n_contributors':len(pool),'n_supported_centers':support,
                             'n_zero_labeled_neighbors':int((~valid).sum()),'reconstruction_error':error,'iterations':iterations})
            summaries[-1].update(n_states_found=n_states_found,smallest_state_centers=smallest_state)
    for left, right in itertools.combinations(results, 2):
        ns1, scale1, rank1, labels1, h1 = left
        ns2, scale2, rank2, labels2, h2 = right
        valid = pd.notna(labels1) & pd.notna(labels2)
        row = {'left_namespace':ns1,'right_namespace':ns2,'n_common_supported_centers':int(valid.sum()),
               'assignment_ari':float(adjusted_rand_score(labels1[valid], labels2[valid])) if valid.any() else None,
               'matched_component_cosine':None}
        if h1 is not None and h2 is not None and rank1 == rank2:
            norm1 = np.linalg.norm(h1, axis=1); norm2 = np.linalg.norm(h2, axis=1)
            cos = h1 @ h2.T / np.maximum(norm1[:,None]*norm2[None,:], 1e-15)
            i,j = linear_sum_assignment(-cos)
            row['matched_component_cosine'] = float(cos[i,j].mean())
        stability.append(row)
    combined = pd.concat(assignments, ignore_index=True)
    combined.to_csv(args.out_dir/'assignments.csv', index=False)
    pd.DataFrame(summaries).to_csv(args.out_dir/'fit_summary.csv', index=False)
    pd.DataFrame(compositions, columns=['namespace','assignment','neighbor_phenotype','mean_proportion','n_centers']).to_csv(args.out_dir/'state_composition.csv', index=False)
    pd.DataFrame(stability, columns=['left_namespace','right_namespace','n_common_supported_centers','assignment_ari','matched_component_cosine']).to_csv(args.out_dir/'stability.csv', index=False)
    if cfg['mode'] == 'nmf':
        pd.DataFrame(loadings, columns=['namespace','component','neighbor_phenotype','loading']).to_csv(args.out_dir/'nmf_loadings.csv', index=False)
        pd.DataFrame(weights, columns=[cfg['cell_key'],'namespace','component','weight']).to_csv(args.out_dir/'nmf_weights.csv', index=False)
    support = combined.groupby(['namespace',cfg['biological_unit_key'],cfg['roi_key'],'status'], dropna=False).size().rename('n_centers').reset_index()
    support.to_csv(args.out_dir/'biological_unit_support.csv', index=False)
    grade = 'PASS' if all(s['status'] == 'PASS' for s in summaries) else 'WARN'
    run.finish({'status':grade,'scope':'spatial_phenotype_composition','self_excluded':True,
                'cross_roi_edges':0,'feature_kind':'neighbor_phenotype'}, evidence_grade='exploratory',
                authority=cfg.get('authority','exploratory'), center_pool_contract={
                    'n_centers':len(centers),'n_contributors':len(pool),'missing_contributor_labels':int(table.iloc[pool][cfg['phenotype_key']].isna().sum())})
    print(f'{grade} phenotype-neighborhood {cfg["mode"]}')


if __name__ == '__main__':
    main()
