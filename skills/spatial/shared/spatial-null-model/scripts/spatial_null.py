"""Fixed-anchor within-ROI/block/stratum target-label nulls; no group permutation."""
from pathlib import Path
import argparse
import json
import sys

shared = next(p/'_shared/spatial' for p in Path(__file__).resolve().parents if (p/'_shared/spatial/run_contract.py').is_file())
sys.path.insert(0, str(shared))
from run_contract import Run, cell_table, mask, neighborhoods

import numpy as np
import pandas as pd


def bh(values, family_size=None):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    keep = np.flatnonzero(np.isfinite(values))
    if len(keep):
        order = keep[np.argsort(values[keep], kind='stable')]
        m = family_size if family_size is not None else len(values)
        q = values[order] * m / np.arange(1, len(order)+1)
        result[order] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
    return result


def permutation_units(table, cfg, pool):
    size = cfg['block_size']
    origin = np.asarray(cfg['block_origin'], float)
    if not np.isfinite(size) or size <= 0 or origin.shape != (2,) or not np.isfinite(origin).all():
        raise ValueError('Finite positive block_size and explicit [x,y] block_origin required')
    blocks = np.floor((table[[cfg['x'],cfg['y']]].to_numpy(float) - origin) / size).astype(np.int64)
    target = table[cfg['target_label_key']]
    strata = cfg.get('strata_keys', [])
    groups, excluded = {}, []
    for i in pool:
        if pd.isna(target.iloc[i]):
            excluded.append({'cell_key':table.iloc[i][cfg['cell_key']],'reason':'missing_target_label_preserved'})
            continue
        if strata and table.iloc[i][strata].isna().any():
            excluded.append({'cell_key':table.iloc[i][cfg['cell_key']],'reason':'missing_stratum_not_testable'})
            continue
        key = (str(table.iloc[i][cfg['roi_key']]), int(blocks[i,0]), int(blocks[i,1]),
               *(str(table.iloc[i][s]) for s in strata))
        groups.setdefault(key, []).append(int(i))
    groups = {key:np.asarray(indices, dtype=int) for key,indices in groups.items()}
    return groups, excluded


def permute_labels(labels, groups, rng):
    shuffled = labels.copy()
    for indices in groups.values():
        shuffled[indices] = rng.permutation(labels[indices])
    return shuffled


def anchor_scores(table, cfg, anchors, pool, graph, labels, hypothesis):
    target = hypothesis['target']
    metric = hypothesis['metric']
    coords = table[[cfg['x'],cfg['y']]].to_numpy(float)
    roi = table[cfg['roi_key']].to_numpy()
    scores = np.full(len(anchors), np.nan)
    for row, anchor in enumerate(anchors):
        # Missing target labels remain excluded evidence, not negatives.
        neighbors = graph[int(anchor)]
        known = neighbors[pd.notna(labels[neighbors])]
        candidate = pool[(roi[pool] == roi[anchor]) & (pool != anchor)]
        candidate = candidate[pd.notna(labels[candidate])]
        if not len(candidate):
            continue
        if metric == 'nearest':
            targets = candidate[labels[candidate] == target]
            scores[row] = min(float(np.linalg.norm(coords[targets]-coords[anchor], axis=1).min()), cfg['radius']) if len(targets) else cfg['radius']
        elif metric == 'proximity':
            scores[row] = float(np.any(labels[known] == target))
        elif metric == 'cooccurrence':
            scores[row] = float(np.count_nonzero(labels[known] == target))
        elif metric == 'local_enrichment':
            if len(known):
                scores[row] = float(np.mean(labels[known] == target) - np.mean(labels[candidate] == target))
        else:
            raise ValueError(f'Unknown metric: {metric}')
    return scores


def aggregate(table, cfg, anchors, scores):
    frame = table.iloc[anchors][[cfg['roi_key'], cfg['biological_unit_key']]].copy()
    frame['score'] = scores
    rois = frame.groupby([cfg['biological_unit_key'], cfg['roi_key']], sort=True)['score'].mean()
    units = rois.groupby(level=0).mean()
    return float(units.mean()) if units.notna().any() else np.nan, units


def validate_hypotheses(cfg):
    hypotheses = cfg['hypotheses']
    ids = [h['id'] for h in hypotheses]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('Distinct predeclared hypothesis IDs required')
    assigned = [x for group in cfg['families'].values() for x in group]
    if sorted(assigned) != sorted(ids):
        raise ValueError('Families must partition exactly the predeclared hypotheses')
    for h in hypotheses:
        if h['metric'] not in {'nearest','proximity','cooccurrence','local_enrichment'} or h['alternative'] not in {'greater','less','two-sided'}:
            raise ValueError('Explicit supported metric and test direction required')
        if h['id'] not in cfg['families'].get(h['family'], []):
            raise ValueError('Hypothesis family conflicts with declared family inventory')
        if h['alternative'] == 'two-sided' and not np.isfinite(h.get('null_center', np.nan)):
            raise ValueError('Two-sided tests require a predeclared finite null_center')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    validate_hypotheses(cfg)
    if not isinstance(cfg['permutations'], int) or cfg['permutations'] < 1:
        raise ValueError('Positive integer permutation count required')
    cfg['phenotype_key'] = cfg['target_label_key']
    table = cell_table(cfg['cells_csv'], cfg, labels=True)
    anchors = np.flatnonzero(mask(table, cfg['anchor_mask_key']))
    pool = np.flatnonzero(mask(table, cfg['candidate_mask_key']))
    if not len(anchors) or not len(pool):
        raise ValueError('Nonempty fixed anchors and complete candidate pool required')
    groups, excluded = permutation_units(table, cfg, pool)
    labels = table[cfg['target_label_key']].to_numpy(object)
    # Cells with missing strata are not eligible evidence for these conditional tests.
    effective_pool = np.asarray(sorted(i for v in groups.values() for i in v), dtype=int)
    graph = neighborhoods(table, cfg, anchors, effective_pool, radius=cfg['radius'])
    run = Run(args.out_dir, cfg, [args.config,cfg['cells_csv']], __file__)
    if run.resumed:
        print('PASS complete-artifact resume')
        return
    diagnostics, margins = [], []
    for key,indices in groups.items():
        counts = pd.Series(labels[indices]).value_counts().to_dict()
        diagnostics.append({'permutation_unit':json.dumps(key),'n_cells':len(indices),
                            'n_labels':len(counts),'status':'exchangeable' if len(counts)>1 else 'NOT_TESTABLE',
                            'reason':'' if len(counts)>1 else 'no_label_variation'})
        for label,count in counts.items():
            margins.append({'permutation_unit':json.dumps(key),'label':label,'n_cells':int(count)})
    observed, supports, null = [], [], np.empty((cfg['permutations'],len(cfg['hypotheses'])))
    for h in cfg['hypotheses']:
        scores = anchor_scores(table,cfg,anchors,effective_pool,graph,labels,h)
        value, units = aggregate(table,cfg,anchors,scores)
        observed.append(value)
        for unit,score in units.items():
            unit_anchor_mask = table.iloc[anchors][cfg['biological_unit_key']].eq(unit).to_numpy()
            supports.append({'hypothesis_id':h['id'],'biological_unit':unit,'observed':score,
                             'n_anchors':int(unit_anchor_mask.sum()),'n_evaluable_anchors':int(np.isfinite(scores[unit_anchor_mask]).sum()),
                             'n_rois':int(table.iloc[anchors[unit_anchor_mask]][cfg['roi_key']].nunique())})
    rng = np.random.default_rng(cfg['seed'])
    for b in range(cfg['permutations']):
        shuffled = permute_labels(labels,groups,rng)
        for indices in groups.values():
            if not pd.Series(shuffled[indices]).value_counts().equals(pd.Series(labels[indices]).value_counts()):
                # value_counts tie ordering is arbitrary; compare counts by label.
                if pd.Series(shuffled[indices]).value_counts().to_dict() != pd.Series(labels[indices]).value_counts().to_dict():
                    raise AssertionError('Permutation margins changed')
        if not np.array_equal(pd.isna(shuffled),pd.isna(labels)):
            raise AssertionError('Missing target labels were moved or filled')
        for j,h in enumerate(cfg['hypotheses']):
            scores = anchor_scores(table,cfg,anchors,effective_pool,graph,shuffled,h)
            null[b,j], _ = aggregate(table,cfg,anchors,scores)
    rows = []
    for j,h in enumerate(cfg['hypotheses']):
        distribution = null[:,j]
        finite = distribution[np.isfinite(distribution)]
        variable_target_groups = [v for v in groups.values() if np.any(labels[v] == h['target']) and np.any(labels[v] != h['target'])]
        status, reason = 'PASS', ''
        if not np.isfinite(observed[j]) or len(finite) != cfg['permutations']:
            status, reason = 'NOT_TESTABLE', 'no_evaluable_anchors'
        elif not variable_target_groups:
            status, reason = 'NOT_TESTABLE', 'no_exchangeable_target_labels'
        elif np.ptp(finite) < 1e-12:
            status, reason = 'NOT_TESTABLE', 'degenerate_statistic_under_declared_null'
        mean = float(finite.mean()) if len(finite) else None
        std = float(finite.std(ddof=1)) if len(finite)>1 else None
        p, extreme = None, None
        if status == 'PASS':
            if h['alternative'] == 'greater':
                extreme = int(np.count_nonzero(finite >= observed[j]))
            elif h['alternative'] == 'less':
                extreme = int(np.count_nonzero(finite <= observed[j]))
            else:
                center = h['null_center']
                extreme = int(np.count_nonzero(np.abs(finite-center) >= abs(observed[j]-center)))
            p = (extreme+1)/(cfg['permutations']+1)
        rows.append({'hypothesis_id':h['id'],'family':h['family'],'metric':h['metric'],'target':h['target'],
            'alternative':h['alternative'],'status':status,'reason':reason,'observed':observed[j],
            'null_mean':mean,'null_sd':std,'null_q025':float(np.quantile(finite,.025)) if len(finite) else None,
            'null_q975':float(np.quantile(finite,.975)) if len(finite) else None,
            'effect_observed_minus_null':observed[j]-mean if mean is not None else None,
            'standardized_effect':(observed[j]-mean)/std if std and std>0 else None,
            'extreme_permutations':extreme,'permutations':cfg['permutations'],'p_empirical':p})
    summary = pd.DataFrame(rows)
    summary['q_bh'] = np.nan
    for family, ids in cfg['families'].items():
        keep = summary.family.eq(family)
        summary.loc[keep,'q_bh'] = bh(summary.loc[keep,'p_empirical'],family_size=len(ids))
    summary.to_csv(args.out_dir/'null_summary.csv',index=False)
    pd.DataFrame(diagnostics,columns=['permutation_unit','n_cells','n_labels','status','reason']).to_csv(args.out_dir/'permutation_units.csv',index=False)
    pd.DataFrame(margins,columns=['permutation_unit','label','n_cells']).to_csv(args.out_dir/'label_margins.csv',index=False)
    pd.DataFrame(excluded,columns=['cell_key','reason']).to_csv(args.out_dir/'excluded_missing.csv',index=False)
    pd.DataFrame(supports).to_csv(args.out_dir/'biological_unit_support.csv',index=False)
    np.savez_compressed(args.out_dir/'null_distributions.npz',null=null,hypothesis_ids=np.asarray([h['id'] for h in cfg['hypotheses']]),observed=np.asarray(observed))
    run.finish({'status':'PASS','scope':'conditional_fixed_anchor_label_randomization',
                'margins_preserved':True,'missing_labels_preserved':True,'cross_roi_edges':0},
               evidence_grade='exploratory',aggregation='equal ROI means within unit, equal biological-unit means',
               group_inference='NOT_IMPLEMENTED',not_testable_hypotheses=int(summary.status.eq('NOT_TESTABLE').sum()))
    print('PASS spatial null; untestable hypotheses retained explicitly')


if __name__ == '__main__':
    main()
