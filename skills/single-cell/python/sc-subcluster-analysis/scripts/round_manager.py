"""Immutable review rounds and exact-list manual or program-supported automatic decisions; no clustering engine."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

HELPER = Path(__file__).resolve().parents[3] / 'workflow/workflow-sc-analysis/scripts/sc_contracts.py'
spec = importlib.util.spec_from_file_location('sc_contracts_rounds', HELPER)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
ACTIONS = {'retain', 'exclude_foreign', 'exclude_mixed', 'review'}
EXCLUDE = {'exclude_foreign', 'exclude_mixed'}
EVIDENCE_RULE = 'foreign_program_v1'
PLATFORMS = {'scrna', 'snrna', 'xenium_segmented', 'segmented_cells'}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def evidence_files(path):
    """Bind portable evidence manifests to the tables they actually describe."""
    path = Path(path).resolve()
    if path.stat().st_size == 0:
        raise ValueError(f'Empty evidence file: {path}')
    result = {str(path): file_hash(path)}
    if path.suffix == '.json':
        manifest = read(path)
        for item, expected in manifest.get('outputs', {}).items():
            if isinstance(expected, str) and re.fullmatch('[0-9a-f]{64}', expected):
                artifact = (path.parent / item).resolve()
                if file_hash(artifact) != expected:
                    raise ValueError(f'Evidence manifest is stale: {artifact}')
                result[str(artifact)] = expected
    return result


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, ensure_ascii=False)
        f.write('\n')


def stamp():
    return datetime.now(timezone.utc).isoformat()


def name(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', value):
        raise ValueError('Use a simple nonempty round/proposal ID (no paths)')
    return value


def frame(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def count_identity(a):
    if a.n_obs == 0 or a.n_vars == 0 or not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError('Nonempty unique cell/gene IDs required')
    if any(str(k) == '' for k in list(a.obs_names) + list(a.var_names)):
        raise ValueError('Empty key')
    c.validate_matrix(a.layers['counts'], integer=True)
    return {'cells':c.key_fingerprint(a.obs_names), 'genes':c.key_fingerprint(a.var_names),
            'counts':c.matrix_fingerprint(a.layers['counts']), 'shape':list(a.shape)}


def same_counts(a, b):
    count_identity(a); count_identity(b)
    if set(a.obs_names) != set(b.obs_names) or set(a.var_names) != set(b.var_names):
        raise ValueError('Cell/gene universe differs from preceding retained counts')
    aligned = b[a.obs_names, a.var_names]
    if c.matrix_fingerprint(a.layers['counts']) != c.matrix_fingerprint(aligned.layers['counts']):
        raise ValueError('Counts changed')


def require_v2(meta):
    if meta.get('schema') != 'subcluster_loop_v2':
        raise ValueError('Legacy v1 runs are read-only; create a new v2 run and fresh approval')


def harmony_batch(a, key):
    if not key or key not in a.obs or a.obs[key].isna().any() or a.obs[key].astype(str).str.strip().eq('').any():
        raise ValueError('Harmony requires an explicit complete nonempty batch_key')
    if a.obs[key].nunique() < 2:
        raise ValueError('Harmony requires at least two retained batches')


def init(run, source, context):
    run, source = Path(run).resolve(), Path(source).resolve()
    if run.exists():
        raise FileExistsError(run)
    if context.get('platform') not in PLATFORMS:
        raise ValueError('Cell exclusion requires sc/snRNA or documented segmented cells; spot/bin -> mixture review')
    for key in ('species', 'tissue', 'target_lineage', 'counts_origin', 'sample_key', 'gene_scope'):
        if not context.get(key):
            raise ValueError(f'Declare {key}')
    if 'unit_key' not in context:
        raise ValueError('Declare unit_key (null explicitly means unavailable)')
    if context.get('integration_method') != 'harmony':
        raise ValueError('Subcluster v2 requires integration_method=harmony')
    context = dict(context)
    context.setdefault('cleanup_mode', 'manual_per_round')
    if context['cleanup_mode'] not in {'manual_per_round', 'auto_supported'}:
        raise ValueError('Invalid cleanup_mode')
    if context['cleanup_mode'] == 'auto_supported':
        request = Path(context.get('cleanup_request_record', '')).resolve()
        if not request.is_file() or not request.read_text().strip():
            raise ValueError('Automatic cleanup requires the initial explicit cleanup request record')
        context['cleanup_request_record'] = str(request)
        context['cleanup_request_sha256'] = file_hash(request)
    a = ad.read_h5ad(source)
    harmony_batch(a, context.get('batch_key'))
    identity = count_identity(a)
    for key in [context['sample_key']] + ([context['unit_key']] if context['unit_key'] else []):
        if key not in a.obs or a.obs[key].isna().any():
            raise ValueError(f'Missing/incomplete {key}')
    run.mkdir(parents=True)
    initial = pd.DataFrame({'obs_name':a.obs_names, 'sample':a.obs[context['sample_key']].astype(str).to_numpy(),
                            'biological_unit':a.obs[context['unit_key']].astype(str).to_numpy() if context['unit_key'] else '',
                            'destination':'retained', 'removed_round':'', 'proposal_id':'', 'reason':'', 'history':'[]'})
    initial.to_csv(run/'initial_ledger.csv', index=False)
    write(run/'run.json', {'schema':'subcluster_loop_v2', 'created_at':stamp(), 'context':context,
        'source':str(source), 'source_sha256':file_hash(source), 'identity':identity,
        'implementation':{str(p):file_hash(p) for p in (Path(__file__).resolve(), HELPER)},
        'initial_ledger_sha256':file_hash(run/'initial_ledger.csv')})
    return read(run/'run.json')


def run_contract(run):
    run = Path(run).resolve(); meta = read(run/'run.json')
    if file_hash(meta['source']) != meta['source_sha256'] or file_hash(run/'initial_ledger.csv') != meta['initial_ledger_sha256']:
        raise ValueError('Initial source/ledger changed')
    if meta['context']['platform'] not in PLATFORMS:
        raise ValueError('Spot/bin cell exclusion forbidden')
    if meta['context'].get('cleanup_mode', 'manual_per_round') == 'auto_supported':
        if file_hash(meta['context']['cleanup_request_record']) != meta['context']['cleanup_request_sha256']:
            raise ValueError('Initial cleanup authorization changed')
    return run, meta


def round_contract(run, round_id):
    run, meta = run_contract(run)
    directory = run/'rounds'/name(round_id)
    frozen = read(directory/'frozen.json')
    token = frozen.pop('proposal_hash')
    if digest(frozen) != token:
        raise ValueError('Frozen proposal changed; old approval invalid')
    frozen['proposal_hash'] = token
    if digest(meta) != frozen['run_hash'] or file_hash(frozen['source']) != frozen['source_sha256']:
        raise ValueError('Source/run contract changed; old approval invalid')
    for artifact, expected in frozen['artifacts'].items():
        if file_hash(directory/artifact) != expected:
            raise ValueError(f'Frozen {artifact} changed')
    for path, expected in frozen['evidence'].items():
        if file_hash(path) != expected:
            raise ValueError(f'Evidence changed: {path}')
    if frozen['previous_round']:
        previous, _, prev = round_contract(run, frozen['previous_round'])
        applied = checked_apply(previous, prev)
        if digest(applied) != frozen['previous_apply_hash']:
            raise ValueError('Previous applied round changed')
    return directory, meta, frozen


def checked_apply(directory, frozen):
    data = read(directory/'applied.json')
    decision, decision_path = checked_decision(directory, frozen)
    if data['proposal_hash'] != frozen['proposal_hash'] or data.get('decision_sha256', data.get('approval_sha256')) != file_hash(decision_path):
        raise ValueError('Applied approval does not match frozen round')
    for artifact, expected in data['outputs'].items():
        if file_hash(directory/artifact) != expected:
            raise ValueError('Applied output changed')
    return data


def require_text(mapping, fields):
    for field in fields:
        if not isinstance(mapping.get(field), str) or not mapping[field].strip():
            raise ValueError(f'Required substantive review field: {field}')


def validate_program(proposal, a):
    """Validate an auditable program contract, not biological truth or a score cutoff."""
    program = proposal.get('foreign_program')
    if not isinstance(program, dict):
        raise ValueError('Exclusion requires foreign_program; two screening markers are insufficient')
    require_text(program, ('id', 'source', 'specificity', 'complementarity',
                          'cell_coordination', 'target_comparison', 'coverage_limitations'))
    genes = program.get('genes')
    if not isinstance(genes, list) or not genes or any(not isinstance(g,str) or not g.strip() for g in genes) or len(set(genes)) != len(genes):
        raise ValueError('Unique program feature IDs required')
    # Program genes are var_names, as in the collector; never guess symbol aliases.
    mapping = program.get('feature_mapping')
    if mapping != {g:g for g in genes}:
        raise ValueError('Program genes must use audited var_names with identity feature_mapping')
    measured = set(genes) & set(a.var_names)
    coverage = program.get('measured_genes')
    if not isinstance(coverage,list) or len(coverage) != len(set(coverage)) or set(coverage) != measured:
        raise ValueError('Program measured coverage does not match source features')
    markers = proposal.get('foreign_markers')
    if not isinstance(markers,list) or not set(markers) <= measured or len(set(markers)) < 2:
        raise ValueError('Screening markers must map to measured genes of this foreign program')
    if program.get('coordinated_program_supported') is not True:
        raise ValueError('A complete coordinated program review is required, not two markers alone')
    if not isinstance(proposal.get('evidence_files'),list) or not proposal['evidence_files']:
        raise ValueError('Program evidence files required')


def validate_auto_record(record, directory, frozen):
    if frozen.get('evidence_rule_version') != EVIDENCE_RULE:
        raise ValueError('Legacy proposal lacks current program evidence; freeze a new run')
    if record.get('round') != frozen['round'] or record.get('proposal_hash') != frozen['proposal_hash']:
        raise ValueError('Automatic review is stale or from a different round')
    require_text(record, ('reviewer',))
    proposals = {p['id']:p for p in read(directory/'proposals.json') if p['action'] in EXCLUDE}
    decisions = record.get('decisions')
    if not isinstance(decisions,list) or len(decisions) != len(proposals) or {d.get('proposal_id') for d in decisions} != set(proposals):
        raise ValueError('Review every exclusion proposal exactly once')
    selected = []
    for d in decisions:
        require_text(d, ('reason', 'alternatives', 'target_loss_assessment'))
        if d.get('decision') not in {'supported', 'review'}:
            raise ValueError('Automatic decision must be supported or review')
        if d['decision'] == 'supported':
            if any(d.get(k) is not False for k in ('target_loss_concern','biological_mixed_state_unresolved','target_subtype_unresolved','partial_cluster_only')):
                raise ValueError('Unresolved target loss/biology or partial-cluster support cannot authorize whole-cluster removal')
            if d.get('removal_basis') not in {'technical_mixture','outside_target_scope'}:
                raise ValueError('Removal needs technical-mixture or outside-target-scope evidence')
            loss = d.get('target_positive')
            n = int(frame(directory/'cells.csv').proposal_id.eq(d['proposal_id']).sum())
            if not isinstance(loss,dict) or type(loss.get('count')) is not int or not 0 <= loss['count'] <= n:
                raise ValueError('Quantify target-positive count and fraction')
            if not isinstance(loss.get('fraction'),(int,float)) or not np.isclose(loss['fraction'],loss['count']/n):
                raise ValueError('Target-positive fraction must match frozen cluster size')
            require_text(loss, ('definition',))
            if proposals[d['proposal_id']]['action'] == 'exclude_mixed' and not proposals[d['proposal_id']].get('target_cell_loss'):
                raise ValueError('Mixed exclusion needs frozen target loss assessment')
            selected.append(d['proposal_id'])
    return sorted(selected)


def auto_review(run, round_id, expected_hash, record):
    directory, meta, frozen = round_contract(run, round_id)
    require_v2(meta)
    if meta['context'].get('cleanup_mode','manual_per_round') != 'auto_supported':
        raise ValueError('Automatic review requires explicit auto_supported cleanup authorization')
    if (Path(run)/'stop.json').exists() or expected_hash != frozen['proposal_hash']:
        raise ValueError('Stopped run or stale automatic review hash')
    record = Path(record).resolve()
    reviewed = read(record)
    selected = validate_auto_record(reviewed, directory, frozen)
    cells = frame(directory/'cells.csv')
    decision = {'kind':'automatic_program_review', 'round':round_id, 'proposal_hash':expected_hash,
        'source_sha256':frozen['source_sha256'], 'record':str(record), 'record_sha256':file_hash(record),
        'approved_proposals':selected, 'approved_cells_sha256':digest(sorted(cells.loc[cells.proposal_id.isin(selected),'obs_name'])),
        'cleanup_request_sha256':meta['context']['cleanup_request_sha256'], 'recorded_at':stamp()}
    path = directory/'auto_review.json'
    if path.exists():
        old = checked_auto_review(directory, frozen)
        if all(old[k] == decision[k] for k in decision if k != 'recorded_at'):
            return old
        raise FileExistsError('Automatic review immutable; changed evidence/decision requires a fresh run')
    write(path, decision)
    return decision


def checked_auto_review(directory, frozen):
    decision = read(directory/'auto_review.json')
    if decision.get('kind') != 'automatic_program_review' or decision.get('proposal_hash') != frozen['proposal_hash'] or decision.get('round') != frozen['round'] or decision.get('source_sha256') != frozen['source_sha256']:
        raise ValueError('Automatic review is stale or from a different round')
    if file_hash(decision['record']) != decision['record_sha256']:
        raise ValueError('Automatic review record changed')
    selected = validate_auto_record(read(decision['record']), directory, frozen)
    cells = frame(directory/'cells.csv')
    if selected != decision['approved_proposals'] or digest(sorted(cells.loc[cells.proposal_id.isin(selected),'obs_name'])) != decision['approved_cells_sha256']:
        raise ValueError('Automatic exact-cell list changed')
    return decision


def checked_decision(directory, frozen):
    meta = read(directory.parent.parent/'run.json')
    if meta['context'].get('cleanup_mode','manual_per_round') == 'auto_supported':
        decision = checked_auto_review(directory, frozen)
        if decision['cleanup_request_sha256'] != meta['context']['cleanup_request_sha256']:
            raise ValueError('Automatic authorization changed')
        if not decision['approved_proposals']:
            raise ValueError('No supported exclusions; retain cells and stop after review')
        return decision, directory/'auto_review.json'
    return checked_approval(directory, frozen), directory/'approval.json'


def freeze(run, round_id, source, cluster_key, proposals, previous=None, recomputed=False, execution_manifest=None):
    run, meta = run_contract(run)
    require_v2(meta)
    if (run/'stop.json').exists():
        raise ValueError('Loop stopped; create a new run to resume with a new contract')
    directory = run/'rounds'/name(round_id)
    if directory.exists():
        raise FileExistsError(directory)
    existing = sorted((run/'rounds').glob('*/frozen.json')) if (run/'rounds').exists() else []
    if bool(existing) != bool(previous):
        raise ValueError('A subsequent round requires the preceding applied round')
    if previous and any(read(p).get('previous_round') == previous for p in existing):
        raise ValueError('Linear rounds only; predecessor already has a successor')
    source = Path(source).resolve()
    a = ad.read_h5ad(source); identity = count_identity(a)
    if cluster_key not in a.obs or a.obs[cluster_key].isna().any():
        raise ValueError('Complete explicit cluster key required')
    if recomputed:
        check_recluster(a, cluster_key, meta['context']['batch_key'])
        execution_manifest = Path(execution_manifest or source.parent/'run_params.json').resolve()
        check_execution(a, cluster_key, meta['context']['batch_key'], read(execution_manifest))
    if previous:
        pdirectory, _, pfrozen = round_contract(run, previous)
        applied = checked_apply(pdirectory, pfrozen)
        expected = ad.read_h5ad(pdirectory/'retained_counts.h5ad')
        if not recomputed:
            raise ValueError('Next round after filtering requires actual expression recomputation')
        ledger = frame(pdirectory/'ledger.csv')
    else:
        applied = None
        expected = ad.read_h5ad(meta['source'])
        ledger = frame(run/'initial_ledger.csv')
    same_counts(a, expected)
    batch_key = meta['context']['batch_key']
    if batch_key not in a.obs or a.obs[batch_key].astype(str).tolist() != expected.obs.loc[a.obs_names, batch_key].astype(str).tolist():
        raise ValueError('Batch assignment changed')
    indexed_ledger = ledger.set_index('obs_name').loc[a.obs_names]
    for key, column in ((meta['context']['sample_key'], 'sample'), (meta['context']['unit_key'], 'biological_unit')):
        if key and (key not in a.obs or a.obs[key].astype(str).tolist() != indexed_ledger[column].tolist()):
            raise ValueError('Sample/biological-unit assignment changed')
    labels = a.obs[cluster_key].astype(str)
    if len(proposals) != len(set(p['id'] for p in proposals)) or len(proposals) != len(set(str(p['cluster']) for p in proposals)):
        raise ValueError('Unique proposal IDs and one action per cluster required')
    if {str(p['cluster']) for p in proposals} != set(labels):
        raise ValueError('Actions must cover every observed cluster exactly once')
    evidence = evidence_files(execution_manifest) if recomputed else {}
    if recomputed:
        execution = read(execution_manifest)
        for path in (execution['output']['path'], execution['input']['input_path']):
            evidence[str(Path(path).resolve())] = file_hash(path)
    for p in proposals:
        name(p['id'])
        if p['action'] not in ACTIONS or not p.get('reason'):
            raise ValueError('Action and reason required')
        if p['action'] in EXCLUDE:
            for key in ('target_program', 'foreign_lineage', 'foreign_markers', 'unit_support', 'qc_review', 'interpretation', 'evidence_files'):
                if not p.get(key):
                    raise ValueError(f'Exclusion proposal requires {key}')
            validate_program(p, a)
            if len(set(p['foreign_markers'])) < 2:
                raise ValueError('At least two coordinated same-foreign-lineage markers required')
            if p['action'] == 'exclude_mixed' and not p.get('target_cell_loss'):
                raise ValueError('Mixed whole-cluster removal requires explicit target-cell loss assessment')
        for path in p.get('evidence_files', []):
            evidence.update(evidence_files(path))
    lookup = {str(p['cluster']):p for p in proposals}
    cells = pd.DataFrame({'obs_name':a.obs_names, 'cluster':labels.to_numpy(),
        'namespace':round_id+':'+cluster_key, 'proposal_id':[lookup[v]['id'] for v in labels],
        'action':[lookup[v]['action'] for v in labels]})
    cells = cells.merge(ledger[['obs_name','sample','biological_unit']], on='obs_name', validate='one_to_one')
    proposed = set(cells.loc[cells.action.isin(EXCLUDE), 'obs_name'])
    impact = impact_table(ledger, proposed)
    by_cell = cells.set_index('obs_name')
    for i, row in ledger.iterrows():
        if row.obs_name in by_cell.index:
            current = by_cell.loc[row.obs_name]
            history = json.loads(row.history)
            history.append({'round': round_id, 'namespace': round_id+':'+cluster_key,
                            'cluster': current.cluster, 'proposal_id': current.proposal_id,
                            'action': current.action})
            ledger.at[i, 'history'] = json.dumps(history, separators=(',', ':'))
    directory.mkdir(parents=True)
    cells.to_csv(directory/'cells.csv', index=False)
    ledger.to_csv(directory/'review_ledger.csv', index=False)
    impact.to_csv(directory/'proposed_impact.csv', index=False)
    write(directory/'proposals.json', proposals)
    frozen = {'schema':'subcluster_round_v2','evidence_rule_version':EVIDENCE_RULE,'round':round_id,'created_at':stamp(),'run_hash':digest(meta),
        'previous_round':previous,'previous_apply_hash':digest(applied) if applied else None,
        'source':str(source),'source_sha256':file_hash(source),'identity':identity,'cluster_key':cluster_key,
        'namespace':round_id+':'+cluster_key,'entry_mode':'expression_recomputed' if recomputed else 'existing_cluster_review',
        'artifacts':{p.name:file_hash(p) for p in directory.iterdir()},'evidence':evidence,
        'approval_status':'pending','scientific_status':'awaiting_review'}
    frozen['proposal_hash'] = digest(frozen)
    frozen['harmony_context'] = {'integration_method':'harmony', 'batch_key':meta['context']['batch_key']}
    frozen['proposal_hash'] = digest({k:v for k,v in frozen.items() if k != 'proposal_hash'})
    write(directory/'frozen.json', frozen)
    return frozen


def impact_table(ledger, selected):
    data = ledger.copy()
    data['newly_removed'] = data.obs_name.isin(selected).astype(int)
    data['previously_removed'] = data.destination.eq('excluded').astype(int)
    result = data.groupby('sample', dropna=False).agg(initial=('obs_name','size'),
        previously_removed=('previously_removed','sum'),newly_removed=('newly_removed','sum')).reset_index()
    result['cumulative_removed'] = result.previously_removed + result.newly_removed
    result['retained'] = result.initial - result.cumulative_removed
    return result


def approve(run, round_id, expected_hash, selected, approver, record, simulated=False):
    directory, meta, frozen = round_contract(run, round_id)
    require_v2(meta)
    if meta['context'].get('cleanup_mode', 'manual_per_round') != 'manual_per_round':
        raise ValueError('Automatic mode uses auto-review, never human/simulated approval')
    if expected_hash != frozen['proposal_hash']:
        raise ValueError('Approval hash differs from frozen proposal')
    if (Path(run)/'stop.json').exists():
        raise ValueError('Loop stopped')
    proposals = read(directory/'proposals.json')
    allowed = {p['id'] for p in proposals if p['action'] in EXCLUDE}
    selected = sorted(set(selected))
    if not selected or not set(selected) <= allowed or not approver.strip():
        raise ValueError('Select explicit exclusion proposal IDs and record the human approver')
    record = Path(record).resolve()
    if not record.read_text().strip():
        raise ValueError('Save the actual user approval text/reference before recording')
    approval = {'round':round_id, 'proposal_hash':expected_hash, 'source_sha256':frozen['source_sha256'],
        'approved_proposals':selected, 'approver':approver, 'record':str(record), 'record_sha256':file_hash(record),
        'recorded_at':stamp(), 'simulated':bool(simulated)}
    cells = frame(directory/'cells.csv')
    approval['approved_cells_sha256'] = digest(sorted(cells.loc[cells.proposal_id.isin(selected), 'obs_name']))
    path = directory/'approval.json'
    if path.exists():
        old = checked_approval(directory, frozen)
        if all(old[k] == approval[k] for k in approval if k != 'recorded_at'):
            return old  # Same approved list needs no repeat question or mutation.
        raise FileExistsError('Approval immutable; freeze a fresh round/run for a changed list')
    write(path, approval)
    return approval


def checked_approval(directory, frozen):
    approval = read(directory/'approval.json')
    if approval['round'] != frozen['round'] or approval['proposal_hash'] != frozen['proposal_hash'] or approval['source_sha256'] != frozen['source_sha256']:
        raise ValueError('Approval is stale or from a different round')
    if file_hash(approval['record']) != approval['record_sha256']:
        raise ValueError('Approval record changed')
    allowed = {p['id'] for p in read(directory/'proposals.json') if p['action'] in EXCLUDE}
    if not approval['approved_proposals'] or not set(approval['approved_proposals']) <= allowed:
        raise ValueError('Approval contains invalid proposals')
    cells = frame(directory/'cells.csv')
    selected = sorted(cells.loc[cells.proposal_id.isin(approval['approved_proposals']), 'obs_name'])
    if digest(selected) != approval['approved_cells_sha256']:
        raise ValueError('Approved exact-cell list changed')
    return approval


def apply(run, round_id, allow_simulated=False):
    directory, meta, frozen = round_contract(run, round_id)
    require_v2(meta)
    approval, decision_path = checked_decision(directory, frozen)
    if approval.get('simulated', False) and not (allow_simulated and meta['context'].get('synthetic') is True):
        raise ValueError('Simulated approval only allowed for explicitly synthetic runs')
    if (directory/'applied.json').exists():
        return checked_apply(directory, frozen)
    if (Path(run)/'stop.json').exists():
        raise ValueError('Loop stopped')
    cells = frame(directory/'cells.csv')
    selected = cells.proposal_id.isin(approval['approved_proposals'])
    excluded = set(cells.loc[selected, 'obs_name'])
    source = ad.read_h5ad(frozen['source'])
    keep = ~source.obs_names.isin(excluded)
    if not keep.any():
        raise ValueError('Cannot remove every cell')
    previous = frozen['previous_round']
    ledger = frame(directory/'review_ledger.csv')
    proposals = {p['id']:p for p in read(directory/'proposals.json')}
    by_cell = cells.set_index('obs_name')
    for i, row in ledger.iterrows():
        if row.obs_name not in by_cell.index:
            continue
        current = by_cell.loc[row.obs_name]
        if row.obs_name in excluded:
            ledger.loc[i,['destination','removed_round','proposal_id','reason']] = [
                'excluded',round_id,current.proposal_id,proposals[current.proposal_id]['reason']]
    impact = impact_table(frame(Path(run)/'rounds'/previous/'ledger.csv' if previous else Path(run)/'initial_ledger.csv'), excluded)
    # Build a fresh counts-only derivative; old obs fields stay historical evidence.
    counts = source.layers['counts'][keep].copy()
    result = ad.AnnData(X=counts.copy(), obs=source.obs.loc[keep].copy(), var=source.var.copy())
    result.layers['counts'] = counts.copy()
    result.obs['history_'+round_id+'_cluster'] = result.obs.pop(frozen['cluster_key']).astype(str)
    result.var.drop(columns=['highly_variable','means','variances','variances_norm','dispersions','dispersions_norm'], errors='ignore', inplace=True)
    if 'spatial' in source.obsm:
        result.obsm['spatial'] = source.obsm['spatial'][keep].copy()
    result.uns['matrix_sources'] = {'counts':{'kind':'raw_counts','path':frozen['source'],
        'sha256':frozen['source_sha256'],'matrix':'layers/counts','proposal_hash':frozen['proposal_hash']}}
    result.uns['subcluster_loop'] = {'round':round_id, 'requires_expression_recluster':True,
        'proposal_hash':frozen['proposal_hash'], 'decision_sha256':file_hash(decision_path), 'decision_kind':approval.get('kind', 'user_approval')}
    with tempfile.TemporaryDirectory(prefix='.apply-', dir=directory) as temp:
        temp = Path(temp)
        result.write_h5ad(temp/'retained_counts.h5ad')
        reopened = ad.read_h5ad(temp/'retained_counts.h5ad')
        same_counts(source[keep], reopened)
        ledger.to_csv(temp/'ledger.csv', index=False)
        impact.to_csv(temp/'approved_impact.csv', index=False)
        cells.loc[selected].to_csv(temp/'excluded_cells.csv', index=False)
        outputs = {p.name:file_hash(p) for p in temp.iterdir()}
        if any((directory/p).exists() for p in outputs):
            raise FileExistsError('Incomplete prior apply output; inspect before recovery')
        for p in temp.iterdir():
            shutil.move(str(p), directory/p.name)
    applied = {'round':round_id,'proposal_hash':frozen['proposal_hash'],'decision_sha256':file_hash(decision_path), 'decision_kind':approval.get('kind', 'user_approval'),
        'technical_status':'PASS','scientific_status':'requires_recluster_and_review',
        'removed':len(excluded),'retained':result.n_obs,'outputs':outputs}
    write(directory/'applied.json', applied)
    try:
        harmony_batch(result, meta['context']['batch_key'])
    except ValueError as exc:
        stop(run, round_id, 'technical_failure', str(exc))
    return applied


def check_recluster(a, cluster_key, batch_key=None):
    count_identity(a)
    provenance = a.uns.get('expression_recluster_provenance', {})
    if not all(bool(provenance.get('recomputed_'+key, False)) for key in ('hvg','pca','neighbors','umap','leiden')):
        raise ValueError('Incomplete recomputation provenance')
    if cluster_key not in a.obs or a.obs[cluster_key].isna().any():
        raise ValueError('Missing clusters')
    if 'X_umap' not in a.obsm or a.obsm['X_umap'].shape != (a.n_obs,2) or not np.isfinite(a.obsm['X_umap']).all():
        raise ValueError('Invalid current 2D UMAP')
    if 'X_pca' not in a.obsm or 'highly_variable' not in a.var or not a.var.highly_variable.any():
        raise ValueError('Missing actual PCA/HVG evidence')
    if batch_key is not None:
        harmony_batch(a, batch_key)
        if provenance.get('integration_method') != 'harmony' or provenance.get('batch_key') != batch_key or not provenance.get('recomputed_harmony'):
            raise ValueError('Missing Harmony recomputation provenance')
        latent = a.obsm.get('X_harmony')
        if latent is None or latent.shape != a.obsm['X_pca'].shape or not np.isfinite(latent).all():
            raise ValueError('Missing/invalid actual X_harmony')
        if provenance.get('harmony_representation_sha256') != c.matrix_fingerprint(latent) or provenance.get('harmony_pca_sha256') != c.matrix_fingerprint(a.obsm['X_pca']) or not provenance.get('harmony_parameters'):
            raise ValueError('Harmony representation/parameters fingerprint mismatch')
        if a.uns.get('neighbors', {}).get('params', {}).get('use_rep') != 'X_harmony':
            raise ValueError('Main graph must use X_harmony')
        if provenance.get('graph_sha256') != c.matrix_fingerprint(a.obsp['connectivities']):
            raise ValueError('Harmony graph fingerprint mismatch')
    c.validate_log(a)
    graph_key = a.uns.get('neighbors', {}).get('connectivities_key', 'connectivities')
    if graph_key not in a.obsp or a.obsp[graph_key].shape != (a.n_obs,a.n_obs):
        raise ValueError('Missing graph')
    graph = sparse.csr_matrix(a.obsp[graph_key]); c.validate_matrix(graph)
    if graph.nnz == 0 or (abs(graph-graph.T) > 1e-6).nnz:
        raise ValueError('Invalid graph')


def check_execution(a, cluster_key, batch_key, execution):
    params = execution.get('parameters', {})
    if execution.get('status') != 'PASS' or params.get('integration') != 'harmony' or params.get('batch_key') != batch_key or params.get('require_harmony') is not True:
        raise ValueError('Execution manifest requires successful mandatory Harmony')
    stages = execution.get('stages', [])
    if {s.get('stage') for s in stages} != {'preprocessed','pca','graph','clustered'} or not all(s.get('contract_sha256') and s.get('action') in {'recomputed','reused'} for s in stages) or not execution.get('input'):
        raise ValueError('Incomplete Harmony execution stages/input')
    output = execution.get('output', {})
    if not output.get('path') or file_hash(output['path']) != output.get('sha256'):
        raise ValueError('Execution output hash is stale/missing')
    input_info = execution['input']
    if not input_info.get('input_path') or file_hash(input_info['input_path']) != input_info.get('input_sha256'):
        raise ValueError('Execution input hash is stale/missing')
    same_counts(a, ad.read_h5ad(input_info['input_path']))
    computed = ad.read_h5ad(output['path'])
    check_recluster(computed, cluster_key, batch_key)
    same_counts(a, computed)
    aligned = computed[a.obs_names, a.var_names]
    for key in ('X_harmony','X_pca','X_umap'):
        if not np.array_equal(a.obsm[key], aligned.obsm[key]):
            raise ValueError(f'{key} differs from execution artifact')
    if c.matrix_fingerprint(a.obsp['connectivities']) != c.matrix_fingerprint(aligned.obsp['connectivities']):
        raise ValueError('Graph differs from execution artifact')
    if a.obs[cluster_key].astype(str).tolist() != aligned.obs[cluster_key].astype(str).tolist():
        raise ValueError('Clusters differ from execution artifact')


def validate(run, round_id, artifacts=None):
    directory, meta, frozen = round_contract(run, round_id)
    legacy = meta.get('schema') != 'subcluster_loop_v2'
    report = {'round':round_id,'technical_status':'PASS','scope':'frozen_proposal_integrity',
              'scientific_status':'awaiting_review','proposal_hash':frozen['proposal_hash']}
    if (directory/'applied.json').exists():
        applied = checked_apply(directory, frozen)
        report.update(scope='approved_filter_integrity',scientific_status=applied['scientific_status'])
    if artifacts is not None:
        try:
            report.update(validate_bundle(frozen, artifacts, cleanup_pdf=not legacy))
        except Exception as exc:
            report.update(technical_status='FAIL', scope='complete_round_artifacts',
                          scientific_status='not_determined', error=f'{type(exc).__name__}: {exc}')
            if legacy:
                raise
            reports = directory/'validation'; reports.mkdir(exist_ok=True)
            write(reports/(str(len(list(reports.glob('*.json')))+1).zfill(4)+'.json'), report)
            raise
    if legacy:
        return report
    # Append immutable validations; a technical PASS never changes the stop reason.
    reports = directory/'validation'; reports.mkdir(exist_ok=True)
    write(reports/(str(len(list(reports.glob('*.json')))+1).zfill(4)+'.json'),report)
    return report


def validate_bundle(frozen, artifacts, cleanup_pdf=True):
    from PIL import Image
    from matplotlib.colors import to_hex, to_rgb
    required = {'annotated_h5ad','hvg_h5ad','all_markers','top20','evidence','sidecar','palette','pngs','execution_manifest'}
    if not required <= artifacts.keys() or not artifacts['pngs']:
        raise ValueError('Incomplete round artifacts')
    paths = [Path(v).resolve() for k,v in artifacts.items() if k in required-{'pngs'}]
    paths += [Path(p).resolve() for p in artifacts['pngs']]
    if any(not p.is_file() or p.stat().st_size == 0 for p in paths):
        raise ValueError('Missing/empty round artifacts')
    a = ad.read_h5ad(artifacts['annotated_h5ad']); source = ad.read_h5ad(frozen['source'])
    same_counts(source,a)
    key = frozen['cluster_key']
    if not a.obs[key].astype(str).equals(source.obs.loc[a.obs_names,key].astype(str)):
        raise ValueError('Annotated clusters changed since freeze')
    if frozen['entry_mode'] == 'expression_recomputed':
        check_recluster(a,key,frozen.get('harmony_context', {}).get('batch_key'))
    for field in ('cell_identity','cell_subtype','cell_state','display_name','annotation_status'):
        if field not in a.obs or a.obs[field].isna().any():
            raise ValueError(f'Explicit {field} required (unknown allowed)')
    hvg = ad.read_h5ad(artifacts['hvg_h5ad'])
    if not a.obs_names.equals(hvg.obs_names) or not set(hvg.var_names) <= set(a.var_names):
        raise ValueError('HVG work not aligned')
    same_counts(a[:,hvg.var_names],hvg)
    sidecar = frame(artifacts['sidecar']).set_index('obs_name')
    sidecar = c.align_frame(sidecar,a.obs_names)
    for field in (key,'cell_identity','cell_subtype','cell_state','display_name','annotation_status'):
        if field not in sidecar or sidecar[field].tolist() != a.obs[field].astype(str).tolist():
            raise ValueError(f'Sidecar/H5AD disagree: {field}')
    palette = read(artifacts['palette'])
    # Namespace cluster labels by round; semantic labels share one union mapping.
    displayed = artifacts.get('displayed_labels')
    if not displayed or not set(displayed) <= palette.keys():
        raise ValueError('Declare every displayed label and cover them with the saved palette')
    fields = artifacts.get('categorical_fields', {})
    if not fields or any(field not in a.obs for field in fields):
        raise ValueError('Declare categorical_fields mapping plotted fields to label prefixes')
    expected_labels = {prefix+str(label) for field, prefix in fields.items() for label in a.obs[field].unique()}
    if not expected_labels <= set(displayed):
        raise ValueError('Displayed labels omit observed categories')
    colors = [to_hex(palette[k]).lower() for k in palette]
    if len(colors) != len(set(colors)) or any(min(to_rgb(v)) > .94 for v in colors):
        raise ValueError('Colors must be unique and readable on white')
    for png in map(Path,artifacts['pngs']):
        if png.suffix.lower() != '.png':
            raise ValueError('PNG only')
        with Image.open(png) as im:
            im.verify()
        if cleanup_pdf and png.with_suffix('.pdf').exists():
            png.with_suffix('.pdf').unlink()
    markers = pd.read_csv(artifacts['all_markers'])
    if not {'gene','score','logfoldchanges','pvals','pvals_adj','pct_in','pct_out'} <= set(markers):
        raise ValueError('Complete marker table columns missing')
    evidence = read(artifacts['evidence'])
    bound = evidence_files(artifacts['evidence'])
    for field, filename in (('all_markers','all_markers.csv.gz'), ('top20','top20.csv')):
        expected_path = str((Path(artifacts['evidence']).resolve().parent/filename).resolve())
        if expected_path not in bound or file_hash(artifacts[field]) != bound[expected_path]:
            raise ValueError(f'{field} differs from the evidence manifest')
    if evidence.get('schema') != 'annotation_evidence_v1' or evidence.get('source_sha256') != frozen['source_sha256']:
        raise ValueError('Evidence manifest must describe this exact frozen source')
    if evidence['config']['cluster_key'] != key:
        raise ValueError('Evidence uses a different cluster field')
    minimum = int(evidence['config'].get('min_detected_cells', 3))
    eligible = set(a.var_names[np.asarray((a.layers['log1p'] > 0).sum(axis=0)).ravel() >= minimum])
    group_column = 'group' if 'group' in markers else key
    sizes = a.obs[key].astype(str).value_counts()
    if len(sizes) >= 2 and sizes.min() >= 2 and eligible:
        if group_column not in markers or set(markers[group_column].astype(str)) != set(sizes.index):
            raise ValueError('Full marker table omits clusters')
        for _, table in markers.groupby(group_column, observed=True):
            if table.gene.duplicated().any() or set(table.gene) != eligible:
                raise ValueError('Full marker table omits/duplicates tested genes')
    elif not evidence.get('marker_skipped'):
        raise ValueError('Inadequate marker test must be explicitly recorded')
    execution = read(artifacts['execution_manifest'])
    if execution.get('status', execution.get('technical_status')) != 'PASS':
        raise ValueError('Execution manifest did not pass')
    if frozen['entry_mode'] == 'expression_recomputed':
        if not execution.get('stages') or not execution.get('input'):
            raise ValueError('Recomputation requires actual input and stage evidence, not flags alone')
        if frozen.get('harmony_context'):
            check_execution(a, key, frozen['harmony_context']['batch_key'], execution)
        output = execution.get('output', {})
        if not output.get('path') or file_hash(output['path']) != output.get('sha256'):
            raise ValueError('Execution output hash is stale/missing')
        computed = ad.read_h5ad(output['path'])
        same_counts(a, computed)
        aligned = computed[a.obs_names, a.var_names]
        if not np.array_equal(a.obsm['X_umap'], aligned.obsm['X_umap']):
            raise ValueError('UMAP differs from execution artifact')
        if a.obs[key].astype(str).tolist() != aligned.obs[key].astype(str).tolist():
            raise ValueError('Clusters differ from execution artifact')
        if c.matrix_fingerprint(a.obsp['connectivities']) != c.matrix_fingerprint(aligned.obsp['connectivities']):
            raise ValueError('Graph differs from execution artifact')
        if not all(stage.get('contract_sha256') and stage.get('action') in {'recomputed','reused'} for stage in execution['stages']):
            raise ValueError('Incomplete execution stage fingerprints')
    return {'scope':'complete_round_artifacts','outputs':{str(p):file_hash(p) for p in paths},
            'palette_and_png_validation':'PASS','visual_inspection':'agent_required',
            'scientific_status':'review_required'}


def stop(run, round_id, reason, note):
    directory, meta, frozen = round_contract(run, round_id)
    require_v2(meta)
    allowed = {'no_supported_exclusions','evidence_insufficient','user_stopped','awaiting_approval','technical_failure'}
    if reason not in allowed or not note.strip():
        raise ValueError('Explicit scientific stop reason and note required')
    if any(read(p).get('previous_round') == round_id for p in (Path(run)/'rounds').glob('*/frozen.json')):
        raise ValueError('Stop must refer to the terminal review round')
    actions = read(directory/'proposals.json')
    if reason == 'evidence_insufficient' and meta['context'].get('cleanup_mode') == 'auto_supported' and any(p['action'] in EXCLUDE for p in actions):
        if checked_auto_review(directory, frozen)['approved_proposals']:
            raise ValueError('Apply all supported proposals and recompute before scientific stopping')
    if reason == 'no_supported_exclusions' and any(p['action'] in EXCLUDE or p['action']=='review' for p in actions):
        raise ValueError('Unresolved/proposed exclusions cannot be called converged')
    if (directory/'applied.json').exists() and reason in {'no_supported_exclusions','evidence_insufficient'}:
        raise ValueError('Must recluster and review retained cells before scientific stopping')
    record = {'round':round_id,'reason':reason,'note':note,'at':stamp(),
              'converged':reason=='no_supported_exclusions','technical_status':'separately_validated'}
    # Approval wait is a pause snapshot, not a terminal run lock.
    write(Path(run)/('pause.json' if reason=='awaiting_approval' else 'stop.json'),record)
    return record


def summary(run):
    run, meta = run_contract(run)
    rounds=[]
    for p in sorted((run/'rounds').glob('*/frozen.json')):
        directory,_,frozen=round_contract(run,p.parent.name)
        rounds.append({'round':frozen['round'],'namespace':frozen['namespace'],
            'entry_mode':frozen['entry_mode'],'proposal_hash':frozen['proposal_hash'],
            'ledger':str(directory/('ledger.csv' if (directory/'applied.json').exists() else 'review_ledger.csv')),
            'approved':(directory/'approval.json').exists(),
            'automatic_review':checked_auto_review(directory,frozen) if (directory/'auto_review.json').exists() else None,
            'applied':checked_apply(directory,frozen) if (directory/'applied.json').exists() else None})
    return {'schema':meta['schema'],'cleanup_mode':meta['context'].get('cleanup_mode','manual_per_round'),
        'cumulative_removed':sum(r['applied']['removed'] for r in rounds if r['applied']), 'initial_cells':meta['identity']['shape'][0],'rounds':rounds,
        'stop':read(run/'stop.json') if (run/'stop.json').exists() else None,
        'scientific_status':'not_converged' if not (run/'stop.json').exists() else read(run/'stop.json')['reason']}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='command',required=True)
    for command in ('init','freeze','approve','auto-review','apply','validate','stop','summary'):
        q=sub.add_parser(command); q.add_argument('--run',type=Path,required=True)
        if command not in {'init','summary'}: q.add_argument('--round',required=True)
        if command in {'init','freeze'}: q.add_argument('--source',type=Path,required=True)
        if command=='init': q.add_argument('--context',type=Path,required=True)
        if command=='freeze':
            q.add_argument('--cluster-key',required=True); q.add_argument('--proposals',type=Path,required=True)
            q.add_argument('--execution-manifest', type=Path); q.add_argument('--previous'); q.add_argument('--recomputed',action='store_true')
        if command=='approve':
            q.add_argument('--expected-hash',required=True); q.add_argument('--select',nargs='+',required=True)
            q.add_argument('--approver',required=True); q.add_argument('--record',type=Path,required=True)
            q.add_argument('--simulated',action='store_true')
        if command=='auto-review':
            q.add_argument('--expected-hash',required=True); q.add_argument('--record',type=Path,required=True)
        if command=='apply': q.add_argument('--allow-simulated',action='store_true')
        if command=='validate': q.add_argument('--artifacts',type=Path)
        if command=='stop': q.add_argument('--reason',required=True); q.add_argument('--note',required=True)
    a=p.parse_args(argv)
    if a.command=='init': result=init(a.run,a.source,read(a.context))
    elif a.command=='freeze': result=freeze(a.run,a.round,a.source,a.cluster_key,read(a.proposals),a.previous,a.recomputed,a.execution_manifest)
    elif a.command=='approve': result=approve(a.run,a.round,a.expected_hash,a.select,a.approver,a.record,a.simulated)
    elif a.command=='auto-review': result=auto_review(a.run,a.round,a.expected_hash,a.record)
    elif a.command=='apply': result=apply(a.run,a.round,a.allow_simulated)
    elif a.command=='validate': result=validate(a.run,a.round,read(a.artifacts) if a.artifacts else None)
    elif a.command=='stop': result=stop(a.run,a.round,a.reason,a.note)
    else: result=summary(a.run)
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
