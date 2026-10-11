"""Audited explicit gene aliases; stable source keys never changed in-place."""
import numpy as np
import pandas as pd


def map_genes(table, stable_key, symbol_key, aliases=None, *, species):
    if not species or table[stable_key].isna().any() or table[stable_key].duplicated().any():
        raise ValueError('Species and unique stable gene IDs required')
    mapping = {}
    if aliases is not None:
        if not {'alias','canonical','species'}.issubset(aliases):
            raise ValueError('Alias table requires alias, canonical, species')
        subset = aliases[aliases.species.eq(species)]
        if subset[['alias','canonical']].isna().any().any() or (subset.groupby('alias').canonical.nunique()>1).any():
            raise ValueError('Ambiguous alias maps to multiple canonical genes')
        mapping = dict(zip(subset.alias.astype(str),subset.canonical.astype(str)))
    rows = []
    for _, row in table.iterrows():
        symbol = row[symbol_key]
        canonical = mapping.get(str(symbol), str(symbol)) if pd.notna(symbol) else None
        rows.append({'stable_gene_id':str(row[stable_key]),'source_symbol':symbol,
                     'canonical_symbol':canonical,'mapping_status':'missing_symbol' if canonical is None else ('explicit_alias' if canonical != symbol else 'unchanged')})
    audit = pd.DataFrame(rows)
    known = audit.canonical_symbol.notna()
    collisions = audit.loc[known,'canonical_symbol'].duplicated(keep=False)
    if collisions.any():
        raise ValueError('Many-to-one gene mapping collision; explicit biological aggregation policy required: ' +
                         ','.join(audit.loc[known].loc[collisions,'canonical_symbol'].unique()))
    return audit


def vote_review(votes, evidence, *, cell_key, min_fraction, min_margin):
    if not (0 <= min_fraction <= 1 and 0 <= min_margin <= 1):
        raise ValueError('Vote thresholds must be in [0,1]')
    if evidence[cell_key].isna().any() or evidence[cell_key].duplicated().any():
        raise ValueError('Unique nonmissing evidence cell keys required')
    if votes[[cell_key,'reference_label','weight']].isna().any().any():
        raise ValueError('Missing vote cell, label or weight')
    if not np.isfinite(votes.weight.to_numpy(float)).all() or (votes.weight < 0).any():
        raise ValueError('Votes require finite nonnegative weights')
    if not votes[cell_key].isin(evidence[cell_key]).all():
        raise ValueError('Vote cells missing from evidence table')
    rows = []
    for _, ev in evidence.iterrows():
        totals = votes[votes[cell_key].eq(ev[cell_key])].groupby('reference_label').weight.sum().sort_values(ascending=False,kind='stable')
        if totals.empty or totals.sum() <= 0:
            prediction, fraction, margin = None, np.nan, np.nan
            review, status = 'unknown', 'no_reference_support'
        else:
            fractions = totals/totals.sum()
            prediction, fraction = str(fractions.index[0]), float(fractions.iloc[0])
            margin = fraction - (float(fractions.iloc[1]) if len(fractions)>1 else 0)
            # Program audit is independently supplied; missing evidence cannot promote.
            accepted = ev.get('program_audit_status') == 'supported'
            conflict = ev.get('program_audit_status') == 'conflict'
            if accepted and fraction >= min_fraction and margin >= min_margin:
                review, status = prediction, 'supported_prediction_review'
            elif pd.notna(ev.get('supported_broad_label')) and not conflict:
                review, status = str(ev['supported_broad_label']), 'broader_supported_review'
            else:
                review, status = 'unknown', 'conflict_review' if conflict else 'insufficient_evidence'
        rows.append({cell_key:ev[cell_key],'predicted_label':prediction,'knn_vote_fraction':fraction,
                     'knn_vote_margin':margin,'calibrated_confidence':np.nan,
                     'calibration_status':'not_calibrated','review_label':review,'review_status':status})
    return pd.DataFrame(rows)
