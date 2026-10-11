"""Audit supplied gene maps and transfer vote/review evidence; no embedding fit."""
from pathlib import Path
import argparse
import json
import sys

shared = next(p/'_shared/spatial' for p in Path(__file__).resolve().parents if (p/'_shared/spatial/run_contract.py').is_file())
sys.path.insert(0,str(shared))
from run_contract import Run
from gene_mapping import map_genes, vote_review
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if bool(cfg.get('votes_csv')) != bool(cfg.get('evidence_csv')):
        raise ValueError('Votes and evidence must be supplied together')
    inputs = [args.config,cfg['reference_genes_csv'],cfg['target_genes_csv']]
    inputs += [cfg[k] for k in ['aliases_csv','votes_csv','evidence_csv'] if cfg.get(k)]
    run = Run(args.out_dir,cfg,inputs,__file__)
    if run.resumed:
        print('PASS complete-artifact resume')
        return
    aliases = pd.read_csv(cfg['aliases_csv']) if cfg.get('aliases_csv') else None
    audits = {}
    for label in ['reference','target']:
        table = pd.read_csv(cfg[f'{label}_genes_csv'])
        audits[label] = map_genes(table,cfg['stable_gene_key'],cfg['symbol_key'],aliases,species=cfg['species'])
        audits[label].to_csv(args.out_dir/f'{label}_gene_mapping.csv',index=False)
    reference = audits['reference'].dropna(subset=['canonical_symbol'])
    target = audits['target'].dropna(subset=['canonical_symbol'])
    reference.merge(target,on='canonical_symbol',suffixes=('_reference','_target'),validate='one_to_one').to_csv(args.out_dir/'shared_gene_mapping.csv',index=False)
    if cfg.get('votes_csv'):
        key = cfg['cell_key']
        evidence = pd.read_csv(cfg['evidence_csv'],dtype={key:str})
        votes = pd.read_csv(cfg['votes_csv'],dtype={key:str})
        reviews = vote_review(votes,evidence,cell_key=key,min_fraction=cfg['min_vote_fraction'],min_margin=cfg['min_vote_margin'])
        # Source annotations/evidence remain present with unambiguous source prefix.
        source = evidence.rename(columns={c:f'source::{c}' for c in evidence if c!=key})
        reviews.merge(source,on=key,validate='one_to_one').to_csv(args.out_dir/'transfer_review.csv',index=False)
    run.finish({'status':'PASS','scope':'mapping_and_supplied_vote_review_only','calibration':'NOT_CALIBRATED',
                'embedding_fit':'NOT_RUN','biological_evidence':'SUPPLIED_NOT_INDEPENDENTLY_VALIDATED'},evidence_grade='prediction_review')
    print('PASS gene/vote audit; curated annotations unchanged')


if __name__ == '__main__':
    main()
