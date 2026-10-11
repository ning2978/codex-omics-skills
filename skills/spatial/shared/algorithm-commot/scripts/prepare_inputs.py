"""Portable explicit-matrix COMMOT ROI inputs and measured-panel complex audit."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

shared = next(p/'_shared/spatial' for p in Path(__file__).resolve().parents if (p/'_shared/spatial/run_contract.py').is_file())
sys.path.insert(0, str(shared))
from run_contract import Run, matrix, write_json
from gene_mapping import map_genes

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def panel_audit(database, genes, expression, delimiter, min_cell_fraction):
    if not {'ligand','receptor','pathway'}.issubset(database) or database[['ligand','receptor','pathway']].isna().any().any():
        raise ValueError('LR database needs nonmissing ligand, receptor, pathway')
    if database.duplicated(['ligand','receptor','pathway']).any():
        raise ValueError('Duplicate database LR/pathway rows')
    if not delimiter or not 0 <= min_cell_fraction <= 1:
        raise ValueError('Explicit complex delimiter and detection threshold required')
    detected = np.asarray((expression > 0).sum(axis=0)).ravel()/expression.shape[0]
    lookup = dict(zip(genes, detected))
    records = []
    for i,row in database.iterrows():
        parts = {key:str(row[key]).split(delimiter) for key in ('ligand','receptor')}
        if any(not gene for group in parts.values() for gene in group):
            raise ValueError('Malformed LR complex identifier')
        missing = sorted({gene for group in parts.values() for gene in group if gene not in lookup})
        unavailable = sorted({gene for group in parts.values() for gene in group if gene in lookup and lookup[gene]<min_cell_fraction})
        records.append({'database_row':int(i),'ligand':row.ligand,'receptor':row.receptor,'pathway':row.pathway,
            'ligand_subunits':json.dumps(parts['ligand']),'receptor_subunits':json.dumps(parts['receptor']),
            'missing_panel_subunits':json.dumps(missing),'below_detection_subunits':json.dumps(unavailable),
            'all_subunits_measured':not missing,'expression_eligible':not missing and not unavailable,
            'subunit_detection_fractions':json.dumps({gene:lookup.get(gene) for group in parts.values() for gene in group})})
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    required = ['species','database_name','database_species','database_version','coordinate_unit','matrix_source',
                'matrix_kind','feature_key','complex_delimiter','distance_threshold','output_transform','min_cell_fraction']
    if any(k not in cfg for k in required) or cfg['species'] != cfg['database_species']:
        raise ValueError('Explicit matrix, species-matched database, units and parameters required')
    if not cfg['coordinate_unit'] or not np.isfinite(cfg['distance_threshold']) or cfg['distance_threshold'] <= 0:
        raise ValueError('Explicit distance unit and positive distance threshold required')
    inputs = [args.config,cfg['input_h5ad'],cfg['lr_database_csv']]
    if cfg.get('aliases_csv'):
        inputs.append(cfg['aliases_csv'])
    run = Run(args.out_dir,cfg,inputs,__file__)
    if run.resumed:
        print('PASS complete-artifact resume')
        return
    obj = ad.read_h5ad(cfg['input_h5ad'])
    if not obj.n_obs or not obj.n_vars:
        raise ValueError('Nonempty cell by measured-gene matrix required')
    expression = matrix(obj,cfg['matrix_source'],cfg['matrix_kind'])
    spatial = np.asarray(obj.obsm[cfg.get('spatial_key','spatial')],dtype=float)
    if spatial.shape != (obj.n_obs,2) or not np.isfinite(spatial).all():
        raise ValueError('Finite two-dimensional cell coordinates required')
    roi_key, unit_key = cfg['roi_key'],cfg['biological_unit_key']
    if obj.obs[[roi_key,unit_key]].isna().any().any() or (obj.obs.groupby(roi_key,observed=True)[unit_key].nunique()!=1).any():
        raise ValueError('Unique ROI to biological-unit nesting required')
    table = pd.DataFrame({'stable_id':obj.var_names.astype(str),
                          'symbol':obj.var_names.astype(str) if cfg['feature_key']=='var_names' else obj.var[cfg['feature_key']].to_numpy()})
    aliases = pd.read_csv(cfg['aliases_csv']) if cfg.get('aliases_csv') else None
    audit = map_genes(table,'stable_id','symbol',aliases,species=cfg['species'])
    if audit.canonical_symbol.isna().any():
        raise ValueError('Every COMMOT feature needs an audited canonical symbol')
    if cfg['matrix_kind'] == 'counts':
        if cfg['output_transform'] != 'normalize_log1p' or not np.isfinite(cfg.get('target_sum',0)) or cfg.get('target_sum',0)<=0:
            raise ValueError('Counts preparation requires explicit normalize_log1p and positive target_sum')
        sums = np.asarray(expression.sum(axis=1)).ravel()
        scale = np.divide(cfg['target_sum'],sums,out=np.zeros_like(sums,dtype=float),where=sums>0)
        if sparse.issparse(expression):
            expression = sparse.diags(scale) @ expression
            expression.data = np.log1p(expression.data)
        else:
            expression = np.log1p(expression * scale[:,None])
    elif cfg['output_transform'] != 'preserve':
        raise ValueError('Normalized/log1p source must be preserved; no implicit second normalization')
    retained_counts = matrix(obj,cfg['counts_source'],'counts') if cfg.get('counts_source') else None
    var = pd.DataFrame({'source_gene_id':obj.var_names.astype(str),'source_symbol':audit.source_symbol.to_numpy()},
                       index=pd.Index(audit.canonical_symbol.to_numpy(),name='canonical_symbol'))
    database = pd.read_csv(cfg['lr_database_csv']).reset_index(drop=True)
    audit.to_csv(args.out_dir/'gene_mapping_audit.csv',index=False)
    panel_audit(database,var.index,expression,cfg['complex_delimiter'],cfg['min_cell_fraction']).to_csv(args.out_dir/'panel_audit_all_cells.csv',index=False)
    roi_dir = args.out_dir/'roi_inputs'
    roi_dir.mkdir()
    summaries = []
    for roi in sorted(obj.obs[roi_key].astype(str).unique()):
        selection = obj.obs[roi_key].astype(str).eq(roi).to_numpy()
        working = ad.AnnData(X=expression[selection].copy(),obs=obj.obs.loc[selection].copy(),var=var.copy())
        working.obsm['spatial'] = spatial[selection].copy()
        if retained_counts is not None:
            working.layers['counts'] = retained_counts[selection].copy()
        lr = panel_audit(database,var.index,working.X,cfg['complex_delimiter'],cfg['min_cell_fraction'])
        token = hashlib.sha256(roi.encode()).hexdigest()[:16]
        working.uns['communication_input_contract'] = {'species':cfg['species'],'database_name':cfg['database_name'],
            'database_version':cfg['database_version'],'coordinate_unit':cfg['coordinate_unit'],
            'distance_threshold':cfg['distance_threshold'],'matrix_source':cfg['matrix_source'],
            'matrix_kind':cfg['matrix_kind'],'output_transform':cfg['output_transform'],
            'status':'PREPARED_MODEL_NOT_RUN','stable_gene_id_column':'source_gene_id'}
        working.write_h5ad(roi_dir/f'{token}.h5ad')
        database.loc[lr.expression_eligible].to_csv(roi_dir/f'{token}_lr.csv',index=False)
        lr.to_csv(roi_dir/f'{token}_panel_audit.csv',index=False)
        summaries.append({'roi':roi,'biological_unit':str(working.obs[unit_key].iloc[0]),'n_cells':working.n_obs,
                          'n_eligible_lr':int(lr.expression_eligible.sum()),'input_h5ad':f'roi_inputs/{token}.h5ad',
                          'lr_csv':f'roi_inputs/{token}_lr.csv','status':'PREPARED' if lr.expression_eligible.any() else 'NO_ELIGIBLE_LR'})
    pd.DataFrame(summaries).to_csv(args.out_dir/'roi_input_manifest.csv',index=False)
    run.finish({'status':'PASS','scope':'explicit_matrix_panel_and_subunit_audit','model_execution':'NOT_RUN'},
               evidence_grade='exploratory',interpretation='modeled_communication_potential_only')
    print('PASS COMMOT inputs prepared; model not run')


if __name__ == '__main__':
    main()
