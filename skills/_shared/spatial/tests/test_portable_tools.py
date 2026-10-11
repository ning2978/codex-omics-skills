"""Behavior tests in isolated skill trees with no project package or registry."""
from pathlib import Path
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'_shared/spatial'))
from run_contract import cell_table, matrix, neighborhoods, Run
from gene_mapping import map_genes, vote_review


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PortableSpatialTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='portable_spatial_')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.tree = self.work/'skills'
        self.project = self.work/'unrelated_project'
        self.project.mkdir()
        self.outputs = self.project/'outputs'
        self.outputs.mkdir()
        for rel in ['_shared/spatial','spatial/high-resolution/xenium/xenium-segmentation-new',
                    'spatial/high-resolution/xenium/xenium-resegmentation',
                    'spatial/high-resolution/xenium/xenium-scrna-label-transfer',
                    'spatial/high-resolution/xenium/xenium-scimap-neighborhood-analysis',
                    'spatial/shared/spatial-null-model','spatial/shared/algorithm-commot']:
            shutil.copytree(ROOT/rel,self.tree/rel,ignore=shutil.ignore_patterns('__pycache__','tests'))
        if os.environ.get('SPATIAL_TEST_ARTIFACT_DIR'):
            self.addCleanup(self.save_artifacts)
        self.cfg = {'cell_key':'obs_name','roi_key':'roi','biological_unit_key':'donor',
                    'x':'x','y':'y','coordinate_unit':'um','phenotype_key':'phenotype',
                    'label_evidence_grade':'synthetic','upstream':[]}

    def save_artifacts(self):
        target = Path(os.environ['SPATIAL_TEST_ARTIFACT_DIR'])/self._testMethodName
        if any(self.outputs.iterdir()):
            shutil.copytree(self.outputs,target)

    def script(self, name):
        mapping = {
            'seg':'spatial/high-resolution/xenium/xenium-segmentation-new/scripts/segmentation.py',
            'neighbor':'spatial/high-resolution/xenium/xenium-scimap-neighborhood-analysis/scripts/neighborhood.py',
            'null':'spatial/shared/spatial-null-model/scripts/spatial_null.py',
            'comm':'spatial/shared/algorithm-commot/scripts/prepare_inputs.py',
            'reseg':'spatial/high-resolution/xenium/xenium-resegmentation/scripts/check_inputs.py',
            'transfer':'spatial/high-resolution/xenium/xenium-scrna-label-transfer/scripts/audit_transfer.py'}
        return self.tree/mapping[name]

    def command(self, command, expected=0):
        env = {k:v for k,v in os.environ.items() if k!='PYTHONPATH'}
        env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        result = subprocess.run(command,cwd=self.project,env=env,text=True,capture_output=True,timeout=60)
        if os.environ.get('SPATIAL_TEST_COMMAND_LOG'):
            with open(os.environ['SPATIAL_TEST_COMMAND_LOG'],'a') as handle:
                handle.write(json.dumps({'test':self._testMethodName,'command':list(map(str,command)),
                    'cwd':str(self.project),'exit_code':result.returncode,'expected_exit':expected,
                    'status':'PASS' if (result.returncode==0)==(expected==0) else 'FAIL',
                    'stdout':result.stdout,'stderr':result.stderr})+'\n')
        if expected==0:
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        else:
            self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
        return result

    def cli(self,name,cfg,out,expected=0):
        path = self.project/(out+'-config.json')
        path.write_text(json.dumps(cfg))
        return self.command([sys.executable,str(self.script(name)),'--config',str(path),
                             '--out-dir',str(self.outputs/out)],expected)

    def cells(self):
        rows=[]
        for roi,donor,other in [('r1','u1','B'),('r2','u1','C'),('r3','u2','B')]:
            for j,(x,label,anchor,candidate) in enumerate([
                (0,'anchor',1,0),(2,'T',0,1),(12,other,0,1),(15,'T',0,1),(18,other,0,1),
                (4,None,0,1),(1000,'anchor',1,0)]):
                rows.append({'obs_name':f'{roi}:{j}','roi':roi,'donor':donor,'x':x,'y':0,
                             'phenotype':label,'anchor':anchor,'candidate':candidate,'stratum':'tissue'})
        frame=pd.DataFrame(rows)
        path=self.project/'cells.csv';frame.to_csv(path,index=False)
        return frame,path

    def null_config(self,path):
        return {**self.cfg,'cells_csv':str(path),'target_label_key':'phenotype','anchor_mask_key':'anchor',
            'candidate_mask_key':'candidate','block_size':200,'block_origin':[0,0],
            'strata_keys':['stratum'],'radius':6,'permutations':79,'seed':42,
            'families':{'primary':['near','count','enrich','nearest','abs','absent']},
            'hypotheses':[{'id':name,'family':'primary','metric':metric,'target':target,'alternative':alt,
                           **({'null_center':0} if alt=='two-sided' else {})} for name,metric,target,alt in [
                ('near','proximity','T','greater'),('count','cooccurrence','T','greater'),
                ('enrich','local_enrichment','T','greater'),('nearest','nearest','T','less'),
                ('abs','local_enrichment','T','two-sided'),('absent','proximity','Z','greater')]],
            'selection_provenance':{'status':'synthetic_predeclared'}}

    def test_no_project_package_or_registry(self):
        self.assertFalse((self.project/'metadata/current_annotation_registry.json').exists())
        self.command([sys.executable,'-c',"import importlib.util; assert importlib.util.find_spec('xenium_downstream') is None"])

    def test_ring_multilabel_unsorted_and_row_groups(self):
        path=self.project/'rings.parquet'
        rows=[]
        for vertex in range(4):
            for cid,label in [('a',0),('b',0),('a',1)]:
                rows.append({'cell_id':cid,'label_id':label,'vertex_x':vertex,'vertex_y':label+vertex})
        frame=pd.DataFrame(rows);pq.write_table(pa.Table.from_pandas(frame),path,row_group_size=2)
        rings=load('rings_test',self.script('seg').with_name('rings.py'))
        loaded=rings.load_rings(path)
        self.assertEqual(set(loaded['a']),{0,1})
        np.testing.assert_array_equal(loaded['a'][1][:,0],[0,1,2,3])
        self.assertEqual(rings.ring_vertex_counts(path).n_vertices.tolist(),[4,4,4])
        with self.assertRaisesRegex(ValueError,'Non-contiguous'):
            list(rings.iter_rings(path))
        ordered=frame.sort_values(['cell_id','label_id'],kind='stable')
        pq.write_table(pa.Table.from_pandas(ordered),path,row_group_size=2)
        self.assertEqual([len(v) for _,_,v in rings.iter_rings(path)],[4,4,4])

    def segmentation_fixture(self):
        path=self.project/'vendor';path.mkdir()
        pd.DataFrame({'cell_id':['a','b'],'cell_area':[100,80],'nucleus_count':[2,1],
                      'segmentation_method':['boundary_stain','interior_stain']}).to_parquet(path/'cells.parquet')
        rows=[]
        for cid,label,offset in [('a',0,0),('a',1,5),('b',0,20)]:
            for x,y in [(0,0),(3,0),(3,3),(0,3)]:
                rows.append({'cell_id':cid,'label_id':label,'vertex_x':x+offset,'vertex_y':y+offset})
        for name in ['cell','nucleus']:
            pq.write_table(pa.Table.from_pandas(pd.DataFrame(rows)),path/f'{name}_boundaries.parquet',row_group_size=3)
        return path

    def test_segmentation_cli_and_y_direction(self):
        vendor=self.segmentation_fixture()
        self.command([sys.executable,str(self.script('seg')),'--input-dir',str(vendor),'--out-dir',str(self.outputs/'seg'),
                      '--bbox','-1','25','-1','25'])
        manifest=json.loads((self.outputs/'seg/segmentation.json').read_text())
        self.assertEqual((manifest['cell_rings'],manifest['nucleus_rings']),(3,3))
        sys.path.insert(0,str(self.script('seg').parent))
        try:
            module=load('seg_test',self.script('seg'))
            recorded=[]
            original=module.plt.close
            def capture(fig):
                if hasattr(fig,'axes'):
                    recorded.append(fig.axes[0].get_ylim())
                original(fig)
            with patch.object(module.plt,'close',side_effect=capture):
                module.render_current(vendor,self.outputs/'y_down.png')
            self.assertTrue(recorded[-1][0]>recorded[-1][1])
            cells=module.read_cells(vendor/'cells.parquet').set_index('cell_id')
            collection,_,_,keys=module.ring_collection(module.load_rings(vendor/'cell_boundaries.parquet'),cells)
            self.assertEqual(len(collection.get_paths()),3)
            self.assertIn(('a',1),keys)
            with self.assertRaisesRegex(ValueError,'PNG only'):
                module.render_current(vendor,self.outputs/'forbidden.pdf')
        finally:
            sys.path.remove(str(self.script('seg').parent))

    def test_invalid_duplicate_cell_keys_and_roi_nesting(self):
        frame,path=self.cells()
        pd.concat([frame,frame.iloc[:1]]).to_csv(path,index=False)
        with self.assertRaisesRegex(ValueError,'unique'):
            cell_table(path,self.cfg,labels=True)
        frame.loc[1,'donor']='other';frame.to_csv(path,index=False)
        with self.assertRaisesRegex(ValueError,'nested'):
            cell_table(path,self.cfg,labels=True)

    def test_roi_isolation_center_pool_and_zero_neighbors(self):
        frame,path=self.cells()
        module=load('neighbor_test',self.script('neighbor'))
        centers=np.flatnonzero(frame.anchor)
        pool=np.flatnonzero(frame.candidate)
        counts,features,graph=module.count_features(frame,self.cfg,centers,pool,radius=6)
        self.assertEqual(counts[0,features.index('T')],1)
        self.assertEqual(counts[0,features.index('C')],0)
        self.assertEqual(counts[1].sum(),0)
        self.assertTrue(all(frame.iloc[v].roi.eq(frame.iloc[c].roi).all() for c,v in graph.items()))
        self.assertFalse(any(frame.iloc[v].anchor.any() for v in graph.values()))
        g=neighborhoods(frame,self.cfg,centers,pool,k=2)
        self.assertTrue(all(len(v)==2 for v in g.values()))

    def test_nmf_rcn_rank_seed_scale_handoffs(self):
        frame,path=self.cells()
        cfg={**self.cfg,'cells_csv':str(path),'mode':'nmf','cluster_method':'kmeans',
            'representation':'proportions','radii':[6,20],'ranks':[1,2],'seeds':[0,1],
            'n_states':2,'min_centers':2,'max_iter':2000,'center_mask_key':'anchor',
            'neighbor_mask_key':'candidate','namespace':'toy_nmf_v1','assignment_id':'toy_nmf_v1'}
        self.cli('neighbor',cfg,'nmf')
        loading=pd.read_csv(self.outputs/'nmf/nmf_loadings.csv')
        self.assertEqual(set(loading.neighbor_phenotype),{'T','B','C'})
        result=pd.read_csv(self.outputs/'nmf/assignments.csv')
        zero=result[result.obs_name.str.endswith(':6')]
        self.assertTrue(zero.assignment.isna().all())
        self.assertTrue(zero.status.eq('zero_labeled_neighbors').all())
        self.assertTrue(result.assignment.dropna().str.contains('toy_nmf_v1/').all())
        self.assertEqual(len(pd.read_csv(self.outputs/'nmf/fit_summary.csv')),8)
        self.assertGreater(len(pd.read_csv(self.outputs/'nmf/stability.csv')),0)
        cfg.update(mode='rcn',namespace='toy_rcn_v1',assignment_id='toy_rcn_v1',radii=[20])
        self.cli('neighbor',cfg,'rcn')
        assignments=pd.read_csv(self.outputs/'rcn/assignments.csv')
        self.assertTrue(assignments.assignment.dropna().str.contains('::RCN').all())

    def test_matrix_rejects_wrong_source_and_scaled_or_fractional_counts(self):
        obj=ad.AnnData(X=np.array([[1.,-1.],[2.,3.]]))
        with self.assertRaisesRegex(ValueError,'nonnegative'):
            matrix(obj,'X','normalized')
        obj.X=np.array([[.5,1],[2,3.]])
        with self.assertRaisesRegex(ValueError,'integer'):
            matrix(obj,'X','counts')
        with self.assertRaisesRegex(ValueError,'existing matrix'):
            matrix(obj,'layers/counts','counts')
        obj.var_names=['same','same']
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            matrix(obj,'X','log1p')

    def test_alias_mapping_conflicts_and_stable_ids(self):
        genes=pd.DataFrame({'id':['ens1','ens2'],'symbol':['OLD','NEW']})
        aliases=pd.DataFrame({'alias':['OLD'],'canonical':['NEW'],'species':['human']})
        with self.assertRaisesRegex(ValueError,'collision'):
            map_genes(genes,'id','symbol',aliases,species='human')
        aliases=pd.DataFrame({'alias':['OLD','OLD'],'canonical':['A','B'],'species':['human','human']})
        with self.assertRaisesRegex(ValueError,'Ambiguous'):
            map_genes(genes,'id','symbol',aliases,species='human')
        aliases=pd.DataFrame({'alias':['OLD'],'canonical':['A'],'species':['human']})
        audit=map_genes(genes,'id','symbol',aliases,species='human')
        self.assertEqual(audit.stable_gene_id.tolist(),['ens1','ens2'])
        self.assertEqual(genes.symbol.tolist(),['OLD','NEW'])

    def test_transfer_unknown_broader_and_uncalibrated_votes(self):
        votes=pd.DataFrame({'obs_name':['a','a','b','c'],'reference_label':['fine','other','fine','fine'],'weight':[.9,.1,1,1]})
        evidence=pd.DataFrame({'obs_name':['a','b','c','d'],'program_audit_status':['conflict','unavailable','supported','unavailable'],
                               'supported_broad_label':[None,'broad',None,None],'curated_identity':['native']*4})
        result=vote_review(votes,evidence,cell_key='obs_name',min_fraction=.7,min_margin=.2)
        self.assertEqual(result.review_label.tolist(),['unknown','broad','fine','unknown'])
        self.assertAlmostEqual(result.knn_vote_fraction.iloc[0],.9)
        self.assertTrue(result.calibrated_confidence.isna().all())
        self.assertEqual(evidence.curated_identity.tolist(),['native']*4)
        genes=pd.DataFrame({'gene_id':['ens1','ens2'],'symbol':['A','B']})
        genes.to_csv(self.project/'ref.csv',index=False);genes.to_csv(self.project/'panel.csv',index=False)
        votes.to_csv(self.project/'votes.csv',index=False);evidence.to_csv(self.project/'evidence.csv',index=False)
        cfg={'species':'human','stable_gene_key':'gene_id','symbol_key':'symbol',
             'reference_genes_csv':str(self.project/'ref.csv'),'target_genes_csv':str(self.project/'panel.csv'),
             'votes_csv':str(self.project/'votes.csv'),'evidence_csv':str(self.project/'evidence.csv'),
             'cell_key':'obs_name','min_vote_fraction':.7,'min_vote_margin':.2,'upstream':[]}
        self.cli('transfer',cfg,'transfer')
        reopened=pd.read_csv(self.outputs/'transfer/transfer_review.csv')
        self.assertIn('source::curated_identity',reopened)

    def test_null_margins_missing_and_reproducible_seed(self):
        frame,path=self.cells();cfg=self.null_config(path)
        module=load('null_test',self.script('null'))
        groups,_=module.permutation_units(frame,cfg,np.flatnonzero(frame.candidate))
        original=frame.phenotype.to_numpy(object)
        shuffled=module.permute_labels(original,groups,np.random.default_rng(42))
        np.testing.assert_array_equal(pd.isna(shuffled),pd.isna(original))
        for idx in groups.values():
            self.assertEqual(pd.Series(shuffled[idx]).value_counts().to_dict(),pd.Series(original[idx]).value_counts().to_dict())
        self.cli('null',cfg,'null1');self.cli('null',cfg,'null2')
        with np.load(self.outputs/'null1/null_distributions.npz') as a,np.load(self.outputs/'null2/null_distributions.npz') as b:
            np.testing.assert_array_equal(a['null'],b['null'])
        summary=pd.read_csv(self.outputs/'null1/null_summary.csv')
        testable=summary.status.eq('PASS')
        self.assertTrue(testable.any())
        np.testing.assert_allclose(summary.loc[testable,'p_empirical'],(summary.loc[testable,'extreme_permutations']+1)/80)
        self.assertEqual(summary.loc[summary.hypothesis_id.eq('absent'),'status'].iloc[0],'NOT_TESTABLE')
        self.assertTrue(summary.loc[~testable,'q_bh'].isna().all())
        support=pd.read_csv(self.outputs/'null1/biological_unit_support.csv')
        self.assertEqual(set(support.biological_unit),{'u1','u2'})
        self.assertEqual(support[support.biological_unit.eq('u1')].n_rois.unique().tolist(),[2])
        np.testing.assert_allclose(module.bh([.01,.04,np.nan]),[.03,.06,np.nan],equal_nan=True)

    def test_null_untestable_strata_and_families(self):
        frame,path=self.cells();cfg=self.null_config(path)
        cfg.update(strata_keys=[],block_size=1)
        self.cli('null',cfg,'degenerate')
        summary=pd.read_csv(self.outputs/'degenerate/null_summary.csv')
        self.assertTrue(summary.status.eq('NOT_TESTABLE').all())
        self.assertTrue(summary.p_empirical.isna().all())
        cfg['families']['primary'].append('not_declared')
        self.cli('null',cfg,'bad_family',expected=1)

    def test_null_missing_strata_and_missing_targets_remain_unavailable(self):
        frame,path=self.cells();cfg=self.null_config(path)
        frame.loc[2,'stratum']=None
        frame.to_csv(path,index=False)
        self.cli('null',cfg,'missing_strata')
        excluded=pd.read_csv(self.outputs/'missing_strata/excluded_missing.csv')
        self.assertIn('missing_stratum_not_testable',set(excluded.reason))
        self.assertEqual(excluded.reason.eq('missing_target_label_preserved').sum(),3)
        self.assertTrue(pd.read_csv(path).phenotype.isna().sum()==3)

    def test_complete_resume_and_corrupt_partial_changed_parameter_rejection(self):
        _,path=self.cells();cfg=self.null_config(path)
        self.cli('null',cfg,'resume')
        manifest=json.loads((self.outputs/'resume/run_manifest.json').read_text())
        self.assertIn('tools',manifest['signature']);self.assertIn('environment',manifest['signature'])
        result=self.cli('null',cfg,'resume')
        self.assertIn('complete-artifact resume',result.stdout)
        # A source change invalidates resume even when the path and params match.
        original_input=path.read_bytes()
        frame=pd.read_csv(path);frame.loc[1,'x']=3;frame.to_csv(path,index=False)
        self.cli('null',cfg,'resume',expected=1)
        path.write_bytes(original_input)
        # Tool and environment provenance are part of the complete signature.
        helper=self.tree/'_shared/spatial/run_contract.py';original_tool=helper.read_bytes()
        helper.write_bytes(original_tool+b'\n# tool-version change\n')
        self.cli('null',cfg,'resume',expected=1)
        helper.write_bytes(original_tool)
        manifest_path=self.outputs/'resume/run_manifest.json'
        original_manifest=manifest_path.read_bytes()
        tampered=json.loads(original_manifest)
        tampered['signature']['environment']['python_version']='unrecorded-version'
        manifest_path.write_text(json.dumps(tampered))
        self.cli('null',cfg,'resume',expected=1)
        manifest_path.write_bytes(original_manifest)
        extra=self.outputs/'resume/extra.txt';extra.write_text('not in inventory')
        self.cli('null',cfg,'resume',expected=1)
        extra.unlink()
        cfg['seed']=43
        self.cli('null',cfg,'resume',expected=1)
        cfg['seed']=42
        (self.outputs/'resume/null_summary.csv').write_text('corrupt\n')
        self.cli('null',cfg,'resume',expected=1)
        partial=self.outputs/'partial';partial.mkdir();(partial/'unfinished.txt').write_text('partial')
        self.cli('null',cfg,'partial',expected=1)

    def test_commot_preparation_matrix_panel_complex_and_roi_isolation(self):
        genes=['L','R','A','B','U']
        counts=sparse.csr_matrix([[2,1,1,1,0],[1,2,1,0,0],[3,2,1,1,0],[1,1,2,2,0]])
        obj=ad.AnnData(X=np.full((4,5),-1.),obs=pd.DataFrame({'roi':['r1','r1','r2','r2'],'donor':['u1','u1','u2','u2']},index=['a','b','c','d']),
            var=pd.DataFrame({'symbol':genes},index=[f'ens{i}' for i in range(5)]))
        obj.layers['counts']=counts;obj.obsm['spatial']=np.array([[0,0],[1,0],[0,0],[1,0]],float)
        path=self.project/'input.h5ad';obj.write_h5ad(path)
        lr=self.project/'lr.csv'
        pd.DataFrame({'ligand':['L','A_B','L','L'],'receptor':['R','R','R_MISSING','U'],'pathway':['p1','p2','p3','p4']}).to_csv(lr,index=False)
        cfg={'input_h5ad':str(path),'matrix_source':'layers/counts','matrix_kind':'counts',
             'output_transform':'normalize_log1p','target_sum':100,'counts_source':'layers/counts',
             'feature_key':'symbol','species':'human','database_species':'human','database_name':'synthetic',
             'database_version':'1','lr_database_csv':str(lr),'complex_delimiter':'_',
             'min_cell_fraction':.1,'roi_key':'roi','biological_unit_key':'donor','coordinate_unit':'um',
             'distance_threshold':2,'upstream':[]}
        self.cli('comm',cfg,'comm')
        panel=pd.read_csv(self.outputs/'comm/panel_audit_all_cells.csv')
        self.assertFalse(panel.loc[2,'all_subunits_measured'])
        self.assertTrue(panel.loc[1,'all_subunits_measured'])
        self.assertFalse(panel.loc[3,'expression_eligible'])
        self.assertEqual(json.loads(panel.loc[3,'below_detection_subunits']),['U'])
        mapping=pd.read_csv(self.outputs/'comm/gene_mapping_audit.csv')
        self.assertEqual(mapping.stable_gene_id.tolist(),obj.var_names.tolist())
        for file in (self.outputs/'comm/roi_inputs').glob('*.h5ad'):
            result=ad.read_h5ad(file)
            self.assertEqual(result.obs.roi.nunique(),1)
            self.assertTrue((result.X.data>=0).all())
            self.assertIn('counts',result.layers)
            self.assertEqual(result.var.source_gene_id.tolist(),obj.var_names.tolist())
            self.assertEqual(result.uns['communication_input_contract']['status'],'PREPARED_MODEL_NOT_RUN')
        original=ad.read_h5ad(path)
        self.assertTrue((original.X==-1).all())
        cfg.update(matrix_source='X',matrix_kind='normalized',output_transform='preserve')
        self.cli('comm',cfg,'bad_matrix',expected=1)

    def test_resegmentation_input_and_unverified_handoff(self):
        import tifffile
        image=self.project/'pilot.ome.tif'
        tifffile.imwrite(image,np.zeros((2,10,10),dtype=np.uint8),metadata={'axes':'CYX'})
        tx=self.project/'transcripts.csv'
        pd.DataFrame({'tid':['t1','t2'],'x':[1,2],'y':[1,3],'gene':['A','B']}).to_csv(tx,index=False)
        cfg={'method_class':'image_segmentation','method':'example_model','image':str(image),'channels':{'DAPI':0,'membrane':1},
            'segmentation_channels':['DAPI','membrane'],'pixel_size':.5,'image_coordinate_frame':'um_frame',
            'coordinate_frame':'um_frame','orientation':'y_down','coordinate_unit':'um','transcripts':str(tx),
            'transcript_key':'tid','transcript_x':'x','transcript_y':'y','gene_key':'gene','roi_id':'r1','biological_unit':'u1',
            'pilot_bbox':[0,5,0,5],'previous_cell_set_version':'old','new_cell_set_version':'new',
            'previous_assignment_version':'old','new_assignment_version':'new','upstream':[]}
        self.cli('reseg',cfg,'reseg')
        result=json.loads((self.outputs/'reseg/input_validation.json').read_text())
        self.assertEqual(result['model_execution'],'NOT_RUN')
        self.assertEqual(result['model_validation'],'UNVERIFIED')
        self.assertIn('annotation',result['required_downstream_rebuild'])
        cfg['new_cell_set_version']='old'
        self.cli('reseg',cfg,'bad_version',expected=1)
        cfg['new_cell_set_version']='new';cfg['image_coordinate_frame']='unregistered'
        self.cli('reseg',cfg,'bad_frame',expected=1)
        cfg['image_coordinate_frame']='um_frame';cfg.pop('image')
        cfg.update(method_class='transcript_assignment',method='example_assignment',prior_geometry_not_required_reason='synthetic prior-free pilot')
        self.cli('reseg',cfg,'transcript_assignment')


if __name__=='__main__':
    unittest.main(verbosity=2)
