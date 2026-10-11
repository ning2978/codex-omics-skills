import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

SHARED = Path(__file__).resolve().parents[1]
SKILLS = SHARED.parents[1]
sys.path.insert(0, str(SHARED))
from annotation_contract import resolve

def entry(modality):
    p = SKILLS/f'render-{modality}-dotplot/scripts/render_{modality}_dotplot.py'
    spec = importlib.util.spec_from_file_location('test_'+modality, p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class Contracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.guide = self.root/'guide.md'
        self.registry = self.root/'registry.json'
        self.source = {'source_id':'source', 'priority':1, 'path':'source.csv',
            'join':{'sidecar_key':'obs_name','target_key':'obs_name','sidecar_key_unique':True,'row_order_join_forbidden':True},
            'label_columns':{'primary_operational':'label'}, 'fallback':'Use label only when include is true.'}
        self.track = {'track_id':'test','status':'current','sources':[self.source]}
        self.persist(self.track)
    def persist(self, track, other=None):
        self.registry.write_text(json.dumps({'updated_at':'different','use_cases':{'xenium_top_marker':{'default_track':track},'epi_subtype':{'default_track':other or track}}}))
        self.guide.write_text('updated_at: old\n<!-- xenium_top_marker_semantics_begin -->\n```json\n'+json.dumps(self.track)+'\n```\n<!-- xenium_top_marker_semantics_end -->')
    def test_timestamp_only_allowed(self):
        track = copy.deepcopy(self.track); track['updated_at']='new'; self.persist(track)
        self.assertEqual(resolve(self.guide,self.registry,'xenium_top_marker')['track_id'],'test')
    def test_unrelated_role_allowed(self):
        self.persist(self.track, {'track_id':'changed','status':'paused'})
        resolve(self.guide,self.registry,'xenium_top_marker')
    def test_each_semantic_conflict_rejected(self):
        for k,v in [('track_id','wrong'),('status','review'),('sources',[])]:
            with self.subTest(k=k):
                track=copy.deepcopy(self.track); track[k]=v; self.persist(track)
                with self.assertRaisesRegex(ValueError,k): resolve(self.guide,self.registry,'xenium_top_marker')
        for k,v in [('path','other.csv'),('priority',2),('source_id','other'),('label_columns',{'primary':'other'}),('fallback','fill all'),('exclusion_gates',['include == false']),('join',{})]:
            with self.subTest(k=k):
                track=copy.deepcopy(self.track); track['sources'][0][k]=v; self.persist(track)
                with self.assertRaisesRegex(ValueError,k): resolve(self.guide,self.registry,'xenium_top_marker')
    def test_missing_summary_rejected(self):
        self.guide.write_text('# no summary')
        with self.assertRaisesRegex(ValueError,'exactly one'): resolve(self.guide,self.registry,'xenium_top_marker')
    def test_duplicate_summary_rejected(self):
        self.guide.write_text(self.guide.read_text()*2)
        with self.assertRaisesRegex(ValueError,'exactly one'): resolve(self.guide,self.registry,'xenium_top_marker')
    def test_join_shuffle_and_errors_both_modalities(self):
        for modality,key in [('xenium','obs_name'),('scrna','cell_id')]:
            m=entry(modality); frame=pd.DataFrame({key:['c','a','b'],'label':['B','A','A']}); obs=pd.Index(['a','b','c'],name='obs_name')
            aligned=m.align_sidecar_to_obs_names(obs,frame,['label'],key)
            pd.testing.assert_frame_equal(aligned,m.align_sidecar_to_obs_names(obs,frame.sample(frac=1,random_state=2),['label'],key))
            for invalid in [pd.concat([frame,frame.iloc[:1]]),frame.iloc[:2],frame.assign(**{key:['x','a','b']})]:
                with self.assertRaises(ValueError): m.align_sidecar_to_obs_names(obs,invalid,['label'],key)
    def setup_source(self):
        m=entry('xenium'); frame=pd.DataFrame({'obs_name':['b','a'],'label':['B','A'],'include':[True,True]}); frame.to_csv(self.root/'source.csv',index=False)
        aligned=m.align_sidecar_to_obs_names(pd.Index(['a','b'],name='obs_name'),frame,['label'],'obs_name')
        config={'_project_root':self.root,'join':{'sidecar_key':'obs_name'},'groupings':{'g':{'field':'label'}},'project':{'_annotation_bindings':{'g':{'source':self.source,'path':str(self.root/'source.csv'),'use_case':'xenium_top_marker','track_id':'test'}},'_resolved_annotation_roles':{'xenium_top_marker':self.track}}}
        return m,frame,aligned,config
    def test_source_shuffle_verified(self):
        m,frame,aligned,config=self.setup_source(); self.assertEqual(m.validate_annotation_sources(config,aligned)[0]['rows_verified'],2)
    def test_label_tampering_rejected(self):
        m,frame,aligned,config=self.setup_source(); aligned.loc['a','label']='B'
        with self.assertRaisesRegex(ValueError,'differ'): m.validate_annotation_sources(config,aligned)
    def test_false_or_missing_gate_rejected(self):
        for value in [False,None,'unknown']:
            m,frame,aligned,config=self.setup_source(); frame['include']=frame['include'].astype(object); frame.loc[0,'include']=value; frame.to_csv(self.root/'source.csv',index=False)
            with self.assertRaisesRegex(ValueError,'gate failed'): m.validate_annotation_sources(config,aligned)
    def test_outside_source_rejected(self):
        m,frame,aligned,config=self.setup_source(); frame.iloc[:1].to_csv(self.root/'source.csv',index=False)
        with self.assertRaisesRegex(ValueError,'outside'): m.validate_annotation_sources(config,aligned)
    def test_cohort_rejected(self):
        m,frame,aligned,config=self.setup_source(); self.source['cohort_scope']='example_control_only'; frame['sample_role']='P'; frame.to_csv(self.root/'source.csv',index=False)
        with self.assertRaisesRegex(ValueError,'cohort'): m.validate_annotation_sources(config,aligned)

    def test_portable_cohort_predicates_and_explicit_legacy_scope(self):
        m,frame,aligned,config=self.setup_source()
        frame['condition']='treated'; frame.to_csv(self.root/'source.csv',index=False)
        self.source['cohort_gate']={'allowed_values':{'condition':['treated']}}
        m.validate_annotation_sources(config,aligned)
        self.source['cohort_gate']['allowed_values']['condition']=['control']
        with self.assertRaisesRegex(ValueError,'cohort gate failed'):
            m.validate_annotation_sources(config,aligned)
        del self.source['cohort_gate']
        self.source['cohort_scope']='example treatment cohort'
        with self.assertRaisesRegex(ValueError,'explicit typed predicates'):
            m.validate_annotation_sources(config,aligned)
        config['project']['cohort_scope_gates']={'example treatment cohort':{'allowed_values':{'condition':['treated']}}}
        m.validate_annotation_sources(config,aligned)
    def test_priority_rejected(self):
        m,frame,aligned,config=self.setup_source(); higher=copy.deepcopy(self.source); higher.update(priority=10,source_id='higher'); self.track['sources'].append(higher)
        with self.assertRaisesRegex(ValueError,'higher-priority'): m.validate_annotation_sources(config,aligned)
    def test_independent_entrypoints(self):
        x=entry('xenium'); s=entry('scrna'); self.assertEqual(x.MODALITY,'xenium'); self.assertEqual(s.MODALITY,'scrna'); x.run_self_test(); s.run_self_test()

    def prepare_fixture(self, modality='xenium'):
        import yaml
        m=entry(modality)
        for f in ['counts.h5ad','source.csv','markers.yaml']:
            (self.root/f).touch()
        config={'analysis_id':'test','project_root':str(self.root),'inputs':{'counts_h5ad':'counts.h5ad','sidecar':'source.csv','marker_yaml':'markers.yaml'},'join':{**self.source['join'],'target_key_unique':True},'groupings':{'g':{'field':'label','order':['A','B']}},'marker_selection':{},'seurat_native':{'rscript':'Rscript','expected_r_version':'4.4.3','expected_seurat_version':'5.5.1','expected_seuratobject_version':'5.4.0'},'project':{'annotation_use_case':'xenium_top_marker','annotation_guide':'guide.md','annotation_registry':'registry.json'},'output':{'root':'out'}}
        return m,config
    def prepare(self,m,config):
        import yaml
        path=self.root/'run.yaml'; path.write_text(yaml.safe_dump(config)); return m.prepare_config(path)
    def test_prepare_rejects_xenium_cell_id(self):
        m,c=self.prepare_fixture(); c['join']['sidecar_key']='cell_id'
        with self.assertRaisesRegex(ValueError,'cell_id joins are forbidden'): self.prepare(m,c)
    def test_prepare_binds_exact_source_join(self):
        m,c=self.prepare_fixture('scrna'); c['join']['sidecar_key']='cell_id'
        other=copy.deepcopy(self.source); other['source_id']='other'; other['label_columns']={'primary':'other_label'}; other['join']['sidecar_key']='cell_id'; self.track['sources'].append(other); self.persist(self.track)
        with self.assertRaisesRegex(ValueError,'join differs'): self.prepare(m,c)
    def test_prepare_separate_roles(self):
        m,c=self.prepare_fixture(); c['groupings']['epi']={'field':'label','order':['A','B'],'annotation_use_case':'epi_subtype'}
        self.guide.write_text(self.guide.read_text()+'\n<!-- epi_subtype_semantics_begin -->\n```json\n'+json.dumps(self.track)+'\n```\n<!-- epi_subtype_semantics_end -->')
        prepared=self.prepare(m,c); self.assertEqual(set(prepared['project']['_resolved_annotation_roles']),{'xenium_top_marker','epi_subtype'})
    def test_prepare_rejects_retired_rcn_assignment(self):
        m,c=self.prepare_fixture(); c['groupings']['g']['assignment_id']='scvi5AllCellRCN1-10'
        with self.assertRaisesRegex(ValueError,'qualified assignment'): self.prepare(m,c)
    def test_prepare_rejects_paused_role(self):
        m,c=self.prepare_fixture(); self.track['status']='paused_pending'; self.persist(self.track)
        with self.assertRaisesRegex(ValueError,'unavailable'): self.prepare(m,c)

    def test_auto_height_and_explicit_style(self):
        for modality in ('xenium', 'scrna'):
            m,c=self.prepare_fixture(modality)
            c['groupings']['g']['order']=['A']
            c['output'].pop('height_inches', None)
            prepared=self.prepare(m,c)
            self.assertEqual(prepared['groupings']['g']['height_inches'],3.2)
            c['groupings']['g']['order']=[str(i) for i in range(20)]
            prepared=self.prepare(m,c)
            self.assertEqual(prepared['groupings']['g']['height_inches'],8.2)
            c['output']['height_inches']=5.5
            c['style']={'gene_font_pt':13,'strip_height_inches':.7}
            prepared=self.prepare(m,c)
            self.assertEqual(prepared['groupings']['g']['height_inches'],5.5)
            self.assertEqual(prepared['style']['gene_font_pt'],13)

    def test_invalid_style_and_dimensions(self):
        for key,value in [('height_inches',0),('height_inches',float('nan')),('width_inches',float('inf'))]:
            m,c=self.prepare_fixture(); c['output'][key]=value
            with self.assertRaises(ValueError): self.prepare(m,c)
        for style in [{'gene_font_pt':0},{'strip_height_inches':float('nan')},{'unknown':11}]:
            m,c=self.prepare_fixture(); c['style']=style
            with self.assertRaises(ValueError): self.prepare(m,c)

    def test_registered_hash_rejected(self):
        m,frame,aligned,config=self.setup_source(); self.source['sha256']='incorrect'
        with self.assertRaisesRegex(ValueError,'sha256 mismatch'): m.validate_annotation_sources(config,aligned)
    def test_legacy_projection_and_complete_contract_both_checked(self):
        full=json.dumps(self.track)
        self.guide.write_text('<!-- xenium_top_marker_semantics_begin -->\n```json\n'+json.dumps({'track_id':'wrong'})+'\n```\n<!-- xenium_top_marker_semantics_end -->\n<!-- xenium_top_marker_complete_semantics_begin -->\n```json\n'+full+'\n```\n<!-- xenium_top_marker_complete_semantics_end -->')
        with self.assertRaisesRegex(ValueError,'legacy_summary.track_id'): resolve(self.guide,self.registry,'xenium_top_marker')

if __name__=='__main__': unittest.main()
