"""Style integration: old/new native statistics, single groups, sizes and NA parity."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
import yaml

p=argparse.ArgumentParser(); p.add_argument('--baseline',type=Path,help='Optional previous renderer for old/new numeric regression'); p.add_argument('--rscript',required=True); p.add_argument('--work-dir',type=Path,required=True); args=p.parse_args()
root=args.work_dir; root.mkdir(parents=True,exist_ok=False)
skills=Path(__file__).resolve().parents[3]
counts=np.array([[5,0,2,1],[2,0,1,3],[0,8,2,1],[0,3,1,2],[3,1,0,2],[1,4,0,3]],dtype=np.int32)
a=ad.AnnData(X=sparse.csr_matrix(counts),obs=pd.DataFrame(index=[f'c{i}' for i in range(6)]),var=pd.DataFrame(index=['HLA-DRA','B','C','D'])); a.layers['counts']=a.X.copy(); a.write_h5ad(root/'counts.h5ad')
pd.DataFrame({'obs_name':a.obs_names,'label':['one','one','two','two','three','three'],'include':True}).sample(frac=1,random_state=8).to_csv(root/'source.csv',index=False)
source={'source_id':'synthetic','priority':1,'path':'source.csv','label_columns':{'primary':'label','cells':'obs_name'},'fallback':'Use label only when include is true.','join':{'sidecar_key':'obs_name','target_key':'obs_name','sidecar_key_unique':True,'row_order_join_forbidden':True}}
track={'track_id':'synthetic_test_only','status':'current','sources':[source]}
(root/'registry.json').write_text(json.dumps({'use_cases':{'scrna_top_marker':{'default_track':track}}}))
(root/'guide.md').write_text('<!-- scrna_top_marker_semantics_begin -->\n```json\n'+json.dumps(track)+'\n```\n<!-- scrna_top_marker_semantics_end -->')
(root/'markers.yaml').write_text(yaml.safe_dump({'programs':[{'name':'Program','category':'identity','genes':['HLA-DRA','B','C','D']}]}))
config={'analysis_id':'synthetic_parity','project_root':str(root),'inputs':{'counts_h5ad':'counts.h5ad','sidecar':'source.csv','marker_yaml':'markers.yaml'},'join':{**source['join'],'target_key_unique':True},'groupings':{'labels':{'field':'label','order':['one','two','three']}},'marker_selection':{},'seurat_native':{'rscript':args.rscript,'expected_r_version':'4.4.3','expected_seurat_version':'5.5.1','expected_seuratobject_version':'5.4.0'},'project':{'annotation_use_case':'scrna_top_marker','annotation_guide':'guide.md','annotation_registry':'registry.json','run_recluster_guard':False},'output':{'root':'unused','dpi':80,'width_inches':8,'height_inches':5}}

import copy

def run(name, cfg, script):
    cfg=copy.deepcopy(cfg); cfg['output']['root']=name
    if name == 'baseline':
        cfg['seurat_native']['script']=str(args.baseline.parent/'seurat_dotplot_renderer.R')
    path=root/f'{name}.yaml'; path.write_text(yaml.safe_dump(cfg))
    result=subprocess.run([sys.executable,str(script),'render','--config',str(path)],capture_output=True,text=True)
    (root/f'{name}.log').write_text(result.stdout+result.stderr)
    if result.returncode: raise RuntimeError(f'{name}: see {root/name}.log\n'+result.stderr[-3000:])
    result=subprocess.run([sys.executable,str(script),'validate','--output',str(root/name),'--verify-inputs'],capture_output=True,text=True)
    if result.returncode: raise RuntimeError(result.stdout+result.stderr)
    for png in (root/name).rglob('*.png'):
        with Image.open(png) as image:
            image.load(); assert image.getpixel((0,0))[:3] == (255,255,255)
    assert not list((root/name).rglob('*.pdf'))
    return root/name

base=run('baseline',config,args.baseline) if args.baseline else None
reports=[]
for modality in ['xenium','scrna']:
    script=skills/f'render-{modality}-dotplot/scripts/render_{modality}_dotplot.py'
    current=run(modality,config,script)
    rel='02_numeric_sources/native_seurat/labels_seurat_native_pdata.csv'
    if base is not None:
        old=pd.read_csv(base/rel); new=pd.read_csv(current/rel)
        pd.testing.assert_frame_equal(old,new,atol=1e-9,rtol=0)
    cfg=copy.deepcopy(config); cfg['output'].pop('height_inches'); cfg['output']['dpi']=240
    cfg['groupings']['single']={'field':'label','order':['All'], 'label_map':{'one':'All','two':'All','three':'All'},'allow_many_to_one_label_map':True, 'title':'Single-group identity'}
    cfg['groupings']['labels']['title']='Long marker labels and single-gene programs'
    cfg['groupings']['labels']['subtitle']='Synthetic validation'
    cfg['groupings']['labels']['label_map']={'one':'Long group label one','two':'Long group label two','three':'Long group label three'}
    cfg['groupings']['labels']['order']=list(cfg['groupings']['labels']['label_map'].values())
    # Deliberately varied heights in one native invocation.
    cfg['groupings']['cells']={'field':'obs_name','order':list(a.obs_names),'title':'Individual cells'}
    (root/'style_markers.yaml').write_text(yaml.safe_dump({'programs':[{'name':'Long program label','short_label':'Long\nprogram','genes':['HLA-DRA'],'color':'#D5E5EF'},{'name':'Remaining markers','short_label':'Other\nmarkers','genes':['B','C','D']}]}))
    cfg['inputs']['marker_yaml']='style_markers.yaml'
    styled=run(modality+'_auto',cfg,script)
    layouts=json.loads((styled/'00_contract/figure_layout_contract.json').read_text())
    for level,n in [('single',1),('labels',3),('cells',6)]:
        img=Image.open(styled/f'01_figures/{level}_seurat_dotplot.png')
        assert img.height == int(max(3.2,2.2+.30*n)*240), (level,img.size)
    single=pd.read_csv(styled/'02_numeric_sources/single_seurat_dotplot_values.csv.gz')
    assert single['avg_exp_zscore_unclipped'].isna().all()
    np.testing.assert_allclose(single['avg.exp.scaled'],single['avg_exp_log1p'],rtol=0,atol=1e-9)
    # Explicit dimensions/fonts must survive all the way to PNG and layout.
    long=a.copy(); long.var_names=['GENE-WITH-LONG-NAME','B','C','D']; long.write_h5ad(root/f'{modality}_long.h5ad')
    (root/f'{modality}_long_markers.yaml').write_text((root/'style_markers.yaml').read_text().replace('HLA-DRA','GENE-WITH-LONG-NAME'))
    cfg['inputs']['counts_h5ad']=f'{modality}_long.h5ad'
    cfg['inputs']['marker_yaml']=f'{modality}_long_markers.yaml'
    cfg['output']['height_inches']=5.5; cfg['output']['width_inches']=10
    cfg['style']={'gene_font_pt':10,'group_font_pt':12,'legend_font_pt':12,'strip_font_pt':12,'title_font_pt':16,'strip_height_inches':.6}
    explicit=run(modality+'_override',cfg,script)
    assert Image.open(explicit/'01_figures/single_seurat_dotplot.png').size == (2400,1320)
    reports.append({'modality':modality,'numeric_baseline_parity':'PASS' if base is not None else 'NOT_REQUESTED','single_group_native_parity':'PASS','auto_height':'PASS','explicit_overrides':'PASS'})
# Missing positions: one all-zero marker has native NaNs across multiple groups.
zero=copy.deepcopy(config); zero['inputs']['counts_h5ad']='zero.h5ad'
b=a.copy(); b.layers['counts']=sparse.csr_matrix(np.column_stack([counts[:,:3],np.zeros(6)])); b.X=b.layers['counts'].copy(); b.write_h5ad(root/'zero.h5ad')
run('constant_gene',zero,skills/'render-scrna-dotplot/scripts/render_scrna_dotplot.py')
report={'status':'PASS','reports':reports,'constant_gene_missing_positions':'PASS','work_dir':str(root)}
(root/'style_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report))
