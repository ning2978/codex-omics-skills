"""Lightweight resegmentation pilot input/handoff check; never loads a model."""
from pathlib import Path
import argparse
import json
import sys

shared = next(p/'_shared/spatial' for p in Path(__file__).resolve().parents if (p/'_shared/spatial/run_contract.py').is_file())
sys.path.insert(0, str(shared))
from run_contract import Run, write_json

import numpy as np
import pandas as pd


def validate(cfg):
    if cfg['method_class'] not in {'image_segmentation','transcript_assignment'}:
        raise ValueError('Distinguish image segmentation from transcript-assignment inference')
    for key in ['method','coordinate_unit','coordinate_frame','orientation','roi_id','biological_unit',
                'previous_cell_set_version','new_cell_set_version','previous_assignment_version','new_assignment_version']:
        if not cfg.get(key):
            raise ValueError(f'Explicit {key} required')
    if cfg['previous_cell_set_version'] == cfg['new_cell_set_version'] or cfg['previous_assignment_version'] == cfg['new_assignment_version']:
        raise ValueError('New cells and transcript assignments must use new version identifiers')
    if cfg['orientation'] not in {'y_down','y_up'}:
        raise ValueError('Explicit orientation must be y_down or y_up')
    bbox = np.asarray(cfg['pilot_bbox'],float)
    if bbox.shape != (4,) or not np.isfinite(bbox).all() or bbox[0]>=bbox[1] or bbox[2]>=bbox[3]:
        raise ValueError('Finite pilot_bbox xmin,xmax,ymin,ymax required')
    records = pd.read_parquet(cfg['transcripts']) if str(cfg['transcripts']).endswith('.parquet') else pd.read_csv(cfg['transcripts'])
    key = cfg['transcript_key']
    required = [key,cfg['transcript_x'],cfg['transcript_y'],cfg['gene_key']]
    if set(required)-set(records) or records[required].isna().any().any() or records[key].duplicated().any():
        raise ValueError('Unique transcript key and complete coordinate/gene records required')
    xy = records[[cfg['transcript_x'],cfg['transcript_y']]].to_numpy(float)
    if not np.isfinite(xy).all():
        raise ValueError('Nonfinite transcript coordinates')
    if cfg.get('image'):
        import tifffile
        with tifffile.TiffFile(cfg['image']) as image:
            shape = image.series[0].shape
            axes = image.series[0].axes
        if not np.isfinite(cfg.get('pixel_size',0)) or cfg.get('pixel_size',0)<=0 or not cfg.get('image_coordinate_frame'):
            raise ValueError('Image pixel calibration and coordinate frame required')
        if cfg['image_coordinate_frame'] != cfg['coordinate_frame']:
            affine = np.asarray(cfg.get('image_to_transcript_affine'),float)
            if affine.shape!=(3,3) or not np.isfinite(affine).all() or not np.allclose(affine[2],[0,0,1]) or abs(np.linalg.det(affine))<1e-12:
                raise ValueError('Differing coordinate frames need a verified invertible 2D affine')
            if cfg.get('registration_status') != 'validated':
                raise ValueError('Image registration must be validated independently')
        channel_axis = 'C' if 'C' in axes else ('S' if 'S' in axes else None)
        n_channels = shape[axes.index(channel_axis)] if channel_axis else 1
        channels = cfg.get('channels',{})
        if not channels or any(not isinstance(i,int) or not 0<=i<n_channels for i in channels.values()):
            raise ValueError('Explicit available stain/channel indices required')
        if cfg['method_class']=='image_segmentation' and not cfg.get('segmentation_channels'):
            raise ValueError('Select segmentation channels explicitly')
        if any(name not in channels for name in cfg.get('segmentation_channels',[])):
            raise ValueError('Segmentation channel missing from image contract')
    elif cfg['method_class']=='image_segmentation':
        raise ValueError('Image segmentation requires a stain image')
    else:
        shape,axes = None,None
        if not cfg.get('prior_geometry') and not cfg.get('prior_geometry_not_required_reason'):
            raise ValueError('Declare prior geometry or why selected assignment method needs none')
    selected = (xy[:,0]>=bbox[0]) & (xy[:,0]<=bbox[1]) & (xy[:,1]>=bbox[2]) & (xy[:,1]<=bbox[3])
    if not selected.any():
        raise ValueError('Pilot ROI contains no transcripts')
    return {'status':'PASS','input_contract':'VALIDATED','model_execution':'NOT_RUN',
            'model_validation':'UNVERIFIED','image_shape':list(shape) if shape else None,'image_axes':axes,
            'n_transcripts':len(records),'n_pilot_transcripts':int(selected.sum()),
            'new_cell_set_version':cfg['new_cell_set_version'],'new_assignment_version':cfg['new_assignment_version'],
            'required_downstream_rebuild':['counts','cell_qc','annotation','spatial_analysis']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    inputs = [args.config,cfg['transcripts']] + [cfg[k] for k in ['image','prior_geometry'] if cfg.get(k)]
    run = Run(args.out_dir,cfg,inputs,__file__)
    if run.resumed:
        print('PASS complete-artifact resume; model remains NOT_RUN')
        return
    result = validate(cfg)
    write_json(args.out_dir/'input_validation.json',result)
    run.finish({'status':'PASS','scope':'pilot_input_and_handoff_only','model_validation':'UNVERIFIED'})
    print('PASS input/handoff check; segmentation model NOT_RUN')


if __name__ == '__main__':
    main()
