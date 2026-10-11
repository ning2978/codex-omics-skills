"""Portable spatial stage validation and complete-artifact resume. No project imports."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(2**20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def assert_new_outputs(paths):
    existing = [str(p) for p in paths if Path(p).exists()]
    if existing:
        raise FileExistsError(f'Existing outputs protected: {existing}')


def environment():
    versions = {}
    for package in ('numpy', 'pandas', 'scipy', 'scikit-learn', 'anndata', 'pyarrow',
                    'matplotlib', 'Pillow', 'commot', 'POT'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return {'python': sys.executable, 'python_version': platform.python_version(),
            'platform': platform.platform(), 'packages': versions}


def reopen(path):
    path = Path(path)
    if path.stat().st_size == 0:
        raise ValueError(f'Empty output: {path}')
    if path.suffix == '.json':
        json.loads(path.read_text())
    elif path.suffix == '.png':
        from PIL import Image
        with Image.open(path) as img:
            img.verify()
    elif path.suffix == '.npz':
        with np.load(path, allow_pickle=False) as data:
            for key in data.files:
                data[key]
    elif path.suffix == '.h5ad':
        import anndata
        obj = anndata.read_h5ad(path, backed='r')
        obj.file.close()
    elif path.name.endswith(('.csv', '.csv.gz')):
        pd.read_csv(path)


class Run:
    """Refuse partial/mismatched results; reuse only complete hash-checked outputs."""

    def __init__(self, out_dir, config, inputs, tool):
        self.out = Path(out_dir)
        self.environment = environment()
        tool = Path(tool).resolve()
        tool_files = sorted(set(tool.parent.glob('*.py')) | set(Path(__file__).parent.glob('*.py')))
        self.signature = {
            'config': config, 'environment': self.environment,
            'inputs': {str(Path(p).resolve()): file_sha256(p) for p in inputs},
            'tools': {str(p): file_sha256(p) for p in tool_files},
        }
        self.fingerprint = hashlib.sha256(json.dumps(self.signature, sort_keys=True).encode()).hexdigest()
        self.resumed = False
        manifest = self.out / 'run_manifest.json'
        if self.out.exists() and any(self.out.iterdir()):
            if not manifest.is_file():
                raise ValueError('Incomplete output directory; choose a new directory')
            old = json.loads(manifest.read_text())
            if not old.get('complete') or old.get('fingerprint') != self.fingerprint or old.get('signature') != self.signature:
                raise ValueError('Incomplete or mismatched run; choose a new output directory')
            if not old.get('outputs'):
                raise ValueError('Manifest has no complete output inventory')
            present = {str(p.relative_to(self.out)) for p in self.out.rglob('*')
                       if p.is_file() and p.name != 'run_manifest.json'}
            if present != set(old['outputs']):
                raise ValueError('Output inventory changed; resume rejected')
            for relative, digest in old['outputs'].items():
                p = self.out / relative
                if not p.is_file() or file_sha256(p) != digest:
                    raise ValueError(f'Missing/corrupt resume artifact: {relative}')
                reopen(p)
            self.resumed = True
        else:
            self.out.mkdir(parents=True, exist_ok=True)

    def finish(self, validation, **metadata):
        outputs = {}
        for p in sorted(self.out.rglob('*')):
            if p.is_file() and p.name != 'run_manifest.json':
                if p.suffix in {'.pdf', '.tmp'}:
                    raise ValueError(f'Unexpected output: {p}')
                reopen(p)
                outputs[str(p.relative_to(self.out))] = file_sha256(p)
        if not outputs:
            raise ValueError('No outputs produced')
        write_json(self.out/'run_manifest.json', {
            'schema_version': 1, 'complete': True, 'fingerprint': self.fingerprint,
            'signature': self.signature, 'validation': validation,
            'executed_at': datetime.now(timezone.utc).isoformat(),
            'outputs': outputs,
            'authority': self.signature['config'].get('authority', 'exploratory'),
            'upstream': self.signature['config'].get('upstream', []), **metadata,
        })


def cell_table(path, config, *, labels=False):
    table = pd.read_csv(path, dtype={config['cell_key']: str, config['roi_key']: str,
                                    config['biological_unit_key']: str})
    required = [config[k] for k in ('cell_key', 'roi_key', 'biological_unit_key', 'x', 'y')]
    if not config.get('coordinate_unit'):
        raise ValueError('Explicit coordinate_unit required')
    if labels:
        required.append(config['phenotype_key'])
        if not config.get('label_evidence_grade'):
            raise ValueError('Explicit label_evidence_grade required')
    if set(required) - set(table):
        raise ValueError(f'Missing cell columns: {sorted(set(required)-set(table))}')
    keys = table[config['cell_key']]
    if keys.isna().any() or keys.eq('').any() or keys.duplicated().any():
        raise ValueError('Cell key must be unique and nonmissing')
    if table[[config['roi_key'], config['biological_unit_key']]].isna().any().any():
        raise ValueError('Missing ROI or biological unit')
    if (table.groupby(config['roi_key'])[config['biological_unit_key']].nunique() != 1).any():
        raise ValueError('Each ROI must be nested in exactly one biological unit')
    if not np.isfinite(table[[config['x'], config['y']]].to_numpy(float)).all():
        raise ValueError('Nonfinite coordinates')
    return table


def mask(table, key):
    if key is None:
        return np.ones(len(table), dtype=bool)
    values = table[key]
    if values.isna().any() or not values.isin([True, False, 0, 1]).all():
        raise ValueError(f'{key}: mask must contain explicit bool/0/1 values')
    return values.to_numpy(bool)


def matrix(adata, source, kind):
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('Duplicate cell or gene keys')
    if source == 'X':
        value = adata.X
    elif source.startswith('layers/') and source[7:] in adata.layers:
        value = adata.layers[source[7:]]
    else:
        raise ValueError('Explicit existing matrix required: X or layers/<name>')
    if kind not in {'counts', 'normalized', 'log1p'}:
        raise ValueError('Matrix kind must be counts, normalized or log1p')
    values = value.data if sparse.issparse(value) else np.asarray(value).ravel()
    if value.shape != adata.shape or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Expression must be finite, nonnegative and cells by genes')
    if kind == 'counts' and not np.allclose(values, np.rint(values), rtol=0, atol=1e-8):
        raise ValueError('Counts matrix is not integer-valued')
    return value.copy()


def neighborhoods(table, config, centers, pool, *, radius=None, k=None):
    """Independent centers/contributors, explicit self exclusion and ROI isolation."""
    from scipy.spatial import cKDTree
    if (radius is None) == (k is None):
        raise ValueError('Choose exactly one radius or KNN definition')
    if radius is not None and (not np.isfinite(radius) or radius <= 0):
        raise ValueError('Radius must be finite and positive')
    if k is not None and (not isinstance(k, int) or k < 1):
        raise ValueError('KNN k must be a positive integer')
    coords = table[[config['x'], config['y']]].to_numpy(float)
    roi = table[config['roi_key']].to_numpy()
    answer = {int(i): np.empty(0, dtype=int) for i in centers}
    for name in pd.unique(roi[centers]):
        ci = np.asarray(centers)[roi[centers] == name]
        pi = np.asarray(pool)[roi[pool] == name]
        if not len(pi):
            continue
        tree = cKDTree(coords[pi])
        if radius is not None:
            hits = tree.query_ball_point(coords[ci], radius)
            for center, local in zip(ci, hits):
                indices = pi[np.asarray(local, dtype=int)]
                answer[int(center)] = indices[indices != center]
        else:
            # Deterministic distance then stable cell-key tie order; bounded by k.
            for center in ci:
                count = min(k + 1, len(pi))
                dist, _ = tree.query(coords[center], k=count)
                cutoff = float(np.max(np.atleast_1d(dist)))
                candidates = pi[np.asarray(tree.query_ball_point(coords[center], np.nextafter(cutoff, np.inf)), dtype=int)]
                candidates = candidates[candidates != center]
                distance = np.linalg.norm(coords[candidates] - coords[center], axis=1)
                keys = table.iloc[candidates][config['cell_key']].to_numpy(str)
                answer[int(center)] = candidates[np.lexsort((keys, distance))][:k]
    return answer
