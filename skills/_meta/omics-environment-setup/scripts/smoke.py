import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
from scipy import sparse
from PIL import Image
out = Path(sys.argv[1])
rng = np.random.default_rng(714)
a = sc.AnnData(sparse.csr_matrix(rng.poisson(2, (120, 300)).astype('float32')))
a.obs_names = [f'c{i}' for i in range(a.n_obs)]
a.var_names = [f'g{i}' for i in range(a.n_vars)]
a.layers['counts'] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=100)
sc.pp.pca(a, n_comps=15)
sc.pp.neighbors(a, n_neighbors=10, n_pcs=15)
sc.tl.leiden(a, flavor='igraph', n_iterations=2, random_state=714)
sc.tl.umap(a, random_state=714)
h = harmonypy.run_harmony(a.obsm['X_pca'], pd.DataFrame({'batch':['a','b']*60}), 'batch', nclust=4, max_iter_harmony=2, verbose=False)
assert np.isfinite(h.Z_corr).all()
a.write_h5ad(out/'synthetic.h5ad'); b = sc.read_h5ad(out/'synthetic.h5ad')
assert b.shape == (120, 300) and b.obsm['X_umap'].shape == (120,2)
assert (b.layers['counts'] != a.layers['counts']).nnz == 0
assert 'connectivities' in b.obsp and 'leiden' in b.obs
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4,4)); ax.scatter(*b.obsm['X_umap'].T, s=8)
ax.set(xlabel='UMAP 1', ylabel='UMAP 2'); fig.tight_layout(); fig.savefig(out/'umap.png', dpi=120)
with Image.open(out/'umap.png') as image: image.load()
print('PASS counts/HVG/PCA/neighbors/Leiden/UMAP/Harmony/H5AD/PNG')
