"""Execute retained-count stages and layout-only behavior on synthetic data."""
import unittest
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import sc_contracts as c


def analyze(counts, names):
    a=ad.AnnData(counts.copy(),obs=pd.DataFrame(index=names))
    a.layers["counts"]=counts.copy()
    a.uns["matrix_sources"]={"counts":{"source_artifact":"synthetic raw count fixture"}}
    c.prepare_log(a)
    attempts=c.select_hvg(a,n_top_genes=80)
    work=a[:,a.var.highly_variable].copy()
    work.X=work.layers["log1p"].copy()
    sc.pp.scale(work)
    sc.pp.pca(work,n_comps=10,random_state=2)
    a.obsm["X_pca"]=work.obsm["X_pca"].copy()
    sc.pp.neighbors(a,n_neighbors=10,use_rep="X_pca",key_added="baseline",random_state=2)
    sc.tl.leiden(a,adjacency=a.obsp["baseline_connectivities"],key_added="cluster",
                 resolution=.5,random_state=2,flavor="igraph",n_iterations=2,directed=False)
    sc.tl.umap(a,neighbors_key="baseline",random_state=2)
    a.uns["stage_fixture"]={"counts":c.matrix_fingerprint(counts),"hvg_attempts":str(attempts)}
    return a


class RerunLayout(unittest.TestCase):
    def test_filter_recomputes_and_layout_preserves_clusters_graph_counts(self):
        rng=np.random.default_rng(9)
        counts=sparse.csr_matrix(rng.poisson(np.linspace(.1,30,300),size=(70,300)).astype(float))
        names=[f"c{i}" for i in range(70)]
        full=analyze(counts,names)
        clean=analyze(counts[:50],[*names[:50]])
        self.assertEqual(clean.obsp["baseline_connectivities"].shape,(50,50))
        self.assertEqual(clean.obsm["X_umap"].shape,(50,2))
        self.assertNotEqual(full.uns["stage_fixture"]["counts"],clean.uns["stage_fixture"]["counts"])
        self.assertFalse(np.allclose(clean.obsm["X_pca"],full.obsm["X_pca"][:50]))
        labels=clean.obs.cluster.copy()
        graph=c.matrix_fingerprint(clean.obsp["baseline_connectivities"])
        raw=c.matrix_fingerprint(clean.layers["counts"])
        cells=c.key_fingerprint(clean.obs_names)
        c.update_umap_layout(clean,neighbors_key="baseline",expected_cells=cells,
                             expected_graph=graph,random_state=19,min_dist=.1)
        pd.testing.assert_series_equal(clean.obs.cluster,labels)
        self.assertEqual(graph,c.matrix_fingerprint(clean.obsp["baseline_connectivities"]))
        self.assertEqual(raw,c.matrix_fingerprint(clean.layers["counts"]))
        self.assertEqual(clean.uns["umap_layout_provenance"]["mode"],"layout_only")
        with self.assertRaises(ValueError):
            c.update_umap_layout(clean[:-1].copy(),neighbors_key="baseline",expected_cells=cells,expected_graph=graph)


if __name__ == "__main__": unittest.main(verbosity=2)
