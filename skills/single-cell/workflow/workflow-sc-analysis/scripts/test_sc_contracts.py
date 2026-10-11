"""Behavioral regressions for matrix scope, keys, design and chunked reads."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import sc_contracts as c

ROOT = Path(__file__).resolve().parents[4]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture():
    counts = sparse.csr_matrix([[3,1,5,0],[0,1,5,0],[3,1,5,1],[0,1,5,1],
                               [0,3,5,0],[0,3,5,0],[0,3,5,1],[0,3,5,1]])
    a = ad.AnnData(counts.astype(float), obs=pd.DataFrame({"cluster":["A"]*4+["B"]*4}, index=[f"c{i}" for i in range(8)]),
                   var=pd.DataFrame({"gene_symbol":["duplicate","duplicate","house","rare"]}, index=["ID1","ID2","ID3","ID4"]))
    a.layers["counts"] = counts.copy()
    a.uns["matrix_sources"] = {"counts":{"kind":"raw_counts","source_artifact":"synthetic fixture"}}
    return a


class Contracts(unittest.TestCase):
    def test_full_markers_ignore_hvg_scaled_x_and_raw(self):
        import scanpy as sc
        a = fixture()
        before = c.matrix_fingerprint(a.layers["counts"])
        c.prepare_log(a)
        a.var["highly_variable"] = [False,True,True,False]
        a.raw = a[:,["ID2"]].copy()
        a.X = a.layers["log1p"].toarray()
        sc.pp.scale(a)
        self.assertEqual(np.count_nonzero(a.X[:4,0])/4,1.)
        result = c.rank_markers(a,"cluster",min_cells=1)
        classic = result[(result.group=="A") & (result.gene=="ID1")].iloc[0]
        self.assertEqual(classic.pct_in,0.5)
        self.assertEqual(classic.pct_out,0.)
        self.assertEqual(set(result.gene),set(a.var_names))
        self.assertEqual(len(result),8)
        self.assertEqual(before,c.matrix_fingerprint(a.layers["counts"]))
        a.layers["log1p"] = a.X.copy()
        with self.assertRaises(ValueError): c.rank_markers(a,"cluster",min_cells=1)

    def test_hvg_fallback_uses_log_retains_batch_and_records_failure(self):
        import scanpy as sc
        rng = np.random.default_rng(17)
        a = ad.AnnData(sparse.csr_matrix(rng.poisson(np.linspace(.05,20,500),size=(100,500)).astype(float)))
        a.layers["counts"] = a.X.copy()
        a.uns["matrix_sources"] = {"counts":{"source_artifact":"synthetic"}}
        a.obs["capture"] = ["a"]*50+["b"]*50
        c.prepare_log(a)
        real = sc.pp.highly_variable_genes
        def hvg(obj, **kwargs):
            if kwargs["flavor"] == "seurat_v3": raise ValueError("LOESS fixture failure")
            self.assertEqual(kwargs["layer"],"log1p")
            self.assertEqual(kwargs["batch_key"],"capture")
            return real(obj,**kwargs)
        with patch.object(sc.pp,"highly_variable_genes",side_effect=hvg):
            attempts = c.select_hvg(a,n_top_genes=40,batch_key="capture")
        self.assertIn("LOESS",attempts[0]["error"])
        self.assertEqual(attempts[1]["status"],"PASS")
        with patch.object(sc.pp,"highly_variable_genes",side_effect=ValueError("failed")) as fn:
            with self.assertRaises(RuntimeError): c.select_hvg(a,batch_key="capture")
            self.assertEqual(fn.call_count,2)

    def test_capture_keys_shuffled_and_duplicate_symbols(self):
        keys = c.cell_keys(["a","b","a:x"],["same","same","x"])
        self.assertTrue(keys.is_unique)
        frame = pd.DataFrame({"value":[3,2,1]},index=keys[::-1])
        self.assertEqual(c.align_frame(frame,keys).value.tolist(),[1,2,3])
        with self.assertRaises(ValueError): c.cell_keys(["a","a"],["same","same"])
        with self.assertRaises(ValueError): c.align_frame(frame.iloc[:2],keys)
        with self.assertRaises(ValueError): c.align_frame(pd.concat([frame,frame.iloc[:1]]),keys)

    def test_replication_and_confounding(self):
        samples = pd.DataFrame({"biological_unit":["a","b","c","d"],"donor":["a","b","c","d"],
                                "condition":["control","control","case","case"],"technical_batch":["x","x","y","y"]})
        report = c.assess_design(samples)
        self.assertEqual(report["inference_status"],"descriptive_only")
        self.assertTrue(any("confounded" in r for r in report["reasons"]))
        samples["technical_batch"] = "x"
        self.assertEqual(c.assess_design(samples)["inference_status"],"eligible_for_model")
        self.assertEqual(c.assess_design(samples.iloc[[0,2]])["inference_status"],"descriptive_only")
        self.assertEqual(c.assess_design(samples.drop(columns="biological_unit"))["inference_status"],"descriptive_only")
        technical = pd.concat([samples.iloc[[0]],samples.iloc[[0]],samples.iloc[[2]],samples.iloc[[2]]])
        self.assertEqual(c.assess_design(technical)["inference_status"],"descriptive_only")
        pairs = pd.DataFrame({"biological_unit":["a0","a1","b0","b1"],"donor":["a","a","b","b"],"condition":["C","D","C","D"]})
        self.assertEqual(c.assess_design(pairs,paired=True)["inference_status"],"eligible_for_model")

    def test_pseudobulk_raw_sums_and_no_invented_zero_samples(self):
        a = fixture()
        a.obs["cell_identity"] = ["T"]*4+["B"]*4
        a.obs["biological_unit"] = ["one"]*2+["two"]*2+["one"]*4
        a.obs["condition"] = "control"
        result,meta,ledger = c.pseudobulk(a,min_cells=2)
        self.assertEqual(result.shape[0],3)
        self.assertFalse(((meta.cell_identity=="B") & (meta.biological_unit=="two")).any())
        np.testing.assert_array_equal(np.asarray(result.sum(axis=0)),np.asarray(a.layers["counts"].sum(axis=0)))
        a.layers["counts"] = a.layers["counts"].astype(float)*.3
        with self.assertRaises(ValueError): c.pseudobulk(a,min_cells=2)

    def test_velocity_axis_alignment_and_mode_prerequisites(self):
        a = fixture()
        v = a[::-1,::-1].copy()
        v.layers["spliced"] = v.layers["counts"].copy()
        v.layers["unspliced"] = v.layers["counts"].copy()*2
        out = c.merge_velocity_layers(a,v)
        np.testing.assert_array_equal(out.layers["unspliced"].toarray(),a.layers["counts"].toarray()*2)
        c.validate_velocity(out,"stochastic")
        with self.assertRaises(ValueError): c.validate_velocity(out,"dynamical")
        out.var["fit_alpha"] = 1.
        out.uns["recover_dynamics"] = {"mode":"synthetic prerequisite test"}
        c.validate_velocity(out,"dynamical")
        with self.assertRaises(ValueError): c.merge_velocity_layers(a,v[:-1])

    def test_direct_h5ad_csr_csc_dense_and_full_scan(self):
        preflight = load(ROOT/"_meta/omics-analysis-runtime-preflight/scripts/h5ad_contract.py","preflight")
        with tempfile.TemporaryDirectory() as tmp:
            for format in ("csr","csc","dense"):
                a = fixture()
                a.layers["counts"] = a.layers["counts"].tocsc() if format=="csc" else a.layers["counts"].toarray() if format=="dense" else a.layers["counts"]
                a.layers["irrelevant_large_layer"] = a.layers["counts"].copy()
                path=Path(tmp)/(format+".h5ad")
                a.write_h5ad(path)
                with h5py.File(path,"r") as handle:
                    node=handle["layers/counts"]
                    chunks=list(preflight.matrix_chunks(node,a.n_obs,a.n_vars,3))
                    self.assertTrue(all(len(chunk)<=3 for chunk in chunks))
                args=["--input",str(path),"--kind","scrna","--count-source","layer:counts","--chunk-values","3"]
                with patch.object(ad,"read_h5ad",side_effect=AssertionError("Eager read forbidden")),contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(preflight.main(args),0)
                a.layers["counts"] = (a.layers["counts"].toarray() if sparse.issparse(a.layers["counts"]) else a.layers["counts"]).astype(float)
                a.layers["counts"][-1,-1] = np.nan
                a.write_h5ad(path)
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(preflight.main(args+["--sample-rows","2","--sample-cols","2"]),0)
                    self.assertEqual(preflight.main(args+["--full-count-scan"]),1)


if __name__ == "__main__": unittest.main(verbosity=2)
