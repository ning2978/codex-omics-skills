"""Actual Python/Matrix Market/Seurat parity, using SC_BRIDGE_RSCRIPT when set."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

HERE=Path(__file__).resolve().parent
RSCRIPT=os.environ.get("SC_BRIDGE_RSCRIPT")


@unittest.skipUnless(RSCRIPT,"Set SC_BRIDGE_RSCRIPT to an exact Seurat environment executable")
class Bridge(unittest.TestCase):
    def test_roundtrip_counts_duplicate_symbols_shuffled_metadata_and_reductions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            values=np.array([[1,2,0,8],[3,0,1,5],[4,2,0,6]])
            a=ad.AnnData(sparse.csr_matrix(values),obs=pd.DataFrame({"barcode":["AA","AA","BB"],"sample_id":["s1","s2","s1"]},index=["02","01","libA:AA"]),
                         var=pd.DataFrame({"gene_ids":["ENSG_01","ENSG-01","ID3","ID4"],"gene_symbol":["same","same","G3","G4"]},index=["ENSG_01","ENSG-01","ID3","ID4"]))
            a.layers["counts"]=a.X.copy()
            a.obsm["X_umap"]=np.array([[.25,-.4],[4.125,3.5],[-1.75,8.25]])
            a.obsm["X_pca"]=np.arange(9).reshape(3,3)*.123
            a.write_h5ad(root/"input.h5ad")
            exp=root/"export"
            subprocess.run([sys.executable,str(HERE/"export_h5ad_sparse_for_r.py"),"--input-h5ad",str(root/"input.h5ad"),"--output-dir",str(exp),"--mode","sc","--counts-layer","counts","--export-reductions"],check=True,capture_output=True,text=True)
            for path in [exp/"metadata.csv",exp/"reductions/X_umap.csv",exp/"reductions/X_pca.csv"]:
                tab=pd.read_csv(path,dtype={"barcode":str})
                tab.iloc[::-1].to_csv(path,index=False)
            out=root/"result.rds"
            cmd=[RSCRIPT,"--vanilla",str(HERE/"read_sparse_export_to_seurat.R"),"--input-dir",str(exp),"--output-rds",str(out)]
            result=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            checker=root/"inspect.R"
            checker.write_text('''suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly=TRUE)
o <- readRDS(args[[1]])
jsonlite::write_json(list(cells=colnames(o),counts=as.matrix(GetAssayData(o,layer="counts")),
                        umap=Embeddings(o[["umap"]]),pca=Embeddings(o[["pca"]]),
                        source_ids=o@misc$h5ad_bridge_feature_mapping$source_feature_id,
                        symbols=o@misc$h5ad_bridge_feature_mapping$gene_symbol),args[[2]],auto_unbox=TRUE,digits=NA)
''')
            subprocess.run([RSCRIPT,"--vanilla",str(checker),str(out),str(root/"inspect.json")],check=True,capture_output=True)
            observed=json.loads((root/"inspect.json").read_text())
            self.assertEqual(observed["cells"],a.obs_names.tolist())
            self.assertEqual(observed["source_ids"],a.var_names.tolist())
            self.assertEqual(observed["symbols"],["same","same","G3","G4"])
            np.testing.assert_array_equal(observed["counts"],values.T)
            np.testing.assert_allclose(observed["umap"],a.obsm["X_umap"],rtol=0,atol=1e-12)
            np.testing.assert_allclose(observed["pca"],a.obsm["X_pca"],rtol=0,atol=1e-12)
            # Reject a reduction row duplicated across cells instead of inventing matches.
            path=exp/"reductions/X_umap.csv"
            tab=pd.read_csv(path,dtype={"barcode":str})
            tab.loc[1,"barcode"]=tab.loc[0,"barcode"]
            tab.to_csv(path,index=False)
            bad=subprocess.run(cmd[:-1]+[str(root/"bad.rds")],capture_output=True,text=True)
            self.assertNotEqual(bad.returncode,0)
            self.assertFalse((root/"bad.rds").exists())


if __name__ == "__main__": unittest.main(verbosity=2)
