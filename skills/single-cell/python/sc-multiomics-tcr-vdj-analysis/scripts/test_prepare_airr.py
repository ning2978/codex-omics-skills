import tempfile
from pathlib import Path
import unittest

import anndata as ad
import awkward as ak
import numpy as np
import pandas as pd
from scipy import sparse
import scirpy as ir

from prepare_airr import attach_airr, export_airr_metadata


class Repertoire(unittest.TestCase):
    def test_duplicate_barcode_multicontig_shuffled_rows_and_missing_receptor(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=[]
            for chain in ["TRA","TRB"]:
                rows.append({"barcode":"AA-1","is_cell":True,"high_confidence":True,
                             "contig_id":"AA-1_"+chain,"chain":chain,"productive":True,
                             "cdr3":"CASSLGQETQYF","cdr3_nt":"TGTGCCAGCAGCCTGGGGCAGGAGACCCAGTACTTC",
                             "umis":10,"reads":100,"v_gene":"TRAV1-1" if chain=="TRA" else "TRBV1",
                             "d_gene":None,"j_gene":"TRAJ1" if chain=="TRA" else "TRBJ1-1",
                             "c_gene":"TRAC" if chain=="TRA" else "TRBC1"})
            for lib in ["libA","libB"]:
                pd.DataFrame(rows[::-1]).to_csv(Path(tmp)/(lib+".csv"),index=False)
            rna=ad.AnnData(sparse.csr_matrix(np.ones((3,2))),obs=pd.DataFrame({
                "capture_library":["libB","libA","libA"],"barcode":["AA-1","AA-1","missing"],
                "sample_id":["sampleB","sampleA","sampleA"]},index=["B:A","A:A","A:missing"]))
            mdata,report=attach_airr(rna,{lib:Path(tmp)/(lib+".csv") for lib in ["libA","libB"]})
            self.assertEqual(mdata.mod["rna"].obs_names.tolist(),rna.obs_names.tolist())
            airr=mdata.mod["airr"]
            self.assertEqual(airr.n_obs,2)
            self.assertEqual(ak.num(airr.obsm["airr"]).to_list(),[2,2])
            ir.pp.index_chains(airr)
            ir.tl.chain_qc(airr)
            ir.pp.ir_dist(airr,metric="identity",sequence="nt")
            ir.tl.define_clonotypes(airr,receptor_arms="all",dual_ir="primary_only",
                                   within_group="receptor_type",n_jobs=1)
            ir.tl.clonal_expansion(airr,expanded_in="sample_id")
            out=export_airr_metadata(mdata)
            self.assertEqual(out.n_obs,3)
            self.assertTrue(pd.isna(out.obs.loc["A:missing","clone_id"]))
            self.assertEqual(out.obs.loc["B:A","sample_id"],"sampleB")
            self.assertEqual(report.matched.sum(),2)
            bad=pd.DataFrame(rows+[rows[0]])
            bad.to_csv(Path(tmp)/"duplicate.csv",index=False)
            with self.assertRaises(ValueError): attach_airr(rna,{"libA":Path(tmp)/"duplicate.csv"})


if __name__ == "__main__": unittest.main(verbosity=2)
