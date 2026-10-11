---
name: xenium-segmentation-new
description: Inspect and render existing Xenium cells and all cell/nucleus boundary rings with portable parquet QC and PNG tools. Use for segmentation integrity or polygon display; model-based resegmentation routes separately.
---

# Inspect existing Xenium segmentation

Inspect existing segmented cells without changing their cell set. Boundary-stain segmentation is accepted; keep segmentation error, compartment mismatch, adjacent transcripts and biological mixed states as distinct review alternatives. Route a requested new segmentation/assignment to `xenium-resegmentation`.

Read [the ring contract](references/rings-and-tools.md) before reading or rendering boundaries. Read [stage handoffs](../../../shared/workflow-st-analysis/references/stage-contract.md) when chaining stages. A geometry-only check needs no expression matrix or annotation registry; consuming labels requires their project adapter.

The bundled tool works outside this repository:

```bash
python <skill-dir>/scripts/segmentation.py --input-dir <vendor-output> \
  --out-dir <new-output> --bbox <xmin> <xmax> <ymin> <ymax>
```

Omit bbox for the full field. Nuclei are drawn by default; `--no-nucleus` explicitly hides the overlay while QC still audits nuclei. `--color-by segmentation_method` colors cell metadata. `--upstream-json` records earlier-stage provenance. Existing complete results are reused only when input, tool, environment, parameters and all output hashes match. Partial directories and mismatches are rejected.

Verify all `(cell_id,label_id)` rings survive, source vertex order is preserved across row groups, unsorted reading works, duplicate cell keys fail, all nuclei are drawn, y increases downward and PNGs decode. Missing/truncated rings are failures or explicit QC warnings; neither PASS nor WARN validates biological segmentation accuracy. The tools load full boundary tables before clipping complete polygons; budget memory for a full slide.

Keep project package commands and registry/cohort rules in the repository adapter. Historical readers and their outputs remain unchanged.
