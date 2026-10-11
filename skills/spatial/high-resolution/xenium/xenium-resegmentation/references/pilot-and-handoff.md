# Select, pilot and compare

Image instance segmentation uses stain morphology to infer masks: Cellpose, StarDist, Mesmer or vendor image segmentation depending on channels/model training. Verify 2D/3D input, resolution, nuclear/membrane/interior channel meaning, TIFF axes, model/version, object diameter and stain-domain applicability. A nuclei-only model is not a validated whole-cell segmentation. Transcript assignment to masks is a separate step with explicit outside-mask/overlap handling.

Transcript-based methods such as Baysor/Proseg infer attribution/cell structure from transcript coordinates and gene information, possibly with a segmentation prior. Record algorithm-specific inputs, prior strength, molecule quality filtering and treatment of unassigned molecules. Transcript-derived boundaries need not equal morphological cell edges. These methods are not image segmenters and should not be described as rerunning stain segmentation.

Input checker config:

```json
{
  "method_class":"image_segmentation", "method":"cellpose",
  "image":"/path/pilot.ome.tif", "channels":{"DAPI":0,"membrane":1},
  "segmentation_channels":["DAPI","membrane"], "pixel_size":0.5,
  "image_coordinate_frame":"registered_um", "coordinate_frame":"registered_um",
  "orientation":"y_down", "coordinate_unit":"um",
  "transcripts":"/path/pilot_transcripts.csv", "transcript_key":"transcript_id",
  "transcript_x":"x", "transcript_y":"y", "gene_key":"gene",
  "roi_id":"roi1", "biological_unit":"donor1", "pilot_bbox":[0,100,0,100],
  "previous_cell_set_version":"vendor_v1", "new_cell_set_version":"pilot_v2",
  "previous_assignment_version":"vendor_v1", "new_assignment_version":"pilot_v2",
  "upstream":[]
}
```

For transcript assignment, use `method_class=transcript_assignment`, the exact method and a `prior_geometry` input or `prior_geometry_not_required_reason`; image is optional. Different image/transcript frames require `image_to_transcript_affine` (3x3), `registration_status=validated` and independent registration evidence. A matrix passing algebra checks does not validate registration scientifically. Pixel size uses the declared coordinate unit per pixel. Select a pilot transcript table first; the lightweight checker reads it in memory, not a full slide. Every transcript ID is unique/nonmissing, genes and coordinates are complete, and the pilot bbox contains transcripts. Preserve quality/unassigned status columns separately.

Choose ROIs from declared tissue/density/stain strata before inspecting desired outcomes. Compare original and new output using cell count, area/nucleus distributions, valid polygon topology, nucleus multiplicity, mask overlap (IoU), split/merge counts, transcript assignment fraction/unassigned fraction and per-cell count/marker consistency. Record manual/blinded reference and uncertainty; no single metric proves improvement. Match old/new cells through explicit geometry overlap, retaining one-to-many/many-to-one cases rather than forcing a bijection. Report model failures and ROI scope.

Handoff includes new mask/rings, versioned cell IDs, molecule-to-cell assignment (including unassigned), cell-by-gene counts with stable gene IDs, old/new relationship audit, model/environment/parameter fingerprints, comparison summaries and PNGs. Rebuild downstream matrices/QC/annotation/spatial features; expression reclustering recomputes counts-derived HVG/PCA/correction/neighbors/UMAP/Leiden. Old annotations cannot be copied by row order. Only tested model paths are marked validated; this round validates input/handoff rules, not model execution or real-data segmentation accuracy.
