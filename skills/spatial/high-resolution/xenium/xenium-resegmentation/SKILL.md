---
name: xenium-resegmentation
description: Plan or pilot explicitly requested Xenium resegmentation or transcript reassignment, distinguishing stain-image segmentation from transcript-based inference and rebuilding versioned downstream cell artifacts.
---

# Requested Xenium resegmentation

Use only when the task requests a new segmentation or transcript-assignment analysis. An existing segmented-cell QC/plot request routes to `xenium-segmentation-new`. Choose the method for the observed failure and available data, not merely because a model environment exists.

Read [method choice and pilot comparison](references/pilot-and-handoff.md) before model work. The bundled `scripts/check_inputs.py --config <json> --out-dir <new>` checks lightweight inputs and versioned handoff only. A passing input check leaves model execution `NOT_RUN` and model validation `UNVERIFIED`.

Declare image/transcript coordinate frames, units/orientation, stains/channels, ROI and biological unit, original geometry/counts and transcript IDs. Audit image registration before combining coordinate frames. Image instance segmentation and transcript-based cell assignment answer different questions; neither produces validated cell identities automatically.

Choose representative small ROIs spanning density, stain quality and plausible errors. Lock comparison metrics and preserve original segmentation. Run the selected model only within the authorized scope and exact verified environment. No install or full model platform is bundled here.

Compare cell/nucleus geometry, split/merge behavior, transcript assignment, morphology and marker consistency with blinded/manual reference where possible. Keep segmentation error, broad-compartment mismatch, adjacent transcripts and mixed biological states distinct. Boundary-stain alone is accepted evidence of a method, not a removal rule.

Resegmentation changes the cell set and transcript attribution version. Export new unique IDs plus explicit overlap/many-to-many old/new links; rebuild counts, cell QC, annotations and spatial results. Preserve all old IDs, sources and products. New results remain candidate/exploratory until independently validated and approved under the project's change-control contract.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
