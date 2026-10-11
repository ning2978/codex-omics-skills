---
name: algorithm-cell2location
description: Use for cell2location-specific reference signature estimation, spatial abundance modeling, training parameter selection, runtime diagnosis, output validation, or review of cell2location deconvolution across Visium-like projects.
---

# Cell2location

## Role

Own cell2location method choices. Pair with `workflow-deconvolution-run` for launch, monitoring, and generic completion handling; pair with `workflow-st-analysis` when spatial inputs are not yet standardized.

## Input Contract

Require:

- a single-cell reference with raw integer counts, trusted label key, and optional batch key;
- a spatial object with raw counts, coordinates, sample/library key, and one gene namespace;
- an explicit shared-gene report after removing duplicates;
- recorded species, reference version, sample scope, and count-layer names.

Do not train from normalized or log-transformed matrices when the model expects counts.

## Model Stages

1. Fit `RegressionModel` on the reference.
2. Export cell-state signatures with stable names and verify every requested label is represented.
3. Set up the spatial `Cell2location` model with a defensible `N_cells_per_location` and `detection_alpha`.
4. Train, export posterior abundance summaries, and preserve model/run metadata.

Keep reference and spatial model artifacts separate so a validated signature model can be reused intentionally.

## Parameter Decisions

Read `references/parameter-profiles.md` whenever selecting concrete epochs, batch size, learning rate, expected cells per location, or detection alpha.

- Use the `recorded-author-workflow` profile only for its documented reproduction target.
- Treat `historical-project-adapted` as provenance, not a universal default.
- Estimate `N_cells_per_location` from registered nuclei/cell segmentation when available; summarize counts per in-tissue spot and use a defensible central value as the pilot center.
- Without image-derived counts, state the assay/tissue prior, test a small grid around it, and label the choice uncertain.
- Keep `detection_alpha` at the recorded profile value unless technical evidence motivates a sensitivity run; do not tune it toward a preferred biological result.
- Record every deviation and the reason; never rely on hidden package defaults for a final run.

Read `references/profile-language.md` when documenting profile provenance in scripts or reports.

## Runtime Guidance

- Run one main training process per GPU unless parallel allocation is deliberate.
- Use a writable `NUMBA_CACHE_DIR` and a non-interactive runner.
- Distinguish slow reference setup or posterior export from a stalled training loop.
- Require finite loss history without divergence and document whether training plateaued; do not invent one universal convergence threshold across package versions and datasets.
- Delegate PID capture, log polling, duplicate-run checks, and generic failure detection to `workflow-deconvolution-run`.

## Output Contract

Preserve at least:

- reference model and exported signature matrix;
- spatial model directory;
- per-sample or combined `spatial.h5ad`;
- named abundance payload, commonly `means_cell_abundance_w_sf` and/or `q05_cell_abundance_w_sf`, while verifying the actual installed-version keys;
- cell-state names aligned to abundance columns;
- run parameters, package versions, input paths, and sample-level summary;
- training and posterior-export logs.

## Validation

- Reopen the result object and verify `n_obs`, spatial coordinates, sample keys, and abundance dimensions.
- Confirm abundance columns exactly match the exported cell-state names.
- Check for missing, all-zero, non-finite, or implausibly dominant abundance columns.
- Record posterior sampling/export arguments and whether every reported table uses posterior mean, q05, or another quantity.
- Verify observed genes and shared genes match the recorded report.
- Confirm the selected parameter profile and all overrides are present in run metadata.
- Review model diagnostics and spatial marker plausibility before biological interpretation.

## Boundaries

- Use `cell2location-result-visualization` for substantial plotting or grouped result summaries.
- Do not treat abundance as direct cell counts without checking the stored posterior quantity.
- Do not reuse one project's expected cells per spot as a cross-project default.

## Figure presentation

Default to PNG with no explanatory text below the figure. Put methods, provenance
and qualifications in documentation/result records; preserve meaningful axes,
ticks, panel labels, program strips, legends and color bars. Choose titles for
the task and inspect exported PNGs for overlap, cropping and excessive margins.
