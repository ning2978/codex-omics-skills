# CellCharter API Notes

## Verified API Surface

Prefer CellCharter installs where these calls are available:

- `cc.gr.aggregate_neighbors`
- `cc.gr.remove_long_links`
- `cc.tl.Cluster`
- `cc.tl.ClusterAutoK`
- `cc.gr.nhood_enrichment`
- `cc.pl.nhood_enrichment`

Wrap clustering and representation building in small helper functions so future CellCharter API changes can be localized.

## Robustness Patterns

Prefer carrying these patterns into related scripts:

- Use `labels_to_ordered_categorical` that tolerates non-integer labels.
- Sanitize categorical NaN values before writing `.h5ad`.
- Save AutoK stability as a long table plus a per-K summary table.
- Pass `spatial_key` explicitly to `squidpy.gr.spatial_neighbors`.
- Log edge counts before and after `remove_long_links`.
- Check cross-sample edges are zero after graph construction.
- Use a separate UMAP neighbors key to avoid overwriting existing neighbors.
- Align cell-type proportion tables by `adata.obs_names`, not raw `.values` position.
- Avoid accidental condition/sample mask leakage in multi-panel spatial plots.

## Feature Transform Notes

For RCTD/cell-type composition workflows, `clr_zscore` is the conservative default because RCTD weights are compositional. Use a data-adaptive pseudocount based on observed positive values when zeros are present.

For expression/scVI workflows, do not use RCTD transforms. Use the scVI latent as the CellCharter input representation.

## Graph Notes

Pixel-space Visium coordinates typically use `coord_type="generic"` with Delaunay graph construction, followed by long-link removal. Canonical Visium array grid coordinates can use `coord_type="grid"`.

When subsetting samples, ensure `uns["spatial"]` does not cause misleading plots by retaining unrelated image keys without matching spots.
