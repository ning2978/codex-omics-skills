# Visium HD H&E Coordinate Contract

10x Visium HD spatial coordinates commonly come from `spatial/tissue_positions.parquet`.

Use these columns as the source of truth:

- `pxl_col_in_fullres`: full-resolution image x/column coordinate.
- `pxl_row_in_fullres`: full-resolution image y/row coordinate.

For plotting on embedded 10x images:

- If using `images["hires"]`, compute:
  - `he_col = pxl_col_in_fullres * tissue_hires_scalef`
  - `he_row = pxl_row_in_fullres * tissue_hires_scalef`
- If using `images["lowres"]`, compute:
  - `he_col = pxl_col_in_fullres * tissue_lowres_scalef`
  - `he_row = pxl_row_in_fullres * tissue_lowres_scalef`

Do not invert the y axis when plotting over `imshow(image)`: image coordinates already use row-down orientation. Invert y only for coordinate-only scatter plots without an image background.

For CellCharter multi-sample outputs, always group by the sample/library key before selecting the corresponding image and scalefactor.
