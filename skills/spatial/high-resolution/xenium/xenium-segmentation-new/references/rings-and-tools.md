# Ring and portable-tool contract

A single vendor output supplies `cells.parquet`, `cell_boundaries.parquet`, `nucleus_boundaries.parquet`. `cells.cell_id` is unique/nonmissing within this sample. Boundary rows require `cell_id,label_id,vertex_x,vertex_y`, coordinates in microns with top-left origin and y down. Boundaries reference existing cells. Across samples, create an explicit sample-prefixed unique cell key; preserve vendor IDs.

`scripts/rings.py` exports:

- `load_rings(path, cell_ids=None)`: cell -> label -> complete vertices. Accepts interleaved/unsorted rows, preserving relative source row order within each ring.
- `iter_rings(path, chunk_row_groups=1)`: yields complete rings across row-group splits; requires contiguous ring keys and rejects a completed ring reappearing. Use the full reader for unsorted data.
- `ring_vertex_counts(path)`: sums every ring across all row groups without reading vertex coordinates.

Never keep just the largest ring, connect separate nuclei, clip individual vertices before grouping, or take a row-group maximum as a ring's total. Vertex number measures completeness, not biological complexity. For unordered vertices within a ring, require a vendor ordering field and explicit upstream repair; do not sort vertices by x/y.

The tool writes ring-count CSVs, segmentation-method breakdown when present, QC PNG, complete-ring overlay PNG, summaries and `run_manifest.json`. Loaded ring counts describe the full tables, not visible crop intersections. No count matrix or model is generated. Every output is reopened; PNG-only output is enforced.

Dependencies are the existing Python environment's numpy, pandas, scipy, pyarrow, matplotlib and Pillow. Distribution requires this skill directory and `<skills-root>/_shared/spatial/*.py` in their relative tree; no `xenium_downstream` or project registry is imported. Shared run contracts fingerprint both tool and shared modules.

The source readers/renderers were extracted from the prior validated project implementation; the project copy remains frozen. Behavior tests validate multi-ring, cross-row-group and unsorted input in an isolated skill tree. Full-slide memory/performance and segmentation accuracy require a real-data task.
