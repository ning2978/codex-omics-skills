# scRNA-seq native-Seurat DotPlot input contract

The renderer consumes one run YAML and one marker YAML. Relative paths resolve
against `project_root`; inputs are read-only.

## Run YAML

```yaml
analysis_id: scrna_native_dotplot_v1
status: exploratory_not_current
project_root: /absolute/project/root

inputs:
  counts_h5ad: results/source/full_gene.h5ad
  sidecar: results/source/labels.csv.gz
  counts_layer: counts
  marker_yaml: configs/marker_programs.yaml
  expected_cells: 1000           # optional
  expected_genes: 20000          # optional
  counts_h5ad_sha256: "..."      # optional
  sidecar_sha256: "..."          # optional
  marker_yaml_sha256: "..."      # optional

join:
  sidecar_key: obs_name           # or registry-declared unique cell_id
  target_key: obs_name
  sidecar_key_unique: true
  target_key_unique: true
  row_order_join_forbidden: true

groupings:
  broad:
    field: cell_type_broad_draft
    order: [Epithelial, Stromal, Myeloid, T_NK]
    label_map:                    # optional one-to-one source -> display mapping
      Epithelial_source: Epithelial
      Stromal_source: Stromal
      Myeloid_source: Myeloid
      T_NK_source: T_NK
    allow_many_to_one_label_map: false  # set true only for an explicit merge
    title: scRNA-seq native Seurat DotPlot
    expected_counts:              # optional
      Epithelial: 100
      Stromal: 80
      Myeloid: 60
      T_NK: 40

marker_selection:
  include_categories: []          # empty means all
  exclude_categories: []
  include_programs: []
  exclude_programs: []
  deduplication:
    strategy: first_occurrence    # error, first_occurrence, move_selected_to_program
  minimum_max_group_pct_exp: 0.0
  filter_groupings: [broad]       # omitted means all

dotplot:
  normalization_scale_factor: 10000.0
  scaling_ddof: 1
  col_min: -2.5
  col_mid: 0.0
  col_max: 2.5
  dot_min: 0.0
  scale_by: radius
  dot_scale: 6.0
  color_low: "#2166AC"
  color_mid: "#FFFFFF"
  color_high: "#B2182B"

seurat_native:
  rscript: /absolute/R/environment/bin/Rscript
  expected_r_version: "4.4.3"
  expected_seurat_version: "5.5.1"
  expected_seuratobject_version: "5.4.0"
  avg_exp_abs_tolerance: 1.0e-9
  avg_exp_log1p_abs_tolerance: 1.0e-9
  zscore_abs_tolerance: 1.0e-9
  scaled_abs_tolerance: 1.0e-9
  pct_exp_abs_tolerance: 1.0e-10

project:
  annotation_use_case: scrna_top_marker
  annotation_guide: docs/CURRENT_ANNOTATIONS.md
  annotation_registry: metadata/current_annotation_registry.json
  expected_annotation_contract:   # optional exact guard
    track_id: example_scrna_track_v1
    status: current
  protected_files:
    - docs/CURRENT_ANNOTATIONS.md
    - metadata/current_annotation_registry.json
  run_recluster_guard: false  # enable only with the project guard installed

output:
  root: results/scrna_native_dotplot_v1
  width_inches: 20.0            # requested minimum width
  auto_expand_width: true       # default; false fails if width is too small
  # height_inches: 6.4         # optional explicit override; omitted/null = automatic
  dpi: 240
```

Each grouping must contain at least one unique label. Without `label_map`,
observed labels must equal its order; with `label_map`, source labels must equal
the map keys and its unique values must equal the order. Missing or extra labels
fail closed. A many-to-one map also fails unless
`allow_many_to_one_label_map: true` explicitly authorizes the merge. The
explicit sidecar key must be unique and have
exact set equality with H5AD `obs_names`; shuffled rows are valid. `target_key`
is always the AnnData `obs_name` index. A sidecar key named `cell_id` is valid
only when the selected registry source declares that unique sample-prefixed key
maps to `obs_name`; it is never joined to an AnnData `cell_id` column.

`layers["counts"]` must use AnnData CSR encoding. Counts must be finite,
non-negative, integer-valued, contain no stored explicit zeros, and yield a
positive complete-library size for every cell.

## Marker YAML

```yaml
programs:
  - name: Epithelial
    category: identity
    short_label: Epithelial
    color: "#D5E5EF"       # optional
    genes:
      - {gene: EPCAM, role: core}
      - {gene: KRT8, role: support}
  - name: Myeloid
    category: identity
    short_label: Myeloid
    genes: [LST1, TYROBP, FCER1G]
```

The grouped `programs: {category: [...]}` form and `core`/`support` lists are
also accepted. Scalar genes inherit their list role. Missing colors use the
deterministic selected-program sequence `#D5E5EF`, `#F0DDCB`, `#DAE8D3`, `#E4DBEF`, `#F1D5D8`,
`#D6E8E5`, `#EAE4C9`, `#DDE0EC`.

For `move_selected_to_program`, declare a `shared_program` with `name`, `genes`,
optional `category`, `short_label`, `color`, `role`, and `position`; set
`remaining_duplicates` to `error` or `first_occurrence`.

Filtering follows selection and deduplication. Retain a measured gene when its
maximum `pct.exp` across `filter_groupings` meets the threshold. Unavailable and
below-threshold genes remain audited but do not enter Seurat or the figure.

## Output contract

Each grouping produces one white PNG, one unmodified native Seurat `p$data`
CSV, and one annotated parity numeric table. The package also contains marker
and strip audits, join and annotation-resolution records, Python-Seurat parity,
protected-input before/after SHA256, artifact reopen checks, completion JSON,
figure-source index, validation summary, and a self-verifying SHA256 manifest.

No PDF, H5AD/H5, sidecar, or UMAP may appear in the package. Run
`validate --verify-inputs` before calling it canonical.

Gene labels default to a 60-degree rotation, `hjust=1`, `vjust=1`, and 9 pt
text. The effective width is
`max(width_inches, 6.0 + marker_count * max(0.42, 0.18 + 0.025 *
longest_gene_length))` when `auto_expand_width` is true. Existing numeric-width
YAML remains valid because auto-expansion defaults to true. With
`auto_expand_width: false`, a configured width below the computed lower bound
fails before Seurat rendering. The layout contract, run contract, and
validation summary record the label geometry, requested width, computed lower
bound, and final width.

## Registry binding (shared engine)

Each grouping may specify `annotation_use_case` (defaults to the existing
project setting), `annotation_source_id` (inferred only if its field resolves
to exactly one source), and `annotation_gate_mode: display` (default) or
`fine_analysis`. RCN assignment fields require `assignment_id` containing an
exact registered ID or full namespace. Different roles stay in separate
fields. Unregistered grouping fields fail closed; prepare registered fields
in a derivative sidecar instead of renaming labels or filling missing values.

The guide must contain exactly one `<!-- ROLE_semantics_begin -->` / JSON
code block / `<!-- ROLE_semantics_end -->` summary equal to that role's full
registry default track. When a legacy consumer requires a partial summary,
preserve it and add `ROLE_complete_semantics_begin/end` for the full track;
the reader checks both blocks. Only `updated_at` is ignored. Existing incomplete
summaries must be reviewed and extended before rendering; this is not an
annotation promotion. The role reader can compare complete files while
returning compact contracts:

```bash
python ${SKILLS_ROOT}/_shared/dotplot/annotation_contract.py --root /project --role xenium_top_marker --role epi_subtype
```

The input sidecar may be a filtered derivative, but its keys must exactly
match the counts H5AD and be a subset of its selected registered source.
Every plotted label must equal the source value by key, pass source include
and cohort gates, and respect higher-priority overlaps. Source hashes and
verified gates are recorded. `fine_analysis` additionally requires an explicit
registry fine-analysis gate; absence fails closed. Current paused roles cannot
supply labels. No source input is mutated.

## Shared style overrides

Read [the shared style contract](../../_shared/dotplot/style.md). Omitted/null
height computes `max(3.2, 2.2 + 0.30 * n_groups)` separately for every figure.
Explicit `output.height_inches` is preserved. Width keeps its existing automatic
minimum rule (20 inches by default; no universal 24-inch minimum).

```yaml
style:                         # all optional, positive finite values
  gene_font_pt: 9
  group_font_pt: 11
  strip_font_pt: 11
  legend_font_pt: 11
  title_font_pt: 14
  strip_height_inches: 0.5
```

Set `groupings.<name>.title` and optional `subtitle` for the actual lineage,
grouping and annotation status. Marker `short_label` accepts explicit newlines
(e.g. `"Luminal\nspecific"`) and recorded abbreviations; full program names
remain in audits. Per-program `color` overrides the shared pastel sequence.
Single-group `avg_exp_scaled` is native log1p mean, without clipping or fabricated
z-score; record the lack of cross-group scaling in result records and documentation. Keep the Expression (unscaled) legend, without a bottom caption. Display
colors squish out-of-range values. Numeric parity includes matching NA positions.

## Portable cohort predicates and optional project guard

Declare cohort restrictions in the registered source as `cohort_gate.sample_ids`
or `cohort_gate.allowed_values` (column name to a nonempty list of allowed values).
A legacy prose-only `cohort_scope` requires an explicit
`project.cohort_scope_gates.<exact_scope>.allowed_values` mapping in the run YAML;
otherwise rendering fails closed. No disease or source-code mapping is inferred.
These predicates add to source-row, hash, priority and applicability checks.

The optional recluster guard invokes the consuming project's
`xenium_downstream.recluster_guard`; it is not bundled here. A project guard
contract auto-enables it when the setting is omitted; disabling a present guard requires a recorded reason and must respect the
consuming project policy. Pure rendering does not recompute clustering.
