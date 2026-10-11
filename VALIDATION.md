# Publication validation — 2026-10-11

Validation used synthetic fixtures in temporary directories outside the release
tree. No patient data, project registry or analysis result is published. The
local operational skills were read as sources and were not edited.

| Check | Result |
|---|---|
| Catalog, taxonomy, relative resources, portability and Python syntax | PASS, 43 skills |
| Skill Creator frontmatter validation | PASS, 43 skills |
| R parse validation | PASS, 5 scripts |
| Temporary installer tests | PASS, 5 tests: full, two individual DotPlots, categories, repeats/force and dry-run/error behavior |
| DotPlot contracts | PASS, 23 tests plus both entrypoint self-tests |
| Taxonomy registry semantics | PASS, 3 tests |
| Portable spatial tools | PASS, 15 tests |
| Single-cell contracts and retained-count/layout execution | PASS, 8 tests |
| AnnData/Matrix Market/Seurat round trip | PASS, 1 test |
| Multi-library AIRR/Scirpy integration | PASS, 1 test |
| edgeR paired quasi-likelihood smoke | PASS, 200 synthetic genes, residual df 3 |
| Native Seurat DotPlot integration | PASS, both modalities, single/multiple groups, constant genes, automatic dimensions and explicit overrides |
| Previous local vs published native DotPlot statistics | PASS, unchanged native numeric tables |
| Installed-copy rendering without an explicit R renderer path | PASS, both modalities; shared renderer resolved from installed tree |
| PNG decoding, white background and PNG-only packages | PASS |
| Visual review | PASS, 15 installed-copy PNGs reviewed in contact sheets, representative PNGs inspected at full size |

The unit suites total 56 tests. DotPlot checks additionally render and validate
complete packages, compare independent Python values to unmodified Seurat
`DotPlot$data`, and require matching missing-value positions. Single-group
z-scores remain missing; expression is native log1p mean. Explanations are in
layout/numeric records and documentation; no explanatory bottom caption is drawn.
Axes, program strips, labels and both native legends remain visible. Explicit
larger test dimensions are respected rather than silently compacted.

The tested R stack was R 4.4.3, Seurat 5.5.1, SeuratObject 5.4.0 and edgeR 4.4.0.
Python tests used Python 3.11.14, NumPy 2.3.5, pandas 2.3.3, SciPy 1.16.3,
AnnData 0.12.7, Scanpy 1.11.5 and Scirpy 0.22.4 in the existing analysis
environment; no packages were installed.
The system Python lacked PyYAML for Skill Creator's validator, so that check was
rerun successfully with the existing analysis interpreter. The repository's
standalone packaging checks require only the Python standard library.

## Scope of portability adaptations

- Required shared resources are installed transitively alongside their consumers.
- R rendering resolves a regular shared file instead of local symbolic links.
- Existing CLI options and skip/nonzero/explicit-force behavior are retained.
- Private paths and historical cohort examples are removed; missing project
  scripts are described as configurable project runners, not bundled tools.
- Source/registry/priority/join/gate validation remains active. Cohort predicates
  are explicit and fail closed when a prose-only scope lacks a mapping.
- Only publication copies are adapted. Historical Git history is not rewritten.

## Runtime work not executed

No full biological dataset analysis or GPU training was run. In particular,
scVI and cell2location CUDA smoke scripts, full COMMOT transport fitting,
resegmentation model fitting, Tangram, CellCharter model fitting, RCTD fitting,
Monocle2 trajectory, pySCENIC, velocity/CellRank and external project renderers
were not runtime-validated for this release. Their Python/R syntax or instruction
contracts were checked where applicable. Spatial preparation/geometry/null
helpers and synthetic expression/layout behavior were exercised as listed above;
this does not certify every external analysis environment or biological claim.

Commands to reproduce the bundled checks are in [README.md](README.md).
