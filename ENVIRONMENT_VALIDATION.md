# Environment setup release validation — 2026-10-13

This additive record covers the 44-skill release. Historical
[VALIDATION.md](VALIDATION.md) remains unchanged. Only synthetic data was used;
original analysis environments were not upgraded or overwritten.

| Check | Result |
|---|---|
| Skill format, catalog, taxonomy, portability, relative resources, syntax | PASS, 44 skills; `_meta` = 4 |
| Deployment and installer regression tests | PASS, 16 tests |
| Miniforge fixed release 26.7.2-0, SHA256, isolated bootstrap | PASS on Linux x86_64 |
| Python fresh environment and full synthetic Scanpy/Harmony/H5AD smoke | PASS |
| R fresh environment and Seurat clustering/UMAP/RDS/PNG smoke | PASS |
| Existing Python–R keyed sparse bridge and reduction round trip | PASS |
| edgeR paired QL fit | PASS, 200 genes and residual df 3 |
| Native scRNA/Xenium DotPlot, numeric parity, single groups, constant genes, sizing | PASS |
| PNG decoding | PASS, 17 outputs per final acceptance run |
| Visual inspection | PASS, 15 native DotPlots plus both basic smoke PNGs |
| Installed local skill using preserved local dependency tree | PASS, full `verify --profile both` |
| macOS arm64 locks | PASS for target macOS 14 cross-platform solve only |
| Native macOS installation and smoke | **PENDING — not run** |
| Remote CI / GitHub push | NOT RUN; manual dispatch workflow provided |

The 16 regression tests cover platform routing, checksum rejection, existing
prefix preservation, failed installation/reporting, path argument boundaries,
inspect/plan without installation, compatible reuse, missing-environment verify,
incomplete bootstrap, concurrent cache rejection, lock coverage and dependency
closure, alongside the existing installer cases. A real `plan` call with spaces
in both output paths left the installation directory absent.

Exact core versions and lock/config hashes are recorded in
[the portable acceptance manifest](skills/_meta/omics-environment-setup/environments/validation.json).
Python: 3.11.17, Scanpy 1.11.5, AnnData 0.12.7, Harmony 0.0.10,
Leidenalg 0.10.2, UMAP 0.5.12. R: 4.4.3, Seurat 5.3.1,
SeuratObject 5.2.0, Matrix 1.7.6, edgeR 4.4.0.
The scikit-misc conda build is 0.5.3 but its Python metadata reports 0.0.0;
both are retained in runtime reports, and an actual LOESS fit passed.

## Issues found and resolved

- Concurrent creation initially shared a package cache and collided on downloads.
  Default creation now isolates caches; explicit cache reuse takes an exclusive
  advisory lock. Initial failed reports/logs remain available outside the release.
- The real sparse bridge required R.utils to read compressed TSV with data.table.
  Both R profiles/locks now include it. Verification used a newly created R
  prefix, retaining the earlier environment.
- Existing host channel notifications caused unnecessary network retries.
  Creation now names the configured channels and disables those notifications.
- The shared DotPlot synthetic acceptance script accepts `--runtime-versions`,
  reading the exact R interpreter versions into its run contract. Its original
  default version contract and optional historical baseline mode are preserved.

Local installation added the new skill, minimally extended the local taxonomy,
and updated its shared synthetic DotPlot test helper with backups. Other local
skills, project adapters, historical outputs and research data were preserved.
The published `install_skills.py` still only copies skill files.

Raw evidence is retained in an external acceptance directory: final
`create-both-final/report.json`, `verify-installed/report.json`, corresponding
logs, PNGs and machine-readable numeric records. No local account paths,
installers, environment directories, synthetic data or logs are bundled.

No GPU model, specialized method, complete Seurat analysis runner or research
analysis is certified by this release. macOS remains pending until the native
arm64 manual CI job (which asserts actual architecture) is executed successfully.
