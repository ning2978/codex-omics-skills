---
name: omics-environment-setup
description: Inspect, plan, create and verify isolated conda environments for basic Scanpy or R/Seurat analysis on Linux x86_64 and Apple Silicon macOS, including verified Miniforge bootstrap. Use for biomedical software setup or rebuilding an incompatible environment, not for running research analyses or GPU model setup.
---

# Omics environment setup

Use [environment-installation-policy](../environment-installation-policy/SKILL.md)
for conda/pip tool selection. This skill installs software; `install_skills.py`
only copies skills. Read [platforms and operations](references/platforms.md)
for bootstrap prerequisites, CLI options, platform acceptance and recovery.

1. Resolve the native OS/architecture and requested profile (`scrna-python`,
   `seurat-r`, or `both`). Use exact prefixes for existing environments; do not
   infer readiness from an environment name or an activated shell.
2. Run `scripts/setup.py inspect` or `plan` with a new `--report-dir`.
   These operations write reports only, never install or change shell startup.
   A successful inspection is not a runtime PASS. Review each `compatible`
   result, configured package versions, and proposed action.
3. When installation is requested, run `create` with `--install-dir` and a new
   report directory. The script uses existing conda, or downloads a fixed
   Miniforge release verified against a bundled SHA256; `--bootstrap` forces
   an independent bootstrap. No sudo or shell initialization is used.
   It installs from the native platform lock, then verifies with synthetic data.
4. Compatible prefixes are reused after smoke verification. Incompatible ones
   are preserved and replaced by independent prefixes under the installation
   directory; existing destinations get unique siblings. Never upgrade an
   existing analysis environment or remove a failed environment automatically.
5. `verify` checks exact interpreters, required package imports/version ranges,
   and synthetic analysis round trips without installing packages. `both`
   additionally tests the existing Python–R bridge, native DotPlot renderer,
   edgeR and PNG decoding. Single-profile verification cannot certify the bridge.
6. Deliver `report.json`, logs, exact interpreter paths and unresolved checks.
   Inspect generated PNGs before calling visual acceptance complete. Keep outputs
   outside the skills/release tree. Failures return nonzero; keep logs and rebuild
   at a new prefix after correcting the cause.

The platform locks pin transitive packages; profile JSON files document the
intended version series and required packages. A solved lock is not evidence of
installation. See [acceptance status](references/validation.md); native macOS
acceptance remains pending until actually run. Do not describe all 44 skills as
runtime-certified by the basic environments.

For dataset contracts, resources, CUDA and real scVI/cell2location pilots, hand
off to [omics-analysis-runtime-preflight](../omics-analysis-runtime-preflight/SKILL.md)
and the method skill. Monocle2 and other specialized methods need their own
compatible environments. This skill does not provide a complete Seurat analysis
runner, choose biological annotations, or run real research data.
