# Platforms and operations

Run the CLI with Python 3.10+ (standard library only). Bootstrap also needs bash,
HTTPS connectivity and the usual OS C runtime; it does not require conda already
installed. On macOS use native arm64 Python/Terminal, not Rosetta. Linux needs
x86_64 and a current glibc distribution (glibc 2.28+ required, Ubuntu 22.04+ recommended); macOS 14+
is required for these conda-forge locks. Package virtual requirements may
reject older systems. Do not override virtual packages to bypass compatibility.

Windows: enable WSL2 using Microsoft's installation instructions, install an
x86_64 Ubuntu distribution, and run this workflow inside Linux, storing prefixes
in its Linux filesystem. Windows-native, Windows ARM, Intel Mac, Linux ARM,
CUDA/GPU models and specialized methods are not automated in this version.

```bash
python3 skills/_meta/omics-environment-setup/scripts/setup.py inspect \
  --profile both --install-dir "$HOME/omics-environments" \
  --report-dir "$HOME/omics-reports/inspect-01"
python3 skills/_meta/omics-environment-setup/scripts/setup.py plan \
  --profile both --report-dir "$HOME/omics-reports/plan-01"
python3 skills/_meta/omics-environment-setup/scripts/setup.py create \
  --profile both --bootstrap --install-dir "$HOME/omics-environments" \
  --report-dir "$HOME/omics-reports/create-01"
python3 skills/_meta/omics-environment-setup/scripts/setup.py verify \
  --python-prefix "$HOME/omics-environments/envs/scrna-python" \
  --r-prefix "$HOME/omics-environments/envs/seurat-r" \
  --report-dir "$HOME/omics-reports/verify-01"
```

All paths are passed as argument arrays, not shell interpolation. Quote paths in
shell invocations. Some Miniforge installers reject spaces in their own prefix;
if bootstrap fails on such a path, retain the log and use a bootstrap installation
directory without spaces. Environment prefixes and report paths accept spaces.
Use `--conda /path/to/bin/conda` to choose an existing installation explicitly;
`--bootstrap` overrides discovery. No `conda init`, sudo or shell edits occur.
Do not activate the new environments implicitly: use reported `bin/python` or
`bin/Rscript`. An `inspect` PASS means the inspection completed, not that the
candidate environments passed smoke tests. Supply the two prefix options to
inspect an existing installation outside the default install directory.

By default each creation uses its own package cache to avoid download collisions.
For a rebuild, `--package-cache /path/to/completed-cache` can reuse downloaded
packages while creating a new environment. The script takes an exclusive advisory
lock and refuses simultaneous setup use. Do not run unrelated conda operations
against that cache at the same time. Cache reuse never updates the old environment.

Use a new report directory for each operation. Logs, synthetic artifacts and
failure metadata remain there. Creation records requested versus actual prefixes
and configuration/lock hashes; version probes record the active package versions.
A failed installation is never reused without passing probes and smoke tests.
If a compatible existing environment fails its smoke test, stop and retain it;
rebuild using a new installation directory without passing that candidate prefix.
Network failure: inspect the preserved log, correct proxy/certificates/connectivity,
then retry in new directories. SHA256 failure: never execute the download.
Dependency/ABI failure: preserve the environment, rebuild from locks, and route
method-specific incompatibility to runtime-preflight instead of upgrading in place.

## Lock maintenance

Each platform has a profile JSON and an explicit conda package URL/build/checksum
lock. These are curated fresh solves, not exports of a personal environment.
For cross-platform lock generation only, Linux maintainers may set
`CONDA_OVERRIDE_OSX=14.0` to declare the macOS solve target; this does not run
macOS software. Unset all `CONDA_OVERRIDE_*` variables before native creation
or verification (the CLI rejects them).
Maintainers use `conda create --dry-run --json --platform PLATFORM` with the
profile's channels/specifications and an empty `CONDA_PKGS_DIRS`; convert using
`scripts/lock_from_solve.py SOLVE_JSON OUTPUT.explicit.txt`. Do not label the
result validated until native `create` and `verify` pass. Re-run acceptance when
locks, probes or smoke tests change. Miniforge release/hash pins are separate in
`environments/bootstrap.json`. Retain failed solve/install logs outside releases.

The manually dispatched `environment-acceptance` GitHub workflow uses Ubuntu and
native Apple Silicon runners. It checks architecture before creation, runs both
profiles, and uploads reports on success or failure. CI execution and publication
are separate user-authorized release actions; adding the workflow does not run it.
