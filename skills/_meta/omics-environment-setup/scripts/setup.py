#!/usr/bin/env python3
"""Inspect, plan, create or verify isolated CPU omics environments (stdlib CLI)."""
from __future__ import annotations
import argparse
try:
    import fcntl
except ImportError:  # Windows is rejected by detect before any installation.
    fcntl = None
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
import uuid
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE.parent / 'environments'
PROFILES = ('scrna-python', 'seurat-r')


def detect(system=None, machine=None):
    system, machine = system or platform.system(), machine or platform.machine()
    if (system, machine) == ('Linux', 'x86_64'):
        return 'linux-64'
    if (system, machine) == ('Darwin', 'arm64'):
        return 'osx-arm64'
    raise ValueError(f'Unsupported native platform {system}/{machine}. Windows: use x86_64 WSL2 Linux; Intel Mac and emulation are outside this release.')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(url, target, expected):
    with urllib.request.urlopen(url, timeout=120) as src, target.open('wb') as dst:
        shutil.copyfileobj(src, dst)
    if digest(target) != expected:
        raise RuntimeError('SHA256 mismatch; installer was not executed. Preserve download for audit; retry with a new report directory and verify the release checksum.')


def clean_env():
    env = os.environ.copy()
    for key in ('PYTHONPATH', 'PYTHONHOME', 'R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE'):
        env.pop(key, None)
    env.update(R_LIBS_USER=os.devnull, R_LIBS_SITE=os.devnull, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
    return env


def run(cmd, directory, label, env=None):
    log = directory / f'{label}.log'
    with log.open('w') as stream:
        result = subprocess.run([str(x) for x in cmd], env=env or clean_env(), stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{label} failed (exit {result.returncode}); see {log}. Check network/channel availability and package/ABI errors; rebuild at a new prefix, never update the original environment.')
    return log


def executable(prefix, profile):
    return prefix / 'bin' / ('python' if profile == 'scrna-python' else 'Rscript')


def probe(prefix, profile):
    exe = executable(prefix, profile)
    if not exe.is_file():
        return {'compatible': False, 'reason': f'Missing {exe}', 'interpreter': str(exe)}
    try:
        proc = subprocess.run([str(exe), *(['-I', '-B'] if profile == 'scrna-python' else ['--vanilla']),
                               str(HERE / ('probe.py' if profile == 'scrna-python' else 'probe.R'))],
                              capture_output=True, text=True, timeout=120, env=clean_env())
        if proc.returncode:
            raise ValueError(proc.stderr[-4000:])
        data = json.loads(proc.stdout)
        if data.get('machine') not in (platform.machine(), 'aarch64' if platform.machine() == 'arm64' else platform.machine()):
            raise ValueError('Interpreter architecture differs from native host')
        cfg = json.loads((CONFIG / f'{detect()}-{profile}.json').read_text())
        versions = data.get('versions', {})
        for spec in cfg['specs']:
            package, _, version = spec.partition('=')
            observed = versions.get(package)
            if observed is None or (version and observed != version and not observed.startswith(version + '.')):
                raise ValueError(f'{spec}: observed {observed!r}')
        return {'compatible': True, 'interpreter': str(exe), **data}
    except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
        return {'compatible': False, 'interpreter': str(exe), 'reason': str(exc)}


def bootstrap(root, report, subdir):
    prefix = root / 'miniforge'
    exe = prefix / 'bin/conda'
    if exe.is_file():
        return exe
    if prefix.exists():
        raise RuntimeError(f'Incomplete bootstrap at {prefix}; preserve it and choose a new --install-dir.')
    meta = json.loads((CONFIG / 'bootstrap.json').read_text())
    asset = meta[subdir]
    installer = report / asset['asset']
    download(f"https://github.com/conda-forge/miniforge/releases/download/{meta['release']}/{asset['asset']}", installer, asset['sha256'])
    run(['bash', installer, '-b', '-p', prefix], report, 'bootstrap')
    return exe


def verify(prefix, profile, report):
    info = probe(prefix, profile)
    if not info['compatible']:
        raise RuntimeError(f'Incompatible {profile}: {info}')
    exe = executable(prefix, profile)
    work = report / profile
    work.mkdir()
    run([exe, *(['-I', '-B'] if profile == 'scrna-python' else ['--vanilla']),
         HERE / ('smoke.py' if profile == 'scrna-python' else 'smoke.R'), work], report, profile + '-smoke')
    if profile == 'seurat-r':
        run([exe, '--vanilla', HERE.parents[2] / 'single-cell/python/sc-sample-level-comparison/scripts/edger_runtime_smoke.R'], report, 'edger')
        info['edgeR'] = 'PASS'
    info['smoke'] = 'PASS'
    info['conda_inventory'] = [
        {k: record.get(k) for k in ('name', 'version', 'build', 'subdir', 'channel', 'sha256', 'md5')}
        for record in (json.loads(p.read_text()) for p in sorted((prefix / 'conda-meta').glob('*.json')))
    ]
    return info


def joint_checks(prefixes, report):
    skills = HERE.parents[2]
    py = executable(prefixes['scrna-python'], 'scrna-python')
    r = executable(prefixes['seurat-r'], 'seurat-r')
    env = clean_env(); env['SC_BRIDGE_RSCRIPT'] = str(r)
    run([py, '-B', '-m', 'unittest', 'discover', '-s', skills / 'single-cell/shared/h5ad-to-r-bridge/scripts', '-p', 'test_bridge_roundtrip.py', '-v'], report, 'bridge', env)
    run([py, skills / '_shared/dotplot/tests/verify_style.py', '--rscript', r, '--work-dir', report / 'dotplot', '--runtime-versions'], report, 'dotplot', env)
    run([py, '-I', '-B', HERE / 'decode_pngs.py', report], report, 'png-decode', env)
    return {'bridge': 'PASS', 'edgeR': 'PASS', 'native_dotplot': 'PASS', 'png_decode': 'PASS'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['inspect', 'plan', 'create', 'verify'])
    parser.add_argument('--profile', choices=[*PROFILES, 'both'], default='both')
    parser.add_argument('--install-dir', type=Path, default=Path.home() / 'omics-environments')
    parser.add_argument('--report-dir', type=Path, required=True, help='New directory, retained even on failure')
    parser.add_argument('--conda', type=Path, help='Explicit existing conda; otherwise discover PATH, then bootstrap')
    parser.add_argument('--package-cache', type=Path, help='Optional reusable conda cache; protected against concurrent setup use')
    parser.add_argument('--bootstrap', action='store_true', help='Use isolated Miniforge even if PATH contains conda (create only)')
    parser.add_argument('--python-prefix', type=Path, help='Existing environment to inspect/reuse, never upgraded')
    parser.add_argument('--r-prefix', type=Path, help='Existing environment to inspect/reuse, never upgraded')
    args = parser.parse_args(argv)
    report = args.report_dir.expanduser().resolve()
    report.mkdir(parents=True, exist_ok=False)
    result = {'operation': args.operation, 'started_at': datetime.now(timezone.utc).isoformat(),
              'system': platform.system(), 'machine': platform.machine(), 'status': 'FAIL', 'environments': {},
              'bootstrap_release': json.loads((CONFIG / 'bootstrap.json').read_text())}
    try:
        subdir = detect(); result['platform'] = subdir
        if args.operation in ('create', 'verify'):
            if subdir == 'osx-arm64' and int(platform.mac_ver()[0].split('.')[0]) < 14:
                raise ValueError('These locks require macOS 14 or newer')
            if subdir == 'linux-64':
                libc, version = platform.libc_ver()
                if libc != 'glibc' or tuple(map(int, version.split('.')[:2])) < (2, 28):
                    raise ValueError('These locks require glibc 2.28 or newer')
            if any(k.startswith('CONDA_OVERRIDE_') for k in os.environ):
                raise ValueError('Unset CONDA_OVERRIDE_* for native installation/verification')
        root = args.install_dir.expanduser().resolve()
        selected = PROFILES if args.profile == 'both' else (args.profile,)
        prefixes = {}
        for name in selected:
            given = args.python_prefix if name == 'scrna-python' else args.r_prefix
            prefix = (given or root / 'envs' / name).expanduser().resolve()
            cfg = CONFIG / f'{subdir}-{name}.json'
            lock = CONFIG / f'{subdir}-{name}.explicit.txt'
            info = probe(prefix, name)
            info['initial_inspection'] = dict(info)
            info.update(prefix=str(prefix), config=json.loads(cfg.read_text()), config_sha256=digest(cfg),
                        lock_sha256=digest(lock) if lock.exists() else None)
            result['environments'][name] = info
            prefixes[name] = prefix
        if args.operation in ('inspect', 'plan'):
            for name, info in result['environments'].items():
                info['next_action'] = 'verify_then_reuse' if info['compatible'] else 'create_independent_environment'
            result['status'] = 'PASS'  # Successful inspection, not runtime certification.
        else:
            conda = None
            if args.operation == 'create':
                # Do not bootstrap at all if all requested environments can be reused.
                if any(not info['compatible'] for info in result['environments'].values()):
                    root.mkdir(parents=True, exist_ok=True)
                    conda = (bootstrap(root, report, subdir) if args.bootstrap else
                             args.conda or shutil.which('conda') or bootstrap(root, report, subdir))
                    result['conda'] = str(Path(conda).resolve())
            for name in selected:
                info = result['environments'][name]
                if args.operation == 'create' and not info['compatible']:
                    prefix = root / 'envs' / name
                    if prefix.exists():
                        prefix = prefix.with_name(name + '-' + uuid.uuid4().hex[:10])
                    prefixes[name] = prefix
                    info['created_prefix'] = str(prefix)
                    lock = CONFIG / f'{subdir}-{name}.explicit.txt'
                    if not lock.exists():
                        raise RuntimeError(f'Missing platform lock: {lock}; do not silently solve a new stack.')
                    cache = (args.package_cache or root / 'pkgs' / (name + '-' + uuid.uuid4().hex[:10])).expanduser().resolve()
                    cache.mkdir(parents=True, exist_ok=True)
                    env = clean_env(); env['CONDA_PKGS_DIRS'] = str(cache)
                    env['CONDA_SUBDIR'] = subdir
                    env['CONDA_NUMBER_CHANNEL_NOTICES'] = '0'
                    with (cache / '.omics-setup.lock').open('a') as cache_guard:
                        try:
                            fcntl.flock(cache_guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        except BlockingIOError as exc:
                            raise RuntimeError('Package cache is in use by another setup; choose a separate cache') from exc
                        run([conda, 'create', '-y', '--override-channels', '-c', 'conda-forge', '-c', 'bioconda', '-p', prefix, '--file', lock], report, name + '-create', env)
                info.update(verify(prefixes[name], name, report))
                info.pop('reason', None)
                info['prefix'] = str(prefixes[name])
            result['joint_checks'] = joint_checks(prefixes, report) if len(selected) == 2 else {'status': 'NOT_RUN', 'reason': 'Select both profiles for bridge and native renderer acceptance'}
            result['status'] = 'PASS'
    except Exception as exc:
        result['error'] = str(exc)
        print(str(exc), file=sys.stderr)
    finally:
        result['finished_at'] = datetime.now(timezone.utc).isoformat()
        (report / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
