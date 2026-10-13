"""Maintainer helper: convert an isolated conda --dry-run --json solve to a lock.

Use an empty CONDA_PKGS_DIRS so every LINK record has a corresponding FETCH.
This records a solve only; native installation and smoke acceptance remain separate.
"""
import argparse
import json
from pathlib import Path
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('solve', type=Path); p.add_argument('output', type=Path)
a = p.parse_args()
d = json.loads(a.solve.read_text())
assert d['success'], d
fetch = d['actions']['FETCH']; link = d['actions']['LINK']
assert len(fetch) == len(link), 'Use a fresh package cache; incomplete FETCH records cannot form a lock'
assert {(x['name'],x['version'],x['build']) for x in fetch} == {(x['name'],x['version'],x['build_string']) for x in link}
rows = []
for item in sorted(fetch,key=lambda x:x['name']):
    assert item['url'].startswith('https://conda.anaconda.org/')
    assert len(item['md5']) == 32
    rows.append(item['url'] + '#' + item['md5'])
a.output.write_text('# Generated from an isolated conda solve; see validation status.\n@EXPLICIT\n' + '\n'.join(rows) + '\n')
