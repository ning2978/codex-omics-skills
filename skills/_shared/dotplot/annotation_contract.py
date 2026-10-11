"""Role-scoped guide/registry comparison. No expression or annotation mutation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re


def semantic(value):
    if isinstance(value, dict):
        return {k: semantic(v) for k, v in value.items() if k != 'updated_at'}
    if isinstance(value, list):
        return [semantic(v) for v in value]
    return value


def differences(a, b, path=''):
    if isinstance(a, dict) and isinstance(b, dict):
        result = []
        for k in sorted(a.keys() | b.keys()):
            p = f'{path}.{k}'
            if k not in a or k not in b:
                result.append(f'{p}: missing from {"guide" if k not in a else "registry"}')
            else:
                result.extend(differences(a[k], b[k], p))
        return result
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f'{path}: guide length={len(a)}, registry length={len(b)}']
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in differences(x, y, f'{path}[{i}]')]
    return [] if type(a) is type(b) and a == b else [f'{path}: guide={a!r}, registry={b!r}']


def resolve(guide: Path, registry: Path, role: str):
    # Read both complete files, return only the requested role.
    text = guide.read_text(encoding='utf-8')
    data = json.loads(registry.read_text(encoding='utf-8'))
    track = data.get('use_cases', {}).get(role, {}).get('default_track')
    if not isinstance(track, dict) or not track.get('track_id'):
        raise ValueError(f'{role}: missing registry default_track')
    def blocks_for(name):
        pattern = (r'<!--\s*' + re.escape(name) + r'_semantics_begin\s*-->\s*```json\s*'
                   r'(.*?)\s*```\s*<!--\s*' + re.escape(name) + r'_semantics_end\s*-->')
        return re.findall(pattern, text, re.S)
    legacy = blocks_for(role)
    complete = blocks_for(role + '_complete')
    blocks = complete or legacy
    if len(blocks) != 1 or len(legacy) > 1:
        raise ValueError(f'{role}: require exactly one complete guide semantic block; found {len(blocks)}')
    delta = differences(semantic(json.loads(blocks[0])), semantic(track), role)
    if complete and legacy:
        # Validate every field of the preserved legacy projection too, without
        # forcing historical exact-projection consumers to migrate.
        def project(actual, template):
            if isinstance(template, dict) and isinstance(actual, dict):
                return {k: project(actual.get(k), v) for k, v in template.items()}
            if isinstance(template, list) and isinstance(actual, list) and len(template) == len(actual):
                return [project(x, y) for x, y in zip(actual, template)]
            return actual
        partial = semantic(json.loads(legacy[0]))
        delta.extend(differences(partial, project(semantic(track), partial), role + '.legacy_summary'))
    if delta:
        raise ValueError('Annotation semantic conflict:\n' + '\n'.join(delta))
    return track


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path.cwd())
    p.add_argument('--role', action='append', required=True)
    p.add_argument('--full', action='store_true')
    a = p.parse_args()
    result = {}
    for role in dict.fromkeys(a.role):
        t = resolve(a.root/'docs/CURRENT_ANNOTATIONS.md', a.root/'metadata/current_annotation_registry.json', role)
        result[role] = t if a.full else {k: t[k] for k in ('track_id', 'status', 'resolution_rule', 'assignment_resolution') if k in t}
        if not a.full:
            result[role]['sources'] = [{k: s[k] for k in ('source_id','priority','path','join','label_columns','exclusion_gates','cohort_gate','fallback') if k in s} for s in t.get('sources', [])]
    print(json.dumps({'semantic_status':'PASS','roles':result}, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
