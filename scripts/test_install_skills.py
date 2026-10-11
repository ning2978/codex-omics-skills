"""Installer behavior in temporary targets; never touches the user's skills."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from install_skills import expand_selection

ROOT = Path(__file__).resolve().parents[1]
CATALOG = json.loads((ROOT / 'catalog.json').read_text())['skills']


class InstallTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='omics-install-')
        self.addCleanup(self.temp.cleanup)
        self.target = Path(self.temp.name) / 'skills'

    def run_install(self, *args, code=0):
        result = subprocess.run([sys.executable, str(ROOT/'scripts/install_skills.py'),
                                 '--target', str(self.target), *args], capture_output=True, text=True)
        self.assertEqual(result.returncode, code, result.stdout + result.stderr)
        return result

    def check_closure(self, selected):
        expanded = expand_selection(CATALOG, selected)
        expected = {r['name'] for r in expanded}
        actual = {p.parent.name for p in self.target.rglob('SKILL.md')}
        self.assertEqual(actual, expected)
        for record in expanded:
            for resource in record['shared_resources']:
                source = ROOT / resource
                dest = self.target / source.relative_to(ROOT/'skills')
                for path in source.rglob('*'):
                    if path.is_file() and '__pycache__' not in path.parts:
                        self.assertEqual(path.read_bytes(), (dest/path.relative_to(source)).read_bytes())
        # All relative Markdown resource links must survive the selected install.
        import re
        for path in self.target.rglob('*.md'):
            for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)', path.read_text()):
                if not link.startswith(('http:', 'https:', '#', 'mailto:')):
                    self.assertTrue((path.parent/link.split('#')[0]).exists(), f'{path}: {link}')

    def test_full_repeat_and_force_preserve_unrelated_files(self):
        self.run_install()
        self.check_closure(CATALOG)
        sentinel = self.target/'unrelated.txt'; sentinel.write_text('preserve')
        skill = self.target/'render-scrna-dotplot/SKILL.md'
        skill.write_text('local edit')
        self.run_install(code=1)
        self.assertEqual(skill.read_text(), 'local edit')
        self.run_install('--force')
        self.assertEqual(skill.read_bytes(), (ROOT/'skills/render-scrna-dotplot/SKILL.md').read_bytes())
        self.assertEqual(sentinel.read_text(), 'preserve')

    def test_single_dotplot_has_engine_without_other_skills(self):
        for modality in ['scrna', 'xenium']:
            self.target = Path(self.temp.name)/modality
            name = f'render-{modality}-dotplot'
            self.run_install('--name', name)
            self.check_closure([r for r in CATALOG if r['name']==name])
            self.assertTrue((self.target/'_shared/dotplot/seurat_dotplot_renderer.R').is_file())

    def test_category_includes_external_resource_owners(self):
        self.run_install('--category', 'single-cell/python')
        self.check_closure([r for r in CATALOG if r['category']=='single-cell/python'])
        self.assertTrue((self.target/'single-cell/workflow/workflow-sc-analysis/scripts/sc_contracts.py').is_file())

    def test_spatial_category(self):
        self.run_install('--category', 'spatial')
        self.check_closure([r for r in CATALOG if r['category'].startswith('spatial/')])
        self.assertTrue((self.target/'_shared/spatial/run_contract.py').is_file())

    def test_dry_run_and_invalid_selection_do_not_write(self):
        self.run_install('--dry-run', '--name', 'render-scrna-dotplot')
        self.assertFalse(self.target.exists())
        self.run_install('--name', 'unknown', code=2)
        self.run_install('--category', 'unknown', code=2)
        self.assertFalse(self.target.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
