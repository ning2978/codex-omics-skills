"""Deployment control-flow regression tests; no network or real installation."""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

PATH = Path(__file__).resolve().parents[1] / 'skills/_meta/omics-environment-setup/scripts/setup.py'
spec = importlib.util.spec_from_file_location('omics_setup', PATH)
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)


class SetupTests(unittest.TestCase):
    def call(self, args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return s.main(args)

    def test_profile_locks_are_complete_native_and_pinned(self):
        for subdir in ['linux-64','osx-arm64']:
            for profile in s.PROFILES:
                cfg=json.loads((s.CONFIG/f'{subdir}-{profile}.json').read_text())
                self.assertEqual(cfg['platform'],subdir)
                validation=json.loads((s.CONFIG/'validation.json').read_text())
                record=next(x for x in validation['platforms'] if x['platform']==subdir and x['profile']==profile)
                self.assertEqual(record['config_sha256'],s.digest(s.CONFIG/f'{subdir}-{profile}.json'))
                self.assertEqual(record['lock_sha256'],s.digest(s.CONFIG/f'{subdir}-{profile}.explicit.txt'))
                lines=(s.CONFIG/f'{subdir}-{profile}.explicit.txt').read_text().splitlines()
                self.assertIn('@EXPLICIT',lines)
                urls=[x for x in lines if x.startswith('https://')]
                self.assertGreater(len(urls),50)
                for url in urls:
                    self.assertTrue(f'/{subdir}/' in url or '/noarch/' in url)
                    self.assertRegex(url,r'#[a-f0-9]{32}$')
                for package in cfg['specs']:
                    name=package.split('=')[0]
                    # python-igraph is a conda metapackage and remains explicitly locked.
                    self.assertTrue(any('/'+name+'-' in url for url in urls),name)

    def test_platforms(self):
        self.assertEqual(s.detect('Linux', 'x86_64'), 'linux-64')
        self.assertEqual(s.detect('Darwin', 'arm64'), 'osx-arm64')
        for pair in [('Windows','AMD64'),('Darwin','x86_64'),('Linux','aarch64')]:
            with self.assertRaises(ValueError): s.detect(*pair)

    def test_checksum_mismatch_never_executes(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(s.urllib.request,'urlopen', return_value=io.BytesIO(b'corrupt')), patch.object(s,'run') as run:
            with self.assertRaisesRegex(RuntimeError,'SHA256'):
                s.download('https://example.invalid/file',Path(tmp)/'installer',hashlib.sha256(b'good').hexdigest())
            run.assert_not_called()
            self.assertEqual((Path(tmp)/'installer').read_bytes(),b'corrupt')

    def test_inspect_plan_only_write_reports(self):
        with tempfile.TemporaryDirectory(prefix='omics space ') as tmp:
            root = Path(tmp)
            for op in ['inspect','plan']:
                with patch.object(s,'bootstrap') as boot, patch.object(s,'run') as run:
                    self.assertEqual(self.call([op,'--install-dir',str(root/'never created'),'--report-dir',str(root/op)]),0)
                    boot.assert_not_called(); run.assert_not_called()
                self.assertFalse((root/'never created').exists())
                d=json.loads((root/op/'report.json').read_text())
                self.assertFalse(d['environments']['scrna-python']['compatible'])
                self.assertNotIn('joint_checks',d)

    def test_existing_conflict_preserved_and_failure_reported(self):
        with tempfile.TemporaryDirectory(prefix='omics space ') as tmp:
            root=Path(tmp); old=root/'envs/scrna-python'; old.mkdir(parents=True)
            sentinel=old/'keep'; sentinel.write_text('original')
            with patch.object(s,'probe',return_value={'compatible':False}), patch.object(s,'run',side_effect=RuntimeError('installation failed')) as run:
                self.assertEqual(self.call(['create','--profile','scrna-python','--conda','/fake conda/bin/conda','--install-dir',str(root),'--report-dir',str(root/'report')]),1)
            cmd=run.call_args.args[0]
            self.assertEqual(str(cmd[0]),'/fake conda/bin/conda')
            created=Path(cmd[cmd.index('-p')+1]); self.assertNotEqual(created,old)
            self.assertEqual(sentinel.read_text(),'original')
            d=json.loads((root/'report/report.json').read_text())
            self.assertEqual(d['status'],'FAIL'); self.assertIn('installation failed',d['error'])

    def test_shared_cache_busy_fails_before_install(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(s,'probe',return_value={'compatible':False}), patch.object(s,'run') as run:
            root=Path(tmp); cache=root/'cache'; cache.mkdir()
            with (cache/'.omics-setup.lock').open('a') as guard:
                s.fcntl.flock(guard,s.fcntl.LOCK_EX|s.fcntl.LOCK_NB)
                code=self.call(['create','--profile','scrna-python','--conda','/fake/conda','--install-dir',str(root/'envroot'),'--package-cache',str(cache),'--report-dir',str(root/'report')])
                self.assertEqual(code,1); run.assert_not_called()
                self.assertIn('cache is in use',json.loads((root/'report/report.json').read_text())['error'])

    def test_compatible_reuse_never_installs(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(s,'probe',return_value={'compatible':True}), patch.object(s,'verify',return_value={'compatible':True,'smoke':'PASS'}), patch.object(s,'bootstrap') as boot, patch.object(s,'run') as run:
            root=Path(tmp)
            self.assertEqual(self.call(['create','--profile','scrna-python','--install-dir',str(root/'not created'),'--python-prefix',str(root/'existing'),'--report-dir',str(root/'report')]),0)
            boot.assert_not_called(); run.assert_not_called()
            self.assertFalse((root/'not created').exists())

    def test_verify_missing_environment_fails_without_install(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(s,'bootstrap') as boot, patch.object(s,'run') as run:
            root=Path(tmp)
            self.assertEqual(self.call(['verify','--install-dir',str(root/'missing'),'--report-dir',str(root/'report')]),1)
            boot.assert_not_called(); run.assert_not_called()

    def test_bootstrap_existing_incomplete_preserved(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(s,'download') as download:
            root=Path(tmp); (root/'miniforge').mkdir()
            with self.assertRaisesRegex(RuntimeError,'Incomplete'):
                s.bootstrap(root,root,'linux-64')
            download.assert_not_called()

    def test_run_failure_log_and_argument_boundaries(self):
        with tempfile.TemporaryDirectory(prefix='space test ') as tmp:
            root=Path(tmp); script=root/'fail script.py'
            script.write_text('import sys; print(sys.argv[1]); sys.exit(17)')
            with self.assertRaisesRegex(RuntimeError,'exit 17'):
                s.run([s.sys.executable,script,'a b $literal'],root,'failure')
            self.assertEqual((root/'failure.log').read_text().strip(),'a b $literal')


if __name__ == '__main__': unittest.main()
