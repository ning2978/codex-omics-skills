import json
from pathlib import Path
import tempfile
import unittest
from audit_bio_skills import discover, semantic_audit

ROOT=Path(__file__).resolve().parents[3]


class RegistryAudit(unittest.TestCase):
    def setUp(self):
        self.skills=discover(ROOT,None)

    def test_ordinary_project_explicitly_not_applicable(self):
        with tempfile.TemporaryDirectory() as tmp:
            errors,notes=semantic_audit(self.skills,Path(tmp),ROOT)
            self.assertEqual(errors,[])
            self.assertTrue(any("annotation_registry: NOT_APPLICABLE" in v for v in notes))

    def test_existing_registry_validates_guide_and_assets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/"metadata").mkdir(); (root/"docs").mkdir()
            track={"track_id":"fixture_v1","sources":[{"source_id":"counts","path":"counts.h5ad"}]}
            registry=root/"metadata/current_annotation_registry.json"
            registry.write_text(json.dumps({"use_cases":{"scrna":{"default_track":track}}}))
            guide=root/"docs/CURRENT_ANNOTATIONS.md"
            guide.write_text('<!-- scrna_semantics_begin -->\n```json\n'+json.dumps(track)+'\n```\n<!-- scrna_semantics_end -->\n')
            errors,_=semantic_audit(self.skills,root,ROOT)
            self.assertTrue(any("missing counts.h5ad" in v for v in errors))
            (root/"counts.h5ad").write_bytes(b"asset existence fixture")
            self.assertEqual(semantic_audit(self.skills,root,ROOT)[0],[])
            guide.write_text(guide.read_text().replace("fixture_v1","different"))
            errors,_=semantic_audit(self.skills,root,ROOT)
            self.assertTrue(any("semantic conflict" in v for v in errors))

    def test_malformed_existing_registry_fails_instead_of_skip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/"metadata").mkdir()
            registry=root/"metadata/current_annotation_registry.json"
            for payload in ('not json','[]','{"use_cases":{}}'):
                registry.write_text(payload)
                self.assertTrue(semantic_audit(self.skills,root,ROOT)[0])


if __name__ == "__main__": unittest.main(verbosity=2)
