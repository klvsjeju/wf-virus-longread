import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import clair3_models as models

BASE = 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0'
CATALOG = ROOT / 'clair3_model_catalog.json'

def checkpoint(path):
    path.mkdir(parents=True)
    for name in ('pileup.index', 'pileup.data-00000-of-00001', 'full_alignment.index', 'full_alignment.data-00000-of-00001'):
        (path / name).write_bytes(b'test weights')
    return path

class ModelTests(unittest.TestCase):
    def test_custom_direct_without_basecaller(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); src=checkpoint(p/'custom')
            info=models.prepare('unknown',CATALOG,p/'output',p/'info.json',src,ROOT/'model_map.tsv')
            self.assertEqual(info['source'],'custom')
            self.assertEqual(len(info['checkpoint_sha256']),4)
            models.validate_checkpoint(p/'output')

    def test_bundled_exact_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); checkpoint(p/'bundled/r1041_e82_400bps_sup_v520')
            with patch.object(models,'download') as download:
                info=models.prepare(BASE,CATALOG,p/'output',p/'info.json',bundled_roots=[p/'bundled'])
            download.assert_not_called()
            self.assertEqual(info['source'],'bundled')

    def test_download_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); src=checkpoint(p/'src')
            def fake_download(url,dest):
                with tarfile.open(dest,'w:gz') as tar:
                    tar.add(src,arcname='r1041_e82_400bps_sup_v520')
                return 'a'*64
            with patch.object(models,'download',side_effect=fake_download):
                info=models.prepare(BASE,CATALOG,p/'output',p/'info.json')
            self.assertEqual(info['source'],'download')
            self.assertEqual(info['archive_sha256'],'a'*64)
            self.assertFalse(info['checksum_verified'])
            models.validate_checkpoint(p/'output')

    def test_unknown_and_empty_checkpoints_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with self.assertRaisesRegex(ValueError,'No exact automatic'):
                models.prepare('unknown',CATALOG,p/'out',p/'info')
            src=checkpoint(p/'empty'); (src/'pileup.index').write_bytes(b'')
            with self.assertRaisesRegex(ValueError,'missing/empty'):
                models.validate_checkpoint(src)

    def test_unsafe_archive_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for name in ('../escaped','/absolute','model/link'):
                with self.subTest(name=name):
                    with tarfile.open(p/'bad.tar.gz','w:gz') as tar:
                        item=tarfile.TarInfo(name)
                        if name.endswith('link'):
                            item.type=tarfile.SYMTYPE; item.linkname='/outside'
                        else:item.size=1
                        tar.addfile(item,io.BytesIO(b'x'))
                    with self.assertRaisesRegex(ValueError,'Unsafe'):
                        models.extract_archive(p/'bad.tar.gz',p/'extract')

    def test_catalog_and_ui(self):
        catalog=json.loads(CATALOG.read_text())
        self.assertEqual(len(catalog['models']),16)
        schema=json.loads((ROOT/'nextflow_schema.json').read_text())
        props={k:v for g in schema['definitions'].values() for k,v in g['properties'].items()}
        self.assertNotIn('model_dir',props)
        self.assertNotIn('model_map',props)
        self.assertIn('override_basecaller_cfg',props)
        for item in schema['allOf']:
            self.assertNotIn('model_dir', item.get('then',{}).get('required',[]))

if __name__ == '__main__':
    unittest.main()
