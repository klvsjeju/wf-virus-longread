import importlib.util
import gzip
import json
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('inputs', Path(__file__).resolve().parents[1] / 'bin/inputs.py')
inputs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inputs)

class InputTests(unittest.TestCase):
    def test_detection_merge_and_override(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            a = p / 'a.fastq'
            b = p / 'b.fq.gz'
            # A model tag after record 1000 must still be detected.
            a.write_text('@read\nAC\n+\nII\n' * 1001 + '@last basecall_model_version_id=model_a\nAC\n+\nII\n')
            with gzip.open(b, 'wt') as f:
                f.write('@other model_version_id=model_b\nAC\n+\nII\n')
            out, meta = p / 'reads.gz', p / 'meta.json'
            inputs.merge([a], out, meta)
            self.assertEqual(json.loads(meta.read_text())['selected'], 'model_a')
            with self.assertRaisesRegex(ValueError, 'Mixed models'):
                inputs.merge([a, b], out, meta)
            inputs.merge([a, b], out, meta, 'manual_model')
            data = json.loads(meta.read_text())
            self.assertEqual(data['selected'], 'manual_model')
            self.assertEqual(data['detected_models'], ['model_a', 'model_b'])
            with gzip.open(out, 'rt') as f:
                self.assertEqual(len(f.readlines()), 1003 * 4)

    def test_missing_tags_and_malformed_fastq(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            a, out, meta = p / 'reads.fq', p / 'out.gz', p / 'meta.json'
            a.write_text('@r\nAC\n+\nII\n')
            inputs.merge([a], out, meta)
            self.assertIsNone(json.loads(meta.read_text())['selected'])
            inputs.merge([a], out, meta, 'manual')
            self.assertEqual(json.loads(meta.read_text())['selected'], 'manual')
            a.write_text('@r\nAC\n+\nI\n')
            with self.assertRaisesRegex(ValueError, 'length mismatch'):
                inputs.merge([a], out, meta)

if __name__ == '__main__':
    unittest.main()
