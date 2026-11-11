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


class HeaderModelTests(unittest.TestCase):
    MODEL = 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0'

    def test_sam_read_group_formats(self):
        for suffix in ('', '_barcode01', '_barcode1234', '_5mCG_5hmCG@v2_barcode01'):
            for tag in ('RG:Z:', 'RG='):
                header = f'@read\t{tag}abc123_{self.MODEL}{suffix}\tqs:f:20\n'
                with self.subTest(tag=tag, suffix=suffix):
                    self.assertEqual(inputs.extract_models(header), ({self.MODEL}, 'read_group'))

    def test_explicit_model_has_precedence(self):
        old = 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0'
        header = f'@read basecall_model_version_id={self.MODEL} RG:Z:run_{old}_barcode01'
        self.assertEqual(inputs.extract_models(header), ({self.MODEL}, 'explicit_model_tag'))

    def test_no_guessing_and_unknown_tags(self):
        for header in (f'@{self.MODEL} barcode=barcode01', '@read RG:Z:run_unknown_barcode01', '@read model_version_id=unknown', f'@read note={self.MODEL}'):
            self.assertEqual(inputs.extract_models(header)[0], set())
        self.assertEqual(inputs.extract_models(f'@read model_version_id=unknown RG:Z:run_{self.MODEL}')[0], {self.MODEL})

    def test_rg_merge_plain_gzip_and_mixed_detection(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); a=p/'a.fastq'; b=p/'b.fastq.gz'; out=p/'out.gz'; meta=p/'meta.json'
            a.write_text(f'@a RG:Z:run_{self.MODEL}_barcode01\nAC\n+\nII\n')
            with gzip.open(b,'wt') as f:
                f.write(f'@b basecall_model_version_id={self.MODEL}\nAC\n+\nII\n')
            inputs.merge([a,b],out,meta)
            data=json.loads(meta.read_text())
            self.assertEqual(data['selected'],self.MODEL)
            self.assertEqual(data['model_detection_sources'],{'read_group':1,'explicit_model_tag':1})
            self.assertEqual(data['records_with_model'],2)
            a.write_text('@a RG:Z:run_dna_r10.4.1_e8.2_400bps_hac@v5.0.0_barcode01\nAC\n+\nII\n')
            with self.assertRaisesRegex(ValueError,'Mixed models'):
                inputs.merge([a,b],out,meta)


class FastcatMetadataTests(unittest.TestCase):
    MODEL = 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0'

    def write_stats(self, directory, model, count):
        stats = directory / 'stats'
        stats.mkdir()
        (stats / 'per-file-stats.tsv').write_text('filename\tn_seqs\nreads.fq\t1\n')
        (stats / 'per-file-basecallers.tsv').write_text(
            f'filename\tbasecaller\tcount\nreads.fq\t{model}\t{count}\n')
        return stats

    def test_complete_fastcat_statistics_skip_header_scan(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats = self.write_stats(root, self.MODEL, 1)
            # A complete report avoids reopening/recompressing original reads.
            inputs.fastcat_metadata(stats, root / 'basecaller.json', [root / 'not_reopened.fq'])
            result = json.loads((root / 'basecaller.json').read_text())
            self.assertEqual(result['selected'], self.MODEL)
            self.assertEqual(result['records_scanned'], 1)

    def test_missing_fastcat_model_uses_previous_rg_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats = self.write_stats(root, '', 1)
            reads = root / 'reads.fq'
            reads.write_text(f'@read RG:Z:run_{self.MODEL}_barcode01\nAC\n+\nII\n')
            inputs.fastcat_metadata(stats, root / 'basecaller.json', [reads])
            result = json.loads((root / 'basecaller.json').read_text())
            self.assertEqual(result['selected'], self.MODEL)
            self.assertEqual(result['model_detection_sources']['read_group'], 1)

    def test_missing_or_empty_inputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats = self.write_stats(root, self.MODEL, 1)
            (stats / 'per-file-stats.tsv').write_text('filename\tn_seqs\nreads.fq\t0\n')
            with self.assertRaisesRegex(ValueError, 'every input'):
                inputs.fastcat_metadata(stats, root / 'basecaller.json', [root / 'empty.fq'])

    def test_filtered_empty_chunk_and_retained_header_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats = self.write_stats(root, '', 1)
            (stats / 'per-file-stats.tsv').write_text('filename\tn_seqs\na.fq\t0\nb.fq\t1\n')
            retained = root / 'retained.fastq.gz'
            with gzip.open(retained, 'wt') as handle:
                handle.write(f'@read RG:Z:run_{self.MODEL}_barcode01\nACGT\n+\nIIII\n')
            # Original chunks need not be reopened; filtered-out models cannot leak in.
            inputs.fastcat_metadata(stats, root / 'basecaller.json',
                                   [root / 'a.fq', root / 'b.fq'], filtered_reads=retained)
            result = json.loads((root / 'basecaller.json').read_text())
            self.assertEqual(result['selected'], self.MODEL)
            self.assertEqual(result['records_scanned'], 1)

if __name__ == '__main__':
    unittest.main()
