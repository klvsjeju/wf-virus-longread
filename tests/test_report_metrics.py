"""Check report statistics against small, independently calculable inputs."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import report


class MetricsTests(unittest.TestCase):
    def test_depth_zero_positions_contigs_and_bounded_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'depth'
            path.write_text('ref\t1\t0\nref\t2\t30\nref\t3\t50\nother\t1\t100\n')
            rows = report.coverage_summary(path, 40, max_points=2)
            self.assertEqual(rows[0]['positions'], 3)
            self.assertAlmostEqual(rows[0]['mean_depth'], 80 / 3)
            self.assertAlmostEqual(rows[0]['covered_percent'], 200 / 3)
            self.assertAlmostEqual(rows[0]['at_50x_percent'], 100 / 3)
            self.assertEqual(rows[0]['below_mask'], 2)
            self.assertEqual(rows[0]['plot'], [{'position': 1.5, 'depth': 15.0}, {'position': 3.0, 'depth': 50.0}])
            self.assertEqual(rows[1]['reference'], 'other')
            self.assertEqual(rows[1]['mean_depth'], 100)

    def test_histogram_weighting_and_no_reads(self):
        stats = report.read_summary({'length': [[100, 101, 1], [200, 201, 3]]}, {})
        self.assertEqual(stats['read_count'], 4)
        self.assertEqual(stats['median_length'], 200)
        self.assertEqual(stats['total_bases'], 700)
        self.assertIsNone(report.read_summary({}, {})['median_length'])

    def test_input_and_retained_read_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            stats = Path(directory)
            (stats / 'per-file-stats.tsv').write_text('n_seqs\tr_record_ok\tr_too_short\n2\t4\t2\n0\t1\t1\n')
            result = report.read_summary({'length': [[1000,1001,2]]}, {}, stats)
            self.assertEqual(result['input_read_count'], 5)
            self.assertEqual(result['read_count'], 2)
            self.assertEqual(result['excluded_short_reads'], 3)

    def test_bcftools_stats_and_annotation_limits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'variants.stats'
            path.write_text('# header\nSN\t0\tnumber of SNPs:\t3\nTSTV\t0\t2\t1\t2.0\nST\t0\tT>C\t2\nIDD\t0\t-2\t1\n')
            data = report.variant_summary(path)
            self.assertEqual(data['counts']['number of SNPs'], 3)
            self.assertEqual(data['ts_tv'], 2)
            self.assertEqual(data['substitutions'], {'T>C': 2})
            self.assertEqual(data['indels'], [{'length': -2, 'count': 1}])
            path.write_text('TSTV\t0\t0\t0\t0\n')
            self.assertIsNone(report.variant_summary(path)['ts_tv'])
            path = root / 'var.tsv'
            path.write_text('POS\tANN[*].IMPACT\tANN[*].EFFECT\n1\tHIGH;HIGH\tstop_gained\n2\tLOW\tsynonymous_variant\n')
            data = report.annotations_summary(path, limit=1)
            self.assertEqual(data['records'], 2)
            self.assertEqual(len(data['preview']), 1)
            self.assertEqual(data['impacts'], {'HIGH': 1, 'LOW': 1})

    @unittest.skipUnless(importlib.util.find_spec('ezcharts'), 'Requires the wf-common container')
    def test_labs_reference_render_and_all_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'clair3').mkdir(); (root / 'alignment').mkdir(); (root / 'annot').mkdir()
            (root / 'clair3/s_consensus.fa').write_text('>s\nACNN\n')
            (root / 'clair3/variants.stats').write_text('SN\t0\tnumber of records:\t2\nTSTV\t0\t1\t1\t1\nST\t0\tA>G\t1\nST\t0\tC>A\t1\nIDD\t0\t-1\t1\nIDD\t0\t2\t1\n')
            (root / 'alignment/s.depth').write_text('ref\t1\t0\nref\t2\t50\n')
            (root / 'annot/var.tsv').write_text('POS\tANN[*].IMPACT\tANN[*].EFFECT\n2\tHIGH\tstop_gained\n')
            meta = root / 'meta.json'; meta.write_text('{"selected":"model<a>"}')
            stats = root / 'stats'; stats.mkdir()
            (stats / 'length.hist').write_text('100\t101\t2\n200\t201\t1\n')
            (stats / 'quality.hist').write_text('10\t10.1\t2\n40\t40.1\t1\n')
            sample = report.sample_summary('s', 'reference', meta, root, stats)
            report.aggregate([sample], {'run_annotation': True, 'workflow_version': '0.3.3'}, root)
            text = (root / 'wf-virus-longread-report.html').read_text()
            for marker in ['Sample summary', 'Variant calling', 'Genome coverage', 'Variant annotations', 'model&lt;a&gt;', 'Software versions', 'Workflow parameters', 'ezcharts']:
                self.assertTrue(marker in text, marker)
            self.assertTrue('Custom wf-virus-longread workflow maintained by klvsjeju.' in text)
            self.assertFalse('workflow provided by Oxford Nanopore Technologies.' in text)
            report.aggregate([], {'missing_samples': [{'barcode':'barcode03', 'alias':'Missing_Patient'}]}, root)
            self.assertIn('Missing_Patient', (root / 'wf-virus-longread-report.html').read_text())
            self.assertEqual(json.loads((root / 'results.json').read_text())['samples'], [])
