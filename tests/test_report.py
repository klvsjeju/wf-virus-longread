import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
HAS_EZCHARTS = importlib.util.find_spec('ezcharts') is not None

spec = importlib.util.spec_from_file_location('report', Path(__file__).resolve().parents[1] / 'bin/report.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

class ReportTests(unittest.TestCase):
    @unittest.skipUnless(HAS_EZCHARTS, 'HTML rendering is tested in the wf-common container')
    def test_denovo_and_html_escaping(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            fasta=p/'contigs.fa'; fasta.write_text('>a\nACGT\n>b\nAA\n')
            meta=p/'basecaller.json'; meta.write_text('{"selected": "model<a>"}')
            row=report.sample_summary('sample1','denovo',meta,fasta)
            self.assertEqual(row['assembly'], {'sequences':2,'total_bases':6,'longest':4})
            self.assertNotIn('Consensus FASTA',row['files'])
            report.aggregate([row],{'flye_opts':'--meta <test>', 'effective_min_read_length':1000},p)
            text=(p/'wf-virus-longread-report.html').read_text()
            self.assertTrue('model&lt;a&gt;' in text)
            self.assertTrue('--meta &lt;test&gt;' in text)
            self.assertTrue('De novo minimum read length: 1000 bp.' in text)
            self.assertTrue('Fastcat histograms and retained read counts describe reads that passed the filter.' in text)
            self.assertEqual(json.loads((p/'results.json').read_text())['samples'][0]['sample'],'sample1')

    def test_reference_optional_mask(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'clair3').mkdir()
            (p/'clair3/s_consensus.fa').write_text('>s\nACNN\n')
            meta=p/'metadata.json'; meta.write_text('{"selected":null}')
            self.assertNotIn('Mask BED',report.sample_summary('s','reference',meta,p)['files'])
            (p/'clair3/low_coverage.bed').write_text('ref\t2\t4\n')
            row=report.sample_summary('s','reference',meta,p)
            self.assertEqual(row['files']['Variants'], 's/clair3/s_norm.pass.vcf.gz')
            self.assertEqual(row['files']['Consensus FASTA'], 's/clair3/s_consensus.fa')
            self.assertEqual(row['files']['QC report'], 's/qc/NanoPlot/s_NanoPlot-report.html')
            self.assertEqual(row['files']['Basecalling metadata'], 's/metadata/basecaller.json')
            self.assertEqual(row['variant_statistics_source'], 's_norm.vcf.gz')
            self.assertIn('Mask BED',row['files'])
            self.assertNotIn('Contigs FASTA',row['files'])
