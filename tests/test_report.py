import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('report', Path(__file__).resolve().parents[1] / 'bin/report.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

class ReportTests(unittest.TestCase):
    def test_denovo_and_html_escaping(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            fasta=p/'contigs.fa'; fasta.write_text('>a\nACGT\n>b\nAA\n')
            meta=p/'basecaller.json'; meta.write_text('{"selected": "model<a>"}')
            row=report.sample_summary('sample1','denovo',meta,fasta)
            self.assertEqual(row['assembly'], {'sequences':2,'total_bases':6,'longest':4})
            self.assertNotIn('Consensus FASTA',row['files'])
            report.aggregate([row],{'flye_opts':'--meta <test>'},p)
            text=(p/'wf-vzv-report.html').read_text()
            self.assertIn('model&lt;a&gt;',text)
            self.assertIn('--meta &lt;test&gt;',text)
            self.assertEqual(json.loads((p/'results.json').read_text())['samples'][0]['sample'],'sample1')

    def test_reference_optional_mask(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'clair3').mkdir()
            (p/'clair3/s_consensus.fa').write_text('>s\nACNN\n')
            meta=p/'metadata.json'; meta.write_text('{"selected":null}')
            self.assertNotIn('Mask BED',report.sample_summary('s','reference',meta,p)['files'])
            (p/'clair3/low_coverage.bed').write_text('ref\t2\t4\n')
            row=report.sample_summary('s','reference',meta,p)
            self.assertIn('Mask BED',row['files'])
            self.assertNotIn('Contigs FASTA',row['files'])
