"""Exercise real report argument parsing without requiring HTML dependencies."""
import base64
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('report_cli', Path(__file__).resolve().parents[1] / 'bin/report.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

class ReportCliTests(unittest.TestCase):
    def test_aggregate_argument_orders_and_sample_counts(self):
        encoded = base64.b64encode(b'{"workflow_version":"0.3.3"}').decode()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); versions = root / 'versions'; versions.mkdir()
            (versions / 'mapping.versions.tsv').write_text('Tool\tVersion\nsamtools\t1.22.1\n')
            paths=[]
            for i in range(2):
                path=root/f'summary{i}.json';path.write_text(json.dumps({'sample':f's{i}'}));paths.append(str(path))
            for count in (0,1,2):
                for order in ('before','between','after'):
                    summaries=paths[:count]; option=['--versions-dir',str(versions)]
                    args={'before':option+[encoded]+summaries,
                          'between':[encoded]+option+summaries,
                          'after':[encoded]+summaries+option}[order]
                    with self.subTest(count=count,order=order), patch.object(sys,'argv',['report.py','aggregate']+args), patch.object(report,'aggregate') as aggregate:
                        report.main()
                        rows,settings,output=aggregate.call_args.args
                        self.assertEqual(len(rows),count)
                        self.assertEqual(settings['software_versions'],{'samtools':'1.22.1'})
                        self.assertEqual(output,Path('.'))
