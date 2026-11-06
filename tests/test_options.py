import base64
import importlib.util
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('run_options', Path(__file__).resolve().parents[1] / 'bin/run_options.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
def encode(s):
    return base64.b64encode(s.encode()).decode()

class OptionsTests(unittest.TestCase):
    def test_options_and_quoted_expression(self):
        self.assertEqual(runner.parse_options('flye', encode('--meta --iterations 2')), ['--meta', '--iterations', '2'])
        self.assertEqual(runner.parse_options('depth', encode('-q 15 -Q 20')), ['-q', '15', '-Q', '20'])
        self.assertEqual(runner.parse_options('consensus', encode('--exclude "INFO/DP < 10"')), ['--exclude', 'INFO/DP < 10'])
        self.assertEqual(runner.parse_options('flye', ''), [])

    def test_owned_arguments_and_invalid_quotes(self):
        for kind, option in [('flye', '--out-dir=elsewhere'), ('flye', '-t8'), ('flye', '--nano-raw'), ('depth', '-oout'), ('depth', '-H'), ('consensus', '--mask=x'), ('consensus', '-sother')]:
            with self.subTest(option=option), self.assertRaises(ValueError):
                runner.parse_options(kind, encode(option))
        with self.assertRaises(ValueError):
            runner.parse_options('flye', encode('--meta "'))

    def test_execution_keeps_values_as_arguments(self):
        argv = ['run_options.py', 'consensus', encode('--exclude "INFO/DP < 10"'), '--', 'bcftools', 'consensus', '-f', 'reference.fa', 'in.vcf.gz']
        with patch.object(sys, 'argv', argv), patch.object(subprocess, 'run') as run:
            runner.main()
        run.assert_called_once_with(['bcftools', 'consensus', '--exclude', 'INFO/DP < 10', '-f', 'reference.fa', 'in.vcf.gz'], check=True)

if __name__ == '__main__':
    unittest.main()
