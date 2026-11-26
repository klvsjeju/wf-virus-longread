"""Exercise the production bcftools selection against threshold boundaries."""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

BCFTOOLS = shutil.which('bcftools') or str(Path.home() / 'miniconda3/envs/vzv_longread/bin/bcftools')

@unittest.skipUnless(Path(BCFTOOLS).is_file(), 'bcftools is required for filter integration')
class VariantFilterTests(unittest.TestCase):
    def test_inclusive_thresholds_missing_values_and_caller_status(self):
        script = (Path(__file__).resolve().parents[1] / 'bin/reference.sh').read_text()
        selection = script[script.index('variant_filter='):script.index('bcftools stats')]
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / 'result/varcalling'
            folder.mkdir(parents=True)
            header = ('##fileformat=VCFv4.2\n##contig=<ID=ref,length=100>\n'
                      '##FILTER=<ID=LowQual,Description="Caller failure">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample depth">\n'
                      '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n')
            cases = [(20,30,0.5,'PASS'),(19.99,30,0.9,'PASS'),(20,29,0.9,'PASS'),
                     (21,31,0.8,'PASS'),('.',30,0.9,'PASS'),(20,'.',0.9,'PASS'),
                     (30,40,0.9,'LowQual'),(30,40,0.4999,'PASS'),(30,40,'.','PASS')]
            text = header + ''.join(f'ref\t{i}\t.\tA\tC\t{q}\t{f}\t.\tGT:DP:AF\t1:{d}:{af}\n'
                                    for i,(q,d,af,f) in enumerate(cases,1))
            source = folder / 's_merge_output.vcf.gz'
            reference = Path(tmp) / 'ref.fa'
            reference.write_text('>ref\n'+'A'*100+'\n')
            subprocess.run([BCFTOOLS,'view','-Oz','-o',str(source)],input=text,text=True,check=True,capture_output=True)
            env = dict(os.environ, SAMPLE='s', REF=str(reference), PATH=str(Path(BCFTOOLS).parent)+os.pathsep+os.environ['PATH'])
            for qual,depth,af,norm,expected in [(20,30,0.5,['1','4','7','8','9'],['1','4']),
                                                (21,31,0.5,['4','7','8','9'],['4']),
                                                (20,30,0.8,['1','4','7','8','9'],['4'])]:
                with self.subTest(qual=qual,depth=depth):
                    env.update(MIN_VARIANT_QUAL=str(qual),MIN_VARIANT_DEPTH=str(depth),MIN_VARIANT_AF=str(af))
                    subprocess.run(['bash','-euo','pipefail','-c',selection],cwd=tmp,env=env,check=True,capture_output=True)
                    result=subprocess.run([BCFTOOLS,'query','-f','%POS\n',str(folder/'s_norm.cons.vcf.gz')],capture_output=True,text=True,check=True)
                    self.assertEqual(result.stdout.splitlines(),expected)
                    normalized=subprocess.run([BCFTOOLS,'query','-f','%POS\n',str(folder/'s_norm.vcf.gz')],capture_output=True,text=True,check=True)
                    self.assertEqual(normalized.stdout.splitlines(),norm)
            count=subprocess.run([BCFTOOLS,'view','-H',str(source)],capture_output=True,text=True,check=True)
            self.assertEqual(len(count.stdout.splitlines()),len(cases))
