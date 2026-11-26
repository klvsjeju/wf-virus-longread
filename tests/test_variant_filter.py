"""Run production normalization and independent selections with real bcftools."""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest
import sys
import base64
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from analysis_options import prepare, resolve

BCFTOOLS=shutil.which('bcftools') or str(Path.home()/'miniconda3/envs/vzv_longread/bin/bcftools')

class AnalysisOptionsTests(unittest.TestCase):
    def test_advanced_controls_are_consolidated_and_mode_labelled(self):
        import json
        root=Path(__file__).resolve().parents[1]
        fields=json.loads((root/'nextflow_schema.json').read_text())['definitions']['advanced']['properties']
        for key in ('min_variant_qual','min_variant_depth','min_variant_af','mask_depth','mask_depth_opts'):
            self.assertNotIn(key,fields)
        for key in ('variant_opts','consensus_opts','mask_opts'):
            self.assertTrue(fields[key]['title'].startswith('[reference-based]'))
        for key in ('flye_genome_size','min_read_length','flye_opts'):
            self.assertTrue(fields[key]['title'].startswith('[de novo]'))
        for kind,key in [('variant','variant_opts'),('consensus','consensus_opts'),('mask','mask_opts')]:
            resolved=resolve(kind,fields[key]['default'])
            expected={'variant':{'min_dp':30,'min_qual':20},'consensus':{'min_dp':30,'min_qual':20,'min_af':0.5},'mask':{'min_depth':0}}[kind]
            for name,value in expected.items():self.assertEqual(resolved[name],value)

    def test_defaults_zero_and_validation(self):
        self.assertEqual(resolve('consensus','')['min_af'],0.5)
        for kind,text in [('variant','--min-qual 0 --min-dp 0'),('consensus','--min-qual 0 --min-dp 0 --min-af 0')]:
            self.assertEqual(resolve(kind,text)['expression'],'')
        self.assertNotIn('FORMAT/DP',resolve('variant','--min-dp 0')['expression'])
        for kind,text in [('variant','--min-dp -1'),('consensus','--min-af 50'),('variant','--min-dp 1.5'),('variant','--min-qual nan')]:
            with self.subTest(text=text),self.assertRaises(SystemExit):resolve(kind,text)
        self.assertEqual(resolve('mask','--min-depth 10 -q 15 -Q 20')['extra'],['-q','15','-Q','20'])
        self.assertEqual(resolve('consensus','--min-af 0 --mask-with N')['extra'],['--mask-with','N'])
        with self.assertRaises(ValueError):resolve('mask','--min-depth 10 -o elsewhere')

@unittest.skipUnless(Path(BCFTOOLS).is_file(),'Requires bcftools')
class VariantFilterTests(unittest.TestCase):
    def test_independent_thresholds_and_zero_disable(self):
        script=(Path(__file__).resolve().parents[1]/'bin/reference.sh').read_text()
        selection=script[script.index('# Both outputs'):script.index('bcftools stats')]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder=root/'result/varcalling';folder.mkdir(parents=True)
            reference=root/'ref.fa';reference.write_text('>ref\n'+'A'*100+'\n')
            header=('##fileformat=VCFv4.2\n##contig=<ID=ref,length=100>\n'
                    '##FILTER=<ID=LowQual,Description="Caller failure">\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n')
            cases=[(20,30,0.5,'PASS'),(19.99,30,0.9,'PASS'),(20,29,0.9,'PASS'),
                   (21,31,0.8,'PASS'),('.',30,0.9,'PASS'),(20,'.',0.9,'PASS'),
                   (30,40,0.9,'LowQual'),(30,40,0.4999,'PASS'),(30,40,'.','PASS')]
            text=header+''.join(f'ref\t{i}\t.\tA\tC\t{q}\t{f}\t.\tGT:DP:AF\t1:{d}:{af}\n' for i,(q,d,af,f) in enumerate(cases,1))
            source=folder/'s_merge_output.vcf.gz'
            subprocess.run([BCFTOOLS,'view','-Oz','-o',str(source)],input=text,text=True,check=True,capture_output=True)
            all_calls=[str(i) for i in range(1,10)]
            tests=[('', '', ['1','4','7','8','9'], ['1','4']),
                   ('--min-qual 100 --min-dp 100','',[],['1','4']),
                   ('--min-qual 0 --min-dp 0','--min-qual 0 --min-dp 0 --min-af 0',all_calls,['1','2','3','4','5','6','8','9']),
                   ('--min-dp 0','--min-dp 0 --min-af 0',['1','3','4','6','7','8','9'],['1','3','4','6','8','9']),
                   ('','--min-af 0.8',['1','4','7','8','9'],['4'])]
            for variant,consensus,want_norm,want_cons in tests:
                with self.subTest(variant=variant,consensus=consensus):
                    prepare(variant,consensus,'--min-depth 10 -q 15',root/'resolved.json')
                    env=dict(os.environ,SAMPLE='s',REF=str(reference),PATH=str(Path(BCFTOOLS).parent)+os.pathsep+os.environ['PATH'])
                    subprocess.run(['bash','-euo','pipefail','-c','source resolved.env\n'+selection],cwd=root,env=env,check=True,capture_output=True)
                    for suffix,expected in [('norm',want_norm),('norm.cons',want_cons)]:
                        got=subprocess.run([BCFTOOLS,'query','-f','%POS\n',str(folder/f's_{suffix}.vcf.gz')],text=True,capture_output=True,check=True)
                        self.assertEqual(got.stdout.splitlines(),expected)
            self.assertEqual(len(subprocess.run([BCFTOOLS,'view','-H',str(source)],text=True,capture_output=True,check=True).stdout.splitlines()),9)
            # Exercise the entire postprocessing script, including option decoding,
            # depth generation, mask propagation and consensus output.
            alignment=root/'result/alignment';alignment.mkdir()
            samtools=str(Path(BCFTOOLS).with_name('samtools'))
            sam='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:100\nread\t0\tref\t1\t60\t100M\t*\t0\t0\t'+'A'*100+'\t'+'I'*100+'\n'
            subprocess.run([samtools,'view','-b','-o',str(alignment/'s_sorted.bam')],input=sam,text=True,check=True,capture_output=True)
            subprocess.run([samtools,'index',str(alignment/'s_sorted.bam')],check=True,capture_output=True)
            helper=Path(__file__).resolve().parents[1]/'bin/run_options.py'
            enc=lambda text:base64.b64encode(text.encode()).decode()
            env.update(OPTIONS_RUNNER=str(helper),VARIANT_OPTS_B64=enc('--min-dp 0 --min-qual 0'),
                       CONSENSUS_OPTS_B64=enc('--min-dp 30 --min-qual 20 --min-af 0.5'),MASK_OPTS_B64=enc('--min-depth 2'))
            subprocess.run(['bash',str(helper.with_name('reference.sh'))],cwd=root,env=env,check=True,capture_output=True)
            self.assertTrue((folder/'s_consensus.fa').exists())
            self.assertEqual(len((folder/'low_coverage.bed').read_text().splitlines()),100)
            import json
            resolved=json.loads((root/'result/analysis_options.json').read_text())
            self.assertEqual(resolved['mask']['min_depth'],2)
            self.assertEqual(resolved['variant']['expression'],'')
