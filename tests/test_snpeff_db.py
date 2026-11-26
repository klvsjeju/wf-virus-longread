import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('snpeff_db',Path(__file__).resolve().parents[1]/'bin/snpeff_db.py')
db=importlib.util.module_from_spec(spec);spec.loader.exec_module(db)
GFF='##gff-version 3\nNC_1.1\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\nNC_1.1\ttest\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1\n'

class SnpeffDatabaseTests(unittest.TestCase):
    def test_repeated_descriptions_preserved_structural_duplicates_rejected(self):
        attrs=db.attributes('ID=g1;note=first%3Bpart;note=second%2Cpart;db_xref=A;db_xref=B;function=one;function=two')
        self.assertEqual(attrs['note'],'first%3Bpart,second%2Cpart')
        self.assertEqual(attrs['db_xref'],'A,B')
        self.assertEqual(attrs['function'],'one,two')
        for key in ('ID','Parent','codon_start','transl_table','Derives_from'):
            with self.assertRaisesRegex(ValueError,'structural attribute'):
                db.attributes(f'{key}=a;{key}=b')

    def test_description_merge_is_recorded_in_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);ref,gff=self.fixture(root)
            gff.write_text(GFF.replace('ID=g1','ID=g1;note=first;note=second'))
            normalized,validation=db.validate_gff(gff,db.fasta(ref))
            self.assertIn('note=first,second',normalized)
            self.assertEqual(len(validation['attribute_normalization']),1)

    def fixture(self,root):
        ref=root/'reference.fa';ref.write_text('>NC_1.1 description\nATGAAATAA\n')
        gff=root/'annotation.gff';gff.write_text(GFF)
        return ref,gff

    def fake_build(self,command,**kwargs):
        folder=Path(command[command.index('-dataDir')+1])/command[-1]
        (folder/'snpEffectPredictor.bin').write_bytes(b'example-binary')
        return subprocess.CompletedProcess(command,0,'Built\n','')

    def test_exact_version_bounds_phase_and_parents(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref,gff=self.fixture(root);seq=db.fasta(ref)
            normalized,report=db.validate_gff(gff,seq)
            self.assertIn('Parent=wfvirus_transcript_1',normalized)
            self.assertEqual(report['features']['mRNA'],1)
            self.assertFalse(report['independent_protein_check'])
            for bad,message in [(GFF.replace('NC_1.1','NC_1.2'),'exactly match'),
                                (GFF.replace('\t1\t9\t','\t1\t10\t'),'exceed reference'),
                                (GFF.replace('\t+\t0\t','\t+\t.\t'),'phase'),
                                (GFF.replace('Parent=g1','Parent=missing'),'missing')]:
                gff.write_text(bad)
                with self.assertRaisesRegex(ValueError,message):db.validate_gff(gff,seq)

    def test_embedded_sequence_and_duplicate_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref,gff=self.fixture(root)
            gff.write_text(GFF+'##FASTA\n>NC_1.1\nATGCCCTAA\n')
            with self.assertRaisesRegex(ValueError,'differs'):db.validate_gff(gff,db.fasta(ref))
            ref.write_text('>NC_1.1\nATG\n>NC_1.1\nAAA\n')
            with self.assertRaisesRegex(ValueError,'duplicate'):db.fasta(ref)

    def test_build_reuse_integrity_and_version(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref,gff=self.fixture(root);cache=root/'cache'
            with patch.object(db.subprocess,'run',side_effect=self.fake_build) as run:
                first=db.prepare(ref,gff,cache,root/'out1',version='SnpEff test')
                second=db.prepare(ref,None,cache,root/'out2',version='SnpEff test')
                self.assertFalse(first['cache_reused']);self.assertTrue(second['cache_reused'])
                self.assertEqual(run.call_count,1)
                self.assertEqual(first['db_id'],second['db_id'])
                with self.assertRaisesRegex(ValueError,'No unique verified'):
                    db.prepare(ref,None,cache,root/'different-version',version='SnpEff other')
                (cache/first['db_id']/'data'/first['db_id']/'snpEffectPredictor.bin').write_bytes(b'changed')
                with self.assertRaisesRegex(ValueError,'No unique verified'):
                    db.prepare(ref,None,cache,root/'corrupted',version='SnpEff test')

    def test_ambiguous_annotations_require_gff(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref,gff=self.fixture(root);cache=root/'cache'
            with patch.object(db.subprocess,'run',side_effect=self.fake_build):
                db.prepare(ref,gff,cache,root/'one',version='SnpEff test')
                gff.write_text(GFF+'# second annotation version\n')
                db.prepare(ref,gff,cache,root/'two',version='SnpEff test')
                with self.assertRaisesRegex(ValueError,'No unique verified'):
                    db.prepare(ref,None,cache,root/'three',version='SnpEff test')
                chosen=db.prepare(ref,gff,cache,root/'four',version='SnpEff test')
                self.assertTrue(chosen['cache_reused'])

    def test_new_ui_preserves_unrelated_user_changes(self):
        import json
        schema=json.loads((Path(__file__).resolve().parents[1]/'nextflow_schema.json').read_text())
        keys=schema['definitions']['annotation']['properties']
        self.assertIn('annotation_gff',keys);self.assertIn('snpeff_cache',keys)
        self.assertNotIn('snpeff_config',keys)
        for entry in schema['allOf']:
            self.assertNotIn('snpeff_config',entry.get('then',{}).get('required',[]))
