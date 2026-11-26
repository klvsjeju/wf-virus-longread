from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from annotation_gene_names import normalize, clean_name

class GeneNameTests(unittest.TestCase):
    def test_only_orf_suffix_and_idempotence(self):
        self.assertEqual(clean_name('ORF9A_gene'),'ORF9A')
        for name in ('ORF9A','other_gene','ORF1_gene_extra','gene_ORF1','ORF1_gene2'):
            self.assertEqual(clean_name(name),name)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);vcf=root/'annot.vcf';tsv=root/'var.tsv'
            original=('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                      'ref\t1\tORF1_gene\tA\tC\t20\tPASS\tNOTE=ORF1_gene;ANN=C|missense_variant|MODERATE|ORF1_gene|ORF1_gene|transcript|t1,C|synonymous_variant|LOW|ORF9A_gene|g9|transcript|t9,C|other|LOW|other_gene|g2|transcript|t2\n')
            vcf.write_text(original)
            tsv.write_text('Sample\tANN[*].GENE\tANN[*].EFFECT\nORF1_gene\tORF1_gene;ORF9A_gene;other_gene\tORF1_gene\n')
            self.assertEqual(normalize(vcf,'vcf'),2)
            result=vcf.read_text()
            self.assertIn('NOTE=ORF1_gene',result)
            self.assertIn('|ORF1|ORF1_gene|',result)
            self.assertIn('|ORF9A|g9|',result)
            self.assertIn('|other_gene|g2|',result)
            self.assertEqual(normalize(tsv,'tsv'),2)
            self.assertIn('ORF1_gene\tORF1;ORF9A;other_gene\tORF1_gene',tsv.read_text())
            before=vcf.read_bytes();self.assertEqual(normalize(vcf,'vcf'),0);self.assertEqual(vcf.read_bytes(),before)
