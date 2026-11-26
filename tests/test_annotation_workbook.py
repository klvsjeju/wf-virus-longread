import sys
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile
from xml.etree import ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from annotation_workbook import write_workbook

class WorkbookTests(unittest.TestCase):
    def test_classification_all_rows_columns_and_literal_strings(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);tsv=root/'var.tsv';out=root/'var.xlsx'
            rows=[['001','A','G','=1+1'],['002','AT','A','gene<&>'],['003','A','AT','.'],
                  ['004','AT','GC','MNP'],['005','A','<DEL>','symbolic']]
            original='Sample\tREF\tALT\tNote\n'+''.join('\t'.join(row)+'\n' for row in rows)
            tsv.write_text(original)
            self.assertEqual(write_workbook(tsv,out),{'SNPs':1,'Indels':2,'Other variants':2})
            ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            recovered=[]
            with ZipFile(out) as z:
                for name in z.namelist():ET.fromstring(z.read(name))
                wb=ET.fromstring(z.read('xl/workbook.xml'))
                self.assertEqual([s.attrib['name'] for s in wb.findall('x:sheets/x:sheet',ns)],['SNPs','Indels','Other variants'])
                for i in range(1,4):
                    sheet=ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml'))
                    self.assertEqual(sheet.findall('.//x:f',ns),[])
                    for row in sheet.findall('x:sheetData/x:row',ns)[1:]:
                        recovered.append([c.find('x:is/x:t',ns).text for c in row])
            self.assertCountEqual(recovered,rows)
            self.assertEqual(tsv.read_text(),original)
    def test_empty_annotation_has_three_header_sheets(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);src=root/'var.tsv';src.write_text('Sample\tREF\tALT\n')
            self.assertEqual(sum(write_workbook(src,root/'var.xlsx').values()),0)
