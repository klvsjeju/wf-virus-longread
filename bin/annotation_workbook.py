#!/usr/bin/env python3
"""Write an annotation XLSX using only the Python standard library."""
import argparse
import csv
import re
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.sax.saxutils import escape
from zipfile import ZipFile, ZIP_DEFLATED

NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
SHEETS = ('SNPs', 'Indels', 'Other variants')

def variant_sheet(ref, alt):
    """Classify simple, biallelic sequence alleles; preserve other records."""
    if not re.fullmatch('[ACGTNacgtn]+', ref or '') or not re.fullmatch('[ACGTNacgtn]+', alt or ''):
        return 'Other variants'
    if ref.upper() == alt.upper():
        return 'Other variants'
    if len(ref) == len(alt) == 1:
        return 'SNPs'
    if len(ref) != len(alt):
        return 'Indels'
    return 'Other variants'

def column(index):
    result = ''
    while index:
        index, rem = divmod(index-1, 26)
        result = chr(65+rem)+result
    return result

def row_xml(values, number):
    cells = []
    for i, value in enumerate(values, 1):
        if len(value) > 32767:
            raise ValueError('Annotation cell exceeds Excel 32767-character limit; use var.tsv for full values')
        if re.search('[\x00-\x08\x0b\x0c\x0e-\x1f]', value):
            raise ValueError('Annotation cell contains an XML-incompatible control character')
        # Explicit strings preserve IDs and prevent formula interpretation.
        cells.append(f'<c r="{column(i)}{number}" t="inlineStr"><is><t xml:space="preserve">{escape(value)}</t></is></c>')
    return f'<row r="{number}">'+''.join(cells)+'</row>'

def write_workbook(source, output):
    counts = dict.fromkeys(SHEETS, 0)
    output = Path(output)
    with TemporaryDirectory(dir=output.parent) as tmp:
        paths = [Path(tmp)/f'sheet{i}.xml' for i in range(1,4)]
        handles = [p.open('w',encoding='utf-8') for p in paths]
        try:
            with Path(source).open(newline='',encoding='utf-8-sig') as inp:
                reader = csv.reader(inp,delimiter='\t')
                header = next(reader)
                ref_i, alt_i = header.index('REF'), header.index('ALT')
                if len(header)>16384: raise ValueError('Too many columns for Excel')
                for h in handles:
                    h.write(f'<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="{NS}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" state="frozen"/></sheetView></sheetViews><sheetData>')
                    h.write(row_xml(header,1))
                for line,values in enumerate(reader,2):
                    if not values: continue
                    if len(values)!=len(header):raise ValueError(f'TSV line {line}: column count differs from header')
                    sheet=variant_sheet(values[ref_i],values[alt_i]);counts[sheet]+=1
                    if counts[sheet]>1048575:raise ValueError('Too many annotation rows for one Excel sheet')
                    handles[SHEETS.index(sheet)].write(row_xml(values,counts[sheet]+1))
                for name,h in zip(SHEETS,handles):
                    h.write(f'</sheetData><autoFilter ref="A1:{column(len(header))}{counts[name]+1}"/></worksheet>')
        finally:
            for h in handles:h.close()
        relns='http://schemas.openxmlformats.org/package/2006/relationships'
        types=['<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>','<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>']
        types += [f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,4)]
        staged=Path(tmp)/'workbook.xlsx'
        with ZipFile(staged,'w',ZIP_DEFLATED) as z:
            z.writestr('[Content_Types].xml','<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'+''.join(types)+'</Types>')
            z.writestr('_rels/.rels',f'<Relationships xmlns="{relns}"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
            z.writestr('xl/workbook.xml',f'<workbook xmlns="{NS}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><bookViews><workbookView/></bookViews><sheets>'+''.join(f'<sheet name="{name}" sheetId="{i}" r:id="rId{i}"/>' for i,name in enumerate(SHEETS,1))+'</sheets></workbook>')
            z.writestr('xl/_rels/workbook.xml.rels',f'<Relationships xmlns="{relns}">'+''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,4))+'</Relationships>')
            for i,p in enumerate(paths,1):z.write(p,f'xl/worksheets/sheet{i}.xml')
        staged.replace(output)
    return counts

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tsv',type=Path);parser.add_argument('xlsx',type=Path)
    args=parser.parse_args()
    print(write_workbook(args.tsv,args.xlsx))
