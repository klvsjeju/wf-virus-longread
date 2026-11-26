#!/usr/bin/env python3
"""Remove the descriptive _gene suffix from ORF gene names in output copies."""
import argparse
import os
from pathlib import Path
import re
import tempfile

ORF_GENE = re.compile(r'^(ORF[0-9]+[A-Za-z]?)_gene$')


def clean_name(name):
    return ORF_GENE.sub(r'\1', name)


def normalize(path, kind):
    path = Path(path)
    changes = 0
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='',
                                     dir=path.parent, delete=False) as output:
        staged = Path(output.name)
        try:
            with path.open(encoding='utf-8', newline='') as source:
                gene_column = None
                for line in source:
                    newline = '\r\n' if line.endswith('\r\n') else '\n' if line.endswith('\n') else ''
                    fields = line.rstrip('\r\n').split('\t')
                    changed = False
                    if kind == 'vcf' and not line.startswith('#') and len(fields) >= 8:
                        info = fields[7].split(';')
                        for i, entry in enumerate(info):
                            if not entry.startswith('ANN='):
                                continue
                            annotations = entry[4:].split(',')
                            for j, annotation in enumerate(annotations):
                                values = annotation.split('|')
                                # ANN index 3 is Gene_Name; Gene_ID and all other
                                # fields are preserved, including multi-transcript entries.
                                if len(values) > 3:
                                    name = clean_name(values[3])
                                    if name != values[3]:
                                        values[3] = name
                                        annotations[j] = '|'.join(values)
                                        changes += 1
                                        changed = True
                            info[i] = 'ANN=' + ','.join(annotations)
                        fields[7] = ';'.join(info)
                    elif kind == 'tsv':
                        if gene_column is None:
                            gene_column = next((i for i, value in enumerate(fields)
                                                if value in ('ANN[*].GENE', 'ANN.GENE')), None)
                            if gene_column is None:
                                raise ValueError('Annotation TSV has no ANN gene-name column')
                        else:
                            values = fields[gene_column].split(';')
                            cleaned = [clean_name(value) for value in values]
                            changes += sum(a != b for a, b in zip(values, cleaned))
                            changed = values != cleaned
                            fields[gene_column] = ';'.join(cleaned)
                    output.write('\t'.join(fields) + newline if changed else line)
        except BaseException:
            staged.unlink(missing_ok=True)
            raise
    os.replace(staged, path)
    return changes


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', required=True, type=Path)
    parser.add_argument('--tsv', required=True, type=Path)
    args = parser.parse_args()
    for kind in ('vcf', 'tsv'):
        print(f'{kind}: normalized {normalize(getattr(args, kind), kind)} ORF gene names')
