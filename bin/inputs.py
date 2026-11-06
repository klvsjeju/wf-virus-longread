#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path


def merge(paths, output, metadata, override=''):
    detected = set()
    with gzip.open(output, 'wb') as out:
        for path in paths:
            opener = gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rb') as src:
                count = 0
                while True:
                    header = src.readline()
                    if not header:
                        break
                    seq, plus, qual = src.readline(), src.readline(), src.readline()
                    if not header.startswith(b'@') or not plus.startswith(b'+') or not seq or not qual:
                        raise ValueError(f'Invalid four-line FASTQ: {path}')
                    if len(seq.rstrip(b'\r\n')) != len(qual.rstrip(b'\r\n')):
                        raise ValueError(f'Sequence/quality length mismatch: {path}')
                    text = header.decode('utf-8', errors='replace')
                    for key in ('basecall_model_version_id', 'basecall_model', 'model_version_id'):
                        m = re.search(r'(?:^|\s)' + key + r'=([^\s]+)', text)
                        if m:
                            detected.add(m.group(1))
                    for line in (header, seq, plus, qual):
                        out.write(line.rstrip(b'\r\n') + b'\n')
                    count += 1
                if not count:
                    raise ValueError(f'Empty FASTQ: {path}')
    if len(detected) > 1 and not override:
        raise ValueError(f'Mixed models in sampled FASTQ headers: {sorted(detected)}')
    inferred = next(iter(detected)) if len(detected) == 1 else None
    if len(detected) > 1:
        print(f'WARNING: manual model {override} overrides mixed detected models: {sorted(detected)}', file=sys.stderr)
    if override and inferred and override != inferred:
        print(f'WARNING: manual model {override} overrides detected {inferred}', file=sys.stderr)
    Path(metadata).write_text(json.dumps({'detected': inferred, 'manual': override or None,
                                         'selected': override or inferred,
                                         'detected_models': sorted(detected),
                                         'header_scan': 'all_records'}, indent=2))


def model(root, mapping, metadata):
    root = Path(root).resolve()
    selected = json.loads(Path(metadata).read_text())['selected']
    direct = (root / 'pileup.index').is_file() and (root / 'full_alignment.index').is_file()
    if direct:
        chosen = root
        print('Using explicitly supplied custom Clair3 checkpoint directory', file=sys.stderr)
    else:
        if not selected:
            raise ValueError('No basecaller model found. Select override_basecaller_cfg or provide a direct Clair3 checkpoint folder.')
        with open(mapping, newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        matches = [r['relative_path'] for r in rows if r['basecaller_cfg'] == selected]
        if len(matches) != 1:
            raise ValueError(f'Need exactly one verified Clair3 mapping for {selected}')
        chosen = (root / matches[0]).resolve()
        if not chosen.is_relative_to(root):
            raise ValueError('Model mapping must stay inside model_dir')
    for prefix in ('pileup', 'full_alignment'):
        if not (chosen / f'{prefix}.index').is_file() or not list(chosen.glob(f'{prefix}.data-*')):
            raise ValueError(f'Missing checkpoint: {chosen}/{prefix}')
    return str(chosen)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='command', required=True)
    m = sub.add_parser('merge')
    m.add_argument('--output', required=True)
    m.add_argument('--metadata', required=True)
    m.add_argument('--override', default='')
    m.add_argument('paths', nargs='+')
    r = sub.add_parser('model')
    r.add_argument('--root', required=True)
    r.add_argument('--mapping', required=True)
    r.add_argument('--metadata', required=True)
    a = p.parse_args()
    try:
        if a.command == 'merge':
            merge(a.paths, a.output, a.metadata, a.override)
        else:
            print(model(a.root, a.mapping, a.metadata))
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as e:
        p.exit(2, str(e) + '\n')


if __name__ == '__main__':
    main()
