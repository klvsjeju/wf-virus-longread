#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path


# FASTQ comments can contain legacy key=value metadata or SAM read-group tags.
_MODEL_TAG = re.compile(r'(?:^|\s)(basecall_model_version_id|basecall_model|model_version_id)=([^\s]+)')
_RG_TAG = re.compile(r'(?:^|\s)RG(?::Z:|=)([^\s]+)')
_VERSIONED_MODEL = re.compile(r'(?:^|_)((?:dna|rna)_[A-Za-z0-9_.-]+@v[0-9]+(?:\.[0-9]+)*)(?=$|[_+])')


def extract_models(header):
    """Prefer explicit model tags; otherwise extract a basecaller from RG metadata.

    Never search read IDs, filenames, sequences or barcode values for models.
    This implements FASTQ metadata detection, not BAM header DS parsing.
    """
    explicit = {m.group(2) for m in _MODEL_TAG.finditer(header)
                if m.group(2).lower() not in ('unknown', 'none', 'null', 'n/a')}
    if explicit:
        return explicit, 'explicit_model_tag'
    found = set()
    for rg in _RG_TAG.finditer(header):
        value = rg.group(1)
        found.update(m.group(1) for m in _VERSIONED_MODEL.finditer(value))
    return found, 'read_group' if found else None


def merge(paths, output, metadata, override=''):
    detected = set()
    detection_sources = {}
    records = 0
    records_with_model = 0
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
                    models, source = extract_models(text)
                    detected.update(models)
                    records += 1
                    if models:
                        records_with_model += 1
                        detection_sources[source] = detection_sources.get(source, 0) + 1
                    for line in (header, seq, plus, qual):
                        out.write(line.rstrip(b'\r\n') + b'\n')
                    count += 1
                if not count:
                    raise ValueError(f'Empty FASTQ: {path}')
    write_metadata(metadata, detected, override, detection_sources, records, records_with_model)


def write_metadata(metadata, detected, override, detection_sources, records, records_with_model, method='all_records'):
    if len(detected) > 1 and not override:
        raise ValueError(f'Mixed models in FASTQ headers: {sorted(detected)}')
    if not detected and not override:
        print('WARNING: No basecalling model found in FASTQ model tags or RG:Z read-group tags; inspect basecaller.json diagnostics.', file=sys.stderr)
    inferred = next(iter(detected)) if len(detected) == 1 else None
    if len(detected) > 1:
        print(f'WARNING: manual model {override} overrides mixed detected models: {sorted(detected)}', file=sys.stderr)
    if override and inferred and override != inferred:
        print(f'WARNING: manual model {override} overrides detected {inferred}', file=sys.stderr)
    Path(metadata).write_text(json.dumps({'detected': inferred, 'manual': override or None,
                                         'selected': override or inferred,
                                         'detected_models': sorted(detected),
                                         'header_scan': method,
                                         'model_detection_sources': detection_sources,
                                         'records_scanned': records,
                                         'records_with_model': records_with_model}, indent=2))



def fastcat_metadata(stats, metadata, paths, override='', filtered_reads=None):
    """Use fastcat stats and scan headers when its model counts are incomplete."""
    with (stats / 'per-file-stats.tsv').open(newline='') as handle:
        file_counts = [int(row['n_seqs']) for row in csv.DictReader(handle, delimiter='\t')]
    if (len(file_counts) != len(paths) or not sum(file_counts)
            or any(count < 0 for count in file_counts)
            or (filtered_reads is None and not all(file_counts))):
        raise ValueError('fastcat did not report reads for every input FASTQ')
    records = sum(file_counts)
    with (stats / 'per-file-basecallers.tsv').open(newline='') as handle:
        models = [(row['basecaller'], int(row['count']))
                  for row in csv.DictReader(handle, delimiter='\t')
                  if row['basecaller'].strip().lower() not in ('', 'unknown', 'none', 'null', 'n/a')]
    detected = {name for name, count in models if count > 0}
    tagged = sum(count for _, count in models)
    if tagged > records:
        raise ValueError('fastcat model counts exceed the read count')
    sources = {'fastcat': tagged} if tagged else {}
    if tagged < records:
        # Keep support for explicit and RG header formats missed by fastcat.
        # Read the original files without recompressing them.
        fallback_sources = {}
        fallback_tagged = 0
        for path in ([filtered_reads] if filtered_reads is not None else paths):
            opener = gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rb') as handle:
                for index, line in enumerate(handle):
                    if index % 4 == 0:
                        found, source = extract_models(line.decode('utf-8', errors='replace'))
                        detected.update(found)
                        if found:
                            fallback_tagged += 1
                            fallback_sources[source] = fallback_sources.get(source, 0) + 1
        # Both readers may detect the same records.
        tagged = max(tagged, fallback_tagged)
        sources.update(fallback_sources)
    write_metadata(metadata, detected, override, sources, records, tagged, 'fastcat_with_header_fallback')


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
        try:
            chosen.relative_to(root)
        except ValueError as exc:
            raise ValueError('Model mapping must stay inside model_dir') from exc
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
    f = sub.add_parser('fastcat-metadata')
    f.add_argument('--stats', type=Path, required=True)
    f.add_argument('--metadata', required=True)
    f.add_argument('--override', default='')
    f.add_argument('--filtered-reads', type=Path)
    f.add_argument('paths', nargs='+')
    r = sub.add_parser('model')
    r.add_argument('--root', required=True)
    r.add_argument('--mapping', required=True)
    r.add_argument('--metadata', required=True)
    a = p.parse_args()
    try:
        if a.command == 'merge':
            merge(a.paths, a.output, a.metadata, a.override)
        elif a.command == 'fastcat-metadata':
            fastcat_metadata(a.stats, a.metadata, a.paths, a.override, a.filtered_reads)
        else:
            print(model(a.root, a.mapping, a.metadata))
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as e:
        p.exit(2, str(e) + '\n')


if __name__ == '__main__':
    main()
