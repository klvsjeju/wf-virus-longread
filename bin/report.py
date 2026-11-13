#!/usr/bin/env python3
"""Collect compact report data; render with the upstream EPI2ME Labs components."""
import argparse
import base64
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def fasta_stats(path):
    lengths = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            lengths.append(0)
        elif line.strip():
            if not lengths:
                raise ValueError(f'Invalid FASTA: {path}')
            lengths[-1] += len(line.strip())
    return {'sequences': len(lengths), 'total_bases': sum(lengths), 'longest': max(lengths, default=0)}


def read_histograms(stats):
    """Keep fastcat's pre-binned statistics instead of rereading input FASTQ."""
    histograms = {}
    if stats is not None:
        for kind in ('length', 'quality'):
            path = stats / (kind + '.hist')
            if path.exists():
                histograms[kind] = [[float(a), float(b), int(c)] for a, b, c in
                                    csv.reader(path.open(), delimiter='\t') if int(c) > 0]
    return histograms


def read_summary(histograms, metadata, stats=None):
    lengths = histograms.get('length', [])
    count = sum(row[2] for row in lengths)
    median = None
    if count:
        cumulative = 0
        for start, end, n in lengths:
            cumulative += n
            if cumulative >= count / 2:
                median = int(start)
                break
    result = {'read_count': count or metadata.get('records_scanned'),
            'median_length': median,
            'total_bases': int(sum(a * c for a, b, c in lengths)) if count else None,
            'histograms': histograms}
    if stats is not None and (stats / 'per-file-stats.tsv').exists():
        with (stats / 'per-file-stats.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        result['input_read_count'] = sum(int(row.get('r_record_ok', row['n_seqs'])) for row in rows)
        result['excluded_short_reads'] = sum(int(row.get('r_too_short', 0)) for row in rows)
    else:
        result['input_read_count'] = result['read_count']
        result['excluded_short_reads'] = 0
    return result


def coverage_summary(path, mask_depth, max_points=500):
    """Stream depth twice; keep every position in statistics and bounded plot bins."""
    if not path.exists():
        return []
    totals = {}
    with path.open() as handle:
        for line in handle:
            chrom, pos, depth = line.split()[:3]
            pos, depth = int(pos), int(depth)
            row = totals.setdefault(chrom, {'reference': chrom, 'length': 0, 'positions': 0,
                                          'depth_sum': 0, 'covered': 0, 'at_30x': 0,
                                          'at_50x': 0, 'below_mask': 0})
            row['length'] = max(row['length'], pos)
            row['positions'] += 1
            row['depth_sum'] += depth
            row['covered'] += depth > 0
            row['at_30x'] += depth >= 30
            row['at_50x'] += depth >= 50
            row['below_mask'] += depth < mask_depth
    bins = {chrom: defaultdict(lambda: [0, 0]) for chrom in totals}
    with path.open() as handle:
        for line in handle:
            chrom, pos, depth = line.split()[:3]
            width = max(1, (totals[chrom]['length'] + max_points - 1) // max_points)
            item = bins[chrom][(int(pos) - 1) // width]
            item[0] += int(depth)
            item[1] += 1
    for chrom, row in totals.items():
        n = row['positions']
        row['mean_depth'] = row.pop('depth_sum') / n
        for key in ('covered', 'at_30x', 'at_50x', 'below_mask'):
            row[key + '_percent'] = 100 * row[key] / n
        width = max(1, (row['length'] + max_points - 1) // max_points)
        row['plot'] = [{'position': (i * width + 1 + min((i + 1) * width, row['length'])) / 2,
                        'depth': total / count} for i, (total, count) in sorted(bins[chrom].items())]
    return list(totals.values())


def coverage_with_strands(directory, sample, mask_depth=0):
    rows = coverage_summary(directory / f'{sample}.depth', mask_depth)
    for strand in ('fwd', 'rev'):
        other = {row['reference']: row for row in
                 coverage_summary(directory / f'{sample}.{strand}.depth', mask_depth)}
        for row in rows:
            points = {p['position']: p['depth'] for p in other.get(row['reference'], {}).get('plot', [])}
            for point in row['plot']:
                if point['position'] in points:
                    point[strand] = points[point['position']]
    return rows


def variant_summary(path):
    if not path.exists():
        return None
    result = {'counts': {}, 'substitutions': {}, 'indels': [], 'ts': 0, 'tv': 0, 'ts_tv': None}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if fields[0] == 'SN':
                result['counts'][fields[2].rstrip(':')] = int(fields[3])
            elif fields[0] == 'TSTV':
                result.update(ts=int(fields[2]), tv=int(fields[3]),
                              ts_tv=float(fields[4]) if int(fields[3]) else None)
            elif fields[0] == 'ST':
                result['substitutions'][fields[2]] = int(fields[3])
            elif fields[0] == 'IDD':
                result['indels'].append({'length': int(fields[2]), 'count': int(fields[3])})
    return result


def annotations_summary(path, limit=500):
    if not path.exists():
        return None
    preview, impacts, effects, count = [], Counter(), Counter(), 0
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            count += 1
            if len(preview) < limit:
                preview.append(row)
            for key, target in [('ANN[*].IMPACT', impacts), ('ANN[*].EFFECT', effects)]:
                # Count each unique value once per annotated VCF record.
                target.update(set(row.get(key, '').split(';')) - {'', '.'})
    return {'records': count, 'preview': preview, 'preview_limit': limit,
            'impacts': dict(impacts), 'effects': dict(effects)}


def sample_summary(sample, mode, metadata, analysis, stats=None, flye=None, mask_depth=0, coverage=None):
    prefix = sample
    links = {'QC report': f'{sample}/qc/NanoPlot/{sample}_NanoPlot-report.html',
             'Basecalling metadata': f'{sample}/metadata/basecaller.json'}
    extra = {}
    if mode == 'reference':
        fasta = analysis / 'clair3' / f'{sample}_consensus.fa'
        links.update({'Consensus FASTA': f'{prefix}/clair3/{sample}_consensus.fa',
                      'Variants': f'{prefix}/clair3/{sample}_norm.pass.vcf.gz',
                      'All normalized variants': f'{prefix}/clair3/{sample}_norm.vcf.gz',
                      'Alignment BAM': f'{prefix}/alignment/{sample}_sorted.bam',
                      'Depth': f'{prefix}/alignment/{sample}.depth',
                      'Clair3 model selection': f'{prefix}/model_selection.json'})
        extra['coverage'] = coverage_with_strands(analysis / 'alignment', sample, mask_depth)
        extra['variants'] = variant_summary(analysis / 'clair3/variants.stats')
        extra['variant_statistics_source'] = f'{sample}_norm.vcf.gz'
        extra['annotations'] = annotations_summary(analysis / 'annot/var.tsv')
        if extra['annotations'] is not None:
            links['Annotations'] = f'{prefix}/annot/var.tsv'
            links['Annotated VCF'] = f'{prefix}/annot/{sample}_annot.vcf'
        if (analysis / 'clair3/variants.stats').exists():
            links['Variant statistics'] = f'{prefix}/clair3/variants.stats'
        if (analysis / 'clair3/low_coverage.bed').exists():
            links['Mask BED'] = f'{prefix}/clair3/low_coverage.bed'
    else:
        fasta = analysis
        extra['coverage'] = coverage_with_strands(coverage, sample) if coverage and coverage.is_dir() else []
        if coverage and coverage.is_dir():
            links['Depth'] = f'{sample}/coverage/{sample}.depth'
        links['Contigs FASTA'] = f'{prefix}/{sample}.contigs.fasta'
        links['Flye statistics'] = f'{prefix}/flye/assembly_info.txt'
        info = flye / 'assembly_info.txt' if flye and flye.is_dir() else None
        extra['contigs'] = list(csv.DictReader(info.open(), delimiter='\t')) if info and info.exists() else []
    if extra.get('coverage'):
        directory = 'alignment' if mode == 'reference' else 'coverage'
        for strand in ('fwd', 'rev'):
            if any(strand in p for row in extra['coverage'] for p in row['plot']):
                links[f'{strand} depth'] = f'{sample}/{directory}/{sample}.{strand}.depth'
    model_info = analysis / 'model_selection.json' if mode == 'reference' else None
    selection = json.loads(model_info.read_text()) if model_info and model_info.exists() else None
    basecaller = json.loads(metadata.read_text())
    return dict(extra, clair3_model=selection, sample=sample, mode=mode, basecaller=basecaller,
                assembly=fasta_stats(fasta), reads=read_summary(read_histograms(stats), basecaller, stats), files=links)


def aggregate(summaries, settings, output):
    from labs_report import render
    samples = sorted(summaries, key=lambda x: x['sample'])
    for row in samples:
        row['metadata'] = settings.get('sample_metadata', {}).get(row['sample'], {})
    results = {'missing_samples': settings.get('missing_samples', []), 'workflow': 'wf-virus-longread',
               'settings': settings, 'samples': samples}
    (output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    render(samples, settings, output)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='command', required=True)
    s = sub.add_parser('sample')
    s.add_argument('sample'); s.add_argument('mode'); s.add_argument('metadata', type=Path); s.add_argument('analysis', type=Path)
    s.add_argument('--stats', type=Path); s.add_argument('--flye', type=Path)
    s.add_argument('--coverage', type=Path)
    s.add_argument('--mask-depth', type=int, default=0)
    a = sub.add_parser('aggregate')
    a.add_argument('settings'); a.add_argument('summaries', nargs='*', type=Path)
    args = p.parse_args()
    if args.command == 'sample':
        row = sample_summary(args.sample, args.mode, args.metadata, args.analysis, args.stats, args.flye, args.mask_depth, args.coverage)
        Path('summary.json').write_text(json.dumps(row, indent=2))
    else:
        aggregate([json.loads(p.read_text()) for p in args.summaries], json.loads(base64.b64decode(args.settings)), Path('.'))

if __name__ == '__main__':
    main()
