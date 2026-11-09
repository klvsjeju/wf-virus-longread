#!/usr/bin/env python3
"""Create portable workflow summary artifacts using staged task outputs."""
import argparse
import base64
import html
import json
from pathlib import Path
from urllib.parse import quote


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


def sample_summary(sample, mode, metadata, analysis):
    prefix = f'reference/{sample}/result' if mode == 'reference' else f'denovo/{sample}'
    links = {'QC report': f'qc/{sample}/NanoPlot/{sample}_NanoPlot-report.html',
             'Basecalling metadata': f'metadata/{sample}/basecaller.json'}
    if mode == 'reference':
        fasta = analysis / 'clair3' / f'{sample}_consensus.fa'
        links.update({'Consensus FASTA': f'{prefix}/clair3/{sample}_consensus.fa',
                      'Variants': f'{prefix}/clair3/norm.pass.vcf.gz',
                      'Alignment BAM': f'{prefix}/alignment/{sample}_sorted.bam',
                      'Depth': f'{prefix}/alignment/{sample}.depth',
                      'Clair3 model selection': f'{prefix}/model_selection.json'})
        if (analysis / 'annot/var.tsv').exists():
            links['Annotations'] = f'{prefix}/annot/var.tsv'
        if (analysis / 'clair3/low_coverage.bed').exists():
            links['Mask BED'] = f'{prefix}/clair3/low_coverage.bed'
    else:
        fasta = analysis
        links['Contigs FASTA'] = f'{prefix}/{sample}.contigs.fasta'
        links['Flye statistics'] = f'{prefix}/flye/assembly_info.txt'
    model_info = analysis / 'model_selection.json' if mode == 'reference' else None
    selection = json.loads(model_info.read_text()) if model_info and model_info.exists() else None
    return {'clair3_model': selection, 'sample': sample, 'mode': mode, 'basecaller': json.loads(metadata.read_text()),
            'assembly': fasta_stats(fasta), 'files': links}


def aggregate(summaries, settings, output):
    samples = sorted(summaries, key=lambda x: x['sample'])
    results = {'workflow': 'wf-virus-longread', 'settings': settings, 'samples': samples}
    (output / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    esc = lambda x: html.escape(str(x), quote=True)
    sections = []
    for row in samples:
        stats = row['assembly']
        links = ''.join(f'<li><a href="{quote(path, safe="/")}">{esc(label)}</a></li>' for label, path in row['files'].items())
        model = row.get('clair3_model') or {}
        model_html = f'<p>Clair3: {esc(model.get("model_id", "Not used"))} · Source: {esc(model.get("source", "Not used"))}</p>'
        sections.append(f'<section><h2>{esc(row["sample"])}</h2><p>Mode: {esc(row["mode"])} · Basecalling model: {esc(row["basecaller"].get("selected") or "Unknown")}</p>' +
                        model_html + f'<p>Sequences: {stats["sequences"]} · Total bases: {stats["total_bases"]} · Longest: {stats["longest"]}</p><ul>{links}</ul></section>')
    settings_html = ''.join(f'<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>' for k, v in sorted(settings.items()))
    document = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>wf-virus-longread report</title>'
    document += '<style>body{font:16px system-ui,sans-serif;max-width:1100px;margin:40px auto;padding:0 24px;color:#172c38}section{border-top:1px solid #ccd;padding:16px 0}table{border-collapse:collapse}th,td{border:1px solid #ccd;padding:8px;text-align:left;overflow-wrap:anywhere}a{color:#006c85}</style>'
    document += '<h1>wf-virus-longread analysis report</h1><p>Completed analysis outputs. De novo sequences are Flye contigs; reference sequences are variant-based consensus.</p>'
    document += ''.join(sections) + '<h2>Run settings</h2><table>' + settings_html + '</table>'
    document += '<h2>Execution</h2><ul><li><a href="execution/report.html">Nextflow execution report</a></li><li><a href="execution/timeline.html">Timeline</a></li><li><a href="execution/trace.tsv">Task trace</a></li></ul><p>Execution reports are finalized when Nextflow exits. Open result files through the EPI2ME Files tab if embedded report links are restricted.</p></html>'
    (output / 'wf-virus-longread-report.html').write_text(document)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='command', required=True)
    s = sub.add_parser('sample')
    s.add_argument('sample'); s.add_argument('mode'); s.add_argument('metadata', type=Path); s.add_argument('analysis', type=Path)
    a = sub.add_parser('aggregate')
    a.add_argument('settings'); a.add_argument('summaries', nargs='+', type=Path)
    args = p.parse_args()
    if args.command == 'sample':
        Path('summary.json').write_text(json.dumps(sample_summary(args.sample, args.mode, args.metadata, args.analysis), indent=2))
    else:
        aggregate([json.loads(p.read_text()) for p in args.summaries], json.loads(base64.b64decode(args.settings)), Path('.'))

if __name__ == '__main__':
    main()
