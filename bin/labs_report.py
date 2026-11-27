"""EPI2ME Labs report layout, following wf-bacterial-genomes v2.0.2.

Uses the same upstream LabsReport, Tabs, DataTable, SeqSummary and EZChart
components and theme. Bacterial typing/AMR sections are replaced by viral
consensus, Clair3 and optional snpEff results.
"""
import csv
import importlib.metadata
import json
import html as escaping
import math
import tempfile
from pathlib import Path
from urllib.parse import quote

from dominate import tags as html
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd
from report_branding import apply_branding
from bokeh.models import Title


# Shared red theme for integrated report charts.
CHART_PRIMARY = '#B00020'
COVERAGE_PALETTE = {'total': CHART_PRIMARY, 'fwd': '#D86A7D', 'rev': '#680015'}
SUBSTITUTION_PALETTE = ['#FFF5F6', '#F4CCD3', '#E79AAA', '#D5647B', CHART_PRIMARY, '#720015']


class HistogramSummary(fastcat.SeqSummary):
    """Retain upstream QC plots, fixing histogram yield N50 for small samples."""

    def _draw_summary_plots(self, ldata, qdata, height):
        with Grid(columns=3):
            EZChart(fastcat.read_quality_plot(qdata, min_qual=min(4, math.floor(qdata['start'].min())),
                                              max_qual=max(30, math.ceil(qdata['end'].max())), color=CHART_PRIMARY),
                    'epi2melabs', height=height)
            EZChart(fastcat.read_length_plot(ldata, color=CHART_PRIMARY), 'epi2melabs', height=height)
            # Upstream base_yield_plot indexes past the final histogram bin when
            # the longest bin contributes >= half of the yield. Calculate the
            # weighted N50 directly rather than expanding histogram counts.
            data = ldata.sort_values('start').copy()
            yields = data['start'] * data['count']
            above = yields.iloc[::-1].cumsum().iloc[::-1]
            n50 = data.loc[above >= above.iloc[0] / 2, 'start'].iloc[-1]
            frame = pd.DataFrame({'Length': data['start'] / 1000, 'Yield': above / 1e9})
            plot = ezc.lineplot(data=frame, x='Length', y='Yield', marker=False, palette=[CHART_PRIMARY])
            plot._fig.xaxis.axis_label = 'Read length / kb'
            plot._fig.yaxis.axis_label = 'Yield above length / Gbases'
            plot._fig.add_layout(Title(text='Base yield above read length', text_font_size='1.5em'), 'above')
            plot._fig.add_layout(Title(text=f'Total yield: {int(above.iloc[0]):,} bp. N50: {int(n50):,} bp.',
                                      text_font_size='0.8em'), 'above')
            EZChart(plot, 'epi2melabs', height=height)


def table(rows):
    if rows:
        frame = pd.DataFrame(rows).fillna('-')
        # DataTable inserts raw HTML cells. Escape external values explicitly.
        frame = frame.applymap(lambda value: escaping.escape(value) if isinstance(value, str) else value)
        frame.columns = [escaping.escape(str(column)) for column in frame.columns]
        DataTable.from_pandas(frame, use_index=False)
    else:
        html.p('No data available.')


def per_sample(samples, draw):
    if not samples:
        html.p('No samples with analysis outputs.')
        return
    tabs = Tabs()
    with tabs.add_dropdown_menu():
        for sample in samples:
            with tabs.add_dropdown_tab(sample['sample']):
                draw(sample)


def variant_plots(sample):
    stats = sample.get('variants')
    if not stats:
        html.p('No variant statistics available.')
        return
    # Pair complementary substitutions, as in the reference workflow report.
    complement = str.maketrans('ACGT', 'TGCA')
    counts = pd.DataFrame(0.0, index=['A', 'C'], columns=list('ACGT'))
    for change, count in stats['substitutions'].items():
        ref, alt = change.split('>')
        if ref not in 'AC' and ref in 'GT':
            ref, alt = ref.translate(complement), alt.translate(complement)
        if ref in counts.index and alt in counts.columns:
            counts.loc[ref, alt] += count
    with Grid(columns=2):
        total = counts.to_numpy().sum()
        if total:
            plot = ezc.heatmap((counts * 100 / total).round(1).T, vmin=0, cmap=SUBSTITUTION_PALETTE)
            plot.title = dict(text='Substitution types (%)')
            plot.series[0].label.formatter = '{@[2]} %'
            EZChart(plot, 'epi2melabs')
        else:
            html.p('No substitutions were found for this sample.')
        if stats['indels']:
            frame = pd.DataFrame(stats['indels'])
            plot = ezc.histplot(data=frame['length'], x='length', discrete=True, weights=frame['count'], title='Insertion and deletion lengths', color=CHART_PRIMARY)
            plot._fig.xaxis.axis_label = 'Indel length / bases (deletions negative)'
            plot._fig.yaxis.axis_label = 'Number of sites'
            EZChart(plot, 'epi2melabs')
        else:
            html.p('No indels were found for this sample.')


def coverage_plots(sample):
    rows = sample.get('coverage', [])
    if not rows:
        html.p('No coverage data for sample, check input data and assembly quality.')
        return
    table([{'Reference': row['reference'], 'Reference span (bp)': row['length'],
            'Positions assessed': row['positions'], 'Mean depth': round(row['mean_depth'], 2),
            'Covered (%)': round(row['covered_percent'], 3),
            'Coverage at 30x (%)': round(row['at_30x_percent'], 3),
            'Coverage at 50x (%)': round(row['at_50x_percent'], 3)} for row in rows])
    for row in sorted(rows, key=lambda item: item['length'], reverse=True)[:20]:
        if row['plot']:
            data = [{'position': p['position'], 'depth': p[key], 'strand': label}
                    for p in row['plot'] for key, label in [('depth','total'),('fwd','fwd'),('rev','rev')] if key in p]
            plot = ezc.lineplot(data=pd.DataFrame(data), x='position', y='depth', hue='strand',
                                title=row['reference'], marker=False, palette=COVERAGE_PALETTE)
            plot._fig.xaxis.axis_label = 'Position along reference / bp'
            plot._fig.yaxis.axis_label = 'Mean sequencing depth'
            EZChart(plot, 'epi2melabs')


def assembly_table(sample):
    contigs = sample.get('contigs', [])
    names = {'#seq_name': 'Contig', 'length': 'Length (bp)', 'cov.': 'Coverage',
             'circ.': 'Circular', 'repeat': 'Repeat', 'mult.': 'Multiplicity'}
    if contigs:
        table([{label: row.get(key, '-') for key, label in names.items()} for row in contigs])
    else:
        html.p('No contig-level assembly data are available for this sample. Check input data and assembly quality.')


def consensus_table(sample):
    assembly = sample['assembly']
    table([{'Consensus sequences': assembly['sequences'], 'Total length (bp)': assembly['total_bases'],
            'Longest sequence (bp)': assembly['longest']}])
    html.p('Consensus sequences are generated from independently DP/QUAL/AF-filtered PASS variants with the selected '
           'bcftools consensus options and optional low-coverage mask. '
           'Report variant statistics describe the independently DP/QUAL-selected normalized VCF. '
           'Consensus uses PASS variants and can be restricted further by consensus options.')
    if 'Mask BED' in sample['files']:
        html.a('Download low-coverage mask BED', href=quote(sample['files']['Mask BED'], safe='/'))


def annotation_table(sample):
    data = sample.get('annotations')
    if data is None:
        html.p('No annotation output is available for this sample.')
        return
    table(data['preview'])


def model_table(sample):
    model = sample.get('clair3_model') or {}
    basecaller = sample['basecaller']
    table([{'Basecalling model': basecaller.get('selected') or 'Unknown',
            'Clair3 model': model.get('model_id', 'Not used'),
            'Model source': model.get('source', 'Not used')}])
    if model:
        html.p('Full checkpoint identity and hashes are recorded in the model selection JSON.')


def file_links(sample):
    with html.ul():
        for label, path in sample['files'].items():
            with html.li():
                html.a(label, href=quote(path, safe='/'))


def render(samples, settings, output):
    """Use the reference workflow's report shell, navigation and chart styling."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        params_path = root / 'params.json'
        params_path.write_text(json.dumps({escaping.escape(str(key)): escaping.escape(str(value))
                                          for key, value in settings.items()
                                          if key not in {'containers', 'container_override', 'software_versions'}}, indent=2))
        versions_path = root / 'versions.csv'
        versions = [('wf-virus-longread', settings.get('workflow_version', 'unknown'))]
        detected = settings.get('software_versions', {})
        for tool in ('Nextflow', 'fastcat', 'bgzip', 'mini_align', 'minimap2', 'samtools',
                     'bcftools', 'NanoPlot', 'Flye', 'Clair3', 'bamCoverage', 'snpEff', 'SnpSift', 'Python', 'Java'):
            used = (tool not in {'Flye','Clair3','bamCoverage','snpEff','SnpSift','Java'} or
                    (tool == 'Flye' and not settings.get('reference_based_assembly')) or
                    (tool in {'Clair3','bamCoverage'} and settings.get('reference_based_assembly')) or
                    (tool in {'snpEff','SnpSift','Java'} and settings.get('run_annotation') and settings.get('reference_based_assembly')))
            version = settings.get('nextflow_version','Not recorded') if tool == 'Nextflow' else detected.get(tool, 'Not detected' if used else 'Not used in this run')
            versions.append((tool, version))
        for package in ('ezcharts', 'pandas', 'bokeh'):

            versions.append((package, importlib.metadata.version(package)))
        with versions_path.open('w') as handle:
            csv.writer(handle).writerows(versions)
        report = labs.LabsReport('Virus Long-read Summary Report', 'wf-virus-longread',
                                 str(params_path), str(versions_path),
                                 settings.get('workflow_version', 'unknown'))
        # Keep the Labs theme while identifying this as a custom workflow.
        subtitle = report.banner.children[0].children[1]
        subtitle.children.clear()
        subtitle.add('Custom wf-virus-longread workflow maintained by klvsjeju. '
                     'Report components by Oxford Nanopore Technologies.')
        with report.add_section('Sample summary', 'Samples'):
            html.p('Analysis summary for each selected sample. Missing inputs are listed separately.')
            table([{'Barcode': row.get('metadata', {}).get('barcode') or '-', 'Alias': row['sample'],
                    'Mode': row['mode'], 'Input read count': row['reads'].get('input_read_count', row['reads']['read_count']),
                    'Retained read count': row['reads']['read_count'],
                    'Reads below minimum length': row['reads'].get('excluded_short_reads', 0),
                    'Median retained read length (bp)': row['reads']['median_length'],
                    ('Contig count' if row['mode']=='denovo' else 'Consensus sequences'): row['assembly']['sequences'],
                    'Total assembly / consensus length (bp)': row['assembly']['total_bases']}
                   for row in samples])
            missing = settings.get('missing_samples', [])
            if missing:
                html.h3('Missing input samples')
                table([{'Barcode': row['barcode'], 'Alias': row['alias'], 'Status': 'Missing input'}
                       for row in missing])
        if any(row['reads']['histograms'].get('length') and row['reads']['histograms'].get('quality') for row in samples):
            with report.add_section('Read summary', 'Reads'):
                html.p('Read quality, read length and base yield from fastcat histograms. '
                       'Use the dropdown menu to select a sample. Median lengths use histogram bins.')
                if settings.get('effective_min_read_length', 0):
                    html.p(f"De novo minimum read length: {settings['effective_min_read_length']} bp. "
                           'Fastcat histograms and retained read counts describe reads that passed the filter. '
                           'Input counts and reads below the minimum length are listed separately. '
                           'NanoPlot QC and Flye use only reads that meet the minimum length.')
                histogram_dirs, names = [], []
                for index, sample in enumerate(samples):
                    hist = sample['reads']['histograms']
                    if not hist.get('length') or not hist.get('quality'):
                        continue
                    directory = root / str(index)
                    directory.mkdir()
                    for kind, rows in hist.items():
                        with (directory / (kind + '.hist')).open('w') as handle:
                            for start, end, count in rows:
                                if kind == 'length':
                                    start, end = int(start), int(end)
                                handle.write(f'{start}\t{end}\t{count}\n')
                    histogram_dirs.append(str(directory))
                    names.append(sample['sample'])
                if names:
                    HistogramSummary(seq_summary=tuple(histogram_dirs), sample_names=tuple(names),
                                     alignment_stats=False)
                else:
                    html.p('No read histograms available. See the per-sample NanoPlot QC reports under Result files.')
        denovo = [row for row in samples if row['mode'] == 'denovo']
        reference = [row for row in samples if row['mode'] == 'reference']
        if denovo:
            with report.add_section('Assembly summary statistics', 'Assembly'):
                html.p('Flye de novo contigs and internal assembly statistics. No downstream consensus polishing is performed.')
                per_sample(denovo, assembly_table)
        if reference:
            with report.add_section('Variant calling', 'Variants'):
                html.p('Statistics from bcftools stats on the QUAL/DP-filtered normalized Clair3 VCF '
                       '(<sample>_norm.vcf.gz), after the configured variant DP/QUAL selection. Consensus DP/QUAL/AF selection is independent. Qualifying PASS and non-PASS records are included. '
                       'Complementary substitutions are grouped, as in the reference workflow.')
                table([dict({'Sample': row['sample']}, **row['variants']['counts'])
                       for row in reference if row.get('variants')])
                html.h3('Transitions and transversions')
                table([{'Sample': row['sample'], 'Transitions': row['variants']['ts'],
                        'Transversions': row['variants']['tv'], 'Ts/Tv': row['variants']['ts_tv']}
                       for row in reference if row.get('variants')])
                per_sample(reference, variant_plots)
            with report.add_section('Consensus and masking', 'Consensus'):
                per_sample(reference, consensus_table)
            if settings.get('run_annotation'):
                with report.add_section('Variant annotations', 'snpEff'):
                    per_sample(reference, annotation_table)
        with report.add_section('Genome coverage', 'Depth'):
            html.p('Read depth after alignment filtering: relative to the reference in reference mode, '
                   'or to the Flye contigs in de novo mode. Percentages use assessed positions; '
                   'Plots show mean depth in at most 500 bins per reference, for up to 20 longest references. '
                   'Total, forward and reverse read depths are shown. Depth files retain full resolution. ')
            per_sample(samples, coverage_plots)
        with report.add_section('Basecalling and Clair3 models', 'Models'):
            html.p('The basecalling model describes the existing reads. This workflow does not run basecalling.')
            per_sample(samples, model_table)
        with report.add_section('Result files', 'Files'):
            per_sample(samples, file_links)
            html.p('Open files through the EPI2ME Files tab if embedded report links are restricted.')
        with report.add_section('Execution', 'Execution'):
            table([{'Field': 'Run name', 'Value': settings.get('run_name', '-')},
                   {'Field': 'Session ID', 'Value': settings.get('session_id', '-')}])
            with html.ul():
                for label, path in [('Nextflow execution report', 'execution/report.html'),
                                    ('Execution timeline', 'execution/timeline.html'), ('Task trace', 'execution/trace.tsv')]:
                    with html.li():
                        html.a(label, href=path)
            html.p('Execution reports are finalized when Nextflow exits.')
        apply_branding(report)
        report.write(str(output / 'wf-virus-longread-report.html'))
