#!/usr/bin/env python3
"""Prepare a reference-specific snpEff database from GFF3; reuse verified caches."""
import argparse
import fcntl
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

POLICY = 1


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fasta(path):
    sequences = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                tokens = line[1:].split()
                if not tokens or tokens[0] in sequences:
                    raise ValueError('FASTA has an empty or duplicate sequence ID')
                name = tokens[0]
                sequences[name] = []
            elif name is None:
                raise ValueError('Sequence data precede the first FASTA header')
            else:
                seq = line.upper()
                if re.search('[^ACGTRYSWKMBDHVN]', seq):
                    raise ValueError(f'Invalid DNA symbols in reference sequence {name}')
                sequences[name].append(seq)
    result = {key: ''.join(value) for key, value in sequences.items()}
    if not result or any(not value for value in result.values()):
        raise ValueError('Reference FASTA contains no sequence or an empty record')
    return result


def reference_hash(sequences):
    digest = hashlib.sha256()
    for key, value in sorted(sequences.items()):
        digest.update((key + '\n' + value + '\n').encode())
    return digest.hexdigest()


def attributes(text):
    result = {}
    if text == '.':
        return result
    for item in text.split(';'):
        if not item:
            continue
        if '=' not in item:
            raise ValueError('GFF3 attributes must use key=value syntax')
        key, value = item.split('=', 1)
        if key in result:
            raise ValueError(f'Duplicate GFF3 attribute {key}')
        result[key] = value
    return result


def validate_gff(path, sequences):
    """Check IDs/coordinates/parentage. Only infer direct gene->CDS transcripts."""
    rows, comments, identities, embedded = [], [], {}, []
    in_fasta = False
    for number, raw in enumerate(Path(path).read_text().splitlines(), 1):
        if raw == '##FASTA':
            in_fasta = True
            continue
        if in_fasta:
            embedded.append(raw)
            continue
        if not raw or raw.startswith('#'):
            if raw.startswith('##gff-version') and raw.split() != ['##gff-version', '3']:
                raise ValueError('Only GFF3 is supported')
            if raw.startswith('##sequence-region'):
                parts = raw.split()
                if len(parts) != 4 or parts[1] not in sequences:
                    raise ValueError(f'GFF sequence-region at line {number} must use an exact reference ID, including its version: found {parts[1:]!r}, reference IDs {list(sequences)}; no automatic renaming is performed')
                if int(parts[2]) < 1 or int(parts[3]) > len(sequences[parts[1]]) or int(parts[2]) > int(parts[3]):
                    raise ValueError(f'Invalid GFF sequence-region bounds at line {number}')
            if raw:
                comments.append(raw)
            continue
        fields = raw.split('\t')
        if len(fields) != 9:
            raise ValueError(f'GFF line {number} must have nine tab-separated columns')
        chrom, source, kind, start, end, score, strand, phase, attrs = fields
        if chrom not in sequences:
            raise ValueError(f'GFF line {number}: sequence ID {chrom!r} does not exactly match reference IDs {list(sequences)}. Accession versions are significant; no renaming is performed.')
        try:
            lo, hi = int(start), int(end)
        except ValueError:
            raise ValueError(f'GFF line {number}: coordinates must be integers')
        if not 1 <= lo <= hi <= len(sequences[chrom]):
            raise ValueError(f'GFF line {number}: coordinates {lo}..{hi} exceed reference {chrom} length {len(sequences[chrom])}')
        if strand not in ('+', '-', '.', '?'):
            raise ValueError(f'GFF line {number}: invalid strand')
        if kind == 'CDS' and (strand not in ('+', '-') or phase not in ('0', '1', '2')):
            raise ValueError(f'GFF line {number}: CDS requires +/- strand and phase 0, 1 or 2')
        attrs = attributes(attrs)
        if kind == 'CDS' and unquote(attrs.get('transl_table','1')) != '1':
            raise ValueError(f'GFF line {number}: automatic mode uses the Standard genetic code; non-standard transl_table requires a separately configured database')
        row = {'fields': fields, 'attrs': attrs, 'start': lo, 'end': hi, 'line': number}
        identity = unquote(attrs['ID']) if 'ID' in attrs else None
        if kind in ('gene', 'mRNA', 'transcript') and not identity:
            raise ValueError(f'GFF line {number}: {kind} requires an ID')
        if identity:
            previous = identities.get(identity)
            # Repeated CDS IDs are valid for discontinuous coding features.
            if previous and (kind != 'CDS' or previous['fields'][2] != 'CDS'
                             or previous['fields'][0] != chrom or previous['fields'][6] != strand
                             or previous['attrs'].get('Parent') != attrs.get('Parent')):
                raise ValueError(f'GFF line {number}: duplicate/inconsistent ID {identity}')
            identities[identity] = row
        rows.append(row)
    if embedded:
        with tempfile.TemporaryDirectory() as tmp:
            temp = Path(tmp) / 'embedded.fa'
            temp.write_text('\n'.join(embedded) + '\n')
            for name, seq in fasta(temp).items():
                if sequences.get(name) != seq:
                    raise ValueError(f'Embedded GFF FASTA sequence {name} differs from the selected reference')
    if not any(row['fields'][2] == 'CDS' for row in rows):
        raise ValueError('GFF has no CDS features; provide a coding gene annotation for protein effect prediction')
    transcript_kinds = ('mRNA', 'transcript')
    synthetic, changes = {}, []
    for row in rows:
        kind = row['fields'][2]
        parents = [unquote(p) for p in row['attrs'].get('Parent', '').split(',') if p]
        if kind in (*transcript_kinds, 'CDS', 'exon') and not parents:
            raise ValueError(f'GFF line {row["line"]}: {kind} requires Parent links')
        for parent in parents:
            if parent not in identities:
                raise ValueError(f'GFF line {row["line"]}: Parent {parent} is missing')
            ancestor = identities[parent]
            if (ancestor['fields'][0] != row['fields'][0] or row['start'] < ancestor['start']
                    or row['end'] > ancestor['end']
                    or (ancestor['fields'][6] in ('+', '-') and row['fields'][6] in ('+', '-')
                        and ancestor['fields'][6] != row['fields'][6])):
                raise ValueError(f'GFF line {row["line"]}: Parent sequence, bounds or strand mismatch')
            if kind in transcript_kinds and ancestor['fields'][2] != 'gene':
                raise ValueError(f'GFF line {row["line"]}: transcript Parent must be a gene')
            if kind == 'exon' and ancestor['fields'][2] not in transcript_kinds:
                raise ValueError(f'GFF line {row["line"]}: exon Parent must be a transcript')
            if kind == 'CDS' and ancestor['fields'][2] not in (*transcript_kinds, 'gene'):
                raise ValueError(f'GFF line {row["line"]}: CDS Parent must be a transcript or gene')
        if kind == 'CDS':
            adjusted = []
            for raw_parent in row['attrs']['Parent'].split(','):
                parent = unquote(raw_parent)
                if identities[parent]['fields'][2] == 'gene':
                    tid = f'wfvirus_transcript_{len(synthetic)+1}'
                    if parent in synthetic:
                        tid = synthetic[parent]['attrs']['ID']
                    else:
                        if tid in identities:
                            raise ValueError(f'Reserved generated transcript ID {tid} already exists')
                        fields = identities[parent]['fields'].copy()
                        fields[1], fields[2], fields[7] = 'wf-virus-longread', 'mRNA', '.'
                        synthetic[parent] = {'fields': fields, 'attrs': {'ID': tid, 'Parent': raw_parent}}
                        changes.append(f'Created transcript {tid} for directly gene-parented CDS of {parent}')
                    adjusted.append(tid)
                else:
                    adjusted.append(raw_parent)
            row['attrs']['Parent'] = ','.join(adjusted)
    ordered = sorted([*synthetic.values(), *rows], key=lambda row: {'gene':0,'mRNA':1,'transcript':1}.get(row['fields'][2],2))
    result = ['##gff-version 3'] + [c for c in comments if not c.startswith('##gff-version')]
    for row in ordered:
        fields = row['fields'].copy()
        fields[8] = ';'.join(f'{key}={value}' for key, value in row['attrs'].items())
        result.append('\t'.join(fields))
    return '\n'.join(result) + '\n', {
        'features': dict(Counter(row['fields'][2] for row in ordered)),
        'unannotated_sequences': sorted(set(sequences) - {row['fields'][0] for row in rows}),
        'normalization': changes,
        'checks': ['exact sequence IDs including versions', 'reference coordinate bounds',
                   'CDS strand and phase', 'feature IDs and Parent links', 'parent bounds and strands'],
        'embedded_fasta_checked': bool(embedded),
        'independent_cds_check': False, 'independent_protein_check': False,
        'identity_note': 'Names and coordinate structure verified; GFF without embedded FASTA cannot independently prove its source genome sequence.',
    }


def tool_version():
    result = subprocess.run(['snpEff', '-version'], capture_output=True, text=True, check=True)
    return next((line.strip() for line in (result.stdout + '\n' + result.stderr).splitlines()
                 if 'SnpEff' in line), 'unknown')


def verify_bundle(bundle, refhash, version):
    try:
        metadata = json.loads((bundle / 'metadata.json').read_text())
        if (metadata['reference_hash'] != refhash or metadata['tool_version'] != version
                or metadata['policy'] != POLICY):
            return None
        for name, digest in metadata['files'].items():
            p = bundle / name
            if not p.resolve().is_relative_to(bundle.resolve()) or not p.is_file() or sha(p) != digest:
                return None
        required = ['snpEff.config', f'data/{metadata["db_id"]}/genes.gff',
                    f'data/{metadata["db_id"]}/sequences.fa', f'data/{metadata["db_id"]}/snpEffectPredictor.bin']
        if any(name not in metadata['files'] for name in required):
            return None
        return metadata
    except (OSError, ValueError, KeyError, TypeError):
        return None


def prepare(reference, gff, cache, output, version=None):
    sequences = fasta(reference)
    refhash = reference_hash(sequences)
    version = version or tool_version()
    if version == 'unknown':
        raise ValueError('Unable to determine the snpEff version')
    cache, output = Path(cache).resolve(), Path(output)
    cache.mkdir(parents=True, exist_ok=True)
    original_hash = None
    if gff:
        normalized, validation = validate_gff(gff, sequences)
        original_hash = sha(gff)
        key = hashlib.sha256(json.dumps([POLICY, refhash, original_hash, version], separators=(',', ':')).encode()).hexdigest()
        dbid = 'wfvirus_' + key[:24]
    else:
        matches = [(p, verify_bundle(p, refhash, version)) for p in cache.glob('wfvirus_*') if p.is_dir()]
        matches = [(p, m) for p, m in matches if m]
        if len(matches) != 1:
            raise ValueError('No unique verified database matches this reference and snpEff version. Select the matching GFF3 file to build or identify a database; never remove accession version suffixes to force a match.')
        dbid = matches[0][0].name
    target = cache / dbid
    with (cache / ('.' + dbid + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        metadata = verify_bundle(target, refhash, version)
        reused = bool(metadata and (not gff or metadata.get('gff_sha256') == original_hash))
        if not reused:
            if not gff:
                raise ValueError('Cached database failed integrity checks; provide the matching GFF3 to rebuild')
            if target.exists():
                raise ValueError(f'Cached database is incomplete or changed: {target}. Move it out of the cache before rebuilding.')
            with tempfile.TemporaryDirectory(prefix='.build-', dir=cache) as tmp:
                bundle = Path(tmp) / 'bundle'
                data = bundle / 'data' / dbid
                data.mkdir(parents=True)
                (bundle / 'snpEff.config').write_text(f'data.dir = ./data\n{dbid}.genome : User reference\n')
                (data / 'genes.gff').write_text(normalized)
                (data / 'sequences.fa').write_text(''.join(f'>{name}\n{seq}\n' for name, seq in sequences.items()))
                shutil.copyfile(gff, bundle / 'source.gff')
                command = ['snpEff', 'build', '-gff3', '-v', '-noCheckCds', '-noCheckProtein',
                           '-c', str(bundle / 'snpEff.config'), '-dataDir', str(bundle / 'data'), dbid]
                result = subprocess.run(command, capture_output=True, text=True)
                log = result.stdout + '\n' + result.stderr
                (bundle / 'build.log').write_text(log)
                binary = data / 'snpEffectPredictor.bin'
                if result.returncode or not binary.is_file() or binary.stat().st_size == 0 or re.search(r'^\s*(ERROR|FATAL)(?:\b|_)', log, re.M):
                    Path('snpeff_build_failed.log').write_text(log)
                    raise ValueError('snpEff database build failed; inspect snpeff_build_failed.log in the task work directory')
                metadata = {'policy': POLICY, 'db_id': dbid, 'reference_hash': refhash,
                            'reference_file_sha256': sha(reference), 'gff_sha256': original_hash,
                            'tool_version': version, 'codon_table': 'Standard',
                            'reference_sequences': {name: {'length': len(seq), 'sha256': hashlib.sha256(seq.encode()).hexdigest()} for name, seq in sequences.items()},
                            'validation': validation,
                            'warnings': [line.strip() for line in log.splitlines() if 'WARNING' in line][:200],
                            'files': {str(p.relative_to(bundle)): sha(p) for p in bundle.rglob('*') if p.is_file()}}
                (bundle / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
                bundle.rename(target)
        shutil.copytree(target, output)
    provenance = dict(metadata, mode='automatic_gff', cache_reused=reused)
    (output / 'metadata.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--gff', type=Path)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=Path('prepared_db'))
    args = parser.parse_args()
    try:
        metadata = prepare(args.reference, args.gff, args.cache, args.output)
        print(f"snpEff DB {metadata['db_id']}: {'reused' if metadata['cache_reused'] else 'built'}")
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, f'Annotation database preparation failed: {error}\n')


if __name__ == '__main__':
    main()
