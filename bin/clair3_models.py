#!/usr/bin/env python3
"""Prepare exact-match Clair3 1.x checkpoints from Rerio or a custom directory."""
import argparse
import hashlib
import json
import shutil
import sys
import tarfile
import tempfile
import time
from pathlib import Path, PurePosixPath
from urllib.error import URLError
from urllib.request import urlopen

from inputs import model as custom_model


def validate_checkpoint(root):
    root = Path(root)
    for prefix in ('pileup', 'full_alignment'):
        index = root / f'{prefix}.index'
        shards = list(root.glob(f'{prefix}.data-*'))
        if not index.is_file() or index.stat().st_size == 0 or not shards or any(not p.is_file() or p.stat().st_size == 0 for p in shards):
            raise ValueError(f'Invalid Clair3 1.x checkpoint: missing/empty {prefix} files in {root}')
    return root


def copy_checkpoint(source, output):
    validate_checkpoint(source)
    shutil.copytree(source, output)
    validate_checkpoint(output)


def extract_archive(archive, destination):
    """Extract ordinary files/directories only, refusing traversal and links."""
    with tarfile.open(archive, 'r:gz') as tar:
        members = tar.getmembers()
        for item in members:
            path = PurePosixPath(item.name)
            if path.is_absolute() or '..' in path.parts or '\\' in item.name or not (item.isfile() or item.isdir()):
                raise ValueError(f'Unsafe model archive member: {item.name}')
        for item in members:
            target = destination / item.name
            if item.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with tar.extractfile(item) as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)


def download(url, destination):
    if not url.startswith('https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/'):
        raise ValueError('Automatic models must use the official HTTPS model CDN')
    for attempt in range(3):
        try:
            digest = hashlib.sha256()
            with urlopen(url, timeout=120) as src, destination.open('wb') as dst:
                if not src.geturl().startswith('https://'):
                    raise ValueError('Refusing a non-HTTPS model redirect')
                for chunk in iter(lambda: src.read(1024 * 1024), b''):
                    digest.update(chunk)
                    dst.write(chunk)
            return digest.hexdigest()
        except (OSError, URLError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def prepare(basecaller, catalog_path, output, provenance, custom_root=None, mapping=None, bundled_roots=None):
    output, provenance = Path(output), Path(provenance)
    if output.exists():
        raise ValueError(f'Refusing to overwrite model output: {output}')
    catalog = json.loads(Path(catalog_path).read_text())
    info = {'basecaller_cfg': basecaller, 'clair3_version': '1.2.0', 'checkpoint_format': 'tensorflow-v1'}
    with tempfile.TemporaryDirectory(prefix='.clair3-', dir=output.parent) as tmp:
        tmp = Path(tmp)
        if custom_root:
            metadata = tmp / 'basecaller.json'
            metadata.write_text(json.dumps({'selected': None if basecaller == 'unknown' else basecaller}))
            chosen = Path(custom_model(custom_root, mapping, metadata))
            info.update(source='custom', model_id=chosen.name)
        else:
            entry = catalog['models'].get(basecaller)
            if not entry:
                raise ValueError(f'No exact automatic Clair3 model for {basecaller}. Select the correct basecalling model or provide a compatible Clair3 1.x model_dir (and model_map for a library). No approximate fallback is used.')
            model_id = entry['model_id']
            info.update(model_id=model_id, catalog_revision=catalog['rerio_revision'], url=entry['url'])
            chosen = None
            for root in bundled_roots or ():
                candidate = Path(root) / model_id
                if candidate.exists():
                    chosen = validate_checkpoint(candidate)
                    info.update(source='bundled')
                    break
            if chosen is None:
                archive = tmp / 'model.tar.gz'
                print(f'Downloading Clair3 model {model_id} for {basecaller}', file=sys.stderr)
                checksum = download(entry['url'], archive)
                if entry.get('sha256') and checksum != entry['sha256']:
                    raise ValueError(f'Model archive checksum mismatch for {model_id}')
                extracted = tmp / 'extracted'
                extracted.mkdir()
                extract_archive(archive, extracted)
                candidates = sorted({p.parent for p in extracted.rglob('pileup.index')})
                if len(candidates) != 1:
                    raise ValueError('Expected exactly one Clair3 checkpoint directory in downloaded archive')
                chosen = validate_checkpoint(candidates[0])
                info.update(source='download', archive_sha256=checksum,
                            checksum_verified=bool(entry.get('sha256')))
        prepared = tmp / 'prepared'
        copy_checkpoint(chosen, prepared)
        info['checkpoint_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in sorted(prepared.iterdir()) if p.is_file()}
        prepared.rename(output)
    provenance.write_text(json.dumps(info, indent=2) + '\n')
    return info


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--basecaller', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--output', default='clair3_model')
    p.add_argument('--provenance', default='model_selection.json')
    p.add_argument('--custom-root', default='')
    p.add_argument('--mapping', required=True)
    args = p.parse_args()
    try:
        prepare(args.basecaller, args.catalog, args.output, args.provenance,
                args.custom_root or None, args.mapping,
                bundled_roots=('/opt/clair3-models', '/opt/models', '/opt/conda/bin/models'))
    except (ValueError, OSError, tarfile.TarError, KeyError) as e:
        p.exit(2, f'Clair3 model preparation failed: {e}\n')


if __name__ == '__main__':
    main()
