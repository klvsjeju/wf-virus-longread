#!/usr/bin/env python3
"""Resolve independent reference filters and combined masking controls."""
import argparse
import base64
import json
import math
import os
from pathlib import Path
import shlex
from run_options import parse_options


def number(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError('Thresholds must be finite non-negative numbers')
    return value


def integer(value):
    value = number(value)
    if not value.is_integer():
        raise argparse.ArgumentTypeError('Depth thresholds must be integers')
    return int(value)


def encode(value):
    return base64.b64encode(value.encode()).decode()


def resolve(kind, text):
    parser = argparse.ArgumentParser(prog=f'{kind} options', allow_abbrev=False)
    if kind == 'mask':
        parser.add_argument('--min-depth', type=integer, default=0)
    else:
        parser.add_argument('--min-dp', type=integer, default=30)
        parser.add_argument('--min-qual', type=number, default=10)
        if kind == 'consensus':
            parser.add_argument('--min-af', type=number, default=0.5)
    args, extra = parser.parse_known_args(shlex.split(text))
    if kind == 'variant' and extra:
        parser.error('Unknown arguments: ' + shlex.join(extra))
    if kind == 'consensus' and args.min_af > 1:
        parser.error('--min-af must be a fraction from 0 to 1 (0.5 = 50%)')
    if kind != 'variant':
        parse_options('depth' if kind == 'mask' else 'consensus', encode(shlex.join(extra)))
    values = vars(args)
    clauses = []
    for key, field in [('min_qual','QUAL'),('min_dp','FORMAT/DP'),('min_af','FORMAT/AF')]:
        value = values.get(key, 0)
        # Disabled means no comparison or missing-value requirement for that field.
        if value > 0:
            clauses += [f'{field}!="."', f'{field}>={value:g}']
    return dict(values, expression=' && '.join(clauses), extra=extra)


def prepare(variant, consensus, mask, output):
    options = {kind: resolve(kind, text) for kind, text in
               [('variant',variant),('consensus',consensus),('mask',mask)]}
    options['consensus']['require_caller_pass'] = True
    output = Path(output)
    output.write_text(json.dumps(options, indent=2)+'\n')
    env = {'VARIANT_FILTER': options['variant']['expression'],
           'CONSENSUS_FILTER': options['consensus']['expression'],
           'MASK_DEPTH': str(options['mask']['min_depth']),
           'MASK_DEPTH_OPTS_B64': encode(shlex.join(options['mask']['extra'])),
           'CONSENSUS_OPTS_B64': encode(shlex.join(options['consensus']['extra']))}
    output.with_suffix('.env').write_text(''.join(f'{key}={shlex.quote(value)}\n' for key,value in env.items()))
    return options


if __name__ == '__main__':
    prepare(base64.b64decode(os.environ.get('VARIANT_OPTS_B64','')).decode(),
            base64.b64decode(os.environ.get('CONSENSUS_OPTS_B64','')).decode(),
            base64.b64decode(os.environ.get('MASK_OPTS_B64','')).decode(),
            'result/analysis_options.json')
