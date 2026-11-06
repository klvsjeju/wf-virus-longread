#!/usr/bin/env python3
"""Execute tool arguments without shell expansion; keep workflow-owned IO fixed."""
import argparse
import base64
import shlex
import subprocess

RESERVED = {
    'flye': {'--nano-hq', '--nano-raw', '--nano-corr', '--pacbio-raw', '--pacbio-corr', '--pacbio-hifi', '--subassemblies', '--genome-size', '-g', '--threads', '-t', '--out-dir', '-o', '--help', '-h', '--version', '-v', '--resume', '--resume-from', '--stop-after', '--polish-target'},
    'depth': {'-a', '-aa', '-H', '-X', '-o', '-f', '--output', '--help', '--version'},
    'consensus': {'-f', '--fasta-ref', '-s', '--samples', '-m', '--mask', '-o', '--output', '-c', '--chain', '--help', '--version'},
}
SHORT = {'flye': {'-g', '-t', '-o'}, 'depth': {'-o', '-f'}, 'consensus': {'-f', '-s', '-m', '-o', '-c'}}

def parse_options(kind, encoded):
    args = shlex.split(base64.b64decode(encoded).decode())
    for arg in args:
        key = arg.split('=', 1)[0]
        if key in RESERVED[kind] or any(arg.startswith(k) and len(arg) > len(k) for k in SHORT[kind]):
            raise ValueError(f'{kind}: {arg} is managed by the workflow; use the dedicated parameter')
        if arg == '--' or not arg.startswith('-') and arg in {'|', ';', '&&', '>'}:
            raise ValueError('Shell operators are not supported')
    return args

def main():
    p = argparse.ArgumentParser()
    p.add_argument('kind', choices=RESERVED)
    p.add_argument('encoded')
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    try:
        extra = parse_options(a.kind, a.encoded)
    except (ValueError, UnicodeError) as e:
        p.error(str(e))
    command = a.command
    if command[:1] == ['--']:
        command = command[1:]
    # Options precede the workflow's positional input and managed flags.
    prefix = 2 if a.kind in {'depth', 'consensus'} else 1
    subprocess.run(command[:prefix] + extra + command[prefix:], check=True)

if __name__ == '__main__':
    main()
