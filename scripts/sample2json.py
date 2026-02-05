#!/usr/bin/env python3
"""
sample2json.py

Scan a directory of FASTQ files and produce a samples JSON in the format
used by the Snakefile: {sample: {sample_type: {assay: [r1,r2]}}}

Expected filename pattern (default):
  <sample>_<sampleType>_<assay>_R1(_001).fastq(.gz)
You can override the regex with --pattern.

Example:
  scripts/sample2json.py -i data/ -o snakemake/samples.json

"""
import argparse
import json
import re
from pathlib import Path


DEFAULT_PATTERN = r'(?P<sample>[^_]+)_(?P<type>[^_]+)_(?P<assay>[^_]+)_(?P<read>R[12])'


def find_fastqs(indir, pattern):
    p = Path(indir)
    rx = re.compile(pattern)
    files = list(p.rglob('*.fastq')) + list(p.rglob('*.fastq.gz')) + list(p.rglob('*.fq')) + list(p.rglob('*.fq.gz'))
    files = [f for f in files if f.is_file()]
    samples = {}
    for f in files:
        m = rx.search(f.name)
        if not m:
            continue
        sample = m.group('sample')
        stype = m.group('type')
        assay = m.group('assay')
        read = m.group('read')
        samples.setdefault(sample, {}).setdefault(stype, {}).setdefault(assay, [])
        samples[sample][stype][assay].append(str(f))
    # sort entries so R1 comes before R2 when possible
    for sample in samples:
        for stype in samples[sample]:
            for assay in samples[sample][stype]:
                samples[sample][stype][assay] = sorted(samples[sample][stype][assay])
    return samples


def main():
    ap = argparse.ArgumentParser(description='Create samples JSON for Snakemake from FASTQ filenames')
    ap.add_argument('-i', '--input-dir', required=True, help='Input directory to scan for FASTQ files')
    ap.add_argument('-o', '--output', required=True, help='Output JSON file path')
    ap.add_argument('--pattern', default=DEFAULT_PATTERN, help='Regex pattern with named groups: sample,type,assay,read')
    ap.add_argument('--relative-to', help='If set, make paths relative to this directory')
    args = ap.parse_args()

    samples = find_fastqs(args.input_dir, args.pattern)
    out = Path(args.output)
    if args.relative_to:
        base = Path(args.relative_to).resolve()
        # make paths relative
        for s in samples:
            for t in samples[s]:
                for a in samples[s][t]:
                    samples[s][t][a] = [str(Path(x).resolve().relative_to(base)) for x in samples[s][t][a]]

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w') as fh:
        json.dump(samples, fh, indent=2)
    print(f'Wrote samples JSON to {out}')


if __name__ == '__main__':
    main()
