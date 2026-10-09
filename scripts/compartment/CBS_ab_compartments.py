#!/usr/bin/env python3
"""Generate A/B compartment versus CTCF affinity and ChIP RPKM tables."""

import argparse
import csv
import importlib
import os
from collections import Counter

common = importlib.import_module('CBS_compartment_common')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cell-type', required=True, help='Example: K562, HUVEC, GM12878')
    parser.add_argument('--eigenvector', required=True, help='Four-column 1-based closed bedGraph/WIG')
    parser.add_argument('--sites-bed', required=True, help='The original group/chrom/start/end/token/EMSA BED')
    parser.add_argument('--peaks', required=True, help='Replicate-count peak TSV')
    parser.add_argument('--reads-rep1', type=int, required=True, help='Replicate 1 mapped reads')
    parser.add_argument('--reads-rep2', type=int, required=True, help='Replicate 2 mapped reads')
    parser.add_argument('--emsa-cutoff', type=float, default=3.1,
                        help='Historical highEMSA cutoff on summed peak EMSA (default: 3.1)')
    parser.add_argument('--outdir', default='13review_ab_results')
    args = parser.parse_args()
    if args.reads_rep1 <= 0 or args.reads_rep2 <= 0:
        parser.error('mapped reads must be positive')
    if not args.cell_type.replace('_', '').replace('-', '').isalnum():
        parser.error('cell type may contain only letters, numbers, underscores and hyphens')

    os.makedirs(args.outdir, exist_ok=True)
    intervals = common.read_intervals(args.eigenvector, one_based=True)
    sites, site_counts = common.read_sites(args.sites_bed)
    output = os.path.join(args.outdir, args.cell_type + '_ab_peak_affinity.tsv')
    summary = os.path.join(args.outdir, args.cell_type + '_ab_summary.tsv')
    counts = Counter()
    fields = ['id', 'cell_type', 'peak_id', 'chrom', 'start', 'end', 'ChIP', 'EMSA',
              'compartments', 'CBScenter', 'EMSAcluster', 'n_sites', 'eigenvector']

    with open(output, 'w', newline='') as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for peak_id, chrom, start, end, raw1, raw2 in common.read_peaks(args.peaks):
            counts['total_peaks'] += 1
            score = common.lookup_interval(intervals, chrom, (start + end) // 2)
            if score is None:
                counts['missing_or_invalid_eigenvector'] += 1
                continue
            if score == 0:
                counts['zero_eigenvector_excluded'] += 1
                continue
            compartment = '1.A_Compartment' if score > 0 else '2.B_Compartment'
            counts[compartment + '_peaks'] += 1
            matched = common.sites_in_peak(sites, chrom, start, end)
            if not matched:
                counts[compartment + '_without_site'] += 1
                continue
            counts[compartment + '_with_site'] += 1
            if len(matched) > 1:
                counts[compartment + '_multiple_sites'] += 1
            affinity = sum(item[1] for item in matched)
            counts['written'] += 1
            writer.writerow({'id': counts['written'], 'cell_type': args.cell_type,
                             'peak_id': peak_id, 'chrom': chrom, 'start': start,
                             'end': end, 'ChIP': common.chip_rpkm(raw1, raw2, args.reads_rep1,
                                                           args.reads_rep2, end - start),
                             'EMSA': affinity, 'compartments': compartment,
                             'CBScenter': '_'.join(str(item[0]) for item in matched) + '_',
                             'EMSAcluster': 'highEMSA' if affinity > args.emsa_cutoff else 'lowEMSA',
                             'n_sites': len(matched), 'eigenvector': score})

    counts.update(site_counts)
    with open(summary, 'w', newline='') as destination:
        writer = csv.writer(destination, delimiter='\t')
        writer.writerow(['cell_type', 'metric', 'value'])
        for metric in sorted(counts):
            writer.writerow([args.cell_type, metric, counts[metric]])
    print('Processed table: %s' % os.path.abspath(output))
    print('Summary table: %s' % os.path.abspath(summary))
    print('Peaks written: %d' % counts['written'])


if __name__ == '__main__':
    main()
