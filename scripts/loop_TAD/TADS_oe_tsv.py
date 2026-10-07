#!/usr/bin/env python3
# coding=utf-8
# Writes TAD enrichment as a TSV file in the current working directory by default.

import argparse
import csv
import re
from pathlib import Path

parser = argparse.ArgumentParser(description='Calculate CTCF affinity-class enrichment at TAD boundaries.')
parser.add_argument('--sites-file', default='12w_hg18_motif.csv')
parser.add_argument('--domains-file', default='41586_2012_BFnature11082_MOESM330_ESM_hg18-domain_IMR90.csv')
parser.add_argument('--output', default='TAD_enrichment.tsv',
                    help='Output TSV path; relative paths use the current working directory.')
args = parser.parse_args()
output_path = Path(args.output).expanduser().resolve()
output_path.parent.mkdir(parents=True, exist_ok=True)


def aftCounter(affinity, aftlist):
    if affinity > 4.07:
        aftlist[0] += 1   # strong
    elif affinity > 2.25:
        aftlist[1] += 1   # middle
    else:
        aftlist[2] += 1   # weak
    return aftlist


# ============================================================
# 1. Read all CBS and count the background for each affinity class.
# ============================================================

all_affinity = [0, 0, 0]
ctcf_sites = {}

count = 0

for eachLine in open(args.sites_file):

    count += 1

    if count > 1:
        each = eachLine.split(',')

        chrom = each[1]
        start = int(each[2])
        end = int(each[3])
        center = int((start + end) / 2)

        score = float(each[6])

        # Preserve the original background: count all CBS before filtering.
        all_affinity = aftCounter(
            score,
            all_affinity
        )

        # Preserve the original valid-sequence filter for site lookup.
        if re.search('None', each[5]) is None:

            ctcf_sites[(chrom, center)] = score

        else:
            print(eachLine)


# ============================================================
# 2. Count CBS in inward windows at TAD boundaries.
# ============================================================

# The actual code scans 25 kb inward from each end of a domain.
BOUNDARY_SIZE = 25000

n_affinity = [0, 0, 0]


for eachLine in open(args.domains_file):

    each = eachLine.split()[0].split(',')

    chrom = each[0]
    start = int(each[1])
    end = int(each[2])


    # --------------------------------------------------------
    # Left boundary: scan 25 kb into the domain from start.
    # [start, start + 50 kb)
    # --------------------------------------------------------

    for i in range(BOUNDARY_SIZE):

        key = (chrom, start + i)

        if key in ctcf_sites:

            score = ctcf_sites[key]

            n_affinity = aftCounter(
                score,
                n_affinity
            )


    # --------------------------------------------------------
    # Right boundary: scan 25 kb into the domain from end.
    # (end - 50 kb, end]
    # --------------------------------------------------------

    for i in range(BOUNDARY_SIZE):

        key = (chrom, end - i)

        if key in ctcf_sites:

            score = ctcf_sites[key]

            n_affinity = aftCounter(
                score,
                n_affinity
            )


# ============================================================
# 3. Calculate the observed fraction within each affinity class.
# ============================================================

observed = [
    n_affinity[i] / all_affinity[i]
    for i in range(3)
]


# ============================================================
# 4. Calculate the expected fraction across all affinity classes.
# ============================================================

used_total = sum(n_affinity)
all_total = sum(all_affinity)

expected = used_total / all_total


# ============================================================
# 5. Calculate observed-to-expected ratios.
# ============================================================

oe = [
    observed[i] / expected
    for i in range(3)
]


# ============================================================
# 6. Write a machine-readable table and print the original summary.
with open(output_path, 'w', newline='') as fjob:
    writer = csv.writer(fjob, delimiter='\t')
    writer.writerow(['affinity_class', 'boundary_window_bp_each_end',
                     'n_boundary_cbs_occurrences', 'n_background_cbs',
                     'observed_fraction', 'expected_fraction',
                     'observed_over_expected'])
    for i, affinity_class in enumerate(('strong', 'middle', 'weak')):
        writer.writerow([affinity_class, BOUNDARY_SIZE, n_affinity[i],
                         all_affinity[i], observed[i], expected, oe[i]])

print('Wrote TSV:', output_path)
# ============================================================

print('n_affinity:')
print(n_affinity)

print('all_affinity:')
print(all_affinity)

print('Observed:')
print(observed)

print('Expected:')
print(expected)

print('O/E:')
print(oe)

print()

print(
    'Strong O/E = %s' % oe[0]
)

print(
    'Middle O/E = %s' % oe[1]
)

print(
    'Weak O/E = %s' % oe[2]
)
