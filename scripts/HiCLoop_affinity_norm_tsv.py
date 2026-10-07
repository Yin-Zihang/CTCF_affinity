#!/usr/bin/env python3
# coding=utf-8
# Writes a Hi-C loop normalization TSV in the current working directory by default.

import argparse
import csv
import re
from pathlib import Path

parser = argparse.ArgumentParser(description='Summarize affinity classes in CTCF Hi-C loops.')
parser.add_argument('--sites-file', default='final_split_to_12_groups_motif.csv')
parser.add_argument('--loops-file', default='K562_loop_5k_filter.txt')
parser.add_argument('--output', default='HiCLoop_affinity_normalization.tsv',
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
# 1. Read CTCF sites.
# ============================================================

ctcf_sites = {}

count = 0
for eachLine in open(args.sites_file):
    count += 1

    if count > 1:
        each = eachLine.split(',')

        l = each[4]
        start = int(each[2])

        # Preserve the original start + 1 rule for 41-bp sequences.
        if l == '41':
            start += 1

        end = int(each[3])
        center = int((start + end) / 2)

        score = float(each[6])
        seq = each[5]
        direction = each[9]

        if re.search('None', each[5]) is None:
            ctcf_sites[(each[1], center)] = (
                start,
                end,
                score,
                seq,
                direction
            )
        else:
            print(eachLine)


# ============================================================
# 2. Hi-C loop-strength bin boundaries.
# ============================================================

tloopStrength = [
    6.741221875,
    9.94278506,
    13.99217675,
    20.039869500000002
]


# ============================================================
# 3. Initialize counters.
# ============================================================

# Five loop-strength groups, each with strong / middle / weak counts.
n_affinity = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]

# Track unique CBS participating in any loop, by affinity class.
a1 = {}
a2 = {}
a3 = {}

nondup_affinity = [a1, a2, a3]

# Number of loops in each loop-strength group.
loopS = [0, 0, 0, 0, 0]


# ============================================================
# 4. Scan Hi-C loops.
# ============================================================

for eachLine in open(args.loops_file):

    each = eachLine.split()

    chromL = each[0]
    startL = int(each[1])
    endL = int(each[2])

    chromR = each[3]
    startR = int(each[4])
    endR = int(each[5])

    loopstrength = float(each[6])

    # --------------------------------------------------------
    # Preserve the original five loop-strength groups.
    # --------------------------------------------------------

    if loopstrength <= tloopStrength[0]:
        strength = 0
    elif loopstrength <= tloopStrength[1]:
        strength = 1
    elif loopstrength <= tloopStrength[2]:
        strength = 2
    elif loopstrength <= tloopStrength[3]:
        strength = 3
    else:
        strength = 4

    # Count loops before filtering interchromosomal interactions, as originally.
    loopS[strength] += 1

    if chromR != chromL:
        continue


    # ========================================================
    # Left anchor: retain only plus-strand CBS.
    # ========================================================

    for i in range(endL - startL):

        key = (chromL, startL + i)

        if key in ctcf_sites:

            CBSstart = ctcf_sites[key][0]
            CBSend = ctcf_sites[key][1]
            score = ctcf_sites[key][2]
            direction = ctcf_sites[key][4]

            # Exclude minus-strand CBS from the left anchor.
            if direction == '-':
                continue

            if score > 4.07:
                flag = 0
            elif score > 2.25:
                flag = 1
            else:
                flag = 2

            # occurrence count
            n_affinity[strength] = aftCounter(
                score,
                n_affinity[strength]
            )

            # unique CBS
            nondup_affinity[flag][
                chromL, CBSstart, CBSend
            ] = 0


    # ========================================================
    # Right anchor: retain only minus-strand CBS.
    # ========================================================

    for i in range(endR - startR):

        key = (chromR, startR + i)

        if key in ctcf_sites:

            CBSstart = ctcf_sites[key][0]
            CBSend = ctcf_sites[key][1]
            score = ctcf_sites[key][2]
            direction = ctcf_sites[key][4]

            # Exclude plus-strand CBS from the right anchor.
            if direction == '+':
                continue

            if score > 4.07:
                flag = 0
            elif score > 2.25:
                flag = 1
            else:
                flag = 2

            # occurrence count
            n_affinity[strength] = aftCounter(
                score,
                n_affinity[strength]
            )

            # Preserve the original key; chromL equals chromR after filtering.
            nondup_affinity[flag][
                chromL, CBSstart, CBSend
            ] = 0


# ============================================================
# 5. Calculate normalized affinity-class occurrence rates.
# ============================================================

affinity_norm = []

for i in range(5):

    affinity_norm.append([
        1000 * n_affinity[i][0]
        / len(nondup_affinity[0])
        / loopS[i],

        1000 * n_affinity[i][1]
        / len(nondup_affinity[1])
        / loopS[i],

        1000 * n_affinity[i][2]
        / len(nondup_affinity[2])
        / loopS[i]
    ])


# ============================================================
# 6. Write a machine-readable table and print the original summary.
# The unique-CBS denominator covers all five loop-strength groups.
# The loop denominator includes interchromosomal interactions, as originally.
with open(output_path, 'w', newline='') as fjob:
    writer = csv.writer(fjob, delimiter='\t')
    writer.writerow(['loop_strength_group', 'affinity_class',
                     'n_cbs_occurrences', 'n_unique_cbs_across_all_groups',
                     'n_loops', 'affinity_norm_per_1000_loops'])
    for i in range(5):
        for j, affinity_class in enumerate(('strong', 'middle', 'weak')):
            writer.writerow([i + 1, affinity_class,
                             n_affinity[i][j], len(nondup_affinity[j]),
                             loopS[i], affinity_norm[i][j]])

print('Wrote TSV:', output_path)
# ============================================================

print('n_affinity:')
print(n_affinity)

print('unique affinity CBS:')
print([
    len(nondup_affinity[0]),
    len(nondup_affinity[1]),
    len(nondup_affinity[2])
])

print('loopS:')
print(loopS)

print('affinity_norm:')
print(affinity_norm)

print(
    'g1=c(%s, %s, %s),' %
    tuple(affinity_norm[0])
)

print(
    'g2=c(%s, %s, %s),' %
    tuple(affinity_norm[1])
)

print(
    'g3=c(%s, %s, %s),' %
    tuple(affinity_norm[2])
)

print(
    'g4=c(%s, %s, %s),' %
    tuple(affinity_norm[3])
)

print(
    'g5=c(%s, %s, %s),' %
    tuple(affinity_norm[4])
)
