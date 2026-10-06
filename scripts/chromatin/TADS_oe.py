#!/usr/bin/env python
# coding=utf-8

import re


def aftCounter(affinity, aftlist):
    if affinity > 4.07:
        aftlist[0] += 1   # strong
    elif affinity > 2.25:
        aftlist[1] += 1   # middle
    else:
        aftlist[2] += 1   # weak
    return aftlist


# ============================================================
# 1. Read all CBS and count the total number of strong, middle, and weak
# ============================================================

all_affinity = [0, 0, 0]
ctcf_sites = {}

count = 0

for eachLine in open('../12w_hg18_motif.csv'):

    count += 1

    if count > 1:
        each = eachLine.split(',')

        chrom = each[1]
        start = int(each[2])
        end = int(each[3])
        center = int((start + end) / 2)

        score = float(each[6])

      
        # All CBS are counted into the overall affinity background
        all_affinity = aftCounter(
            score,
            all_affinity
        )

        # Keep the original valid sequence selection
        if re.search('None', each[5]) is None:

            ctcf_sites[(chrom, center)] = score

        else:
            print(eachLine)


# ============================================================
# 2. Search for CBS near the TAD boundary
# ============================================================


BOUNDARY_SIZE = 25000

n_affinity = [0, 0, 0]


for eachLine in open(
    '41586_2012_BFnature11082_MOESM330_ESM_hg18-domain_IMR90.csv'
):

    each = eachLine.split()[0].split(',')

    chrom = each[0]
    start = int(each[1])
    end = int(each[2])


    # --------------------------------------------------------
    # left TAD boundary
    # Scan 25 kb into the TAD from the start
    # [start, start + 25 kb)
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
    # right TAD boundary
    # Scan 50 kb from the end into the TAD
    # (end - 25 kb, end]
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
# 3. Calculate observed proportion
# ============================================================

observed = [
    n_affinity[i] / all_affinity[i]
    for i in range(3)
]


# ============================================================
# 4. Calculate the overall expected proportion
# ============================================================

used_total = sum(n_affinity)
all_total = sum(all_affinity)

expected = used_total / all_total


# ============================================================
# 5. Calculate O/E
# ============================================================

oe = [
    observed[i] / expected
    for i in range(3)
]


# ============================================================
# 6. output
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
