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
# 1. read CTCF sites
# ============================================================

ctcf_sites = {}

count = 0
for eachLine in open(
    '/run/media/guoya/diska/CTCF_IMP/11_split_to_12_groups/fimo/final_split_to_12_groups_motif.csv'
):
    count += 1

    if count > 1:
        each = eachLine.split(',')

        l = each[4]
        start = int(each[2])

        # 保持原代码逻辑
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
# 2. Hi-C loop strength group threshold
# ============================================================

tloopStrength = [
    6.741221875,
    9.94278506,
    13.99217675,
    20.039869500000002
]


# ============================================================
# 3. initialize
# ============================================================

# 5 loop strength groups
# calculate by group: strong / middle / weak
n_affinity = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]

#Record the unique CBS that have participated in the loop in each affinity group
a2 = {}
a3 = {}

nondup_affinity = [a1, a2, a3]

# The number of loops in each loop strength group
loopS = [0, 0, 0, 0, 0]


# ============================================================
# 4. scan Hi-C loops
# ============================================================

for eachLine in open('K562_loop_5k_filter.txt'):

    each = eachLine.split()

    chromL = each[0]
    startL = int(each[1])
    endL = int(each[2])

    chromR = each[3]
    startR = int(each[4])
    endR = int(each[5])

    loopstrength = float(each[6])

    # --------------------------------------------------------
    # Divide into 5 strength groups based on the original code
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

    # First count the loops, then filter out the inter-chromosomal loops
    loopS[strength] += 1

    if chromR != chromL:
        continue


    # ========================================================
    # left anchor
    # keep only + strand CBS
    # ========================================================

    for i in range(endL - startL):

        key = (chromL, startL + i)

        if key in ctcf_sites:

            CBSstart = ctcf_sites[key][0]
            CBSend = ctcf_sites[key][1]
            score = ctcf_sites[key][2]
            direction = ctcf_sites[key][4]

            
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
    # right anchor
    # keep only - strand CBS
    # ========================================================

    for i in range(endR - startR):

        key = (chromR, startR + i)

        if key in ctcf_sites:

            CBSstart = ctcf_sites[key][0]
            CBSend = ctcf_sites[key][1]
            score = ctcf_sites[key][2]
            direction = ctcf_sites[key][4]

            # 原代码：右 anchor 排除 +
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

            nondup_affinity[flag][
                chromL, CBSstart, CBSend
            ] = 0


# ============================================================
# 5.  affinity_norm
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
# 6. output
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
