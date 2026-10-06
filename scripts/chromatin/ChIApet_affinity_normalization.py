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
# 1. Read all CTCF sites
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

        
        # 41-bp sequence:  start + 1
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
# 2. Initialize
# ============================================================

# 4 loop-strength groups:
# 200, 300, 400, >=500
n_affinity = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]

# Record the unique CBS that participated in loops in each affinity class
a1 = {}
a2 = {}
a3 = {}

nondup_affinity = [a1, a2, a3]

# loopS[0] = strength 200
# loopS[1] = strength 300
# loopS[2] = strength 400
# loopS[3] = strength >=500
loopS = [0, 0, 0, 0, 0]


# ============================================================
# 3. Scan ChIA-PET loops
# ============================================================

for eachLine in open('wgEncodeGisChiaPetK562CtcfInteractionsRep1_filter.bed'):

    each = eachLine.split()

    chromL = each[3].split('-')[0].split(':')[0]
    startL = int(each[3].split('-')[0].split(':')[1].split('..')[0])
    endL = int(each[3].split('-')[0].split(':')[1].split('..')[1])

    chromR = each[3].split('-')[1].split(':')[0]
    startR = int(each[3].split('-')[1].split(':')[1].split('..')[0])
    endR = int(
        each[3].split('-')[1].split(':')[1]
        .split('..')[1]
        .split(',')[0]
    )

    loopstrength = int(each[4])

    # All loops >=500 are grouped into the 500 category
    if loopstrength >= 500:
        loopstrength = 500

    loop_group = loopstrength // 100 - 2
    loopS[loop_group] += 1

    if chromR != chromL:
        continue


    # --------------------------------------------------------
    # Left anchor: keep only strand CBS
    # --------------------------------------------------------

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

            n_affinity[loop_group] = aftCounter(
                score,
                n_affinity[loop_group]
            )

            # unique CBS
            nondup_affinity[flag][
                chromL, CBSstart, CBSend
            ] = 0


    # --------------------------------------------------------
    # Right anchor: keep only - strand CBS
    # --------------------------------------------------------

    for i in range(endR - startR):

        key = (chromR, startR + i)

        if key in ctcf_sites:

            CBSstart = ctcf_sites[key][0]
            CBSend = ctcf_sites[key][1]
            score = ctcf_sites[key][2]
            direction = ctcf_sites[key][4]

            if direction == '+':
                continue

            if score > 4.07:
                flag = 0
            elif score > 2.25:
                flag = 1
            else:
                flag = 2

            n_affinity[loop_group] = aftCounter(
                score,
                n_affinity[loop_group]
            )

            nondup_affinity[flag][
                chromL, CBSstart, CBSend
            ] = 0


# ============================================================
# 4. Calculate affinity_norm
# ============================================================

affinity_norm = []

for i in range(4):

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
# 5. output
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
print(loopS[:4])

print('affinity_norm:')
print(affinity_norm)

print(
    'g200=c(%s, %s, %s),' %
    tuple(affinity_norm[0])
)

print(
    'g300=c(%s, %s, %s),' %
    tuple(affinity_norm[1])
)

print(
    'g400=c(%s, %s, %s),' %
    tuple(affinity_norm[2])
)

print(
    'g500=c(%s, %s, %s),' %
    tuple(affinity_norm[3])
)
