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
# 1. 读取 CTCF sites
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
# 2. Hi-C loop strength 分组阈值
# ============================================================

tloopStrength = [
    6.741221875,
    9.94278506,
    13.99217675,
    20.039869500000002
]


# ============================================================
# 3. 初始化
# ============================================================

# 5 个 loop strength groups
# 每组依次统计 strong / middle / weak
n_affinity = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
]

# 记录每个 affinity group 中参与过 loop 的 unique CBS
a1 = {}
a2 = {}
a3 = {}

nondup_affinity = [a1, a2, a3]

# 每个 loop strength group 的 loop 数
loopS = [0, 0, 0, 0, 0]


# ============================================================
# 4. 扫描 Hi-C loops
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
    # 根据原代码划分 5 个 strength groups
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

    # 保持原代码逻辑：
    # 先统计 loopS，再过滤跨染色体 loop
    loopS[strength] += 1

    if chromR != chromL:
        continue


    # ========================================================
    # 左 anchor
    # 只保留 + strand CBS
    # ========================================================

    for i in range(endL - startL):

        key = (chromL, startL + i)

        if key in ctcf_sites:

            CBSstart = ctcf_sites[key][0]
            CBSend = ctcf_sites[key][1]
            score = ctcf_sites[key][2]
            direction = ctcf_sites[key][4]

            # 原代码：左 anchor 排除 -
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
    # 右 anchor
    # 只保留 - strand CBS
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

            # 保持原代码写法
            # 因为已经过滤 chromR != chromL，
            # 所以这里 chromL == chromR
            nondup_affinity[flag][
                chromL, CBSstart, CBSend
            ] = 0


# ============================================================
# 5. 计算 affinity_norm
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
# 6. 输出
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
