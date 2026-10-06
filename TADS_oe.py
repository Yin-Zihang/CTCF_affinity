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
# 1. 读取全部 CBS，并统计全体 strong / middle / weak 数量
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

        # 保持原代码：
        # 所有 CBS 都计入全体 affinity background
        all_affinity = aftCounter(
            score,
            all_affinity
        )

        # 保持原来的有效序列筛选
        if re.search('None', each[5]) is None:

            ctcf_sites[(chrom, center)] = score

        else:
            print(eachLine)


# ============================================================
# 2. 搜索 TAD boundary 附近的 CBS
# ============================================================

# 原代码使用 50 kb
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
    # 左侧 TAD boundary
    # 从 start 向 TAD 内部扫描 50 kb
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
    # 右侧 TAD boundary
    # 从 end 向 TAD 内部扫描 50 kb
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
# 3. 计算 observed proportion
# ============================================================

observed = [
    n_affinity[i] / all_affinity[i]
    for i in range(3)
]


# ============================================================
# 4. 计算总体 expected proportion
# ============================================================

used_total = sum(n_affinity)
all_total = sum(all_affinity)

expected = used_total / all_total


# ============================================================
# 5. 计算 O/E
# ============================================================

oe = [
    observed[i] / expected
    for i in range(3)
]


# ============================================================
# 6. 输出
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
