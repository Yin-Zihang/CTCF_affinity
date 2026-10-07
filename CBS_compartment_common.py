#!/usr/bin/env python3
"""Shared input parsing for CTCF compartment analyses (Python 3, standard library)."""

import bisect
import csv
import math
from collections import defaultdict


def read_intervals(path, label_column=3, one_based=False, labels=None):
    """Read sorted, non-overlapping intervals into a chromosome index."""
    intervals = defaultdict(list)
    with open(path) as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.strip().split()
            if len(fields) <= label_column:
                raise ValueError('%s:%d: expected at least %d columns' %
                                 (path, line_number, label_column + 1))
            try:
                start = int(fields[1]) - int(one_based)
                end = int(fields[2])
            except ValueError:
                if line_number == 1:
                    continue
                raise ValueError('%s:%d: invalid interval' % (path, line_number))
            if start < 0 or end <= start:
                raise ValueError('%s:%d: invalid coordinates' % (path, line_number))
            label = fields[label_column]
            if labels is not None:
                label = labels.get(label)
            elif label.upper() not in ('NA', 'NAN', 'NONE', '.'):
                try:
                    label = float(label)
                except ValueError:
                    label = None
                if label is not None and not math.isfinite(label):
                    label = None
            else:
                label = None
            intervals[fields[0]].append((start, end, label))
    index = {}
    for chrom, entries in intervals.items():
        entries.sort()
        for previous, current in zip(entries, entries[1:]):
            if current[0] < previous[1]:
                raise ValueError('%s: overlapping intervals on %s' % (path, chrom))
        index[chrom] = ([entry[0] for entry in entries], entries)
    if not index:
        raise ValueError('%s: no intervals found' % path)
    return index


def lookup_interval(index, chrom, position):
    item = index.get(chrom)
    if item is None:
        return None
    starts, entries = item
    location = bisect.bisect_right(starts, position) - 1
    if location < 0 or position >= entries[location][1]:
        return None
    return entries[location][2]


def affinity_from_pair(value, path, line_number):
    scores = value.strip().split('_')
    if len(scores) != 2:
        raise ValueError('%s:%d: affinity column needs two values' %
                         (path, line_number))
    valid = [float(score) for score in scores if score != 'None']
    if not valid:
        return None
    if not all(math.isfinite(score) for score in valid):
        raise ValueError('%s:%d: non-finite affinity' % (path, line_number))
    return sum(valid) / len(valid)


def read_sites(path):
    """The original site BED uses group, chrom, start, end, token, score, strand."""
    sites = defaultdict(dict)
    counts = defaultdict(int)
    with open(path) as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.strip().split()
            if len(fields) < 6:
                raise ValueError('%s:%d: site BED needs six columns' %
                                 (path, line_number))
            try:
                start, end = int(fields[2]), int(fields[3])
            except ValueError:
                if line_number == 1:
                    continue
                raise ValueError('%s:%d: invalid site coordinates' %
                                 (path, line_number))
            if end <= start:
                raise ValueError('%s:%d: invalid site length' % (path, line_number))
            affinity = affinity_from_pair(fields[5], path, line_number)
            if affinity is None:
                counts['missing_affinity'] += 1
                continue
            center = (start + end) // 2
            if center in sites[fields[1]]:
                counts['duplicate_site_centers_overwritten'] += 1
            sites[fields[1]][center] = affinity
            counts['site_rows_loaded'] += 1
    ordered = {}
    for chrom, values in sites.items():
        centers = sorted(values)
        ordered[chrom] = (centers, values)
    return ordered, counts


def sites_in_peak(sites, chrom, start, end):
    item = sites.get(chrom)
    if item is None:
        return []
    centers, values = item
    left = bisect.bisect_left(centers, start)
    right = bisect.bisect_left(centers, end)
    return [(center, values[center]) for center in centers[left:right]]


def read_peaks(path):
    with open(path, newline='') as source:
        table = csv.DictReader(source, delimiter='\t')
        required = ('ID', 'chrom', 'start', 'end',
                    'Raw_Counts_CTCF_rep1', 'Raw_Counts_CTCF_rep2')
        if not table.fieldnames or not all(name in table.fieldnames for name in required):
            raise ValueError('%s: expected columns: %s' % (path, ', '.join(required)))
        for row in table:
            start, end = int(row['start']), int(row['end'])
            if start < 0 or end <= start:
                raise ValueError('%s: invalid peak %s' % (path, row['ID']))
            yield (row['ID'], row['chrom'], start, end,
                   float(row['Raw_Counts_CTCF_rep1']),
                   float(row['Raw_Counts_CTCF_rep2']))


def chip_rpkm(count1, count2, depth1, depth2, length):
    return ((count1 / depth1 + count2 / depth2) / 2) * 1000000000 / length
