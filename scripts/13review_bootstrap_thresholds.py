#!/usr/bin/env python3
# 13review_bootstrap_thresholds.py
# python 13review_bootstrap_thresholds.py --count-file 12_count_dup_cell_new.txt --n-cells 61 --bootstrap 1000
# python 13review_bootstrap_thresholds.py --remap-existing

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import expit, logit


# ============================================================
# 0. parameters & output
# ============================================================
parser = argparse.ArgumentParser()
parser.add_argument('--work-dir', type=Path, default=Path(__file__).resolve().parent)
parser.add_argument('--count-file', type=Path, default=Path('12_count_dup_cell_new.txt'))
parser.add_argument('--group-dir', type=Path, default=Path('11_split_to_270_groups'))
parser.add_argument('--score-file', type=Path, default=None)
parser.add_argument('--outdir', type=Path, default=Path('13review_bootstrap_results'))
parser.add_argument('--prefix', default='13review')
parser.add_argument('--n-cells', type=int, default=61)
parser.add_argument('--n-groups', type=int, default=280)
parser.add_argument('--bootstrap', type=int, default=1000)
parser.add_argument('--seed', type=int, default=20260927)
parser.add_argument('--rank-to-line-scale', type=int, default=10000)
parser.add_argument('--remap-existing', action='store_true')
parser.add_argument('--plot-existing', action='store_true',
                    help='Plot the graph just based on the generated CSV, without rereading the original data or refitting.')
parser.add_argument('--diagnose-counts-only', action='store_true')
parser.add_argument('--skip-site-table', action='store_true')
parser.add_argument('--expected-low', type=float, default=2.25)
parser.add_argument('--expected-high', type=float, default=4.07)
parser.add_argument('--reference-tolerance', type=float, default=0.15)
parser.add_argument('--strict-reference', action='store_true')

FIT_COLUMNS = [
    'model', 'n_bins', 'rss', 'r2', 'params',
    'left_rank', 'inflection_rank', 'right_rank',
    'left_emsa', 'inflection_emsa', 'right_emsa',
    'left_score_file_line', 'inflection_score_file_line', 'right_score_file_line',
    'left_score_file_group', 'inflection_score_file_group', 'right_score_file_group',
    'low_rank', 'high_rank', 'low_emsa', 'high_emsa',
    'low_score_file_line', 'high_score_file_line'
]


def save_csv(filename, columns, data):
    with open(filename, 'w', newline='') as fjob:
        writer = csv.DictWriter(fjob, fieldnames=columns)
        writer.writeheader()
        for each in data:
            writer.writerow(each)


def read_csv(filename):
    with open(filename, newline='') as f:
        return list(csv.DictReader(f))


# ============================================================
# 1. Sequence direction handling consistent with the original code
# ============================================================
def reverse(seq):
    new = ''
    l = len(seq)
    for i in range(l):
        dna = seq[l - 1 - i]
        if dna == 'A':
            new += 'T'
        elif dna == 'T':
            new += 'A'
        elif dna == 'C':
            new += 'G'
        elif dna == 'G':
            new += 'C'
    return new


# ============================================================
# 2. Calculate a point for each group: the x-axis is the group number divided by 10, and the y-axis is the median occupancy divided by the number of cells
# ============================================================
def make_group_table(group_count, group_emsa, n_cells, merge=1):
    if len(group_count) % merge != 0:
        raise ValueError('The number of groups can't be evenly divided by the number to merge')
    table = []
    for start in range(0, len(group_count), merge):
        all_count = []
        all_emsa = []
        for j in range(start, start + merge):
            all_count.extend(group_count[j])
            all_emsa.extend(group_emsa[j])
        first_group = start + 1
        last_group = start + merge
        median_count = float(np.median(all_count))
        table.append({
            'group_first': first_group,
            'group_last': last_group,
            'rank_x': (first_group + last_group) / 20,
            'n_records': len(all_count),
            'median_emsa': float(np.median(all_emsa)),  
            'median_count': median_count,
            'frequency_pct': median_count * 100 / n_cells
        })
    return table


# ============================================================
# 3. Fitting; when scal is negative, the curve decreases with the ranking
# ============================================================
def model_value(x, *param, model):
    if model == 'original_logistic_logx':
        Asym, xmid, scal = param
        return 3 + Asym * expit((np.log(x) - xmid) / scal)
    if model == 'free_floor_logistic_logx':
        floor, Asym, xmid, scal = param
        return floor + Asym * expit((np.log(x) - xmid) / scal)
    if model == 'original_floor_logistic_x':
        Asym, xmid, scal = param
        return 3 + Asym * expit((x - xmid) / scal)
    raise ValueError('unknown model：' + model)


def fit_model(table, model):
    x = np.array([row['rank_x'] for row in table], dtype=float)
    y = np.array([row['frequency_pct'] for row in table], dtype=float)
    if np.any(x <= 0) or np.any(~np.isfinite(y)):
        raise ValueError('The order must be greater than 0, and the occupancy frequency must be a finite value')

    Asym0 = max(1.0, float(max(y) - 3))
    half_index = int(np.argmin(np.abs(y - (3 + Asym0 / 2))))
    if model == 'original_floor_logistic_x':
        xmid0 = float(x[half_index])
    else:
        xmid0 = float(np.log(x[half_index]))

    best_param = None
    best_rss = float('inf')
    for sign in (-1, 1):
        if model == 'original_logistic_logx':
            start = [Asym0, xmid0, sign * 0.3]
            lower = [0, -20, -20] if sign < 0 else [0, -20, 0.01]
            upper = [200, 20, -0.01] if sign < 0 else [200, 20, 20]
        elif model == 'free_floor_logistic_logx':
            floor0 = max(0.0, min(99.0, float(min(y))))
            start = [floor0, max(1.0, float(max(y) - floor0)), xmid0, sign * 0.3]
            lower = [0, 0, -20, -20] if sign < 0 else [0, 0, -20, 0.01]
            upper = [100, 200, 20, -0.01] if sign < 0 else [100, 200, 20, 20]
        else:
            start = [Asym0, xmid0, sign * 1.0]
            lower = [0, -100, -100] if sign < 0 else [0, -100, 0.01]
            upper = [200, 100, -0.01] if sign < 0 else [200, 100, 100]
        try:
            fitted_param, _ = curve_fit(
                lambda xx, *pp: model_value(xx, *pp, model=model),
                x, y, p0=start, bounds=(lower, upper), maxfev=20000
            )
            fitted_y = model_value(x, *fitted_param, model=model)
            rss = float(np.sum((y - fitted_y) ** 2))
            if rss < best_rss:
                best_rss = rss
                best_param = fitted_param
        except (RuntimeError, ValueError, OverflowError):
            continue
    if best_param is None:
        raise RuntimeError('Both forward and backward fitting didn’t work：' + model)
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1 - best_rss / tss if tss > 0 else float('nan')
    return best_param, best_rss, r2


# ============================================================
# 4. Find the three positions on the horizontal axis: left f'''=0, middle f''=0, right f'''=0
#    The threshold only takes the two positions on the left and right; the middle inflection point is used as a reference output.
# ============================================================
def get_three_ranks(param, model):
    if model == 'free_floor_logistic_logx':
        xmid = param[2]
        scal = param[3]
    else:
        xmid = param[1]
        scal = param[2]

    if model == 'original_floor_logistic_x':
        s1 = (3 - math.sqrt(3)) / 6
        s3 = (3 + math.sqrt(3)) / 6
        left = xmid + scal * logit(s1)
        right = xmid + scal * logit(s3)
        return min(left, right), float(xmid), max(left, right)

    p = 1 / scal
    if abs(p) <= 1:
        raise ValueError('The second derivative of scal doesn't have any valid zeros')
    center = (p - 1) / (2 * p)
    delta = math.sqrt(3 * (p * p - 1)) / (6 * abs(p))
    s1 = center - delta
    s3 = center + delta
    if s1 <= 0 or s3 >= 1:
        raise ValueError('The zero of the third derivative is out of the valid range')

    root1 = math.exp(xmid + scal * logit(s1))
    root3 = math.exp(xmid + scal * logit(s3))
    middle = math.exp(xmid + scal * logit(center))
    return min(root1, root3), middle, max(root1, root3)


# ============================================================
# 5. Check affinity by the line number of the fully sorted file
#    Example: x=5.0150 -> Line number round(5.0150*10000)=50150 -> Column 5=2.2505
#    The first line of the file is the header; line numbers refer to physical lines, not the median of adjacent groups.
# ============================================================
def find_affinity(rank, rank_rows, scale):
    line_number = math.floor(rank * scale + 0.5)
    if line_number not in rank_rows:
        raise ValueError('Rank %.8f corresponds to line %s, but the affinity file doesn't have that line' %
                         (rank, line_number))
    group_id, affinity = rank_rows[line_number]
    expected_group = math.ceil(rank * 10)
    if abs(int(group_id) - expected_group) > 1:
        raise ValueError('The group number %s and order %.8f on line %s do not match, please check the order file' %
                         (line_number, group_id, rank))
    return affinity, line_number, group_id


def add_affinity_columns(row, rank_rows, scale):
    # The row already has three sort positions; write the corresponding file line numbers and affinity into the results as is
    left_rank = float(row['left_rank'])
    middle_rank = float(row['inflection_rank'])
    right_rank = float(row['right_rank'])
    left_affinity, left_line, left_group = find_affinity(left_rank, rank_rows, scale)
    middle_affinity, middle_line, middle_group = find_affinity(middle_rank, rank_rows, scale)
    right_affinity, right_line, right_group = find_affinity(right_rank, rank_rows, scale)

    row['left_emsa'] = left_affinity
    row['inflection_emsa'] = middle_affinity
    row['right_emsa'] = right_affinity
    row['left_score_file_line'] = left_line
    row['inflection_score_file_line'] = middle_line
    row['right_score_file_line'] = right_line
    row['left_score_file_group'] = left_group
    row['inflection_score_file_group'] = middle_group
    row['right_score_file_group'] = right_group

    if left_affinity <= right_affinity:
        row['low_rank'] = left_rank
        row['high_rank'] = right_rank
        row['low_emsa'] = left_affinity
        row['high_emsa'] = right_affinity
        row['low_score_file_line'] = left_line
        row['high_score_file_line'] = right_line
    else:
        row['low_rank'] = right_rank
        row['high_rank'] = left_rank
        row['low_emsa'] = right_affinity
        row['high_emsa'] = left_affinity
        row['low_score_file_line'] = right_line
        row['high_score_file_line'] = left_line
    return row


def fit_one_table(table, model, rank_rows, scale):
    param, rss, r2 = fit_model(table, model)
    left, middle, right = get_three_ranks(param, model)
    if left < table[0]['rank_x'] or right > table[-1]['rank_x']:
        raise ValueError('The threshold is beyond the range of the fitted data')
    row = {
        'model': model, 'n_bins': len(table),
        'rss': rss, 'r2': r2,
        'params': ';'.join('%.12g' % value for value in param),
        'left_rank': left, 'inflection_rank': middle, 'right_rank': right
    }
    return add_affinity_columns(row, rank_rows, scale)


# ============================================================
# 6. Calculate the 95% CI from 1000 results
# ============================================================
def write_ci(filename, baseline, bootstrap_rows, bootstrap_total):
    success = []
    for row in bootstrap_rows:
        if str(row['success']) == '1':
            success.append(row)
    if len(success) < max(20, bootstrap_total * 0.8):
        raise RuntimeError('Successfully fitted %s/%s, did not reach the preset number' %
                           (len(success), bootstrap_total))
    summary = []
    for name in ('low_rank', 'high_rank', 'low_emsa', 'high_emsa'):
        numbers = []
        for row in success:
            numbers.append(float(row[name]))
        ci_low, ci_high = np.percentile(numbers, [2.5, 97.5])
        summary.append({
            'parameter': name,
            'baseline': baseline[name],
            'bootstrap_median': float(np.median(numbers)),
            'ci_2_5': float(ci_low),
            'ci_97_5': float(ci_high),
            'n_success': len(success),
            'n_failed': bootstrap_total - len(success)
        })
    columns = ['parameter', 'baseline', 'bootstrap_median', 'ci_2_5',
               'ci_97_5', 'n_success', 'n_failed']
    save_csv(filename, columns, summary)


# ============================================================
# 7. Plotting: All curves use the fitting parameters already saved in the table
#    plot-existing read-only CSV; won’t resample or refit
# ============================================================
def save_figure(fig, plot_dir, name):
    fig.savefig(plot_dir / (name + '.png'), dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir / (name + '.svg'), bbox_inches='tight')

def read_plot_binning_points(site_file, original_bins):
    # Merge adjacent groups based on the occupancy frequency of each site; don't take the median of the old group and then median it again.
    n_groups = len(original_bins)
    values_by_group = [[] for _ in range(n_groups)]
    with open(site_file, newline='') as f:
        for row in csv.DictReader(f):
            group_number = int(row['group_number'])
            if group_number < 1 or group_number > n_groups:
                raise ValueError('Invalid group number found in site_records：%s' % group_number)
            values_by_group[group_number - 1].append(float(row['occupancy_percent']))
    points_by_merge = {}
    for merge in (1, 2, 4):
        if n_groups % merge:
            raise ValueError('The number of groups %s cannot evenly divide %s' % (n_groups, merge))
        points = []
        for start in range(0, n_groups, merge):
            values = []
            for i in range(start, start + merge):
                values.extend(values_by_group[i])
            if not values:
                raise ValueError('Merged group %s–%s has no matching sites' % (start + 1, start + merge))
            points.append(((2 * start + merge + 1) / 20,
                           float(np.median(values))))
        points_by_merge[merge] = points
    for row, point in zip(original_bins, points_by_merge[1]):
        if not np.isclose(float(row['rank_x']), point[0], rtol=0, atol=1e-8):
            raise ValueError('The order of site_records doesn't match baseline_bins')
        if not np.isclose(float(row['frequency_pct']), point[1], rtol=0, atol=1e-7):
            raise ValueError('The occupancy frequency of site_records doesn’t match that of baseline_bins')
    return points_by_merge


def draw_plots(result_dir, prefix, bin_dir=None):
    try:
        import matplotlib
        matplotlib.use('Agg') 
        import matplotlib.pyplot as plt
    except ImportError:
        print('WARNING: matplotlib isn't installed, the data table has been generated; once installed, you can use --plot-existing to plot',
              file=sys.stderr)
        return

    if bin_dir is None:
        bin_dir = result_dir
    bins = read_csv(bin_dir / ('%s_baseline_bins.csv' % prefix))
    baseline = read_csv(result_dir / ('%s_baseline_fit.csv' % prefix))[0]
    sensitivity = read_csv(result_dir / ('%s_sensitivity.csv' % prefix))
    replicate_file = result_dir / ('%s_bootstrap_replicates.csv' % prefix)
    ci_file = result_dir / ('%s_bootstrap_ci.csv' % prefix)

    x = np.array([float(row['rank_x']) for row in bins])
    y = np.array([float(row['frequency_pct']) for row in bins])
    site_file = bin_dir / ('%s_site_records.csv' % prefix)
    binning_points = None
    if site_file.is_file():
        binning_points = read_plot_binning_points(site_file, bins)
    else:
        print('WARNING:Missing %s, the binning chart will only show the fitted curve, not the erroneous scatter points' % site_file,
              file=sys.stderr)
    grid = np.linspace(max(0.01, x[0]), x[-1], 1200)
    plot_dir = result_dir / ('%s_plots' % prefix)
    plot_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})

    # Figure 1: Original 280 scatter points, original model curve, two affinity boundaries, and the middle turning point
    fig, ax = plt.subplots(figsize=(9, 5.5))
    param = [float(value) for value in baseline['params'].split(';')]
    fitted_y = model_value(grid, *param, model=baseline['model'])
    ax.scatter(x, y, s=12, color='#999999', alpha=0.7, label='Group median occupancy')
    ax.plot(grid, fitted_y, color='#d62728', lw=2.2, label='Original logistic fit')
    ax.axvline(float(baseline['high_rank']), color='#7b3294', ls='--', lw=1.6,
               label='High affinity = %.4f' % float(baseline['high_emsa']))
    ax.axvline(float(baseline['low_rank']), color='#008837', ls='--', lw=1.6,
               label='Low affinity = %.4f' % float(baseline['low_emsa']))
    ax.axvline(float(baseline['inflection_rank']), color='#777777', ls=':', lw=1.2,
               label='Middle inflection (rank %.3f)' % float(baseline['inflection_rank']))
    ax.set(xlim=(0, x[-1]), ylim=(-3, 106), xlabel='Rank coordinate (group / 10)',
           ylabel='Occupancy frequency (%)', title='Original CTCF occupancy fit')
    ax.legend(loc='upper right', frameon=False, fontsize=8)
    fig.tight_layout()
    save_figure(fig, plot_dir, '%s_original_fit' % prefix)
    plt.close(fig)

    # Figure 2: The fitting curve and the first, second, and third derivatives 
    Asym, xmid, scal = param
    p = 1 / scal
    s = expit((np.log(grid) - xmid) / scal)
    first_derivative = Asym * p * s * (1 - s) / grid
    second_derivative = first_derivative / grid * (p * (1 - 2 * s) - 1)
    third_derivative = first_derivative / (grid * grid) * (
        (p * (1 - 2 * s) - 2) * (p * (1 - 2 * s) - 1)
        - 2 * p * p * s * (1 - s))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.scatter(x, y, s=10, color='#aaaaaa', alpha=0.6, label='Group median occupancy')
    ax.plot(grid, fitted_y, color='#d62728', lw=2.0, label='Fitted curve')
    ax.plot(grid, first_derivative, color='#62c462', lw=1.3, label="First derivative")
    ax.plot(grid, second_derivative, color='#76bedb', lw=1.3, label="Second derivative")
    ax.plot(grid, third_derivative, color='#dbbf32', lw=1.3, label="Third derivative")
    ax.axvline(float(baseline['high_rank']), color='#7b3294', ls='--', lw=1.2)
    ax.axvline(float(baseline['inflection_rank']), color='#777777', ls=':', lw=1.2)
    ax.axvline(float(baseline['low_rank']), color='#008837', ls='--', lw=1.2)
    ax.axhline(0, color='#888888', lw=0.8)
    ax.set(xlim=(0, x[-1]), ylim=(-30, 106),
           xlabel='Rank coordinate (group / 10)', ylabel='Frequency / derivative',
           title='Original logistic fit and derivatives')
    ax.legend(loc='upper right', frameon=False, fontsize=8)
    fig.tight_layout()
    save_figure(fig, plot_dir, '%s_derivatives' % prefix)
    plt.close(fig)

    # Figure 3: Three models, with the full image on top of each column and a zoomed-in transformation area below; only 280 sets of data are used to compare the models.
    model_order = ['original_logistic_logx', 'free_floor_logistic_logx',
                   'original_floor_logistic_x']
    model_title = ['Fixed floor 3, logistic(log rank)',
                   'Free floor, logistic(log rank)',
                   'Fixed floor 3, logistic(rank)']
    model_color = ['#d62728', '#1f77b4', '#ff7f0e']
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharey=True)
    for k in range(3):
        result = None
        for row in sensitivity:
            if row['merge_adjacent_groups'] == '1' and row['model'] == model_order[k] and row['success'] == '1':
                result = row
                break
        if result is None:
            for j in range(2):
                axes[j, k].text(0.5, 0.5, 'Fit failed', ha='center', transform=axes[j, k].transAxes)
            continue
        param = [float(value) for value in result['params'].split(';')]
        fitted_y = model_value(grid, *param, model=result['model'])
        for j in range(2):
            ax = axes[j, k]
            ax.scatter(x, y, s=9, color='#999999', alpha=0.55)
            ax.plot(grid, fitted_y, color=model_color[k], lw=2)
            ax.axvline(float(result['high_rank']), color='#7b3294', ls='--', lw=1.2)
            ax.axvline(float(result['low_rank']), color='#008837', ls='--', lw=1.2)
            ax.axhline(100, color='#dddddd', lw=0.8)
            ax.set_xlim(0, x[-1] if j == 0 else min(9, x[-1]))
            ax.set_ylim(-3, 110)
            ax.set_xlabel('Rank coordinate (group / 10)')
            ax.text(0.97, 0.95, 'High: %.4f\nLow: %.4f\nR²: %.4f' %
                    (float(result['high_emsa']), float(result['low_emsa']), float(result['r2'])),
                    ha='right', va='top', transform=ax.transAxes, fontsize=8,
                    bbox={'facecolor': 'white', 'alpha': 0.8, 'edgecolor': 'none'})
        axes[0, k].set_title(model_title[k])
    axes[0, 0].set_ylabel('Occupancy frequency (%)')
    axes[1, 0].set_ylabel('Occupancy frequency (%)')
    fig.suptitle('Model sensitivity: full range (top), threshold area (bottom); purple = high, green = low')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, plot_dir, '%s_model_sensitivity' % prefix)
    plt.close(fig)

    # Figure 4: The original values are shown on the left/middle with different line styles and staggered markers.；
    #      The right picture shows the differences relative to the 280 curves, zooming in on small differences without changing the fitting results.
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    colors = {1: '#d62728', 2: '#1f77b4', 4: '#008837'}
    line_styles = {1: '-', 2: '--', 4: ':'}
    markers = {1: 'o', 2: 's', 4: '^'}
    marker_start = {1: 35, 2: 85, 4: 135}
    curves = {}
    labels = {}
    for row in sensitivity:
        if row['model'] != 'original_logistic_logx' or row['success'] != '1':
            continue
        merge = int(row['merge_adjacent_groups'])
        param = [float(value) for value in row['params'].split(';')]
        curves[merge] = model_value(grid, *param, model=row['model'])
        labels[merge] = '%s bins: low %.3f, high %.3f' % (
            row['n_bins'], float(row['low_emsa']), float(row['high_emsa']))

    for panel in (0, 1):
        ax = axes[panel]
        if binning_points is not None:
            for merge in (1, 2, 4):
                point_x = [point[0] for point in binning_points[merge]]
                point_y = [point[1] for point in binning_points[merge]]
                ax.scatter(point_x, point_y, s=9, color=colors[merge],
                           marker=markers[merge], alpha=0.30, zorder=1)
        
        for merge in (4, 2, 1):
            if merge not in curves:
                continue
            ax.plot(grid, curves[merge], color=colors[merge],
                    ls=line_styles[merge], lw=1.8, marker=markers[merge],
                    markevery=(marker_start[merge], 210), ms=4,
                    mfc='white', mec=colors[merge], zorder=5 - merge / 4,
                    label=labels[merge])
        ax.set(xlim=(0, x[-1] if panel == 0 else min(9, x[-1])),
               ylim=(-3, 106), xlabel='Rank coordinate (group / 10)')
    axes[0].set_ylabel('Occupancy frequency (%)')
    axes[0].set_title('Full rank range')
    axes[1].set_title('Threshold area (rank 0–9, enlarged x axis)')
    axes[1].legend(loc='upper right', frameon=False, fontsize=8)

    # The 280 curves are defined as 0; the blue/green lines represent the frequency difference, in percentage points.
    ax = axes[2]
    if 1 in curves:
        for merge in (4, 2, 1):
            if merge not in curves:
                continue
            difference = curves[merge] - curves[1]
            ax.plot(grid, difference, color=colors[merge],
                    ls=line_styles[merge], lw=1.8, marker=markers[merge],
                    markevery=(marker_start[merge], 210), ms=4,
                    mfc='white', mec=colors[merge], label='%s bins' % (280 // merge))
    ax.set(xlim=(0, min(9, x[-1])),
           xlabel='Rank coordinate (group / 10)',
           ylabel='Difference vs 280-bin curve (percentage points)',
           title='Small differences magnified')
    ax.axhline(0, color='#777777', lw=0.8)
    ax.legend(frameon=False, fontsize=8)
    if binning_points is None:
        fig.suptitle('Binning sensitivity: fitted curves only (site records unavailable)')
    else:
        fig.suptitle('Binning sensitivity: 280/140/70 medians and fitted curves')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(fig, plot_dir, '%s_binning_sensitivity' % prefix)
    plt.close(fig)

    # Figure 5: Distribution of threshold affinity over 1000 runs and the 2.5%/97.5% percentiles
    if replicate_file.is_file() and ci_file.is_file():
        replicates = read_csv(replicate_file)
        ci_rows = read_csv(ci_file)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for k, name in enumerate(('low_emsa', 'high_emsa')):
            values = []
            for row in replicates:
                if row['success'] == '1':
                    values.append(float(row[name]))
            summary = None
            for row in ci_rows:
                if row['parameter'] == name:
                    summary = row
                    break
            if not values or summary is None:
                continue
            ax = axes[k]
            ax.hist(values, bins=30, color=('#008837' if k == 0 else '#7b3294'),
                    alpha=0.8, edgecolor='white')
            ax.axvline(float(summary['baseline']), color='black', lw=1.5,
                       label='Original estimate: %.4f' % float(summary['baseline']))
            ax.axvline(float(summary['ci_2_5']), color='#d62728', ls='--', lw=1.3,
                       label='95%% CI: %.3f–%.3f' %
                       (float(summary['ci_2_5']), float(summary['ci_97_5'])))
            ax.axvline(float(summary['ci_97_5']), color='#d62728', ls='--', lw=1.3)
            ax.set(xlabel='MpEMSA affinity', ylabel='Bootstrap replicate count',
                   title=('Low affinity threshold' if k == 0 else 'High affinity threshold'))
            ax.legend(frameon=False, fontsize=8)
        fig.suptitle('Bootstrap threshold distributions (%s successful fits)' %
                     ci_rows[0]['n_success'])
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        save_figure(fig, plot_dir, '%s_bootstrap_distributions' % prefix)
        plt.close(fig)
    print('Plots written to', plot_dir.resolve())


# ============================================================
# 8. Once bootstrap is done, just recheck affinity, no need to fit again
# ============================================================
def remap_existing(args, rank_rows, score_file):
    corrected_dir = args.outdir / 'rank_lookup'
    corrected_dir.mkdir(parents=True, exist_ok=True)
    baseline = None
    replicates = None

    for name in ('baseline_fit', 'sensitivity', 'bootstrap_replicates'):
        old_file = args.outdir / ('%s_%s.csv' % (args.prefix, name))
        data = read_csv(old_file)
        if len(data) == 0:
            raise ValueError('The result file is empty：%s' % old_file)

        for row in data:
            if row.get('success', '1') == '1':
                add_affinity_columns(row, rank_rows, args.rank_to_line_scale)
        if name == 'baseline_fit':
            baseline = data[0]
            baseline['reference_match'] = (
                abs(baseline['low_emsa'] - args.expected_low) <= args.reference_tolerance
                and abs(baseline['high_emsa'] - args.expected_high) <= args.reference_tolerance
            )
        if name == 'bootstrap_replicates':
            replicates = data

        columns = []
        for row in data:
            for key in row:
                if key not in columns:
                    columns.append(key)
        save_csv(corrected_dir / old_file.name, columns, data)

    write_ci(corrected_dir / ('%s_bootstrap_ci.csv' % args.prefix),
             baseline, replicates, len(replicates))
    audit_rows = [
        {'item': 'affinity_mapping', 'value': 'nearest physical line in ranked score TSV'},
        {'item': 'rank_to_file_line_scale', 'value': args.rank_to_line_scale},
        {'item': 'score_file', 'value': str(score_file)},
        {'item': 'original_fit_csvs', 'value': str(args.outdir)},
        {'item': 'output_csvs', 'value': str(corrected_dir)}
    ]
    save_csv(corrected_dir / ('%s_mapping_audit.csv' % args.prefix),
             ['item', 'value'], audit_rows)
    print('Original fitted rank roots preserved:',
          baseline['left_rank'], baseline['inflection_rank'], baseline['right_rank'])
    print('Direct ranked-file affinity (low, high):',
          baseline['low_emsa'], baseline['high_emsa'])
    print('Source TSV physical lines (low, high):',
          baseline['low_score_file_line'], baseline['high_score_file_line'])
    print('Remapped CSVs written to', corrected_dir.resolve())
    draw_plots(corrected_dir, args.prefix, args.outdir)


def main():
    args = parser.parse_args()
    if args.n_cells < 1 or args.n_groups < 8 or args.bootstrap < 0 or args.rank_to_line_scale < 1:
        raise ValueError('n-cells > 0, n-groups >= 8, bootstrap >= 0, rank-to-line-scale > 0')

    args.work_dir = args.work_dir.resolve()
    if not args.work_dir.is_dir():
        raise FileNotFoundError(args.work_dir)
    if not args.count_file.is_absolute():
        args.count_file = args.work_dir / args.count_file
    if not args.group_dir.is_absolute():
        args.group_dir = args.work_dir / args.group_dir
    if not args.outdir.is_absolute():
        args.outdir = args.work_dir / args.outdir
    if args.outdir.resolve() == args.work_dir:
        raise ValueError('Please put the results in a subfolder under the original directory')
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.score_file is None:
        args.score_file = args.group_dir / 'final_split_to_270_groups.tsv'
    elif not args.score_file.is_absolute():
        args.score_file = args.work_dir / args.score_file

    if args.plot_existing:
        draw_plots(args.outdir, args.prefix)
        return

    # First, read the affinity file. Use each physical line number as the key; for repeated sequences, the affinity from the original code is used, with the latter overwriting the former.
    score_by_seq = {}
    rank_rows = {}
    score_conflicts = 0
    for line_number, eachLine in enumerate(open(args.score_file), 1):
        if line_number == 1:
            continue  
        if not eachLine.strip():
            continue
        each = eachLine.split()
        if len(each) < 5:
            raise ValueError('%s:%s affinity table has insufficient entries' % (args.score_file, line_number))
        seq = each[2]
        affinity = float(each[4])
        rank_rows[line_number] = (each[0], affinity)
        if seq in score_by_seq and score_by_seq[seq] != affinity:
            score_conflicts += 1
        score_by_seq[seq] = affinity

    if args.remap_existing:
        remap_existing(args, rank_rows, args.score_file)
        return

    # ========================================================
    # 9. Read the occupancy count, handle 41 bp, negative strand, and duplicate keys just like the original code
    # ========================================================
    count_by_seq = {}
    audit = {'count_input_rows': 0, 'count_duplicate_rows': 0,
             'count_conflicting_overwrites': 0, 'count_rows_with_non_acgt': 0,
             'count_non_acgt_bases_dropped_on_reverse': 0,
             'count_empty_sequence_rows': 0}
    non_acgt_examples = []
    diagnostic_file = args.outdir / ('%s_non_acgt_rows.tsv' % args.prefix)
    with open(diagnostic_file, 'w', newline='') as fjob:
        writer = csv.writer(fjob, delimiter='\t')
        writer.writerow(['line_number', 'raw_token', 'sequence_before_reverse',
                         'strand', 'invalid_position_base_codepoint', 'count'])
        for line_number, eachLine in enumerate(open(args.count_file), 1):
            if not eachLine.strip():
                continue
            each = eachLine.split()
            if len(each) < 6:
                raise ValueError('%s:%s Not enough input columns' % (args.count_file, line_number))
            token = each[3].split('_')
            if len(token) < 2:
                raise ValueError('%s:%s sequence field is missing an underscore' % (args.count_file, line_number))
            length = token[0]
            seq = token[1]  
            count = int(each[4])
            orient = each[5]
            if length == '41':
                seq = seq[1:]
            if seq == '':
                audit['count_empty_sequence_rows'] += 1
            invalid = []
            for pos, dna in enumerate(seq, 1):
                if dna not in 'ACGT':
                    invalid.append((pos, dna, ord(dna)))
            if invalid:
                audit['count_rows_with_non_acgt'] += 1
                writer.writerow([line_number, each[3], seq, orient, repr(invalid), count])
                if len(non_acgt_examples) < 5:
                    non_acgt_examples.append('line %s: %s' % (line_number, invalid))
                if args.diagnose_counts_only and audit['count_rows_with_non_acgt'] <= 10:
                    print('NON-ACGT %s:%s raw_token=%r sequence_before_reverse=%r strand=%r invalid_position_base_codepoint=%s' %
                          (args.count_file, line_number, each[3], seq, orient, invalid),
                          file=sys.stderr)
            if orient == '-':
                audit['count_non_acgt_bases_dropped_on_reverse'] += len(invalid)
                seq = reverse(seq)
            elif orient != '+':
                raise ValueError('%s:%s Unknown chain direction %s' % (args.count_file, line_number, orient))
            if count < 0 or count > args.n_cells:
                raise ValueError('%s:%s count=%s out of range' % (args.count_file, line_number, count))

            audit['count_input_rows'] += 1
            if seq in count_by_seq:
                audit['count_duplicate_rows'] += 1
                if count_by_seq[seq] != count:
                    audit['count_conflicting_overwrites'] += 1
            count_by_seq[seq] = count  

    audit['count_unique_sequences'] = len(count_by_seq)
    audit['count_explicit_zeros'] = sum(value == 0 for value in count_by_seq.values())
    audit['count_non_acgt_examples'] = '; '.join(non_acgt_examples)

    if args.diagnose_counts_only:
        rows = []
        for key, value in audit.items():
            rows.append({'item': key, 'value': value})
        save_csv(args.outdir / ('%s_count_diagnostic_summary.csv' % args.prefix),
                 ['item', 'value'], rows)
        print('Input diagnostics:', audit['count_rows_with_non_acgt'],
              'non-ACGT rows; full list in', diagnostic_file)
        return

    # ========================================================
    # 10. Read the 41/42 files group by group; skip positions that aren’t in the files.
    # ========================================================
    group_count = []
    group_emsa = []
    group_rows = 0
    missing_count = 0
    repeated_rows = 0
    seen = set()
    site_file = args.outdir / ('%s_site_records.csv' % args.prefix)
    site_handle = None
    if not args.skip_site_table:
        site_handle = open(site_file, 'w', newline='')
        site_writer = csv.writer(site_handle)
        site_writer.writerow(['sequence', 'group_number', 'rank_x', 'EMSA',
                              'occupancy_count', 'occupancy_percent', 'source_length'])
    try:
        for i in range(1, args.n_groups + 1):
            this_count = []
            this_emsa = []
            for length in (41, 42):
                group_file = args.group_dir / ('%s_group%s.txt' % (length, i))
                for eachLine in open(group_file):
                    if not eachLine.strip():
                        continue
                    seq = eachLine.split()[0]
                    group_rows += 1
                    if seq not in count_by_seq:
                        missing_count += 1
                        continue
                    if seq not in score_by_seq:
                        raise ValueError('Affinity not found for sequence %s of %s' % (group_file, seq))
                    if seq in seen:
                        repeated_rows += 1
                    seen.add(seq)
                    count = count_by_seq[seq]
                    affinity = score_by_seq[seq]
                    this_count.append(count)
                    this_emsa.append(affinity)
                    if site_handle is not None:
                        site_writer.writerow([seq, i, i / 10, affinity,
                                              count, count * 100 / args.n_cells, length])
            if len(this_count) == 0:
                raise ValueError('No matching occupancy records for group %s' % i)
            group_count.append(np.asarray(this_count, dtype=np.int16))
            group_emsa.append(np.asarray(this_emsa, dtype=float))
    finally:
        if site_handle is not None:
            site_handle.close()

    audit['score_conflicting_overwrites'] = score_conflicts
    audit['group_file_rows'] = group_rows
    audit['rows_missing_from_occupancy'] = missing_count
    audit['matched_group_rows'] = group_rows - missing_count
    audit['repeated_normalized_sequences_across_group_rows'] = repeated_rows
    audit['distinct_matched_sequences'] = len(seen)
    audit['count_file'] = str(args.count_file)
    audit['score_file'] = str(args.score_file)
    audit['group_dir'] = str(args.group_dir)
    audit['n_cells'] = args.n_cells
    audit['n_groups_fitted'] = args.n_groups
    audit['bootstrap_iterations'] = args.bootstrap
    audit['seed'] = args.seed
    audit['zero_values_for_absent_sequences'] = 'NO'
    audit['affinity_mapping'] = 'nearest physical line in ranked score TSV'
    audit['rank_to_file_line_scale'] = args.rank_to_line_scale
    audit_rows = []
    for key, value in audit.items():
        audit_rows.append({'item': key, 'value': value})
    save_csv(args.outdir / ('%s_audit.csv' % args.prefix), ['item', 'value'], audit_rows)

    if audit['count_conflicting_overwrites'] or score_conflicts:
        print('WARNING: conflicting duplicate keys were overwritten as in original code; inspect audit.',
              file=sys.stderr)
    if audit['count_rows_with_non_acgt']:
        print('WARNING: %s occupancy rows contain non-ACGT characters; see audit.' %
              audit['count_rows_with_non_acgt'], file=sys.stderr)
    if missing_count:
        print('WARNING: group rows absent from the count file remain excluded, as in original code.',
              file=sys.stderr)

    # ========================================================
    # 11. Fitted the original 280 sets; check if x=5.0150 corresponds to line 50150 in the file
    # ========================================================
    original_table = make_group_table(group_count, group_emsa, args.n_cells)
    save_csv(args.outdir / ('%s_baseline_bins.csv' % args.prefix),
             list(original_table[0]), original_table)
    baseline = fit_one_table(original_table, 'original_logistic_logx',
                             rank_rows, args.rank_to_line_scale)
    baseline['reference_match'] = (
        abs(baseline['low_emsa'] - args.expected_low) <= args.reference_tolerance
        and abs(baseline['high_emsa'] - args.expected_high) <= args.reference_tolerance
    )
    save_csv(args.outdir / ('%s_baseline_fit.csv' % args.prefix),
             list(baseline), [baseline])
    print("Baseline rank roots (left f''', f'', right f'''):",
          baseline['left_rank'], baseline['inflection_rank'], baseline['right_rank'])
    print('Original fit parameters (Asym; xmid; scal):',
          baseline['params'], 'R2:', baseline['r2'])
    print('First/last bins (occupancy %; median EMSA):',
          (original_table[0]['frequency_pct'], original_table[0]['median_emsa']),
          (original_table[-1]['frequency_pct'], original_table[-1]['median_emsa']))
    print('Ranked-file EMSA thresholds (low affinity, high affinity):',
          baseline['low_emsa'], baseline['high_emsa'])
    print('Source TSV physical lines (low, high):',
          baseline['low_score_file_line'], baseline['high_score_file_line'])
    if not baseline['reference_match']:
        warning = 'Baseline differs from reference (%.3f, %.3f); check data version and mapping' % (
            args.expected_low, args.expected_high)
        if args.strict_reference:
            raise RuntimeError(warning)
        print('WARNING:', warning, file=sys.stderr)

    # ========================================================
    # 12. Sensitivity: Combine adjacent groups 1/2/4; fit three curves for each grouping
    # ========================================================
    sensitivity = []
    for merge in (1, 2, 4):
        if args.n_groups % merge != 0:
            continue
        table = make_group_table(group_count, group_emsa, args.n_cells, merge)
        for model in ('original_logistic_logx',
                      'free_floor_logistic_logx',
                      'original_floor_logistic_x'):
            try:
                result = fit_one_table(table, model, rank_rows, args.rank_to_line_scale)
                result['merge_adjacent_groups'] = merge
                result['success'] = 1
                result['error'] = ''
            except (RuntimeError, ValueError, OverflowError) as exc:
                result = {'merge_adjacent_groups': merge, 'success': 0, 'error': str(exc)}
            sensitivity.append(result)
    save_csv(args.outdir / ('%s_sensitivity.csv' % args.prefix),
             ['merge_adjacent_groups', 'success', 'error'] + FIT_COLUMNS, sensitivity)

    # ========================================================
    # 13. Bootstrap: For each group, draw the same number of CBS with replacement, calculate the median 1000 times, and fit each time.
    #     Every time you get the sorting boundary, check the affinity from the original sorted file
    # ========================================================
    if args.bootstrap > 0:
        rng = np.random.default_rng(args.seed)
        bootstrap_rows = []
        for iteration in range(1, args.bootstrap + 1):
            sample_count = []
            sample_emsa = []
            for i in range(args.n_groups):
                n = len(group_count[i])
                index = rng.integers(0, n, size=n)
                sample_count.append(group_count[i][index])
                sample_emsa.append(group_emsa[i][index])
            try:
                sample_table = make_group_table(sample_count, sample_emsa, args.n_cells)
                result = fit_one_table(sample_table, 'original_logistic_logx',
                                       rank_rows, args.rank_to_line_scale)
                result['iteration'] = iteration
                result['success'] = 1
                result['error'] = ''
            except (RuntimeError, ValueError, OverflowError) as exc:
                result = {'iteration': iteration, 'success': 0, 'error': str(exc)}
            bootstrap_rows.append(result)
            if iteration % 100 == 0:
                print('Bootstrap %s/%s' % (iteration, args.bootstrap), file=sys.stderr)
        save_csv(args.outdir / ('%s_bootstrap_replicates.csv' % args.prefix),
                 ['iteration', 'success', 'error'] + FIT_COLUMNS, bootstrap_rows)
        write_ci(args.outdir / ('%s_bootstrap_ci.csv' % args.prefix),
                 baseline, bootstrap_rows, args.bootstrap)
    draw_plots(args.outdir, args.prefix)
    print('Results written to', args.outdir.resolve())


if __name__ == '__main__':
    main()
