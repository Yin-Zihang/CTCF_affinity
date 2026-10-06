#!/usr/bin/env python3
# 13review_bootstrap_thresholds_plot_all.py
# Just make the plot based on the existing CSV, without resampling or refitting.
# Place this file in the CTCF_IMP directory and run it directly:
# python 13review_bootstrap_thresholds_plot_all.py

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit

def read_csv(filename):
    with open(filename, newline='') as f:
        return list(csv.DictReader(f))

def get_binning_points(site_file, original_bins, plot_dir, prefix):
    if not site_file.is_file():
        raise FileNotFoundError('Missing per-base file %s; can't rebuild scatter plot for 140/70 bins' % site_file)
    n_groups = len(original_bins)
    group_values = [[] for _ in range(n_groups)]
    with open(site_file, newline='') as f:
        for row in csv.DictReader(f):
            group_number = int(row['group_number'])
            if group_number < 1 or group_number > n_groups:
                raise ValueError('Invalid group number found in site_records：%s' % group_number)
            group_values[group_number - 1].append(float(row['occupancy_percent']))

    points_by_merge = {}
    all_points = []
    for merge in (1, 2, 4):
        if n_groups % merge:
            raise ValueError('The number of groups %s cannot evenly divide %s' % (n_groups, merge))
        points = []
        for start in range(0, n_groups, merge):
            values = []
            for group_index in range(start, start + merge):
                values.extend(group_values[group_index])
            if not values:
                raise ValueError('Merged group %s–%s has no matching sites' % (start + 1, start + merge))
            point = {'n_bins': n_groups // merge,
                     'merge_adjacent_groups': merge,
                     'group_first': start + 1,
                     'group_last': start + merge,
                     'rank_x': (2 * start + merge + 1) / 20,
                     'n_records': len(values),
                     'frequency_pct': float(np.median(values))}
            points.append(point)
            all_points.append(point)
        points_by_merge[merge] = points

    # Check whether the per-site files match the 280 input sets for this fit to avoid mixing up result directories.
    for source, recomputed in zip(original_bins, points_by_merge[1]):
        if not np.isclose(float(source['rank_x']), recomputed['rank_x'],
                          rtol=0, atol=1e-8):
            raise ValueError('The group numbers of site_records and baseline_bins don't match')
        if not np.isclose(float(source['frequency_pct']),
                          recomputed['frequency_pct'], rtol=0, atol=1e-7):
            raise ValueError('The occupancy frequency of site_records doesn’t match that of baseline_bins')

    output_file = plot_dir / ('%s_binning_plot_points.csv' % prefix)
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(all_points[0]))
        writer.writeheader()
        writer.writerows(all_points)
    print('NEW CSV:', output_file.resolve(), flush=True)
    return points_by_merge

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
    raise ValueError('Unknown model：' + model)

def save_figure(fig, plot_dir, name):
    png = plot_dir / (name + '.png')
    svg = plot_dir / (name + '.svg')
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(svg, bbox_inches='tight')
    print('NEW PNG:', png.resolve(), flush=True)
    print('NEW SVG:', svg.resolve(), flush=True)

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
    grid = np.linspace(max(0.01, x[0]), x[-1], 1200)
    # The new catalog clearly separates the old version that mistakenly used 280 scatter plots.
    plot_dir = result_dir / ('%s_plots_true_bins' % prefix)
    plot_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
    binning_points = get_binning_points(result_dir / ('%s_site_records.csv' % prefix),
                                        bins, plot_dir, prefix)

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
    save_figure(fig, plot_dir, '%s_original_fit_all' % prefix)
    plt.close(fig)

    # Figure 2: The fitting curve and the first, second, and third derivatives corresponding to the old R graph
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
    save_figure(fig, plot_dir, '%s_derivatives_all' % prefix)
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
    save_figure(fig, plot_dir, '%s_model_sensitivity_all' % prefix)
    plt.close(fig)

    # Figure 4: Using four panels from v2: red, blue, and green are drawn respectively; bottom right corner shows enlarged curve differences.
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    colors = {1: '#d62728', 2: '#1f77b4', 4: '#008837'}
    titles = {1: '280 bins (original)', 2: '140 bins (merge pairs)',
              4: '70 bins (merge four)'}
    grid_short = np.linspace(x[0], min(9, x[-1]), 1200)
    curves = {}
    fit_by_merge = {}
    for row in sensitivity:
        if row['model'] != 'original_logistic_logx' or row['success'] != '1':
            continue
        merge = int(row['merge_adjacent_groups'])
        if merge in (1, 2, 4):
            fit_by_merge[merge] = row
    for merge in (1, 2, 4):
        if merge not in fit_by_merge:
            raise ValueError('sensitivity.csv 缺少 %s 组合并的成功拟合结果' % merge)
        row = fit_by_merge[merge]
        param = [float(value) for value in row['params'].split(';')]
        curves[merge] = model_value(grid_short, *param, model=row['model'])
        ax = {1: axes[0, 0], 2: axes[0, 1], 4: axes[1, 0]}[merge]
        point_x = [point['rank_x'] for point in binning_points[merge]]
        point_y = [point['frequency_pct'] for point in binning_points[merge]]
        ax.scatter(point_x, point_y, color='#aaaaaa', s=10, alpha=0.6,
                   label='%s bin medians (n=%s)' % (280 // merge, len(point_x)), zorder=1)
        ax.plot(grid_short, curves[merge], color=colors[merge], lw=2.5,
                label=titles[merge], zorder=3)
        ax.axvline(float(row['high_rank']), color='#666666', ls=':', lw=1)
        ax.axvline(float(row['low_rank']), color='#666666', ls=':', lw=1)
        ax.text(0.97, 0.96, 'High affinity = %.4f\nLow affinity = %.4f' %
                (float(row['high_emsa']), float(row['low_emsa'])),
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                bbox={'facecolor': 'white', 'alpha': 0.85, 'edgecolor': 'none'})
        ax.set(xlim=(0, grid_short[-1]), ylim=(0, 105),
               xlabel='Rank coordinate (group / 10)',
               ylabel='Occupancy frequency (%)')
        ax.set_title(titles[merge], color=colors[merge])
        ax.legend(loc='center right', frameon=False, fontsize=8)

    ax = axes[1, 1]
    for merge in (1, 2, 4):
        ax.plot(grid_short, curves[merge] - curves[1], color=colors[merge],
                lw=2, label='%s bins' % (280 // merge))
    ax.axhline(0, color='#777777', lw=0.8)
    ax.set(xlim=(0, grid_short[-1]),
           xlabel='Rank coordinate (group / 10)',
           ylabel='Difference from 280-bin fit (percentage points)',
           title='Magnified differences between fitted curves')
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Binning sensitivity: transition region (rank 0–9)', fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, plot_dir, '%s_binning_sensitivity_separate_all' % prefix)
    plt.close(fig)

    # Figure 5: Full sorting range. Each binning has its own subplot, with a zoomed-in curve interpolation in the bottom right corner.
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    full_curves = {}
    for merge in (1, 2, 4):
        row = fit_by_merge[merge]
        param = [float(value) for value in row['params'].split(';')]
        full_curves[merge] = model_value(grid, *param, model=row['model'])
        ax = {1: axes[0, 0], 2: axes[0, 1], 4: axes[1, 0]}[merge]
        point_x = [point['rank_x'] for point in binning_points[merge]]
        point_y = [point['frequency_pct'] for point in binning_points[merge]]
        ax.scatter(point_x, point_y, color='#aaaaaa', s=10, alpha=0.6,
                   label='%s bin medians (n=%s)' % (280 // merge, len(point_x)), zorder=1)
        ax.plot(grid, full_curves[merge], color=colors[merge], lw=2.5,
                label=titles[merge], zorder=3)
        ax.axvline(float(row['high_rank']), color='#666666', ls=':', lw=1)
        ax.axvline(float(row['low_rank']), color='#666666', ls=':', lw=1)
        ax.text(0.97, 0.96, 'High affinity = %.4f\nLow affinity = %.4f' %
                (float(row['high_emsa']), float(row['low_emsa'])),
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                bbox={'facecolor': 'white', 'alpha': 0.85, 'edgecolor': 'none'})
        ax.set(xlim=(0, x[-1]), ylim=(0, 105),
               xlabel='Rank coordinate (group / 10)',
               ylabel='Occupancy frequency (%)')
        ax.set_title(titles[merge], color=colors[merge])
        ax.legend(loc='center right', frameon=False, fontsize=8)

    ax = axes[1, 1]
    for merge in (1, 2, 4):
        ax.plot(grid, full_curves[merge] - full_curves[1],
                color=colors[merge], lw=2, label='%s bins' % (280 // merge))
    ax.axhline(0, color='#777777', lw=0.8)
    ax.set(xlim=(0, x[-1]), xlabel='Rank coordinate (group / 10)',
           ylabel='Difference from 280-bin fit (percentage points)',
           title='Magnified differences between fitted curves')
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Binning sensitivity: full rank range (all 280 groups)', fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, plot_dir, '%s_binning_sensitivity_fullrange_all' % prefix)
    plt.close(fig)

    # Figure 6: Distribution of threshold affinity over 1000 runs and the 2.5%/97.5% percentiles
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
        save_figure(fig, plot_dir, '%s_bootstrap_distributions_all' % prefix)
        plt.close(fig)
    print('Plots written to', plot_dir.resolve())

# ============================================================
# Read the existing CSV, output 6 PNGs, 6 SVGs, and 1 scatter check CSV
# ============================================================
def main():
    parser = argparse.ArgumentParser(description='Plot from the saved fit and bootstrap CSV; don’t recalculate the thresholds')
    parser.add_argument('--work-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--outdir', type=Path, default=Path('13review_bootstrap_results_rerun'))
    parser.add_argument('--prefix', default='13review')
    parser.add_argument('--plot-existing', action='store_true',
                        help='Backward compatible: This script always only reads existing CSV files for plotting')
    args = parser.parse_args()

    work_dir = args.work_dir.resolve()
    result_dir = args.outdir if args.outdir.is_absolute() else work_dir / args.outdir
    if not result_dir.is_dir():
        raise FileNotFoundError('Can't find the results folder：%s' % result_dir)
    for name in ('baseline_bins', 'baseline_fit', 'sensitivity',
                 'bootstrap_replicates', 'bootstrap_ci', 'site_records'):
        filename = result_dir / ('%s_%s.csv' % (args.prefix, name))
        if not filename.is_file():
            raise FileNotFoundError('Missing the CSV needed for the output：%s' % filename)
    draw_plots(result_dir, args.prefix)


if __name__ == '__main__':
    main()
