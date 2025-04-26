# 02_filteredh5ad.py

import os
import numpy as np
import scanpy as sc
import scipy
import matplotlib.pyplot as plt
from kneed import KneeLocator
from scipy.signal import find_peaks, savgol_filter
import seaborn as sns
import pandas as pd
from pathlib import Path

# --- Setup directories ---
data_dir = Path('data')
figures_dir = Path('figures/kneeplot')
figures_dir.mkdir(parents=True, exist_ok=True)

# --- Load raw PBMC3k dataset ---
raw_file = data_dir / 'pbmc3k_raw.h5ad'
adata = sc.read_h5ad(raw_file)
print(f"Loaded {raw_file} with shape {adata.shape}")

# --- Calculate basic metrics if missing ---
if 'total_counts' not in adata.obs.columns:
    adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten() if scipy.sparse.issparse(adata.X) else adata.X.sum(axis=1)
if 'n_genes_by_counts' not in adata.obs.columns:
    adata.obs['n_genes_by_counts'] = np.array((adata.X > 0).sum(axis=1)).flatten() if scipy.sparse.issparse(adata.X) else (adata.X > 0).sum(axis=1)

# --- Utility functions ---
def downsample(values, max_points=10000):
    n = len(values)
    step = max(1, n // max_points)
    return values[::step], np.arange(1, n + 1)[::step]

def find_inflection_points(x, y, smoothing_window=101, poly_order=3, prominence=0.05):
    log_x = np.log10(x)
    log_y = np.log10(y + 1)
    y_smooth = savgol_filter(log_y, smoothing_window, poly_order)
    dy = np.gradient(y_smooth, log_x)
    dy_smooth = savgol_filter(dy, smoothing_window, poly_order)
    peaks, _ = find_peaks(-dy_smooth, prominence=prominence)
    valleys, _ = find_peaks(dy_smooth, prominence=prominence)
    all_points = np.sort(np.concatenate([peaks, valleys]))
    return x[all_points], y[all_points]

def find_knee(values):
    sorted_vals = np.sort(values)[::-1]
    down_vals, x = downsample(sorted_vals)
    kl = KneeLocator(x, down_vals, curve='convex', direction='decreasing', online=True, interp_method='interp1d')
    return kl.knee, down_vals, x

def summarize_stats(adata, label):
    return {
        'Label': label,
        'Cells': adata.n_obs,
        'Mean Genes': np.mean(adata.obs['n_genes_by_counts']),
        'Median Genes': np.median(adata.obs['n_genes_by_counts']),
        'Mean UMIs': np.mean(adata.obs['total_counts']),
        'Median UMIs': np.median(adata.obs['total_counts'])
    }

def save_boxplots(adata, label_prefix, plot_dir):
    plot_dir.mkdir(parents=True, exist_ok=True)
    for metric in ['total_counts', 'n_genes_by_counts']:
        plt.figure(figsize=(6, 6))
        sns.boxplot(y=adata.obs[metric])
        plt.yscale('log')
        plt.title(f'{metric} ({label_prefix})')
        plt.tight_layout()
        plt.savefig(plot_dir / f'boxplot_{metric}.png')
        plt.close()

def filter_by_threshold(adata, metric, threshold_value):
    mask = adata.obs[metric] >= threshold_value
    return adata[mask].copy()

# --- Find knees and inflections ---
knee_results = {}
manual_results = {}

for metric in ['total_counts', 'n_genes_by_counts']:
    sorted_vals = np.sort(adata.obs[metric])[::-1]
    down_vals, x = downsample(sorted_vals)
    infl_x, infl_y = find_inflection_points(x, down_vals)
    kl = KneeLocator(x, down_vals, curve='convex', direction='decreasing', online=True, interp_method='interp1d')
    knee_x = kl.knee

    knee_results[metric] = (knee_x, down_vals, x)
    manual_results[metric] = (infl_x, infl_y)

    # Save knee plots
    for method, points, suffix in [('inflection', (infl_x, infl_y), 'manual'), ('kneed', (np.array([knee_x]), None), 'automatic')]:
        method_dir = figures_dir / method
        method_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.loglog(x, down_vals, 'b-')
        if points[0] is not None:
            ax.loglog(points[0], down_vals[points[0].astype(int)-1], 'ro')
        ax.set_title(f'{metric.capitalize()} ({suffix} detection)')
        ax.set_xlabel('Ranked Cells')
        ax.set_ylabel(metric)
        plt.tight_layout()
        plt.savefig(method_dir / f'kneeplot_{metric}.png')
        plt.close()

# --- Filtering ---
threshold_kneed = adata.obs['total_counts'].sort_values(ascending=False).iloc[int(knee_results['total_counts'][0]) - 1]
threshold_inflection = manual_results['total_counts'][1][0] if len(manual_results['total_counts'][1]) > 0 else None

adata_kneed = filter_by_threshold(adata, 'total_counts', threshold_kneed)
adata_inflection = filter_by_threshold(adata, 'total_counts', threshold_inflection)

adata_standard = adata.copy()
sc.pp.filter_cells(adata_standard, min_genes=200)
sc.pp.filter_genes(adata_standard, min_cells=3)
sc.pp.normalize_total(adata_standard, target_sum=1e4)
sc.pp.log1p(adata_standard)

# --- Save datasets ---
adata_standard.write(data_dir / 'pbmc3k_standard.h5ad')
adata_kneed.write(data_dir / 'pbmc3k_kneed.h5ad')
adata_inflection.write(data_dir / 'pbmc3k_inflection.h5ad')

# --- Save boxplots ---
save_boxplots(adata_standard, 'standard', figures_dir / 'standard')
save_boxplots(adata_kneed, 'kneed', figures_dir / 'kneed')
save_boxplots(adata_inflection, 'inflection', figures_dir / 'inflection')

# --- Summarize ---
summary = []
summary.append(summarize_stats(adata_standard, 'Standard'))
summary.append(summarize_stats(adata_kneed, 'Kneed'))
summary.append(summarize_stats(adata_inflection, 'Inflection'))
summary_df = pd.DataFrame(summary)
print("\nSummary Table:")
print(summary_df)

summary_df.to_csv(figures_dir / 'summary_metrics.csv', index=False)

print("\nAll steps completed successfully!")