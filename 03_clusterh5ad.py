# 03_clusterh5ad.py

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import json
from pathlib import Path

# --- Setup ---
sc.settings.seed = 0
np.random.seed(0)

data_dir = Path('data')
figures_dir = Path('figures/UMAP')
hvg_dir = Path('hvg')

for directory in [figures_dir, hvg_dir]:
    directory.mkdir(parents=True, exist_ok=True)

# --- Datasets and Methods ---
methods = ['standard', 'kneed', 'inflection']

# --- Processing function ---
def process_dataset(method_name):
    print(f"\nProcessing {method_name} dataset...")

    input_path = data_dir / f'pbmc3k_{method_name}.h5ad'
    output_fig_dir = figures_dir / method_name
    output_hvg_dir = hvg_dir / method_name

    output_fig_dir.mkdir(parents=True, exist_ok=True)
    output_hvg_dir.mkdir(parents=True, exist_ok=True)

    # Load dataset
    adata = sc.read_h5ad(input_path)

    # Preprocessing steps
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.pct_counts_mt < 5, :]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # Plot highly variable genes
    plt.figure(figsize=(10, 5))
    sc.pl.highly_variable_genes(adata, save=None)
    plt.savefig(output_fig_dir / 'highly_variable_genes.png', dpi=300)
    plt.close()

    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)

    # Determine number of PCs dynamically
    n_cells, n_genes = adata.shape
    n_pcs = min(40, n_cells - 1, n_genes - 1)

    # Preprocessing: PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, save=None)
    plt.savefig(output_fig_dir / 'pca_variance_ratio.png', dpi=300)
    plt.close()

    # Compute neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata, resolution=0.8, flavor='igraph', directed=False)

    # Save UMAP plot
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', frameon=False, save=None)
    plt.savefig(output_fig_dir / 'umap_clusters.png', dpi=300)
    plt.close()

    # Save HVGs per cluster
    def get_hvg_per_cluster(adata, cluster_key='leiden', n_top_genes=50):
        hvg_per_cluster = {}
        for cluster in adata.obs[cluster_key].unique():
            cluster_cells = adata[adata.obs[cluster_key] == cluster]
            mean_expr = np.asarray(cluster_cells.X.mean(axis=0)).flatten()
            gene_names = adata.var_names.tolist()
            top_genes_idx = np.argsort(mean_expr)[-n_top_genes:]
            top_genes = [gene_names[i] for i in top_genes_idx]
            hvg_per_cluster[f"cluster_{cluster}"] = top_genes
        return hvg_per_cluster

    hvg_dict = get_hvg_per_cluster(adata)
    with open(output_hvg_dir / 'hvg_per_cluster.json', 'w') as f:
        json.dump(hvg_dict, f, indent=2)

    # Save processed AnnData
    adata.write(data_dir / f'pbmc3k_{method_name}_clustered.h5ad')

    print(f"Completed {method_name}! Figures and HVG saved.")

# --- Main processing loop ---
for method in methods:
    process_dataset(method)

print("\nAll datasets processed and saved!")
