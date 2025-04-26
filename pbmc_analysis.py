# PBMC Data Analysis with Scanpy
# This script processes PBMC data and generates UMAP visualizations with cluster analysis.

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import json

sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(10, 10))
sc.settings.verbosity = 3

# Load the PBMC dataset
adata = sc.datasets.pbmc3k()
print(f"Dataset shape: {adata.shape}")

# Preprocess the data
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 5, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata, save='_hvg.pdf')
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# Run PCA and elbow plot
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, save='_elbow.pdf')

# Neighbors, UMAP, Clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)
sc.pl.umap(adata, color=['leiden'], legend_loc='on data', frameon=False, save='_clusters.pdf')

# Highly Variable Genes Per Cluster
def get_hvg_per_cluster(adata, cluster_key='leiden', n_top_genes=50):
    hvg_per_cluster = {}
    for cluster in adata.obs[cluster_key].unique():
        cluster_cells = adata[adata.obs[cluster_key] == cluster]
        mean_expr = cluster_cells.X.mean(axis=0)
        gene_names = adata.var_names.tolist()
        top_genes_idx = np.argsort(mean_expr)[-n_top_genes:]
        top_genes = [gene_names[i] for i in top_genes_idx]
        hvg_per_cluster[f"cluster_{cluster}"] = top_genes
    return hvg_per_cluster

hvg_dict = get_hvg_per_cluster(adata)
with open('hvg_per_cluster.json', 'w') as f:
    json.dump(hvg_dict, f, indent=2)

for cluster, genes in hvg_dict.items():
    print(f"\n{cluster} top 10 genes:")
    print(", ".join(genes[:10]))

# Save processed dataset
adata.write('pbmc_processed.h5ad')

print("\nAnalysis complete! Files saved:")
print("1. figures/umap_clusters.pdf - UMAP visualization with clusters")
print("2. figures/pca_elbow.pdf - PCA elbow plot")
print("3. figures/highly_variable_genes_hvg.pdf - Highly variable genes plot")
print("4. hvg_per_cluster.json - Highly variable genes per cluster")
print("5. pbmc_processed.h5ad - Processed data in h5ad format")