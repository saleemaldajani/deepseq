# 05_labelumap.py 

import json
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import seaborn as sns
from pathlib import Path

# --- Setup ---
data_dir = Path('data')
figures_dir = Path('figures/UMAP')
output_dir = figures_dir / 'labeled'
output_dir.mkdir(parents=True, exist_ok=True)

# --- Load predictions ---
with open('processed_clusters.json', 'r') as f:
    predictions = json.load(f)

# Organize predictions by method and cluster
pred_dict = {}
for item in predictions:
    method = item['method']
    cluster = item['cluster']
    label = item['cell_type_prediction']
    if method not in pred_dict:
        pred_dict[method] = {}
    pred_dict[method][cluster] = label

# --- Methods ---
methods = ['standard', 'kneed', 'inflection']

# --- Plot function ---
def plot_labeled_umap(adata, labels, method_name):
    umap = adata.obsm['X_umap']
    cluster_labels = adata.obs['leiden']
    unique_clusters = sorted(cluster_labels.unique(), key=int)

    # Set color palette
    palette = sns.color_palette("Set2", len(unique_clusters))
    cluster_to_color = {str(cluster): palette[i] for i, cluster in enumerate(unique_clusters)}

    plt.figure(figsize=(10, 8))

    for cluster in unique_clusters:
        idx = cluster_labels == cluster
        cluster_key = f"cluster_{cluster}"
        label = labels.get(cluster_key, cluster_key)
        pretty_label = f"{label} ({cluster})"
        plt.scatter(
            umap[idx, 0], umap[idx, 1],
            c=[cluster_to_color[str(cluster)]],
            label=pretty_label,
            s=10, edgecolor='none', alpha=0.9
        )

    plt.title(f"UMAP with Cell Type Labels ({method_name})", fontsize=14)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.gca().set_facecolor("#f9f9f9")
    plt.legend(loc='best', fontsize=8, frameon=True, framealpha=0.8, borderpad=0.5)
    plt.tight_layout()
    plt.savefig(output_dir / f'umap_labeled_{method_name}.png', dpi=300)
    plt.close()

# --- Process each method ---
for method in methods:
    input_path = data_dir / f'pbmc3k_{method}_clustered.h5ad'
    if not input_path.exists():
        print(f"Missing clustered file for {method}, skipping...")
        continue

    adata = sc.read_h5ad(input_path)

    if method not in pred_dict:
        print(f"Missing predictions for {method}, skipping...")
        continue

    plot_labeled_umap(adata, pred_dict[method], method)

print("\nLabeled UMAP plots saved with enhanced styling!")
