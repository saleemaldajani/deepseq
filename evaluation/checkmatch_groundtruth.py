# ======================================
# Imports and Setup
# ======================================

import subprocess
import sys
import os
import json

# Install rapidfuzz if not available
try:
    from rapidfuzz import fuzz
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rapidfuzz"])
    from rapidfuzz import fuzz

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Create plots/groundtruth/ directory if it doesn't exist
os.makedirs("plots/groundtruth", exist_ok=True)

# ======================================
# Marker Definitions (real ground truth)
# ======================================

expanded_markers = {
    'Naive CD4 T cells': ['IL7R', 'CCR7', 'CD3D', 'CD3E', 'SELL', 'TCF7'],
    'Memory CD4 T cells': ['IL7R', 'CD4', 'CD27', 'GZMK'],
    'Effector CD8 T cells': ['CD8A', 'CD8B', 'GZMB', 'GZMA', 'PRF1', 'NKG7'],
    'Regulatory T cells': ['FOXP3', 'CTLA4', 'IL2RA', 'PDCD1'],
    'NK cells': ['GNLY', 'NKG7', 'KLRD1', 'PRF1', 'GZMB', 'ID2', 'FCGR3A'],
    'B cells': ['MS4A1', 'CD19', 'CD79A', 'CD79B', 'CD22', 'POU2AF1'],
    'Plasma cells': ['MZB1', 'IGHG1', 'JCHAIN', 'XBP1'],
    'Classical monocytes': ['CD14', 'LYZ', 'VCAN', 'S100A8', 'S100A9'],
    'Non-classical monocytes': ['FCGR3A', 'MS4A7', 'S100A7'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A', 'LAMP3'],
    'Plasmacytoid DCs': ['IL3RA', 'CLEC4C', 'IRF7', 'IRF8'],
    'Platelets': ['PPBP', 'PF4', 'ITGA2B', 'GP9', 'SPARC'],
    'Neutrophils': ['S100A8', 'S100A9', 'FCGR3B', 'CEACAM8'],
    'Endothelial cells': ['PECAM1', 'VWF', 'ESAM', 'CD34']
}

hierarchical_mapping = {
    'Naive CD4 T cells': 'T cells',
    'Memory CD4 T cells': 'T cells',
    'Effector CD8 T cells': 'T cells',
    'Regulatory T cells': 'T cells',
    'NK cells': 'NK cells',
    'B cells': 'B cells',
    'Plasma cells': 'B cells',
    'Classical monocytes': 'Monocytes',
    'Non-classical monocytes': 'Monocytes',
    'Dendritic cells': 'Dendritic cells',
    'Plasmacytoid DCs': 'Dendritic cells',
    'Platelets': 'Platelets',
    'Neutrophils': 'Neutrophils',
    'Endothelial cells': 'Endothelial cells'
}

# ======================================
# Helper Functions
# ======================================

def normalize_gene(gene):
    if not gene:
        return ""
    gene = gene.upper().strip()
    if '-' in gene:
        gene = gene.split('-')[0]
    if 'P' in gene and gene[-1].isdigit():
        p_index = gene.find('P')
        gene = gene[:p_index]
    return gene

def simple_fuzzy_match(gene1, gene2):
    if not gene1 or not gene2:
        return False
    gene1 = normalize_gene(gene1)
    gene2 = normalize_gene(gene2)
    overlap = set(gene1) & set(gene2)
    score = len(overlap) / max(len(gene1), len(gene2))
    return score >= 0.9

def assign_ground_truth_fuzzy_pure(top_genes):
    scores = {}
    normalized_top_genes = [normalize_gene(gene) for gene in top_genes]
    
    for fine_type, markers in expanded_markers.items():
        normalized_markers = [normalize_gene(marker) for marker in markers]
        match_count = 0
        for gene in normalized_top_genes:
            for marker in normalized_markers:
                if gene == marker or simple_fuzzy_match(gene, marker):
                    match_count += 1
                    break
        scores[fine_type] = match_count

    best_type = max(scores, key=scores.get)
    best_score = scores[best_type]

    if best_score == 0:
        return "Unknown"

    broad_type = hierarchical_mapping.get(best_type, best_type)
    return broad_type

def compute_similarity(pred1, pred2):
    if pred1 is None or pred2 is None:
        return 0.0
    return fuzz.token_set_ratio(pred1, pred2) / 100

# ======================================
# Building Predictions and Ground Truth
# ======================================

def add_predictions(model_name, data_list, method_model_predictions, ground_truth_predictions):
    for item in data_list:
        method = item.get("method", "").strip()
        cluster = item.get("cluster", "").strip()
        prediction = item.get("cell_type_prediction", "").strip()
        top_genes = item.get("top_genes", [])
        
        if method and cluster:
            if model_name not in method_model_predictions[method]:
                method_model_predictions[method][model_name] = {}
            method_model_predictions[method][model_name][cluster] = prediction

            if method not in ground_truth_predictions:
                ground_truth_predictions[method] = {}
            if cluster not in ground_truth_predictions[method]:
                ground_truth_predictions[method][cluster] = assign_ground_truth_fuzzy_pure(top_genes)

def build_supertable(method_model_predictions, ground_truth_predictions, match_type="hard"):
    super_records = []

    for method, models in method_model_predictions.items():
        model_names = list(models.keys()) + ["Ground Truth"]
        for i, model1 in enumerate(model_names):
            for j, model2 in enumerate(model_names):
                if i <= j:
                    if model1 == "Ground Truth":
                        model1_predictions = ground_truth_predictions[method]
                    else:
                        model1_predictions = models.get(model1, {})
                    
                    if model2 == "Ground Truth":
                        model2_predictions = ground_truth_predictions[method]
                    else:
                        model2_predictions = models.get(model2, {})

                    all_clusters = set(model1_predictions.keys()).union(model2_predictions.keys())
                    total_clusters = len(all_clusters)

                    if match_type == "hard":
                        matches = sum(
                            1 if model1_predictions.get(cluster, None) == model2_predictions.get(cluster, None) else 0
                            for cluster in all_clusters
                        )
                    elif match_type == "soft":
                        matches = sum(
                            compute_similarity(model1_predictions.get(cluster, ""), model2_predictions.get(cluster, ""))
                            for cluster in all_clusters
                        )
                    else:
                        raise ValueError("Invalid match_type")
                    
                    match_rate = matches / total_clusters if total_clusters > 0 else 0
                    super_records.append({
                        "Method": method,
                        "Model 1": model1,
                        "Model 2": model2,
                        "Matches": round(matches, 2),
                        "Total Clusters": total_clusters,
                        "Match Rate (%)": round(100 * match_rate, 2)
                    })

    return pd.DataFrame(super_records)

def plot_heatmaps(supertable, title_suffix):
    fig, axes = plt.subplots(1, 3, figsize=(24,6))

    for ax, method in zip(axes, ["standard", "kneed", "inflection"]):
        subset = supertable[supertable["Method"] == method]
        models = sorted(set(subset["Model 1"]).union(subset["Model 2"]))

        # Move Ground Truth last
        models = [m for m in models if m != "Ground Truth"] + ["Ground Truth"]

        match_matrix = pd.DataFrame(0.0, index=models, columns=models)

        for _, row in subset.iterrows():
            match_matrix.loc[row["Model 1"], row["Model 2"]] = row["Match Rate (%)"]
            match_matrix.loc[row["Model 2"], row["Model 1"]] = row["Match Rate (%)"]

        sns.heatmap(match_matrix, annot=True, fmt=".1f", cmap="YlGnBu", ax=ax, cbar=(ax==axes[-1]))
        ax.set_title(f"{method.capitalize()} method {title_suffix}")
        ax.set_xlabel("Model")
        ax.set_ylabel("Model")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    plt.suptitle(f"Cluster Prediction Agreement Across Models - {title_suffix}", fontsize=18)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    filename = f"plots/groundtruth/match_heatmap_{title_suffix.replace(' ', '_').replace('(', '').replace(')', '')}.png"
    plt.savefig(filename, dpi=300)
    print(f"Saved plot to {filename}")

    plt.show()

# ======================================
# Main Execution
# ======================================

files_to_load = {
    "gpt-4o": "processed_clusters_gpt-4o.json",
    "gpt-4o-2024-05-13": "processed_clusters_gpt-4o-2024-05-13.json",
    "gpt-3.5-turbo": "processed_clusters_gpt35turbo.json",
    "llama3.2-1b": "processed_clusters_llama3.2-1b.json"
}

method_model_predictions = {
    "standard": {},
    "kneed": {},
    "inflection": {}
}
ground_truth_predictions = {}

for model_name, filepath in files_to_load.items():
    with open(filepath, 'r') as f:
        data = json.load(f)
    add_predictions(model_name, data, method_model_predictions, ground_truth_predictions)

# Build and save hard supertable
hard_supertable = build_supertable(method_model_predictions, ground_truth_predictions, match_type="hard")
hard_supertable.to_csv("plots/groundtruth/hard_supertable_with_groundtruth_fresh.csv", index=False)
plot_heatmaps(hard_supertable, "(Hard Match Real Ground Truth)")

# Build and save soft supertable
soft_supertable = build_supertable(method_model_predictions, ground_truth_predictions, match_type="soft")
soft_supertable.to_csv("plots/groundtruth/soft_supertable_with_groundtruth_fresh.csv", index=False)
plot_heatmaps(soft_supertable, "(Soft Fuzzy Match Real Ground Truth)")
