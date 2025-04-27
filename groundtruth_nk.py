# ======================================
# Script adapted from ChatGPT discussion:
# https://chatgpt.com/share/680e04cf-9380-8008-b786-7c658359afe6
# ======================================

# --- Setup ---

import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Create output directory
os.makedirs("outputs", exist_ok=True)

# Load your uploaded JSON files
uploaded_files = {
    'gpt-4o': 'processed_clusters_gpt-4o.json',
    'gpt-4o-2024-05-13': 'processed_clusters_gpt-4o-2024-05-13.json',
    'gpt-3.5-turbo': 'processed_clusters_gpt35turbo.json',
    'llama3.2-1b': 'processed_clusters_llama3.2-1b.json'
}

# --- Marker Definitions ---

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

# --- Helper Functions ---

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
        return "Unknown", 0, "Unknown"

    broad_type = hierarchical_mapping.get(best_type, best_type)
    return best_type, best_score, broad_type

def better_normalize_label(label):
    if not label:
        return 'unknown'
    label = label.lower().replace('cells', '').replace('cell', '').replace(' ', '').replace('+', '').replace('(', '').replace(')', '')
    if 'natural' in label and 'nk' in label:
        return 'nkcell'
    if 'b' in label and 't' not in label:
        return 'bcell'
    if 'cd4' in label or 'cd8' in label or 't' in label:
        return 'tcell'
    if 'nk' in label:
        return 'nkcell'
    if 'monocyte' in label:
        return 'monocyte'
    if 'dendritic' in label or 'dc' in label:
        return 'dendriticcell'
    if 'platelet' in label:
        return 'platelet'
    if 'plasma' in label:
        return 'bcell'
    if 'neutrophil' in label:
        return 'neutrophil'
    if 'endothelial' in label:
        return 'endothelial'
    if 'lymphocyte' in label:
        return 'tcell'
    return label

# --- Evaluation ---

fuzzy_summary_pure = []
full_cluster_records = []

for model_name, filepath in uploaded_files.items():
    with open(filepath, 'r') as f:
        cluster_data = json.load(f)

    methods_dict = {}

    for cluster in cluster_data:
        method = cluster['method']
        if method not in methods_dict:
            methods_dict[method] = {'ground_truth': [], 'predicted': [], 'records': []}

        top_genes = cluster.get('top_genes', [])
        predicted_label = cluster.get('cell_type_prediction', 'unknown')
        fine_gt, match_score, broad_gt = assign_ground_truth_fuzzy_pure(top_genes)

        predicted_label_norm = better_normalize_label(predicted_label)
        broad_gt_norm = better_normalize_label(broad_gt)

        is_correct = (predicted_label_norm == broad_gt_norm)

        methods_dict[method]['ground_truth'].append(broad_gt_norm)
        methods_dict[method]['predicted'].append(predicted_label_norm)
        methods_dict[method]['records'].append({
            'Cluster': cluster['cluster'],
            'Predicted': predicted_label,
            'Ground Truth (Fine)': fine_gt,
            'Ground Truth (Broad)': broad_gt,
            'Correct': is_correct,
            'Matching Markers': match_score
        })

        full_cluster_records.append({
            'Model': model_name,
            'Method': method,
            'Cluster': cluster['cluster'],
            'Top Genes': top_genes,
            'Predicted': predicted_label,
            'Ground Truth (Fine)': fine_gt,
            'Ground Truth (Broad)': broad_gt,
            'Matching Markers': match_score,
            'Correct': is_correct
        })

    for method, data in methods_dict.items():
        total = len(data['records'])
        correct = sum(record['Correct'] for record in data['records'])
        accuracy = correct / total if total > 0 else 0
        fuzzy_summary_pure.append({
            'Model': model_name,
            'Method': method,
            'Total Clusters': total,
            'Correct Clusters': correct,
            'Accuracy (%)': round(accuracy * 100, 2)
        })

# --- Save Outputs ---

summary_df = pd.DataFrame(fuzzy_summary_pure)
summary_df.to_csv('outputs/fuzzy_matching_summary.csv', index=False)

full_cluster_df = pd.DataFrame(full_cluster_records)
full_cluster_df.to_csv('outputs/full_cluster_groundtruth_summary.csv', index=False)

# --- Plot 1: Full sorted accuracy barplot (all methods, no bar labels) ---

sorted_summary = summary_df.sort_values(by="Accuracy (%)", ascending=False)

plt.figure(figsize=(12,6))
barplot = sns.barplot(
    data=sorted_summary,
    x="Model",
    y="Accuracy (%)",
    hue="Method",
    palette="viridis",
    dodge=True
)

# Remove percentage labels on bars

plt.ylim(0, 100)
plt.title("Model and Method Accuracy Comparison (All Methods)", fontsize=16)
plt.ylabel("Accuracy (%)", fontsize=14)
plt.xlabel("")
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title="Method", title_fontsize=12, fontsize=10)
plt.tight_layout()
plt.savefig("outputs/model_accuracy_comparison_polished_sorted_nolabels.png", dpi=300)
plt.show()

# --- Plot 2: Sorted accuracy without Inflection ---

filtered_summary = summary_df[~summary_df['Method'].str.contains('inflection', case=False)]
sorted_filtered = filtered_summary.sort_values(by="Accuracy (%)", ascending=False)

plt.figure(figsize=(12,6))
barplot = sns.barplot(
    data=sorted_filtered,
    x="Model",
    y="Accuracy (%)",
    hue="Method",
    palette="viridis",
    dodge=True
)

# Remove percentage labels on bars

plt.ylim(0, 100)
plt.title("Model and Method Accuracy Comparison (No Inflection)", fontsize=16)
plt.ylabel("Accuracy (%)", fontsize=14)
plt.xlabel("")
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title="Method", title_fontsize=12, fontsize=10)
plt.tight_layout()
plt.savefig("outputs/model_accuracy_comparison_polished_sorted_no_inflection_nolabels.png", dpi=300)
plt.show()

print("\n✅ Full evaluation completed and all results saved in /outputs/")
