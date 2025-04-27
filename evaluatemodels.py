# ======================================
# Imports and Setup
# Script adapted from ChatGPT discussion:
# https://chatgpt.com/share/680e0671-ad28-8008-b828-ddee21e5761b
# ======================================

import os
import json
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
from sklearn.metrics import confusion_matrix, adjusted_rand_score, normalized_mutual_info_score
from scipy.optimize import linear_sum_assignment

# ======================================
# Helper Functions
# ======================================

def simple_similarity(pred1, pred2):
    if not pred1 or not pred2:
        return 0.0
    set1 = set(pred1.lower().replace('+', '').replace('-', '').replace(',', '').split())
    set2 = set(pred2.lower().replace('+', '').replace('-', '').replace(',', '').split())
    if not set1 or not set2:
        return 0.0
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union)

def build_confusion_fuzzy(y_true, y_pred, classes_true, classes_pred, threshold=0.6):
    """Build fuzzy confusion matrix based on word overlap."""
    conf_mat = np.zeros((len(classes_true), len(classes_pred)), dtype=int)
    for t, p in zip(y_true, y_pred):
        best_score = 0.0
        best_j = -1
        for j, pred_class in enumerate(classes_pred):
            score = simple_similarity(t, pred_class)
            if score > best_score:
                best_score = score
                best_j = j
        i = list(classes_true).index(t)
        if best_score >= threshold and best_j != -1:
            conf_mat[i, best_j] += 1
    return conf_mat

def evaluate_model(y_true, y_pred, classes_true, classes_pred, model_name, method_name,
                   save_dir, use_fuzzy=False, threshold=0.6):
    """Evaluate model and save confusion matrix + metrics."""
    if use_fuzzy:
        confusion = build_confusion_fuzzy(y_true, y_pred, classes_true, classes_pred, threshold)
    else:
        confusion = confusion_matrix(y_true, y_pred, labels=classes_true)

    cost_matrix = -confusion
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    confusion_aligned = confusion[:, col_ind]
    matched_pred_classes = classes_pred[col_ind]

    # Remap predictions
    label_mapping = {p: t for p, t in zip(matched_pred_classes, classes_true)}
    mapped_preds = [label_mapping.get(p, p) for p in y_pred]

    # Metrics
    ari = adjusted_rand_score(y_true, mapped_preds)
    nmi = normalized_mutual_info_score(y_true, mapped_preds)

    # Save confusion matrix
    confusion_df = pd.DataFrame(confusion_aligned, index=classes_true, columns=matched_pred_classes)
    os.makedirs(save_dir, exist_ok=True)
    confusion_df.to_csv(f"{save_dir}/confusion_matrix_{model_name}_{method_name}.csv")

    # Save metrics
    with open(f"{save_dir}/metrics_{model_name}_{method_name}.txt", "w") as f:
        f.write(f"Adjusted Rand Index (ARI): {ari:.3f}\n")
        f.write(f"Normalized Mutual Information (NMI): {nmi:.3f}\n")

    # Save heatmap
    plt.figure(figsize=(10,8))
    sns.heatmap(confusion_aligned, annot=True, fmt="d", cmap="Blues",
                xticklabels=matched_pred_classes, yticklabels=classes_true)
    plt.xlabel("Predicted Labels (matched)")
    plt.ylabel("Ground Truth Labels")
    match_mode = "Fuzzy" if use_fuzzy else "Hard"
    plt.title(f"Confusion Matrix ({match_mode})\n{model_name} - {method_name}")
    plt.tight_layout()
    plt.savefig(f"{save_dir}/confusion_matrix_{model_name}_{method_name}.png", dpi=300)
    plt.close()

    return ari, nmi

def build_summary_table(metrics_folder, output_csv):
    """Build a summary table from metrics files."""
    records = []
    for filepath in glob.glob(f"{metrics_folder}/metrics_*.txt"):
        filename = os.path.basename(filepath)
        parts = filename.replace(".txt", "").split("_")
        model_name = parts[1]
        method_name = parts[2]

        with open(filepath, "r") as f:
            lines = f.readlines()
            ari = float(lines[0].split(":")[1].strip())
            nmi = float(lines[1].split(":")[1].strip())

        records.append({
            "Model": model_name,
            "Method": method_name,
            "ARI": ari,
            "NMI": nmi
        })

    summary_df = pd.DataFrame(records)
    summary_df.to_csv(output_csv, index=False)
    return summary_df

# ======================================
# Step 1: Load pbmc3k and Preprocess
# ======================================

print("Loading pbmc3k...")
adata = sc.datasets.pbmc3k()

sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.5)
adata.obs['leiden'] = adata.obs['leiden'].astype(int)

# Ground truth assignment
cluster_to_label = {
    0: "CD14+ Monocyte",
    1: "CD4+ T cell",
    2: "NK cell",
    3: "B cell",
    4: "CD8+ T cell",
    5: "Dendritic cell",
    6: "Platelet",
    7: "FCGR3A+ Monocyte",
}
adata.obs['manual_labels'] = adata.obs['leiden'].map(cluster_to_label)

true_labels = adata.obs['manual_labels'].astype(str).values
classes_true = np.unique(true_labels)

# ======================================
# Step 2: Load Your Models
# ======================================

print("Loading model predictions...")

files_to_load = {
    "gpt-4o": "processed_clusters_gpt-4o.json",
    "gpt-4o-2024-05-13": "processed_clusters_gpt-4o-2024-05-13.json",
    "gpt-3.5-turbo": "processed_clusters_gpt35turbo.json",
    "llama3-1b": "processed_clusters_llama31-1b.json"
}

model_method_predictions = {}

for model_name, filepath in files_to_load.items():
    with open(filepath, 'r') as f:
        data = json.load(f)
    for item in data:
        method = item.get("method", "").strip()
        cluster = item.get("cluster", "").strip()
        prediction = item.get("cell_type_prediction", "").strip()
        if method and cluster:
            key = (model_name, method)
            if key not in model_method_predictions:
                model_method_predictions[key] = {}
            model_method_predictions[key][cluster] = prediction

# ======================================
# Step 3: Run Both Hard and Fuzzy Evaluations
# ======================================

print("Evaluating all models (hard match)...")

for (model_name, method_name), preds in model_method_predictions.items():
    pred_labels = adata.obs['leiden'].apply(lambda x: preds.get(f"cluster_{x}", "Unknown")).astype(str).values
    classes_pred = np.unique(pred_labels)
    evaluate_model(true_labels, pred_labels, classes_true, classes_pred,
                   model_name, method_name, save_dir="plots/hard", use_fuzzy=False)

print("\nEvaluating all models (fuzzy match)...")

for (model_name, method_name), preds in model_method_predictions.items():
    pred_labels = adata.obs['leiden'].apply(lambda x: preds.get(f"cluster_{x}", "Unknown")).astype(str).values
    classes_pred = np.unique(pred_labels)
    evaluate_model(true_labels, pred_labels, classes_true, classes_pred,
                   model_name, method_name, save_dir="plots/fuzzy", use_fuzzy=True)

# ======================================
# Step 4: Summarize and Calculate Boost
# ======================================

hard_summary = build_summary_table("plots/hard", "plots/summary_metrics_hard.csv")
fuzzy_summary = build_summary_table("plots/fuzzy", "plots/summary_metrics_fuzzy.csv")

# Merge summaries
boost_df = pd.merge(hard_summary, fuzzy_summary, on=["Model", "Method"], suffixes=("_hard", "_fuzzy"))
boost_df["ARI Boost"] = boost_df["ARI_fuzzy"] - boost_df["ARI_hard"]
boost_df["NMI Boost"] = boost_df["NMI_fuzzy"] - boost_df["NMI_hard"]
boost_df.to_csv("plots/boost_summary.csv", index=False)

print("\n✅ Finished: Hard vs Fuzzy matching comparison complete!")
print("Summary saved in: plots/boost_summary.csv")

# Optional: quick view of top improvements
print("\nTop ARI Boost:")
print(boost_df.sort_values(by="ARI Boost", ascending=False)[["Model", "Method", "ARI Boost"]].head())
