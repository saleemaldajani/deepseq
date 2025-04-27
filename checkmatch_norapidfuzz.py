# ======================================
# Imports and Setup
# ======================================

import os
import json
import subprocess
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Create plots/ directory if it doesn't exist
os.makedirs("plots", exist_ok=True)

# ======================================
# Helper Functions
# ======================================

# Simple fuzzy match based on word overlap
def simple_fuzzy_similarity(pred1, pred2):
    if not pred1 or not pred2:
        return 0.0
    set1 = set(pred1.lower().replace('(', '').replace(')', '').replace(',', '').split())
    set2 = set(pred2.lower().replace('(', '').replace(')', '').split())
    if not set1 or not set2:
        return 0.0
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union)

def build_supertable(method_model_predictions, match_type="hard"):
    super_records = []

    for method, models in method_model_predictions.items():
        model_names = list(models.keys())
        for i, model1 in enumerate(model_names):
            for j, model2 in enumerate(model_names):
                if i <= j:  # only upper triangle + diagonal
                    total_clusters = len(set(models[model1].keys()).union(models[model2].keys()))
                    if match_type == "hard":
                        matches = sum(
                            1 if models[model1].get(cluster, None) == models[model2].get(cluster, None) else 0
                            for cluster in set(models[model1].keys()).union(models[model2].keys())
                        )
                    elif match_type == "soft":
                        matches = sum(
                            simple_fuzzy_similarity(models[model1].get(cluster, ""), models[model2].get(cluster, ""))
                            for cluster in set(models[model1].keys()).union(models[model2].keys())
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
    fig, axes = plt.subplots(1, 3, figsize=(18,6))

    for ax, method in zip(axes, ["standard", "kneed", "inflection"]):
        subset = supertable[supertable["Method"] == method]
        models = sorted(set(subset["Model 1"]).union(subset["Model 2"]))
        match_matrix = pd.DataFrame(0.0, index=models, columns=models)

        for _, row in subset.iterrows():
            match_matrix.loc[row["Model 1"], row["Model 2"]] = row["Match Rate (%)"]
            match_matrix.loc[row["Model 2"], row["Model 1"]] = row["Match Rate (%)"]

        sns.heatmap(match_matrix, annot=True, fmt=".1f", cmap="YlGnBu", ax=ax, cbar=(ax==axes[-1]))
        ax.set_title(f"{method.capitalize()} method {title_suffix}")
        ax.set_xlabel("Model")
        ax.set_ylabel("Model")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)  # <- Horizontal y-axis labels

    plt.suptitle(f"Cluster Prediction Agreement Across Models - {title_suffix}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    filename = f"plots/match_heatmap_{title_suffix.replace(' ', '_').replace('(', '').replace(')', '')}.png"
    plt.savefig(filename, dpi=300)
    print(f"Saved plot to {filename}")

    plt.show()

# ======================================
# Load your JSON files from /home
# ======================================

def add_predictions(model_name, data_list, method_model_predictions):
    for item in data_list:
        method = item.get("method", "").strip()
        cluster = item.get("cluster", "").strip()
        prediction = item.get("cell_type_prediction", "").strip()
        if method and cluster:
            if model_name not in method_model_predictions[method]:
                method_model_predictions[method][model_name] = {}
            method_model_predictions[method][model_name][cluster] = prediction

# Initialize storage
method_model_predictions = {
    "standard": {},
    "kneed": {},
    "inflection": {}
}

# Load all four files
files_to_load = {
    "gpt-4o": "processed_clusters_gpt-4o.json",
    "gpt-4o-2024-05-13": "processed_clusters_gpt-4o-2024-05-13.json",
    "gpt-3.5-turbo": "processed_clusters_gpt35turbo.json",
    "llama3-1b": "processed_clusters_llama31-1b.json"
}

for model_name, filepath in files_to_load.items():
    with open(filepath, 'r') as f:
        data = json.load(f)
    add_predictions(model_name, data, method_model_predictions)

# ======================================
# Run Supertable and Plotting
# ======================================

# Hard matching (exact string match)
hard_supertable = build_supertable(method_model_predictions, match_type="hard")
hard_supertable.to_csv("plots/hard_supertable_no_rapidfuzz.csv", index=False)
plot_heatmaps(hard_supertable, "(Hard Match No RapidFuzz)")

# Soft matching (simple word-overlap)
soft_supertable = build_supertable(method_model_predictions, match_type="soft")
soft_supertable.to_csv("plots/soft_supertable_no_rapidfuzz.csv", index=False)
plot_heatmaps(soft_supertable, "(Soft Match No RapidFuzz)")
