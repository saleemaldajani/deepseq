# 04_classifyclusters.py (Reads dataset names dynamically - Fixed method extraction)

import json
from pathlib import Path
from langchain_ollama import OllamaLLM
import scanpy as sc

# Define paths
hvg_dir = Path("hvg")
data_dir = Path("data")
output_file = Path("processed_clusters.json")

# Initialize Ollama model
model = OllamaLLM(model="llama3.2:1b")

# Template for the classification prompt
def generate_prompt(dataset_name, cluster_id, genes):
    return f"""
You are given the top highly expressed genes in a "{dataset_name}" cluster {cluster_id}.
Top genes: {', '.join(genes)}

Based ONLY on this list of genes, what is the **cell type** most likely represented by this cluster?

ONLY respond with the cell type name.
Do not add explanations, sentences, or commentary. Return JUST the cell type.
"""

# Gather and classify HVGs
processed = []

# Find all datasets dynamically
h5ad_files = list(data_dir.glob("*.h5ad"))

for h5ad_path in h5ad_files:
    method = h5ad_path.stem
    if method.endswith("_clustered"):
        method = method.replace("_clustered", "")
    if method.startswith("pbmc3k_"):
        method = method.replace("pbmc3k_", "")

    json_path = hvg_dir / method / "hvg_per_cluster.json"

    if not json_path.exists():
        print(f"Missing HVG file for method: {method}")
        continue

    # Read dataset name from .h5ad
    adata = sc.read_h5ad(h5ad_path)
    dataset_name = adata.uns.get("dataset_name", method)

    with open(json_path, 'r') as f:
        cluster_genes = json.load(f)

    for cluster_id, genes in cluster_genes.items():
        prompt = generate_prompt(dataset_name, cluster_id, genes[:10])
        try:
            result = model.invoke(input=prompt).strip()
        except Exception as e:
            result = f"[ERROR: {e}]"

        processed.append({
            "method": method,
            "cluster": cluster_id,
            "top_genes": genes[:10],
            "cell_type_prediction": result
        })

# Save all classifications
with open(output_file, 'w') as f:
    json.dump(processed, f, indent=2)

print(f"Classification complete. Results saved to {output_file}")