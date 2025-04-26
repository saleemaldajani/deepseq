# 04_classifyclusters.py (fixed version)

import json
from pathlib import Path
from langchain_ollama import OllamaLLM

# Define paths
hvg_dir = Path("hvg")
output_file = Path("processed_clusters.json")

# Initialize Ollama model
model = OllamaLLM(model="llama3.2:1b")

# Template for the classification prompt
def generate_prompt(cluster_id, genes):
    return f"""
You are given the top highly expressed genes in PBMC cluster {cluster_id}.
Top genes: {', '.join(genes)}

Based ONLY on this list of genes, what is the **cell type** most likely represented by this cluster?

ONLY respond with the cell type name.
Examples of valid answers: "T cells (activated/proliferating)", "NK cells", "Inflammatory monocytes", "B cells", "Monocytes or macrophages", "Plasmacytoid dendritic cells (pDCs)", "Cycling cells / progenitors", "Endothelial cells or progenitors", "Unknown cell type".

Do not add explanations, sentences, or commentary. Return JUST the cell type.
"""

# Gather and classify HVGs
processed = []

for method in ["standard", "kneed", "inflection"]:
    json_path = hvg_dir / method / "hvg_per_cluster.json"
    if not json_path.exists():
        print(f"Missing HVG file for method: {method}")
        continue

    with open(json_path, 'r') as f:
        cluster_genes = json.load(f)

    for cluster_id, genes in cluster_genes.items():
        prompt = generate_prompt(cluster_id, genes[:10])
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
