import json
import os
import asyncio
from pathlib import Path

from agents import Agent, Runner, WebSearchTool, trace, set_default_openai_key

# Set up OpenAI API key
api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise RuntimeError("Please set the OPENAI_API_KEY environment variable.")
set_default_openai_key(api_key)

# Define paths
base_dir = Path(__file__).parent
hvg_dir = base_dir / "hvg"
output_file = base_dir / "processed_clusters_gpt35turbo.json"

# Ensure hvg directory exists
if not hvg_dir.exists() or not hvg_dir.is_dir():
    print(f"Error: HVG directory '{hvg_dir}' not found.")
    exit(1)

# Create classification agent
agent = Agent(
    name="CellTypeClassifier",
    instructions=(
        """
        You are a knowledgeable biology assistant. When given a list of top marker genes
        for a PBMC cluster, use web search to identify which peripheral blood mononuclear cell type
        (e.g., T cell, B cell, monocyte, NK cell) is most likely represented by these markers.
        Provide a concise answer naming the cell type.
        """
    ),
    tools=[WebSearchTool(user_location={"type": "approximate", "city": "New York"})],
)

# Helper to generate search prompt

def generate_prompt(cluster_id, genes):
    gene_list = ", ".join(genes)
    return (
        f"Cluster {cluster_id} top genes: {gene_list}. "
        "Search the web to determine the most likely PBMC cell type represented by these markers. "
        "Reply with only the cell type name."
    )

async def main():
    processed = []
    methods = ["standard", "kneed", "inflection"]

    for method in methods:
        json_path = hvg_dir / method / "hvg_per_cluster.json"
        if not json_path.exists():
            print(f"Missing HVG file for method: {method} at {json_path}")
            continue
        print(f"\nProcessing method '{method}'...")

        try:
            cluster_genes = json.loads(json_path.read_text())
        except Exception as e:
            print(f"Failed to load JSON for {method}: {e}")
            continue

        if not cluster_genes:
            print(f"No clusters found in {json_path}")
            continue

        for cluster_id, genes in cluster_genes.items():
            top_genes = genes[:10]
            print(f"Classifying {method} cluster {cluster_id}: {top_genes}")
            prompt = generate_prompt(cluster_id, top_genes)
            try:
                with trace(f"classify_{method}_{cluster_id}"):
                    result = await Runner.run(agent, prompt)
                cell_type = result.final_output.strip()
                print(f"-> Predicted: {cell_type}")
            except Exception as e:
                cell_type = f"[ERROR: {e}]"
                print(f"Error for {method} cluster {cluster_id}: {e}")

            processed.append({
                "method": method,
                "cluster": cluster_id,
                "top_genes": top_genes,
                "cell_type_prediction": cell_type
            })

    # Save results
    try:
        output_file.write_text(json.dumps(processed, indent=2))
        print(f"\nClassification complete. {len(processed)} entries saved to {output_file}")
    except Exception as e:
        print(f"Failed to write output file: {e}")

if __name__ == "__main__":
    asyncio.run(main())
