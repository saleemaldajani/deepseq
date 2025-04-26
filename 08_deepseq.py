# 07_deepseq.py (final version with explicit port)

import dash
from dash import dcc, html
import plotly.express as px
import scanpy as sc
import json
from pathlib import Path
import pandas as pd

# --- Setup ---
data_dir = Path("data")
methods = ["standard", "kneed", "inflection"]

# Load processed cluster cell-type predictions
with open("processed_clusters.json", "r") as f:
    predictions = json.load(f)

# Organize predictions
pred_dict = {}
for item in predictions:
    method = item['method']
    cluster = item['cluster']
    label = item['cell_type_prediction']
    if method not in pred_dict:
        pred_dict[method] = {}
    pred_dict[method][cluster] = label

# Build plots
plots = []
for method in methods:
    h5ad_path = data_dir / f"pbmc3k_{method}_clustered.h5ad"
    if not h5ad_path.exists():
        continue

    adata = sc.read_h5ad(h5ad_path)
    umap = adata.obsm['X_umap']
    clusters = adata.obs['leiden']

    # Prepare dataframe
    df = pd.DataFrame({
        "UMAP1": umap[:, 0],
        "UMAP2": umap[:, 1],
        "Cluster": clusters
    })

    df["Cell Type"] = df["Cluster"].apply(
        lambda c: pred_dict.get(method, {}).get(f"cluster_{c}", f"cluster_{c}")
    )

    # Make interactive scatter
    fig = px.scatter(
        df, x="UMAP1", y="UMAP2", color="Cell Type",
        hover_data={"Cluster": True, "Cell Type": True},
        title=f"Labeled UMAP - {method.capitalize()}",
        width=800, height=600
    )

    plots.append(dcc.Graph(figure=fig))

# --- Build Dash app ---
app = dash.Dash(__name__)
app.layout = html.Div([
    html.H1("DeepSeq Interactive Labeled UMAPs", style={"textAlign": "center", "marginBottom": "40px"}),
    *plots
])

if __name__ == "__main__":
    app.run(
        debug=True,
        host="0.0.0.0",   # ← so it works from external access too if needed
        port=8050         # ← explicit default port
    )
