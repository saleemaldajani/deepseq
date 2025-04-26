# 07_deepseq.py

import dash
from dash import dcc, html
from pathlib import Path

# Initialize app
app = dash.Dash(__name__, suppress_callback_exceptions=True)
server = app.server

# --- Paths ---
assets_dir = Path("assets")
methods = ["standard", "kneed", "inflection"]

# --- Layout ---
app.layout = html.Div(
    [html.H1("DeepSeq - Labeled UMAPs", style={"textAlign": "center", "marginBottom": "40px"})]
    +
    [
        html.Div([
            html.H2(f"Labeled UMAP - {method.capitalize()}", style={"textAlign": "center"}),
            html.Img(src=f"/assets/UMAP/labeled/umap_labeled_{method}.png", style={"width": "60%", "display": "block", "margin": "auto", "marginBottom": "60px"})
        ]) for method in methods
    ]
)

if __name__ == '__main__':
    app.run(debug=True)
