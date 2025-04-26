# 06_generatereport.py (Extended Version)

import json
from pathlib import Path
import pandas as pd
from fpdf import FPDF
from jinja2 import Environment, FileSystemLoader
import shutil

# --- Paths ---
figures_dir = Path("figures")
kneeplot_dir = figures_dir / "kneeplot"
umap_dir = figures_dir / "UMAP"
labeled_dir = umap_dir / "labeled"
hvg_dir = Path("hvg")
output_dir = Path("reports")
output_dir.mkdir(parents=True, exist_ok=True)

# --- Load processed clusters ---
with open("processed_clusters.json", "r") as f:
    cluster_data = json.load(f)

# --- Group by method ---
method_groups = {}
for item in cluster_data:
    method = item['method']
    if method not in method_groups:
        method_groups[method] = []
    method_groups[method].append(item)

# --- Load summary metrics ---
summary_metrics_path = kneeplot_dir / "summary_metrics.csv"
if summary_metrics_path.exists():
    summary_df = pd.read_csv(summary_metrics_path)
else:
    summary_df = pd.DataFrame()

# --- Setup Jinja2 Template ---
env = Environment(loader=FileSystemLoader("templates"))
html_template = env.get_template("report_template.html")
combined_template = env.get_template("combined_report_template.html")

# --- Generate reports per method ---
for method, clusters in method_groups.items():
    # Load figures
    kneeplot_path = kneeplot_dir / method / "kneeplot_summary.png"
    umap_unlabeled_path = umap_dir / f"umap_{method}.png"
    umap_labeled_path = labeled_dir / f"umap_labeled_{method}.png"

    # Create cluster summary table
    table_data = []
    for cluster in clusters:
        table_data.append({
            "Cluster": cluster['cluster'],
            "Top Genes": ", ".join(cluster['top_genes']),
            "Cell Type": cluster['cell_type_prediction']
        })
    df = pd.DataFrame(table_data)

    # Render HTML
    html_content = html_template.render(
        method=method.capitalize(),
        kneeplot_img=kneeplot_path,
        umap_unlabeled_img=umap_unlabeled_path,
        umap_labeled_img=umap_labeled_path,
        table=df.to_html(index=False, escape=False)
    )

    # Save HTML
    html_file = output_dir / f"report_{method}.html"
    with open(html_file, "w") as f:
        f.write(html_content)

    # Generate PDF (basic text capture)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt=f"Report: {method.capitalize()}", ln=True, align='C')
    pdf.ln(10)
    pdf.multi_cell(0, 10, df.to_string(index=False))
    pdf.output(output_dir / f"report_{method}.pdf")

# --- Generate Combined Report ---
comparison_table = []
all_clusters = sorted(set([item['cluster'] for item in cluster_data]))

for cluster in all_clusters:
    row = {"Cluster": cluster}
    for method in method_groups:
        label = next((item['cell_type_prediction'] for item in method_groups[method] if item['cluster'] == cluster), "-")
        row[method.capitalize()] = label
    comparison_table.append(row)

comparison_df = pd.DataFrame(comparison_table)

# Collect all images side by side (UMAPs and labeled UMAPs)
umap_imgs = {method: f"figures/UMAP/umap_{method}.png" for method in method_groups}
labeled_imgs = {method: f"figures/UMAP/labeled/umap_labeled_{method}.png" for method in method_groups}
kneeplot_imgs = {method: f"figures/kneeplot/{method}/kneeplot_summary.png" for method in method_groups}

combined_html = combined_template.render(
    summary_table=summary_df.to_html(index=False, escape=False),
    comparison_table=comparison_df.to_html(index=False, escape=False),
    umap_imgs=umap_imgs,
    labeled_imgs=labeled_imgs,
    kneeplot_imgs=kneeplot_imgs
)

with open(output_dir / "combined_report.html", "w") as f:
    f.write(combined_html)

print("\nAll reports (individual and combined) generated!")

# --- Note ---
# You need two templates:
# 1. report_template.html (for individual reports)
# 2. combined_report_template.html (for the full combined report)
# Both should be placed inside 'templates/' folder.
