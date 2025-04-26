# prepare_assets.py

from pathlib import Path
import shutil

# Define paths
figures_dir = Path("figures")
assets_dir = Path("assets")

# Create required asset directories
(assets_dir / "kneeplot").mkdir(parents=True, exist_ok=True)
(assets_dir / "UMAP").mkdir(parents=True, exist_ok=True)
(assets_dir / "UMAP" / "labeled").mkdir(parents=True, exist_ok=True)

# Copy kneeplots
kneeplot_dir = figures_dir / "kneeplot"
for file in kneeplot_dir.glob("*.png"):
    shutil.copy(file, assets_dir / "kneeplot")

# Copy UMAPs
umap_dir = figures_dir / "UMAP"
for file in umap_dir.glob("*.png"):
    shutil.copy(file, assets_dir / "UMAP")

# Copy labeled UMAPs
labeled_dir = figures_dir / "UMAP" / "labeled"
for file in labeled_dir.glob("*.png"):
    shutil.copy(file, assets_dir / "UMAP" / "labeled")

print("✅ All figures copied into assets/ folder!")