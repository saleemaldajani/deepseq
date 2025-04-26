import scanpy as sc
from pathlib import Path

# Create the data directory if it doesn't exist
data_dir = Path('data')
data_dir.mkdir(exist_ok=True)

# Load the standard PBMC3k dataset
adata = sc.datasets.pbmc3k()
print(f"Loaded PBMC dataset with shape: {adata.shape}")

# Save the dataset to the 'data' folder
save_path = data_dir / 'pbmc3k.h5ad'
adata.write(save_path)

print(f"Saved PBMC dataset to {save_path}")