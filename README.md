
# deepseq

DeepSeq is a complete pipeline for analyzing PBMC (Peripheral Blood Mononuclear Cells) single-cell RNA sequencing data.  
It automates data filtering, clustering, cell-type classification using LLMs (Ollama), UMAP visualization, report generation, and launches a Dash web application for interactive exploration.

---

## 📋 Pipeline Overview

| Step | Script | Description |
|:----|:-----|:------------------------------------------------|
| 0 | `00_installrequirements.py` | Install Python packages and pull Ollama model |
| 1 | `01_geth5ad.py` | Download standard PBMC3K dataset |
| 2 | `02_filteredh5ad.py` | Filter cells based on UMI/gene thresholds |
| 3 | `03_clusterh5ad.py` | Cluster cells using Leiden algorithm |
| 4 | `04_classifyclusters.py` | Classify clusters using LLM (Ollama with `llama3.2:1b`) |
| 5 | `05_labelumap.py` | Generate labeled UMAP plots |
| 6 | `06_generatereport.py` | Generate summary HTML and PDF reports |
| 7 | `07_prepareassets.py` | Move figures into Dash assets directory |
| 8 | `08_deepseq.py` | Launch Dash web application |

---

## 🚀 How to Run

Just run:

```bash
python 09_runpipeline.py
```

This will:
- Install requirements
- Pull Ollama `llama3.2:1b`
- Process the dataset
- Classify clusters
- Generate plots
- Launch your dashboard automatically at:

```bash
http://127.0.0.1:8050/
```

---

## ✨ Example Output Snippets

### 00_installrequirements.py
```
🚀 Installing Python requirements from requirements.txt...
✅ Python packages installed successfully!

🚀 Pulling Ollama model 'llama3.2:1b'...
✅ Ollama model 'llama3.2:1b' pulled successfully!
```

---

### 02_filteredh5ad.py
```
Loaded data/pbmc3k_raw.h5ad with shape (2700, 32738)

Summary Table:
        Label  Cells   Mean Genes  Median Genes    Mean UMIs  Median UMIs
0    Standard   2700         847           817          2367        2197
1       Kneed   2695         848           817          2370        2200
2  Inflection     19        2223          1997          8584        7171
```

---

### 04_classifyclusters.py
```
Classifying cluster_0 in method 'standard'...
Classifying cluster_1 in method 'standard'...
...

✅ Classification complete. Results saved to processed_clusters.json
```

---

### 06_generatereport.py
```
✅ Generated summary HTML and PDF reports in 'reports/' folder
```

---

### 08_deepseq.py
```
Dash is running on http://127.0.0.1:8050/
 * Serving Flask app '08_deepseq'
 * Debug mode: on
```

---

## 📈 Final Products

After running the pipeline:
- **Interactive Web App**: Dash visualization of clusters
- **HTML and PDF Reports**: Summary metrics and labeled UMAPs
- **CSV Exports**: UMAP coordinates and cluster predictions
- **Labeled Figures**: Saved under `figures/UMAP/labeled/`

---

## 📦 Project Structure

```
deepseq/
├── 00_installrequirements.py
├── 01_geth5ad.py
├── 02_filteredh5ad.py
├── 03_clusterh5ad.py
├── 04_classifyclusters.py
├── 05_labelumap.py
├── 06_generatereport.py
├── 07_prepareassets.py
├── 08_deepseq.py
├── 09_run_pipeline.py
├── requirements.txt
├── data/
├── figures/
├── reports/
├── assets/
```

---

## 🧠 Notes

- If Ollama is not running, `04_classifyclusters.py` will fail — make sure Ollama server is active.
- You can modify clustering parameters (resolution, neighbors) in `03_clusterh5ad.py`.
- You can switch LLMs (e.g., `llama3:instruct`) easily in `04_classifyclusters.py`.

---

## 👏 Credits

- Built with ❤️ using [Scanpy](https://scanpy.readthedocs.io/), [Dash](https://dash.plotly.com/), [Ollama](https://ollama.com/), and [Plotly](https://plotly.com/python/).

--- 

## Evaluation 

checkmatch_groundtruth.py
groundtruth_nk.py