# run_pipeline.py

import subprocess
import webbrowser
import time
import os
import sys

# --- Scripts to run ---
scripts = [
    "00_installrequirements.py",
    "01_geth5ad.py",
    "02_filteredh5ad.py",
    "03_clusterh5ad.py",
    "04_classifyclusters.py",
    "05_labelumap.py",
    "06_generatereport.py",
    "07_prepareassets.py"
]

# --- Run all preprocessing scripts ---
for script in scripts:
    print(f"\n🚀 Running {script}...")
    if not os.path.exists(script):
        print(f"❌ {script} not found! Skipping.")
        continue

    result = subprocess.run(["python", script], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ Error running {script}:")
        print(result.stderr)
        sys.exit(1)
    else:
        print(f"✅ Completed {script}")

# --- Launch Dash app in background ---
print("\n🚀 Launching DeepSeq Dash app...")
dash_process = subprocess.Popen(["python", "08_deepseq.py"])

# --- Open browser to Dash app URL ---
time.sleep(3)  # Give the server a few seconds to start
webbrowser.open("http://127.0.0.1:8050")

print("\n✅ Pipeline completed. Browser should open to DeepSeq dashboard!")
