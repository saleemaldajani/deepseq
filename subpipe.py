# run_pipeline.py

import subprocess
import webbrowser
import time
import os
import sys


# Launch DeepSeq Dash app (this will stay running)
print("🚀 Launching DeepSeq Dash app...\n")
subprocess.run(["python", "08_deepseq.py"], check=True)

# (No more code here — it will stay serving the app)