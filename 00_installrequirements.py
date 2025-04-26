# 00_installrequirements.py

import subprocess
import os
import sys

requirements_file = "requirements.txt"
ollama_model = "llama3.2:1b"

# --- Install Python packages ---
if not os.path.exists(requirements_file):
    print(f"❌ {requirements_file} not found.")
    sys.exit(1)

print(f"🚀 Installing Python requirements from {requirements_file}...\n")

try:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", requirements_file])
    print("\n✅ Python packages installed successfully!")
except subprocess.CalledProcessError as e:
    print("\n❌ Error during pip install!")
    print(e)
    sys.exit(1)

# --- Pull Ollama model ---
print(f"\n🚀 Pulling Ollama model '{ollama_model}'...\n")

try:
    subprocess.check_call(["ollama", "pull", ollama_model])
    print(f"\n✅ Ollama model '{ollama_model}' pulled successfully!")
except subprocess.CalledProcessError as e:
    print("\n❌ Error during Ollama pull!")
    print("Make sure Ollama server is running.")
    print(e)
    sys.exit(1)
except FileNotFoundError:
    print("\n❌ Ollama CLI not found! Please install Ollama and try again.")
    sys.exit(1)
