from openai import OpenAI
import os

# Create a v1 client
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# List all accessible models
models = client.models.list()

# Print their IDs
for m in models:
    print(m.id)
