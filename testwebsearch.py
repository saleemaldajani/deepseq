import os
import asyncio
import logging

from agents import Agent, Runner, WebSearchTool, trace, set_default_openai_key

# === ENABLE DEBUG LOGGING ===
logging.basicConfig(level=logging.DEBUG)

# === SETUP ===
api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise RuntimeError("Set OPENAI_API_KEY in your environment")
set_default_openai_key(api_key)

# === AGENT DEFINITION ===
agent = Agent(
    name="TestSearcher",
    model="gpt-4o",
    instructions="""
You are a knowledgeable biology assistant. When given a list of top marker genes
for a PBMC cluster, use web search to identify which peripheral blood mononuclear cell type
(e.g., T cell, B cell, monocyte, NK cell) is most likely represented by these markers.
Provide a concise answer naming the cell type on the first line.
Then, under a “Sources:” heading, list the URLs you consulted—one per line.
Do not include any other commentary.
""",
    tools=[WebSearchTool(user_location={"type": "approximate", "city": "New York"})],
)

# === SINGLE QUERY ===
async def main():
    genes = ["CD3D", "CD3E", "CD8A", "CD4", "IL7R"]  # example T-cell markers
    prompt = (
        f"Cluster TEST top genes: {', '.join(genes)}.\n"
        "Search the web to identify the PBMC cell type. "
        "Reply first with the cell type name, then list your sources under “Sources:”."
    )

    with trace("test_cluster_search_with_sources"):
        response = await Runner.run(agent, prompt)

    print("\n=== AGENT OUTPUT ===")
    print(response.final_output)

if __name__ == "__main__":
    asyncio.run(main())
