import asyncio
from agents import Agent, Runner, set_default_openai_key
import os

api_key = os.environ.get("OPENAI_API_KEY")
print(f"OPENAI_API_KEY {api_key}")
set_default_openai_key(api_key)

# get OPEANAI_API_KEY from codespaces secrets

    # This is the default and can be omitted



async def main():
    agent = Agent(
        name="Assistant",
        instructions="You only respond in haikus.",
    )

    result = await Runner.run(agent, "Tell me about recursion in programming.")
    print(result.final_output)
    # Function calls itself,
    # Looping in smaller pieces,
    # Endless by design.


if __name__ == "__main__":
    asyncio.run(main())
