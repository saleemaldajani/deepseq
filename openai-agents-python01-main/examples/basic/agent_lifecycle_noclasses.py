import os
import random
import asyncio
from types import SimpleNamespace
from agents import Agent, Runner, function_tool, set_default_openai_key
from pydantic import create_model

# Set and print the OpenAI API key.
api_key = os.environ.get("OPENAI_API_KEY")
print(f"OPENAI_API_KEY {api_key}")
set_default_openai_key(api_key)

def create_hooks(display_name: str):
    """
    Creates an object containing asynchronous hook methods with a shared event counter.
    The returned object has attributes: on_start, on_end, on_handoff, on_tool_start, on_tool_end.
    """
    counter = [0]  # mutable counter stored as a list

    async def on_start(context, agent):
        counter[0] += 1
        print(f"### ({display_name}) {counter[0]}: Agent {agent.name} started")

    async def on_end(context, agent, output):
        counter[0] += 1
        print(f"### ({display_name}) {counter[0]}: Agent {agent.name} ended with output {output}")

    async def on_handoff(context, agent, source):
        counter[0] += 1
        print(f"### ({display_name}) {counter[0]}: Agent {source.name} handed off to {agent.name}")

    async def on_tool_start(context, agent, tool):
        counter[0] += 1
        print(f"### ({display_name}) {counter[0]}: Agent {agent.name} started tool {tool.name}")

    async def on_tool_end(context, agent, tool, result):
        counter[0] += 1
        print(f"### ({display_name}) {counter[0]}: Agent {agent.name} ended tool {tool.name} with result {result}")

    return SimpleNamespace(
        on_start=on_start,
        on_end=on_end,
        on_handoff=on_handoff,
        on_tool_start=on_tool_start,
        on_tool_end=on_tool_end,
    )

@function_tool
def random_number(max: int) -> int:
    """
    Generate a random number up to the provided maximum.
    """
    return random.randint(0, max)

@function_tool
def multiply_by_two(x: int) -> int:
    """Simple multiplication by two."""
    return x * 2

# Dynamically create a Pydantic model for the final result.
FinalResult = create_model("FinalResult", number=(int, ...))

multiply_agent = Agent(
    name="Multiply Agent",
    instructions="Multiply the number by 2 and then return the final result.",
    tools=[multiply_by_two],
    output_type=FinalResult,
    hooks=create_hooks("Multiply Agent"),
)

start_agent = Agent(
    name="Start Agent",
    instructions="Generate a random number. If it's even, stop. If it's odd, hand off to the multiplier agent.",
    tools=[random_number],
    output_type=FinalResult,
    handoffs=[multiply_agent],
    hooks=create_hooks("Start Agent"),
)

async def main() -> None:
    user_input = input("Enter a max number: ")
    await Runner.run(
        start_agent,
        input=f"Generate a random number between 0 and {user_input}.",
    )
    print("Done!")

if __name__ == "__main__":
    asyncio.run(main())
