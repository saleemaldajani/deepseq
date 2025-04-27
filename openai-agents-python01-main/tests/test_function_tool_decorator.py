import asyncio
import json
from typing import Any

import pytest

from agents import function_tool
from agents.run_context import RunContextWrapper


class DummyContext:
    def __init__(self):
        self.data = "something"


def ctx_wrapper() -> RunContextWrapper[DummyContext]:
    return RunContextWrapper(DummyContext())


@function_tool
def sync_no_context_no_args() -> str:
    return "test_1"


@pytest.mark.asyncio
async def test_sync_no_context_no_args_invocation():
    tool = sync_no_context_no_args
    output = await tool.on_invoke_tool(ctx_wrapper(), "")
    assert output == "test_1"


@function_tool
def sync_no_context_with_args(a: int, b: int) -> int:
    return a + b


@pytest.mark.asyncio
async def test_sync_no_context_with_args_invocation():
    tool = sync_no_context_with_args
    input_data = {"a": 5, "b": 7}
    output = await tool.on_invoke_tool(ctx_wrapper(), json.dumps(input_data))
    assert int(output) == 12


@function_tool
def sync_with_context(ctx: RunContextWrapper[DummyContext], name: str) -> str:
    return f"{name}_{ctx.context.data}"


@pytest.mark.asyncio
async def test_sync_with_context_invocation():
    tool = sync_with_context
    input_data = {"name": "Alice"}
    output = await tool.on_invoke_tool(ctx_wrapper(), json.dumps(input_data))
    assert output == "Alice_something"


@function_tool
async def async_no_context(a: int, b: int) -> int:
    await asyncio.sleep(0)  # Just to illustrate async
    return a * b


@pytest.mark.asyncio
async def test_async_no_context_invocation():
    tool = async_no_context
    input_data = {"a": 3, "b": 4}
    output = await tool.on_invoke_tool(ctx_wrapper(), json.dumps(input_data))
    assert int(output) == 12


@function_tool
async def async_with_context(ctx: RunContextWrapper[DummyContext], prefix: str, num: int) -> str:
    await asyncio.sleep(0)
    return f"{prefix}-{num}-{ctx.context.data}"


@pytest.mark.asyncio
async def test_async_with_context_invocation():
    tool = async_with_context
    input_data = {"prefix": "Value", "num": 42}
    output = await tool.on_invoke_tool(ctx_wrapper(), json.dumps(input_data))
    assert output == "Value-42-something"


@function_tool(name_override="my_custom_tool", description_override="custom desc")
def sync_no_context_override() -> str:
    return "override_result"


@pytest.mark.asyncio
async def test_sync_no_context_override_invocation():
    tool = sync_no_context_override
    assert tool.name == "my_custom_tool"
    assert tool.description == "custom desc"
    output = await tool.on_invoke_tool(ctx_wrapper(), "")
    assert output == "override_result"


@function_tool(failure_error_function=None)
def will_fail_on_bad_json(x: int) -> int:
    return x * 2  # pragma: no cover


@pytest.mark.asyncio
async def test_error_on_invalid_json():
    tool = will_fail_on_bad_json
    # Passing an invalid JSON string
    with pytest.raises(Exception) as exc_info:
        await tool.on_invoke_tool(ctx_wrapper(), "{not valid json}")
    assert "Invalid JSON input for tool" in str(exc_info.value)


def sync_error_handler(ctx: RunContextWrapper[Any], error: Exception) -> str:
    return f"error_{error.__class__.__name__}"


@function_tool(failure_error_function=sync_error_handler)
def will_not_fail_on_bad_json(x: int) -> int:
    return x * 2  # pragma: no cover


@pytest.mark.asyncio
async def test_no_error_on_invalid_json():
    tool = will_not_fail_on_bad_json
    # Passing an invalid JSON string
    result = await tool.on_invoke_tool(ctx_wrapper(), "{not valid json}")
    assert result == "error_ModelBehaviorError"


def async_error_handler(ctx: RunContextWrapper[Any], error: Exception) -> str:
    return f"error_{error.__class__.__name__}"


@function_tool(failure_error_function=sync_error_handler)
def will_not_fail_on_bad_json_async(x: int) -> int:
    return x * 2  # pragma: no cover


@pytest.mark.asyncio
async def test_no_error_on_invalid_json_async():
    tool = will_not_fail_on_bad_json_async
    result = await tool.on_invoke_tool(ctx_wrapper(), "{not valid json}")
    assert result == "error_ModelBehaviorError"
