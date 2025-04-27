import json

import pytest
from pydantic import BaseModel
from typing_extensions import TypedDict

from agents import Agent, AgentOutputSchema, ModelBehaviorError, Runner, UserError, _utils
from agents.agent_output import _WRAPPER_DICT_KEY


def test_plain_text_output():
    agent = Agent(name="test")
    output_schema = Runner._get_output_schema(agent)
    assert not output_schema, "Shouldn't have an output tool config without an output type"

    agent = Agent(name="test", output_type=str)
    assert not output_schema, "Shouldn't have an output tool config with str output type"


class Foo(BaseModel):
    bar: str


def test_structured_output_pydantic():
    agent = Agent(name="test", output_type=Foo)
    output_schema = Runner._get_output_schema(agent)
    assert output_schema, "Should have an output tool config with a structured output type"

    assert output_schema.output_type == Foo, "Should have the correct output type"
    assert not output_schema._is_wrapped, "Pydantic objects should not be wrapped"
    for key, value in Foo.model_json_schema().items():
        assert output_schema.json_schema()[key] == value

    json_str = Foo(bar="baz").model_dump_json()
    validated = output_schema.validate_json(json_str)
    assert validated == Foo(bar="baz")


class Bar(TypedDict):
    bar: str


def test_structured_output_typed_dict():
    agent = Agent(name="test", output_type=Bar)
    output_schema = Runner._get_output_schema(agent)
    assert output_schema, "Should have an output tool config with a structured output type"
    assert output_schema.output_type == Bar, "Should have the correct output type"
    assert not output_schema._is_wrapped, "TypedDicts should not be wrapped"

    json_str = json.dumps(Bar(bar="baz"))
    validated = output_schema.validate_json(json_str)
    assert validated == Bar(bar="baz")


def test_structured_output_list():
    agent = Agent(name="test", output_type=list[str])
    output_schema = Runner._get_output_schema(agent)
    assert output_schema, "Should have an output tool config with a structured output type"
    assert output_schema.output_type == list[str], "Should have the correct output type"
    assert output_schema._is_wrapped, "Lists should be wrapped"

    # This is testing implementation details, but it's useful  to make sure this doesn't break
    json_str = json.dumps({_WRAPPER_DICT_KEY: ["foo", "bar"]})
    validated = output_schema.validate_json(json_str)
    assert validated == ["foo", "bar"]


def test_bad_json_raises_error(mocker):
    agent = Agent(name="test", output_type=Foo)
    output_schema = Runner._get_output_schema(agent)
    assert output_schema, "Should have an output tool config with a structured output type"

    with pytest.raises(ModelBehaviorError):
        output_schema.validate_json("not valid json")

    agent = Agent(name="test", output_type=list[str])
    output_schema = Runner._get_output_schema(agent)
    assert output_schema, "Should have an output tool config with a structured output type"

    mock_validate_json = mocker.patch.object(_utils, "validate_json")
    mock_validate_json.return_value = ["foo"]

    with pytest.raises(ModelBehaviorError):
        output_schema.validate_json(json.dumps(["foo"]))

    mock_validate_json.return_value = {"value": "foo"}

    with pytest.raises(ModelBehaviorError):
        output_schema.validate_json(json.dumps(["foo"]))


def test_plain_text_obj_doesnt_produce_schema():
    output_wrapper = AgentOutputSchema(output_type=str)
    with pytest.raises(UserError):
        output_wrapper.json_schema()


def test_structured_output_is_strict():
    output_wrapper = AgentOutputSchema(output_type=Foo)
    assert output_wrapper.strict_json_schema
    for key, value in Foo.model_json_schema().items():
        assert output_wrapper.json_schema()[key] == value

    assert (
        "additionalProperties" in output_wrapper.json_schema()
        and not output_wrapper.json_schema()["additionalProperties"]
    )


def test_setting_strict_false_works():
    output_wrapper = AgentOutputSchema(output_type=Foo, strict_json_schema=False)
    assert not output_wrapper.strict_json_schema
    assert output_wrapper.json_schema() == Foo.model_json_schema()
