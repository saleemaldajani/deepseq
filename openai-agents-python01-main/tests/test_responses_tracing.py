import pytest
from openai import AsyncOpenAI
from openai.types.responses import ResponseCompletedEvent

from agents import ModelSettings, ModelTracing, OpenAIResponsesModel, trace
from agents.tracing.span_data import ResponseSpanData
from tests import fake_model

from .testing_processor import fetch_ordered_spans


class DummyTracing:
    def is_disabled(self):
        return False


class DummyUsage:
    def __init__(self, input_tokens=1, output_tokens=1, total_tokens=2):
        self.input_tokens = input_tokens
        self.output_tokens = output_tokens
        self.total_tokens = total_tokens


class DummyResponse:
    def __init__(self):
        self.id = "dummy-id"
        self.output = []
        self.usage = DummyUsage()

    def __aiter__(self):
        yield ResponseCompletedEvent(
            type="response.completed",
            response=fake_model.get_response_obj(self.output),
        )


@pytest.mark.allow_call_model_methods
@pytest.mark.asyncio
async def test_get_response_creates_trace(monkeypatch):
    with trace(workflow_name="test"):
        # Create an instance of the model
        model = OpenAIResponsesModel(model="test-model", openai_client=AsyncOpenAI(api_key="test"))

        # Mock _fetch_response to return a dummy response with a known id
        async def dummy_fetch_response(
            system_instructions, input, model_settings, tools, output_schema, handoffs, stream
        ):
            return DummyResponse()

        monkeypatch.setattr(model, "_fetch_response", dummy_fetch_response)

        # Call get_response
        await model.get_response(
            "instr", "input", ModelSettings(), [], None, [], ModelTracing.ENABLED
        )

    spans = fetch_ordered_spans()
    assert len(spans) == 1

    assert isinstance(spans[0].span_data, ResponseSpanData)
    assert spans[0].span_data.response is not None
    assert spans[0].span_data.response.id == "dummy-id"


@pytest.mark.allow_call_model_methods
@pytest.mark.asyncio
async def test_non_data_tracing_doesnt_set_response_id(monkeypatch):
    with trace(workflow_name="test"):
        # Create an instance of the model
        model = OpenAIResponsesModel(model="test-model", openai_client=AsyncOpenAI(api_key="test"))

        # Mock _fetch_response to return a dummy response with a known id
        async def dummy_fetch_response(
            system_instructions, input, model_settings, tools, output_schema, handoffs, stream
        ):
            return DummyResponse()

        monkeypatch.setattr(model, "_fetch_response", dummy_fetch_response)

        # Call get_response
        await model.get_response(
            "instr", "input", ModelSettings(), [], None, [], ModelTracing.ENABLED_WITHOUT_DATA
        )

    spans = fetch_ordered_spans()
    assert len(spans) == 1
    assert spans[0].span_data.response is None


@pytest.mark.allow_call_model_methods
@pytest.mark.asyncio
async def test_disable_tracing_does_not_create_span(monkeypatch):
    with trace(workflow_name="test"):
        # Create an instance of the model
        model = OpenAIResponsesModel(model="test-model", openai_client=AsyncOpenAI(api_key="test"))

        # Mock _fetch_response to return a dummy response with a known id
        async def dummy_fetch_response(
            system_instructions, input, model_settings, tools, output_schema, handoffs, stream
        ):
            return DummyResponse()

        monkeypatch.setattr(model, "_fetch_response", dummy_fetch_response)

        # Call get_response
        await model.get_response(
            "instr", "input", ModelSettings(), [], None, [], ModelTracing.DISABLED
        )

    spans = fetch_ordered_spans()
    assert len(spans) == 0


@pytest.mark.allow_call_model_methods
@pytest.mark.asyncio
async def test_stream_response_creates_trace(monkeypatch):
    with trace(workflow_name="test"):
        # Create an instance of the model
        model = OpenAIResponsesModel(model="test-model", openai_client=AsyncOpenAI(api_key="test"))

        # Define a dummy fetch function that returns an async stream with a dummy response
        async def dummy_fetch_response(
            system_instructions, input, model_settings, tools, output_schema, handoffs, stream
        ):
            class DummyStream:
                async def __aiter__(self):
                    yield ResponseCompletedEvent(
                        type="response.completed",
                        response=fake_model.get_response_obj([], "dummy-id-123"),
                    )

            return DummyStream()

        monkeypatch.setattr(model, "_fetch_response", dummy_fetch_response)

        # Consume the stream to trigger processing of the final response
        async for _ in model.stream_response(
            "instr", "input", ModelSettings(), [], None, [], ModelTracing.ENABLED
        ):
            pass

    spans = fetch_ordered_spans()
    assert len(spans) == 1
    assert isinstance(spans[0].span_data, ResponseSpanData)
    assert spans[0].span_data.response is not None
    assert spans[0].span_data.response.id == "dummy-id-123"


@pytest.mark.allow_call_model_methods
@pytest.mark.asyncio
async def test_stream_non_data_tracing_doesnt_set_response_id(monkeypatch):
    with trace(workflow_name="test"):
        # Create an instance of the model
        model = OpenAIResponsesModel(model="test-model", openai_client=AsyncOpenAI(api_key="test"))

        # Define a dummy fetch function that returns an async stream with a dummy response
        async def dummy_fetch_response(
            system_instructions, input, model_settings, tools, output_schema, handoffs, stream
        ):
            class DummyStream:
                async def __aiter__(self):
                    yield ResponseCompletedEvent(
                        type="response.completed",
                        response=fake_model.get_response_obj([], "dummy-id-123"),
                    )

            return DummyStream()

        monkeypatch.setattr(model, "_fetch_response", dummy_fetch_response)

        # Consume the stream to trigger processing of the final response
        async for _ in model.stream_response(
            "instr", "input", ModelSettings(), [], None, [], ModelTracing.ENABLED_WITHOUT_DATA
        ):
            pass

    spans = fetch_ordered_spans()
    assert len(spans) == 1
    assert isinstance(spans[0].span_data, ResponseSpanData)
    assert spans[0].span_data.response is None


@pytest.mark.allow_call_model_methods
@pytest.mark.asyncio
async def test_stream_disabled_tracing_doesnt_create_span(monkeypatch):
    with trace(workflow_name="test"):
        # Create an instance of the model
        model = OpenAIResponsesModel(model="test-model", openai_client=AsyncOpenAI(api_key="test"))

        # Define a dummy fetch function that returns an async stream with a dummy response
        async def dummy_fetch_response(
            system_instructions, input, model_settings, tools, output_schema, handoffs, stream
        ):
            class DummyStream:
                async def __aiter__(self):
                    yield ResponseCompletedEvent(
                        type="response.completed",
                        response=fake_model.get_response_obj([], "dummy-id-123"),
                    )

            return DummyStream()

        monkeypatch.setattr(model, "_fetch_response", dummy_fetch_response)

        # Consume the stream to trigger processing of the final response
        async for _ in model.stream_response(
            "instr", "input", ModelSettings(), [], None, [], ModelTracing.DISABLED
        ):
            pass

    spans = fetch_ordered_spans()
    assert len(spans) == 0
