"""
Tests for the oimodeler.oimModel module.
"""

import copy
import json
from pathlib import Path

import pytest

from oimodeler.oimBasicFourierComponents import oimIRing, oimPt
from oimodeler.oimModel import oimModel

from .helpers import assert_model_equal


class TestSerialisation:
    """Test serialisation of oimModel."""

    @pytest.fixture(scope="module")
    def model(self) -> oimModel:
        """A simple model."""
        return oimModel([oimPt(), oimIRing(d=2)])

    def test_roundtrip(self, model: oimModel, cycles: int = 5) -> None:
        """Tests serialisation and deserialisation of the same object
        for multiple cycles."""
        current = model
        for _ in range(cycles):
            serialised = current.serialize()
            assert_model_equal(current, oimModel.deserialize(serialised))

    def test_deepcopy(self, model: oimModel) -> None:
        """Tests serialisation of a deepcopy."""
        serialised = model.serialize()
        serialised["components"][0][1]["params"]["x"]["value"] = 999
        assert model.components[0].x.value != 999

    def test_shallow_copy(self, model: oimModel) -> None:
        """Tests serialisation of a shallow copy."""
        tmp_model = copy.deepcopy(model)
        serialised = tmp_model.serialize(skip_copy=True)
        serialised["components"][0][1]["params"]["x"]["value"] = 999
        assert tmp_model.components[0].x.value == 999

    def test_json_serialisation(
        self, model: oimModel, tmp_path_factory: Path
    ) -> None:
        """Tests the JSON serialisation."""
        tmp_dir = tmp_path_factory.mktemp("data")
        json_path = tmp_dir / "serialised.json"

        serialised = model.serialize()
        with open(json_path, "w") as f:
            json.dump(serialised, f)

        with open(json_path, "r") as f:
            restored = json.load(f)

        assert restored == serialised
        assert_model_equal(model, oimModel.deserialize(restored))
