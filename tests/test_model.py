"""
Tests for the oimodeler.oimModel module.
"""

import pytest

from oimodeler.oimBasicFourierComponents import oimIRing, oimPt
from oimodeler.oimModel import oimModel

from .helpers import assert_model_equal


class TestSerialisation:
    """Test serialisation of oimModel."""

    @pytest.fixture
    def model(self) -> oimModel:
        return oimModel([oimPt(), oimIRing(d=2)])

    def test_roundtrip(self, model: oimModel, cycles: int = 5) -> None:
        current = model
        for _ in range(cycles):
            serialised = current.serialize()
            assert_model_equal(current, oimModel.deserialize(serialised))

    def test_deepcopy(self, model: oimModel) -> None:
        serialised = model.serialize()
        serialised["components"][0][1]["params"]["x"]["value"] = 999
        assert model.components[0].x.value != 999

    def test_shallow_copy(self, model: oimModel) -> None:
        serialised = model.serialize(skip_copy=True)
        serialised["components"][0][1]["params"]["x"]["value"] = 999
        assert model.components[0].x.value == 999
