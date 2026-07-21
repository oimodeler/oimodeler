import json
import warnings
from typing import List

import pytest

from oimodeler.oimBasicFourierComponents import oimIRing, oimPt
from oimodeler.oimComponent import oimComponent, oimComponentFourier
from oimodeler.oimModel import oimModel
from oimodeler.oimParam import oimInterp, oimParam, oimParamInterpolator


@pytest.fixture
def intp_keys() -> List[str]:
    return [
        "name",
        "description",
        "unit",
        "param0",
        "interparams",
        "interpdescription",
        "class",
    ]


@pytest.fixture
def param_keys() -> List[str]:
    return [
        "name",
        "value",
        "min",
        "max",
        "description",
        "unit",
        "free",
        "error",
    ]


@pytest.fixture
def param(**kwargs) -> oimParam:
    return oimParam(base="default")


@pytest.fixture
def comp(**kwargs) -> oimComponent:
    return oimComponentFourier(elliptic=True, **kwargs)


@pytest.fixture
def model() -> oimModel:
    return oimModel([oimPt(), oimIRing(d=2)])


class TestOimParam:
    """Test serialisation of oimParam."""

    def test_roundtrip(self, param: oimParam) -> None:
        restored = oimParam.deserialize(param.serialize())

        # TODO: Implement oimParam.__eq__ and use it here
        assert restored.name == param.name
        assert restored.value == pytest.approx(param.value)
        assert restored.min == param.min
        assert restored.max == param.max
        assert restored.free == param.free
        assert restored.error == pytest.approx(param.error)
        assert restored.description == param.description
        assert restored.unit.to_string() == param.unit.to_string()

    def test_serialise_deepcopy(self, param: oimParam) -> None:
        ser = param.serialize()
        ser["value"] = 10
        assert param.value == 0

    def test_serialise_shallow_copy(self, param: oimParam) -> None:
        ser = param.serialize(skip_copy=True)
        ser["value"] = 10
        assert param.value == 10

    def test_idempotency(self, param: oimParam) -> None:
        ser1 = param.serialize()
        ser2 = param.serialize()

        assert json.dumps(ser1, sort_keys=True, default=str) == json.dumps(
            ser2, sort_keys=True, default=str
        )

        restored = oimParam.deserialize(ser1)
        ser2 = restored.serialize()

        assert set(ser1.keys()) == set(ser2.keys())


class TestOimParamLinker:

    def test_roundtrip(self): ...

    def test_serialise_deepcopy(self): ...

    def test_serialise_shallow_copy(self): ...


class TestOimParamNorm:

    def test_roundtrip(self): ...

    def test_serialise_deepcopy(self): ...

    def test_serialise_shallow_copy(self): ...


class TestOimParamInterpolator:
    """Test serialisation of oimParam."""

    def test_roundtrip(self): ...

    def test_serialise_deepcopy(self): ...

    def test_serialise_shallow_copy(self): ...


@pytest.mark.parametrize("skip_copy", [True, False])
def test_serialise_param_interpolator(
    skip_copy: bool, intp_keys: List[str]
) -> None:
    starWl_kwargs = {"T": 6500, "R": 3.46, "L": 10**1.35, "dist": 159.3}
    flux = oimInterp("starWl", **starWl_kwargs)
    flux = flux.type(oimParam(base="f"), **flux.kwargs)

    ser = flux.serialize(skip_copy=skip_copy)

    for p_key, value in starWl_kwargs.items():
        assert ser[p_key]["value"] == value
        if skip_copy:
            assert isinstance(flux.__dict__[p_key], dict)
        else:
            assert isinstance(flux.__dict__[p_key], oimParam)

    for intp_key in ["compute_radius", *intp_keys]:
        assert intp_key in ser

    with warnings.catch_warnings(record=True) as w:
        f = oimParamInterpolator.deserialize(ser)

        assert issubclass(w[-1].category, UserWarning)

    # TODO: Implement oimParamInterpolator.__eq__ and use it here


class TestOimComponent:
    """Test serialisation of oimComponent."""

    def test_roundtrip(self): ...

    def test_serialise_deepcopy(self): ...

    def test_serialise_shallow_copy(self): ...


@pytest.mark.parametrize("skip_copy", [True, False])
def test_serialise_component(
    skip_copy: bool,
    param_keys: List[str],
    comp: oimComponentFourier,
) -> None:
    """Test serialisation of oimComponent(Fourier)."""
    ser = comp.serialize(skip_copy=skip_copy)

    assert isinstance(ser, dict)
    assert all(x in ser for x in ["params", "other"])

    for key in comp.params.keys():
        assert (
            oimParam(base=key).serialize(skip_copy=skip_copy)
            == ser["params"][key]
        )

    prev_value = comp.x.value
    ser["params"]["x"]["value"] = 10
    if skip_copy:
        assert comp.x.value == 10
        comp.x.value = prev_value
    else:
        assert comp.x.value != 10

    ser["params"]["x"]["value"] = prev_value

    assert "elliptic" in ser["other"]
    assert ser["other"]["elliptic"]

    c = oimComponent.deserialize(ser)

    # TODO: Implement oimComponent.__eq__ and use it here
    for name, param in comp.params.items():
        assert all(
            getattr(param, k) == getattr(c.params[name], k) for k in param_keys
        )


class TestOimModel:
    """Test serialisation of oimModel."""

    def test_roundtrip(self): ...

    def test_serialise_deepcopy(self): ...

    def test_serialise_shallow_copy(self): ...


@pytest.mark.parametrize("skip_copy", [True, False])
def test_serialise_model(
    skip_copy: bool, model: oimModel, param_keys: List[str]
) -> None:
    """Test serialisation of oimComponent."""
    ser = model.serialize(skip_copy=skip_copy)

    assert isinstance(ser, dict)
    assert all(x in ser for x in ["components", "other"])
    assert "extParams" in ser["other"]

    for i, comp in enumerate(model.components):
        assert comp.__class__.__name__ == ser["components"][i][0]
        assert comp.serialize(skip_copy=skip_copy) == ser["components"][i][1]

    prev_value = model.components[0].x.value
    ser["components"][0][1]["params"]["x"]["value"] = 10
    if skip_copy:
        assert model.components[0].x.value == 10
        model.components[0].x.value = prev_value
    else:
        assert model.components[0].x.value != 10

    ser["components"][0][1]["params"]["x"]["value"] = prev_value

    m = oimModel.deserialize(ser)

    # TODO: Implement oimModel.__eq__ and use it here
    for comp, c in zip(model.components, m.components):
        assert comp.__class__.__name__ == c.__class__.__name__

        for pname, param in comp.params.items():
            assert all(
                getattr(param, k) == getattr(c.params[pname], k)
                for k in param_keys
            )
