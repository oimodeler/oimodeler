import copy
import warnings
from typing import Dict

import pytest

from oimodeler.oimBasicFourierComponents import oimIRing, oimPt
from oimodeler.oimComponent import oimComponent, oimComponentFourier
from oimodeler.oimModel import oimModel
from oimodeler.oimParam import oimInterp, oimParam, oimParamInterpolator

# TODO: Make more tests for the serialised dictinaries itself?


# TODO: Should __eq__ for oimParam, oimComponent, etc. be implemented?
def param_equal(original: oimParam, restored: oimParam) -> bool:
    """Compares an original to a restored oimParam."""
    result = []
    for key, value_org in original.__dict__.items():
        value_res = getattr(restored, key)

        if key in ["value", "error"]:
            value_res = pytest.approx(value_res)
        elif key == "unit":
            value_org, value_res = (
                value_org.to_string(),
                value_res.to_string(),
            )

        result.append(value_org == value_res)

    return all(result)


def intp_equal(
    original: oimParamInterpolator, restored: oimParamInterpolator
) -> bool:
    """Compares an original to a restored oimParamInterpolator."""
    original_dict = copy.deepcopy(original.__dict__)
    results = []
    for key, value in original_dict.items():
        value_res = getattr(restored, key)
        if isinstance(value, oimParam):
            res = param_equal(value, value_res)
        else:
            res = value == value_res

        results.append(res)

    return all(results)


def component_equal(original: oimComponent, restored: oimComponent) -> bool:
    """Compares an original to a restored oimComponent."""
    original_dict = copy.deepcopy(original.__dict__)
    del original_dict["params"]

    results = []
    for key, value in original_dict.items():
        # HACK: Combating python adding the class from which it inherited
        # in front of private attributes
        key = key.replace("_oimComponent_", "")

        value_res = getattr(restored, key)
        if isinstance(value, oimParam):
            res = param_equal(value, value_res)
        else:
            res = value == value_res

        results.append(res)

    for name, param in original.params.items():
        results.append(
            param_equal(param, restored.params.get(name, oimParam()))
        )

    return all(results)


def model_equal(original: oimModel, restored: oimModel) -> bool:
    """Compares an original to a restored oimModel."""
    results = []
    for comp_org, comp_res in zip(original.components, restored.components):
        results.append(component_equal(comp_org, comp_res))

    for param_org, param_res in zip(original.extParams, restored.extParams):
        results.append(param_org == param_res)

    return all(results)


# TODO: Make more tests for oimParam in general
class TestOimParam:
    """Test serialisation of oimParam."""

    @pytest.fixture
    def param(self) -> oimParam:
        return oimParam(base="default")

    def test_roundtrip(self, param: oimParam, cycles: int = 5) -> None:
        current = param
        for _ in range(cycles):
            serialised = current.serialize()
            assert isinstance(serialised, dict)

            current = oimParam.deserialize(serialised)

        assert param_equal(param, current)

    def test_deepcopy(self, param: oimParam) -> None:
        serialised = param.serialize()
        serialised["name"] = "modified"
        assert param.name != "modified"

    def test_shallow_copy(self, param: oimParam) -> None:
        serialised = param.serialize(skip_copy=True)
        serialised["name"] = "modified"
        assert param.name == "modified"


class TestOimParamLinker:

    def test_roundtrip(self) -> None: ...

    def test_deepcopy(self) -> None: ...

    def test_shallow_copy(self) -> None: ...


class TestOimParamNorm:

    def test_roundtrip(self) -> None: ...

    def test_deepcopy(self) -> None: ...

    def test_shallow_copy(self) -> None: ...


class TestOimParamInterpolator:
    """Test serialisation of oimParam."""

    @pytest.fixture
    def star_kwargs(self) -> Dict[str, float]:
        return {"T": 6500, "R": 3.46, "L": 10**1.35, "dist": 159.3}

    @pytest.fixture
    def star_intp(self, star_kwargs: Dict[str, float]) -> oimParamInterpolator:
        intp = oimInterp("starWl", **star_kwargs)
        return intp.type(oimParam(base="f"), **intp.kwargs)

    @pytest.fixture
    def wl_intp(self) -> oimParam:
        intp = oimInterp("wl", wls=[], values=[])
        return intp.type(oimParam(base="f"), **intp.kwargs)

    def test_roundtrip(
        self,
        wl_intp: oimParamInterpolator,
        star_intp: oimParamInterpolator,
        cycles: int = 5,
    ) -> None:
        for intp in [wl_intp, star_intp]:
            current = intp
            for _ in range(cycles):
                with warnings.catch_warnings(record=True):
                    serialised = current.serialize()
                    current = oimParamInterpolator.deserialize(serialised)

            assert intp_equal(intp, current)

    def test_deepcopy(
        self, star_kwargs: Dict[str, float], star_intp: oimParamInterpolator
    ) -> None:
        star_intp = copy.deepcopy(star_intp)
        star_intp.serialize()
        for p_key, value in star_kwargs.items():
            assert isinstance(star_intp.__dict__[p_key], oimParam)

    def test_shallow_copy(
        self, star_kwargs: Dict[str, float], star_intp: oimParamInterpolator
    ) -> None:
        star_intp = copy.deepcopy(star_intp)
        star_intp.serialize(skip_copy=True)
        for p_key, value in star_kwargs.items():
            assert isinstance(star_intp.__dict__[p_key], dict)


class TestOimComponent:
    """Test serialisation of oimComponent."""

    @pytest.fixture
    def component(self) -> oimComponent:
        return oimComponentFourier(elliptic=True)

    def test_roundtrip(self, component: oimComponent, cycles: int = 5) -> None:
        current = component
        for _ in range(cycles):
            serialised = current.serialize()
            assert component_equal(
                current, oimComponent.deserialize(serialised)
            )

    def test_deepcopy(self, component: oimComponent) -> None:
        serialised = component.serialize()

        serialised["params"]["x"]["value"] = 999
        assert component.x.value != 999

        # NOTE: Currently "other" contains only shallow copies
        serialised["other"]["elliptic"] = False
        assert component.elliptic

    def test_shallow_copy(self, component: oimComponent) -> None:
        serialised = component.serialize(skip_copy=True)

        serialised["params"]["x"]["value"] = 999
        assert component.x.value == 999

        # NOTE: Currently "other" contains only shallow copies
        serialised["other"]["elliptic"] = False
        assert component.elliptic


class TestOimModel:
    """Test serialisation of oimModel."""

    @pytest.fixture
    def model(self) -> oimModel:
        return oimModel([oimPt(), oimIRing(d=2)])

    def test_roundtrip(self, model: oimModel, cycles: int = 5) -> None:
        current = model
        for _ in range(cycles):
            serialised = current.serialize()
            assert model_equal(current, oimModel.deserialize(serialised))

    def test_deepcopy(self, model: oimModel) -> None:
        serialised = model.serialize()
        serialised["components"][0][1]["params"]["x"]["value"] = 999
        assert model.components[0].x.value != 999

    def test_shallow_copy(self, model: oimModel) -> None:
        serialised = model.serialize(skip_copy=True)
        serialised["components"][0][1]["params"]["x"]["value"] = 999
        assert model.components[0].x.value == 999
