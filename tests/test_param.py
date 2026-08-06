import copy
import warnings

import astropy.units as u
import numpy as np
import pytest

from oimodeler.oimParam import oimInterp, oimParam, oimParamInterpolator

from .helpers import assert_intp_equal, assert_param_equal


class TestOimParam:
    """Test the oimParam class."""

    @pytest.fixture
    def param(self) -> oimParam:
        return oimParam(value=5 * u.deg, base="din")

    class TestInit:
        """Tests the __init__ method of oimParam."""

        def test_min_max(self) -> None:
            """Tests read-in of mini -> min and maxi -> max."""
            param = oimParam(mini=0, maxi=100)
            assert param.min == 0
            assert param.max == 100
            assert not hasattr(param, "mini")
            assert not hasattr(param, "maxi")

        def test_default(self) -> None:
            """Tests Initialisation without args/kwargs."""
            param = oimParam()
            assert param.name == ""
            assert param.value == 0
            assert param.min == -np.inf
            assert param.max == np.inf
            assert param.description == ""
            assert param.unit == u.one
            assert param.free
            assert param.error == 0

        def test_base(self) -> None:
            """Tests the base kwarg."""
            param = oimParam(base="din")
            assert param.name == "din"
            assert param.value == 0
            assert param.min == 0
            assert param.max == np.inf
            assert param.description == "Inner Diameter"
            assert param.unit == u.mas
            assert param.free
            assert param.error == 0

        def test_override(self) -> None:
            """Tests override of base kwarg."""
            param = oimParam(value=999, base="din")
            assert param.value == 999

        def test_quantity(self) -> None:
            """Initialisation with quantity as value."""
            param = oimParam(value=5 * u.mas)
            assert param.value == 5
            assert param.unit == u.mas

    class TestCall:
        """Tests the __call__ method of oimParam."""

        def test_return(self, param: oimParam) -> None:
            """Tests __call__'s return."""
            assert param() == 5

        def test_static(self, param: oimParam) -> None:
            """Tests wavelength- and time-independent __call__."""
            assert param(wl=3.5e-6, t=59221.14025707) == 5

    class TestQuantity:
        """Tests the __call__ method of oimParam."""

        def test_return(self, param: oimParam) -> None:
            quantity = param.quantity()
            assert isinstance(quantity, u.Quantity)
            assert quantity.value == 5
            assert quantity.unit == u.deg

        def test_alias(self, param: oimParam) -> None:
            assert param.qty() == param.quantity()

    class TestSerialisation:
        """Test serialisation of oimParam."""

        def test_roundtrip(self, param: oimParam, cycles: int = 5) -> None:
            current = param
            for _ in range(cycles):
                serialised = current.serialize()
                assert isinstance(serialised, dict)
                current = oimParam.deserialize(serialised)

            assert_param_equal(param, current)

        def test_deepcopy(self, param: oimParam) -> None:
            serialised = param.serialize()
            serialised["name"] = "modified"
            assert param.name != "modified"

        def test_shallow_copy(self, param: oimParam) -> None:
            serialised = param.serialize(skip_copy=True)
            serialised["name"] = "modified"
            assert param.name == "modified"


class TestOimParamLinker:

    class TestSerialisation:
        def test_roundtrip(self) -> None: ...

        def test_deepcopy(self) -> None: ...

        def test_shallow_copy(self) -> None: ...


class TestOimParamNorm:

    class TestSerialisation:
        def test_roundtrip(self) -> None: ...

        def test_deepcopy(self) -> None: ...

        def test_shallow_copy(self) -> None: ...


class TestOimParamInterpolator:
    """Test of oimParamInterpolator."""

    @pytest.fixture
    def star_kwargs(self) -> dict[str, float]:
        return {"T": 6500, "R": 3.46, "L": 10**1.35, "dist": 159.3}

    @pytest.fixture
    def star_intp(self, star_kwargs: dict[str, float]) -> oimParamInterpolator:
        intp = oimInterp("starWl", **star_kwargs)
        return intp.type(oimParam(base="f"), **intp.kwargs)

    @pytest.fixture
    def wl_intp(self) -> oimParam:
        intp = oimInterp(
            "wl", wl=[3e-6, 3.5e-6, 4e-6], values=[0.9, 0.85, 0.8]
        )
        return intp.type(oimParam(base="f"), **intp.kwargs)

    class TestSerialisation:
        """Test serialisation of oimParamInterpolator."""

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

                assert_intp_equal(intp, current)

        def test_deepcopy(
            self,
            star_kwargs: dict[str, float],
            star_intp: oimParamInterpolator,
        ) -> None:
            star_intp = copy.deepcopy(star_intp)
            star_intp.serialize()
            for p_key, value in star_kwargs.items():
                assert isinstance(star_intp.__dict__[p_key], oimParam)

        def test_shallow_copy(
            self,
            star_kwargs: dict[str, float],
            star_intp: oimParamInterpolator,
        ) -> None:
            star_intp = copy.deepcopy(star_intp)
            star_intp.serialize(skip_copy=True)
            for p_key, value in star_kwargs.items():
                assert isinstance(star_intp.__dict__[p_key], dict)
