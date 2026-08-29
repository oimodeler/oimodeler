"""
Tests for the oimodeler.oimComponent module.
"""

import pytest

from oimodeler.oimComponent import oimComponent, oimComponentFourier
from oimodeler.oimParam import oimInterp, oimParam, oimParamInterpolator

from .helpers import assert_component_equal


@pytest.mark.skip(reason="Test not yet implemented.")
def test_getFourierComponents(): ...


class TestOimComponent:

    @pytest.fixture(scope="module")
    def component(self) -> oimComponent:
        """oimComponent with some values."""
        return oimComponent(x=5, y=10, f=0.5)

    @pytest.fixture(scope="module")
    def fourier_component(self) -> oimComponentFourier:
        """oimFourierComponent with oimParam adding kwarg."""
        return oimComponentFourier(elliptic=True, extincted=True, A_V=0.9)

    class TestEval:
        """Test the oimComponent's __init__ and _eval method."""

        def test_base(self, component: oimComponent) -> None:
            """Test simple __init__/_eval."""
            assert "x" in component.params
            assert "y" in component.params
            assert "f" in component.params
            assert component.params["f"].free

            assert isinstance(component.params["x"], oimParam)
            assert isinstance(component.params["y"], oimParam)
            assert component.params["x"].value == 5
            assert component.params["y"].value == 10
            assert component.params["f"].value == 0.5

            assert component.x == component.params["x"]
            assert component.y == component.params["y"]
            assert component.f == component.params["f"]

            component.x.value = 10
            assert component.params["x"].value == 10

            component.params["x"].value = 5
            assert component.x.value == 5

        def test_param_setting_kwargs(
            self, fourier_component: oimComponentFourier
        ) -> None:
            """Test __init__/_eval with oimFourierComponent for an oimParam setting kwarg."""
            assert fourier_component.elliptic
            assert "pa" in fourier_component.params
            assert "elong" in fourier_component.params

            assert fourier_component.extincted
            assert hasattr(fourier_component, "extargs")
            assert hasattr(fourier_component, "extlaw")
            assert "A_V" in fourier_component.params

        def test_interpolator(self) -> None:
            """Test __init__/_eval with oimInterp for an oimParam value."""
            keyframes, keyvalues = [3e-6, 3.5e-6, 4e-6], [0.9, 0.85, 0.8]
            intp = oimInterp("wl", wl=keyframes, values=keyvalues)
            assert isinstance(intp, oimInterp)
            assert not isinstance(intp, oimParamInterpolator)

            component = oimComponent(f=intp)
            assert component.f == component.params["f"]
            assert isinstance(component.f, oimParamInterpolator)

            for expected, restored in zip(keyframes, component.f.keyframes):
                assert expected == restored.value

            for expected, restored in zip(keyvalues, component.f.keyvalues):
                assert expected == restored.value

    def test_paramstr(self, component: oimComponent) -> None:
        """Test oimComponent's paramstr function."""
        expected = f"x={component.x.value:.2f} y={component.y.value:.2f} f={component.f.value:.2f}"
        assert component._paramstr() == expected

    # NOTE: Is this test needed?
    def test_str(self, component: oimComponent) -> None:
        """Test oimComponent's string representation."""
        expected = f"{component.name}: {component._paramstr()}"
        assert component.__str__() == expected

    # NOTE: Is this test needed?
    def test_repr(self, component: oimComponent) -> None:
        """Test oimComponent's console representation."""
        expected = f"{component.__class__.__name__} at {hex(id(component))}: {component._paramstr()}"
        assert component.__repr__() == expected

    def test_directTranslate(self, component: oimComponent) -> None:
        """Test oimComponent's image space spatial translation."""
        assert component._directTranslate(10, 10, wl=None, t=None) == (5, 0)
        assert component._directTranslate(0, 0, wl=None, t=None) == (-5, -10)

    @pytest.mark.skip(reason="Test not yet implemented.")
    def test_ftTranslateFactor(self, component: oimComponent) -> None:
        """Test oimComponent's fourier space spatial translation."""

    @pytest.mark.skip(reason="Test not yet implemented.")
    def test_getComplexCoherentFlux(self, component: oimComponent) -> None:
        """Test oimComponent's complex coherent flux calculation."""

    @pytest.mark.skip(reason="Test not yet implemented.")
    def test_getImage(self, component: oimComponent) -> None:
        """Test oimComponent class."""

    class TestSerialisation:
        """Test serialisation of oimComponent."""

        def test_roundtrip(
            self, fourier_component: oimComponentFourier, cycles: int = 5
        ) -> None:
            """Tests serialisation and deserialisation of the same object
            for multiple cycles."""
            current = fourier_component
            for _ in range(cycles):
                serialised = current.serialize()
                current = oimComponent.deserialize(serialised)
                assert_component_equal(fourier_component, current)

        def test_deepcopy(self, fourier_component: oimComponent) -> None:
            """Tests serialisation of a deepcopy."""
            serialised = fourier_component.serialize()

            serialised["params"]["x"]["value"] = 999
            assert fourier_component.x.value != 999

            # NOTE: Currently "other" contains only shallow copies
            serialised["other"]["elliptic"] = False
            assert fourier_component.elliptic

        def test_shallow_copy(self, fourier_component: oimComponent) -> None:
            """Tests serialisation of a shallow copy."""
            serialised = fourier_component.serialize(skip_copy=True)

            serialised["params"]["x"]["value"] = 999
            assert fourier_component.x.value == 999

            # NOTE: Currently "other" contains only shallow copies
            serialised["other"]["elliptic"] = False
            assert fourier_component.elliptic
