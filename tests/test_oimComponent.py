from typing import Any, Dict
import numpy as np
import pytest

from oimodeler.oimComponent import oimComponent
from oimodeler.oimParam import oimInterp, oimParam, oimParamGaussianWl


@pytest.fixture()
def params() -> Dict[str, Any]:
    return {
        "x": 5,
        "y": 10,
        "f": oimInterp("GaussWl", val0=2, value=4, x0=2.1656e-6, fwhm=1e-8),
    }


@pytest.fixture
def component(params):
    """Return an instance of oimComponent."""
    return oimComponent(**params)


def test_getFourierComponents(): ...


def test_oimComponent_eval(
    params: Dict[str, Any], component: oimComponent
) -> None:
    """Test the oimComponent's initialization."""
    for name in params.keys():
        assert name in component.params
        assert hasattr(component, name)

        comp_param = component.params[name]
        if name != "f":
            assert isinstance(comp_param, oimParam)
        else:
            assert isinstance(comp_param, oimParamGaussianWl)


def test_get_params(component: oimComponent) -> None:
    """Test oimComponent's get params function."""
    assert all(param in ["x", "y", "f", "dim"] for param in component.params)
    assert component.params["f"].value.free


def test_oimComponent_paramstr(component: oimComponent) -> None:
    """Test oimComponent's paramstr function."""
    assert all(param in component._paramstr() for param in ["x", "y", "f"])


def test_oimComponent_str(component: oimComponent) -> None:
    """Test oimComponent's string representation."""
    string, param_str = component.__str__(), component._paramstr()
    assert param_str in string and component.name in string


def test_oimComponent_repr(component: oimComponent) -> None:
    """Test oimComponent's console representation."""
    representation, param_str = component.__repr__(), component._paramstr()
    assert param_str in representation
    assert component.__class__.__name__ in representation


def test_oimComponent_directTranslate(component: oimComponent) -> None:
    """Test oimComponent's image space spatial translation."""
    assert component._directTranslate(10, 10, wl=None, t=None) == (5, 0)
    assert component._directTranslate(0, 0, wl=None, t=None) == (-5, -10)


# TODO: Do here the aspro tests
def test_oimComponent_ftTranslateFactor(component: oimComponent) -> None:
    """Test oimComponent's fourier space spatial translation."""
    ...


def test_oimComponent_getComplexCoherentFlux(component: oimComponent) -> None:
    """Test oimComponent's complex coherent flux calculation."""
    ucoord = np.linspace(-500, 500)
    complex_coherent_flux = component.getComplexCoherentFlux(
        u=ucoord, v=ucoord.copy()
    )
    assert complex_coherent_flux.shape == ucoord.shape
    assert np.array_equal(complex_coherent_flux, np.zeros(ucoord.shape))


def test_oimComponent_getImage(component: oimComponent) -> None:
    """Test oimComponent class."""
    image = component.getImage(dim=512, pixSize=0.1)
    assert image.size == 512**2
    assert np.array_equal(image, np.zeros((512, 512)))


def test_oimComponent_serialize(component: oimComponent) -> None:
    """Test oimComponent class' serialization."""
    ser = component.serialize()
    breakpoint()


def test_oimComponent_deserialize(): ...
