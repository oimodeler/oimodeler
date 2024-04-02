import numpy as np
import pytest

from oimodeler.oimComponent import oimComponent
from oimodeler.oimParam import oimParam, oimInterp, oimParamGaussianWl


@pytest.fixture
def component():
    """Return an instance of oimComponent."""
    f = oimInterp("GaussWl", val0=2, value=4, x0=2.1656e-6, fwhm=1e-8)
    return oimComponent(x=5, y=10, f=f)


def test_getFourierComponents():
    ...


# NOTE: The shift might be a bit counterintuitive
def test_oimComponent_eval(component: oimComponent) -> None:
    """Test the oimComponent's initialization."""
    params = component.params
    assert isinstance(params["x"], oimParam)
    assert isinstance(params["y"], oimParam)
    assert isinstance(params["f"], oimParamGaussianWl)
    assert params["x"].value == 5
    assert params["y"].value == 10
    assert params["f"].value.value == 4


def test_get_params(component: oimComponent) -> None:
    """Test oimComponent's get params function."""
    assert all(param in ["x", "y", "dim", "f"] for param in component.params)
    assert component.params["f"].value.free


def test_oimComponent_paramstr(component: oimComponent) -> None:
    """Test oimComponent's paramstr function."""
    assert all(param in component._paramstr() for param in ["x", "y", "dim", "f"])


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
            u=ucoord, v=ucoord.copy())
    assert complex_coherent_flux.shape == ucoord.shape
    assert np.array_equal(complex_coherent_flux, np.zeros(ucoord.shape))


def test_oimComponent_getImage(component: oimComponent) -> None:
    """Test oimComponent class."""
    image = component.getImage(dim=512, pixSize=0.1)
    assert image.size == 512**2
    assert np.array_equal(image, np.zeros((512, 512)))
