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


def test_oimComponent_eval(component: oimComponent) -> None:
    """Test the oimComponent's initialization."""
    assert isinstance(component.x, oimParam)
    assert isinstance(component.y, oimParam)
    assert isinstance(component.f, oimParamGaussianWl)
    assert component.x.value == 5
    assert component.y.value == 10
    assert component.f.value.value == 4


def test_oimComponent_paramstr():
    """Test oimComponent class."""
    ...


def test_oimComponent_str():
    """Test oimComponent class."""
    ...


def test_oimComponent_repr():
    """Test oimComponent class."""
    ...


def test_oimComponent_getComplexCoherentFlux():
    """Test oimComponent class."""
    ...


def test_oimComponent_getImage():
    """Test oimComponent class."""
    ...


def test_oimComponent_directTranslate():
    """Test oimComponent class."""
    ...


def test_oimComponent_ftTranslateFactor():
    """Test oimComponent class."""
    ...
