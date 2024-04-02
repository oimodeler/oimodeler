import numpy as np
import pytest

from oimodeler import oimBasicFourierComponents as oimFComp


@pytest.fixture
def dim() -> int:
    """Return the dimension of the grid."""
    return 50


@pytest.fixture
def ucoord(dim: int) -> np.ndarray:
    """Create a meshgrid of uv coordinates."""
    ucoord = np.linspace(-500, 500, dim)
    return np.meshgrid(ucoord, ucoord)[0]


@pytest.fixture
def vcoord(dim: int) -> np.ndarray:
    """Create a meshgrid of uv coordinates."""
    ucoord = np.linspace(-500, 500, dim)
    return np.meshgrid(ucoord, ucoord)[1]


@pytest.fixture
def baselines(ucoord: np.ndarray, vcoord: np.ndarray) -> np.ndarray:
    """Create a baseline grid."""
    return np.hypot(ucoord, vcoord)


@pytest.fixture
def xcoord(dim: int) -> np.ndarray:
    """Create a meshgrid of xy coordinates."""
    xcoord = np.linspace(-100, 100, dim)
    return np.meshgrid(xcoord, xcoord)[0]


@pytest.fixture
def ycoord(dim: int) -> np.ndarray:
    """Create a meshgrid of xy coordinates."""
    xcoord = np.linspace(-100, 100, dim)
    return np.meshgrid(xcoord, xcoord)[1]


@pytest.fixture
def radius(xcoord: np.ndarray, ycoord: np.ndarray) -> np.ndarray:
    """Create a meshgrid of radii."""
    return np.hypot(xcoord, ycoord)


def test_oimPt_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimPt class."""
    pt = oimFComp.oimPt()
    assert pt._visFunction(ucoord, vcoord, baselines, None, None) == 1


def test_oimPt_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimPt class."""
    pt = oimFComp.oimPt()
    image = pt._imageFunction(xcoord, ycoord, None, None)
    assert image.shape == (dim, dim)
    assert image.max() == 1
    # assert image.argmax() == ((dim**2)//2-1)


# TODO: Finish this test
def test_oimBackground_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimBackground class."""
    bg = oimFComp.oimBackground()
    visibility = bg._visFunction(ucoord, vcoord, baselines, None, None)


def test_oimBackground_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimBackground class."""
    bg = oimFComp.oimBackground()
    image = bg._imageFunction(xcoord, ycoord, None, None)
    assert image.shape == (dim, dim)
    assert image.max() == 1
    assert np.all(image == 1)


# TODO: Finish this test
def test_oimUD_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimUD class."""
    ...
    # ud = oimFComp.oimUD()
    # visibility = ud._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimUd_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimUD class."""
    ...
    # ud = oimFComp.oimUD()
    # image = ud._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimEllipse_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimEllipse class."""
    ...
    # ellipse = oimFComp.oimEllipse()
    # visibility = ellipse._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimEllipse_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimEllipse class."""
    ...
    # ellipse = oimFComp.oimEllipse()
    # image = ellipse._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimGauss_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimGaussian class."""
    ...
    # gaussian = oimFComp.oimGauss()
    # visibility = gaussian._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimGauss_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimGaussian class."""
    ...
    # gaussian = oimFComp.oimGauss()
    # image = gaussian._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimEGauss_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimEGaussian class."""
    ...
    # egaussian = oimFComp.oimEGauss()
    # visibility = egaussian._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimIRing_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimIRing class."""
    ...
    # ir = oimFComp.oimIRing()
    # visibility = ir._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimIRing_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimIRing class."""
    ...
    # ir = oimFComp.oimIRing()
    # image = ir._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimEIring_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimEIring class."""
    ...
    # eir = oimFComp.oimEIring()
    # visibility = eir._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimEIring_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimEIring class."""
    ...
    # eir = oimFComp.oimEIring()
    # image = eir._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimRing_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimRing class."""
    ...
    # ring = oimFComp.oimRing()
    # visibility = ring._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimRing_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimRing class."""
    ...
    # ring = oimFComp.oimRing()
    # image = ring._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimRing2_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimRing2 class."""
    ...
    # ring2 = oimFComp.oimRing2()
    # visibility = ring2._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimRing2_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimRing2 class."""
    ...
    # ring2 = oimFComp.oimRing2()
    # image = ring2._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimERing_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimERing class."""
    ...
    # ering = oimFComp.oimERing()
    # visibility = ering._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimERing_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimERing class."""
    ...
    # ering = oimFComp.oimERing()
    # image = ering._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimERing2_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimERing2 class."""
    ...
    # ering2 = oimFComp.oimERing2()
    # visibility = ering2._visFunction(vis, radius, None, None)


# TODO: Finish this test
def test_oimERing2_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimERing2 class."""
    ...
    # ering2 = oimFComp.oimERing2()
    # image = ering2._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimESKIRing_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimESKIRing class."""
    ...
    # esk = oimFComp.oimESKIRing()
    # visibility = esk._visFunction(vis, radius, None, None)


# TODO: Finish this test
def test_oimESKIRing_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimESKIRing class."""
    ...
    # esk = oimFComp.oimESKIRing()
    # image = esk._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimESKRing_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimESKRing class."""
    ...
    # esk = oimFComp.oimESKRing()
    # visibility = esk._visFunction(vis, radius, None, None)


# TODO: Finish this test
def test_oimESKRing_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimESKRing class."""
    ...
    # esk = oimFComp.oimESKRing()
    # image = esk._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimLorentz_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimLorentz class."""
    ...
    # lorentz = oimFComp.oimLorentz()
    # visibility = lorentz._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimLorentz_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimLorentz class."""
    ...
    # lorentz = oimFComp.oimLorentz()
    # image = lorentz._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimELorentz_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimELorentz class."""
    ...
    # elorentz = oimFComp.oimELorentz()
    # visibility = elorentz._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimELorentz_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimELorentz class."""
    ...
    # elorentz = oimFComp.oimELorentz()
    # image = elorentz._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimLinearLDD_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimLinearLDD class."""
    ...
    # ldd = oimFComp.oimLinearLDD()
    # visibility = ldd._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimLinearLDD_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimLinearLDD class."""
    ...
    # ldd = oimFComp.oimLinearLDD()
    # image = ldd._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimQuadLDD_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimQuadLDD class."""
    ...
    # qld = oimFComp.oimQuadLDD()
    # visibility = qld._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimQuadLDD_imageFunction(dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimQuadLDD class."""
    ...
    # qld = oimFComp.oimQuadLDD()
    # image = qld._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
def test_oimPowerLawLDD_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimPowerLawLDD class."""
    ...
    # pldd = oimFComp.oimPowerLawLDD()
    # visibility = pldd._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimPowerLawLDD_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimPowerLawLDD class."""
    ...
    # pld = oimFComp.oimPowerLawLDD()
    # image = pld._imageFunction(xcoord, ycoord, None, None)


# TODO: Finish this test
def test_oimSqrtLDD_visFunction(
        ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimSqrtLDD class."""
    ...
    # sld = oimFComp.oimSqrtLDD()
    # visibility = sld._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
def test_oimSqrtLDD_imageFunction(
        dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
    """Test the image function of the oimSqrtLDD class."""
    ...
    # sld = oimFComp.oimSqrtLDD()
    # image = sld._imageFunction(xcoord, ycoord, None, None)
    # assert image.shape == (dim, dim)


# TODO: Finish this test
# def test_oimConvolutor_visFunction(
#         ucoord: np.ndarray, vcoord: np.ndarray, baselines: np.ndarray) -> None:
#     """Test the visFunction of the oimConvolutor class."""
#     conv = oimFComp.oimConvolutor()
#     visibility = conv._visFunction(ucoord, vcoord, baselines, None, None)


# TODO: Finish this test
# def test_oimConvolutor_imageFunction(
#         dim: int, xcoord: np.ndarray, ycoord: np.ndarray) -> None:
#     """Test the image function of the oimConvolutor class."""
#     conv = oimFComp.oimConvolutor()
#     image = conv._imageFunction(xcoord, ycoord, None, None)
#     assert image.shape == (dim, dim)
