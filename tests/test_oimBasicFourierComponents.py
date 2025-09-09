import numpy as np
import pytest
from numpy.typing import ArrayLike

from oimodeler import oimBasicFourierComponents as oimFComp


@pytest.fixture
def uvcoord() -> ArrayLike:
    """Create a meshgrid of uv coordinates."""
    ucoord = np.linspace(0, 100, 25)
    return np.meshgrid(ucoord, ucoord)


@pytest.fixture
def baselines(uvcoord: ArrayLike) -> ArrayLike:
    """Create a baseline grid."""
    return np.hypot(*uvcoord)


def test_oimPt_visFunction(uvcoord: ArrayLike, baselines: ArrayLike) -> None:
    """Test the visFunction of the oimPt class."""
    pt = oimFComp.oimPt()
    assert pt._visFunction(*uvcoord, baselines, None, None) == 1


# TODO: Finish this test
def test_oimBackground_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimBackground class."""
    bg = oimFComp.oimBackground()


# TODO: Finish this test
def test_oimUD_visFunction(uvcoord: ArrayLike, baselines: ArrayLike) -> None:
    """Test the visFunction of the oimUD class."""
    ud = oimFComp.oimUD()


# TODO: Finish this test
def test_oimEllipse_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimEllipse class."""
    ellipse = oimFComp.oimEllipse()


# TODO: Finish this test
def test_oimGauss_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimGaussian class."""
    gaussian = oimFComp.oimGauss()


# TODO: Finish this test
def test_oimEGauss_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimEGaussian class."""
    egaussian = oimFComp.oimEGauss()


# TODO: Finish this test
def test_oimIRing_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimIRing class."""
    ir = oimFComp.oimIRing()


# TODO: Finish this test
def test_oimEIring_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimEIring class."""
    eir = oimFComp.oimEIRing()


# TODO: Finish this test
def test_oimRing_visFunction(uvcoord: ArrayLike, baselines: ArrayLike) -> None:
    """Test the visFunction of the oimRing class."""
    ring = oimFComp.oimRing()


# TODO: Finish this test
def test_oimRing2_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimRing2 class."""
    ring2 = oimFComp.oimRing2()


# TODO: Finish this test
def test_oimERing_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimERing class."""
    ering = oimFComp.oimERing()


# TODO: Finish this test
def test_oimERing2_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimERing2 class."""
    ...
    # ering2 = oimFComp.oimERing2()
    # visibility = ering2._visFunction(vis, radius, None, None)


# TODO: Finish this test
def test_oimESKIRing_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimESKIRing class."""
    esk = oimFComp.oimESKIRing()


# TODO: Finish this test
def test_oimESKRing_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimESKRing class."""
    esk = oimFComp.oimESKRing()


# TODO: Finish this test
def test_oimLorentz_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimLorentz class."""
    lorentz = oimFComp.oimLorentz()


# TODO: Finish this test
def test_oimELorentz_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimELorentz class."""
    elorentz = oimFComp.oimELorentz()


# TODO: Finish this test
def test_oimLinearLDD_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimLinearLDD class."""
    ldd = oimFComp.oimLinearLDD()


# TODO: Finish this test
def test_oimQuadLDD_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimQuadLDD class."""
    qld = oimFComp.oimQuadLDD()


# TODO: Finish this test
def test_oimPowerLawLDD_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimPowerLawLDD class."""
    pldd = oimFComp.oimPowerLawLDD()


# TODO: Finish this test
def test_oimSqrtLDD_visFunction(
    uvcoord: ArrayLike, baselines: ArrayLike
) -> None:
    """Test the visFunction of the oimSqrtLDD class."""
    sld = oimFComp.oimSqrtLDD()


@pytest.mark.parametrize(
    "pa1, elong1, pa2, elong2",
    [(0, 1, 0, 1), (33, 2, 45, 1.5)],
)
def test_oimConvolutor_visFunction(
    uvcoord: ArrayLike,
    pa1: int,
    elong1: float,
    pa2: int,
    elong2: float,
) -> None:
    """Test the visFunction of the oimConvolutor class."""
    spfu, spfv = np.array(uvcoord) / 3.5e-6
    ring = oimFComp.oimEIRing(pa=pa1, elong=elong1, d=4)
    gauss = oimFComp.oimEGauss(pa=pa2, elong=elong2, fwhm=2)
    ring_vis = ring.getComplexCoherentFlux(spfu, spfv)
    gauss_vis = gauss.getComplexCoherentFlux(spfu, spfv)
    manual_conv_vis = ring_vis * gauss_vis

    conv = oimFComp.oimConvolutor(ring, gauss)
    conv_vis = conv.getComplexCoherentFlux(spfu, spfv)
    assert np.array_equal(conv_vis, manual_conv_vis)
