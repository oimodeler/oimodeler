"""
Tests for the oimodeler.oimBasicFourierComponents module.
"""

import numpy as np
import pytest

from oimodeler import oimBasicFourierComponents as oimFComp


@pytest.fixture
def uvcoord() -> np.ndarray:
    """Create a meshgrid of uv coordinates."""
    ucoord = np.linspace(0, 100, 25)
    return np.meshgrid(ucoord, ucoord)


@pytest.fixture
def baselines(uvcoord: np.ndarray) -> np.ndarray:
    """Create a baseline grid."""
    return np.hypot(*uvcoord)


def test_oimPt_visFunction(uvcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimPt class."""
    assert oimFComp.oimPt()._visFunction(*uvcoord, baselines, None, None) == 1


@pytest.mark.skip(reason="Test not implemented.")
def test_oimBackground_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimBackground class."""
    oimFComp.oimBackground()


@pytest.mark.skip(reason="Test not implemented.")
def test_oimUD_visFunction(uvcoord: np.ndarray, baselines: np.ndarray) -> None:
    """Test the visFunction of the oimUD class."""
    oimFComp.oimUD()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimEllipse_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimEllipse class."""
    oimFComp.oimEllipse()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimGauss_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimGaussian class."""
    oimFComp.oimGauss()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimEGauss_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimEGaussian class."""
    oimFComp.oimEGauss()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimIRing_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimIRing class."""
    oimFComp.oimIRing()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimEIring_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimEIring class."""
    oimFComp.oimEIRing()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimRing_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimRing class."""
    oimFComp.oimRing()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimRing2_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimRing2 class."""
    oimFComp.oimRing2()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimERing_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimERing class."""
    oimFComp.oimERing()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimERing2_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimERing2 class."""
    oimFComp.oimERing2()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimESKIRing_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimESKIRing class."""
    oimFComp.oimESKIRing()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimESKRing_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimESKRing class."""
    oimFComp.oimESKRing()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimLorentz_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimLorentz class."""
    oimFComp.oimLorentz()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimELorentz_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimELorentz class."""
    oimFComp.oimELorentz()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimLinearLDD_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimLinearLDD class."""
    oimFComp.oimLinearLDD()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimQuadLDD_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimQuadLDD class."""
    oimFComp.oimQuadLDD()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimPowerLawLDD_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimPowerLawLDD class."""
    oimFComp.oimPowerLawLDD()


@pytest.mark.skip(reason="Test not yet implemented.")
def test_oimSqrtLDD_visFunction(
    uvcoord: np.ndarray, baselines: np.ndarray
) -> None:
    """Test the visFunction of the oimSqrtLDD class."""
    oimFComp.oimSqrtLDD()


@pytest.mark.parametrize("pa1", (0, 33))
@pytest.mark.parametrize("elong1", (1, 2))
@pytest.mark.parametrize("pa2", (0, 45))
@pytest.mark.parametrize("elong2", (1, 1.5))
def test_oimConvolutor_visFunction(
    uvcoord: np.ndarray,
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
