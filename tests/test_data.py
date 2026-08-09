"""
Tests for the oimodeler.oimData module.
"""

from pathlib import Path

import oimodeler as oim

# TODO: Add more tests for oimodeler.oimData.oimData here


def test_load_gravity(real_data_dir: Path) -> None:
    """Test loading a GRAVITY FITS file."""
    data = oim.oimData(real_data_dir / "GRAVITY" / "gravity_real.fits")
    data.prepareData()
    wl = data.struct_wl[0][0]

    # NOTE: Number of spectral channels for FT camera in GRAVITY = 42
    assert len(wl) == 42
