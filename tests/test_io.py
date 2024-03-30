from pathlib import Path

import oimodeler as oim


def test_load_gravity(global_datadir: Path) -> None:
    """Test loading a GRAVITY FITS file."""
    fits_file = str(global_datadir / "test_gravity_real.fits")
    data = oim.oimData(str(fits_file))
    data.prepareData()
    wl = data.struct_wl[0][0]
    nwl = len(wl)
    # Number of spectral channels for FT camera in GRAVITY = 42
    assert nwl == 42
