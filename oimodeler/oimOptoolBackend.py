from pathlib import Path
from typing import Optional

import astropy.units as u
import optool as op
from scipy import interpolate


# TODO: Check if optool executable is in bin and if not give error?
def load_opacity(file: Path,
                 scat: Optional[bool] = False,
                 wavelength_solution: Optional[u.um] = None) -> u.cm**2/u.g:
    """Reads opacity from a (.dat)-file created with optool.

    If wavelength solution is given, the wavelength range from the optool-file is
    interpolated to the provided wavelength solution

    Parameters
    ----------
    file: Path
        Input (.dat)-file, created with optool, containing the opacity distribution
    scat: bool, optional
        If toggled will assume that the input (.dat)-file contains a scattering matrix
    wavelength_solution: u.um, optional
        A wavelength solution to interpolate the opacity distribution to

    Returns
    -------
    kappa: u.cm**2/u.g
        The opacity distribution. If a wavelength solution has been provided, then it is
        interpolated to it
    """
    _, wl, kappa, *_ = op.readoutputfile(str(file), scat=scat)
    if wavelength_solution is not None:
        if isinstance(wavelength_solution, u.Quantity):
            wavelength_solution = wavelength_solution.value
        kappa = interpolate.interp1d(wl, kappa)(wavelength_solution)
    return kappa
