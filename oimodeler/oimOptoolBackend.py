import pkg_resources
from distutils.spawn import find_executable
from pathlib import Path
from typing import Optional, Any, List, Dict

import astropy.units as u
import numpy as np
import optool as op
from scipy import interpolate


def interpolate_opacity(wl: np.ndarray, opacity: np.ndarray,
                        wavelength_solution: Optional[np.ndarray | u.Quantity] = None
                        ) -> np.ndarray:
    """Interpolates the opacity to a given wavelength solution if provided, otherwise
    returns the opacity as is

    Parameters
    ----------
    wl: np.ndarray
        The wavelength from a optool calculated (.dat)-file
    opacity: np.ndarray
        The opacity from a optool calculated (.dat)-file
    wavelength_solution: np.ndarray | u.um, optional
        A wavelength solution to interpolate the opacity distribution to

    Returns
    -------
    opactiy: np.ndarray
        Either the opacity or interpolated opacity if a wavelength solution is provided
    """
    if wavelength_solution is not None:
        if isinstance(wavelength_solution, u.Quantity):
            wavelength_solution = wavelength_solution.value
    return interpolate.interp1d(wl, opacity)(wavelength_solution)


def load_opacity(file: Path = None,
                 scat: Optional[bool] = False,
                 wavelength_solution: Optional[np.ndarray | u.Quantity] = None) -> np.ndarray:
    """Reads the absorption opacity from a (.dat)-file created with optool.

    If a wavelength solution is given, the wavelength range from the optool-file is
    interpolated to the provided wavelength solution

    Parameters
    ----------
    file: Path, optional
        Input (.dat)-file, created with optool, containing the opacity distribution
    scat: bool, optional
        If toggled will assume that the input (.dat)-file contains a scattering matrix
    wavelength_solution: np.ndarray | u.um, optional
        A wavelength solution to interpolate the opacity distribution to

    Returns
    -------
    opacity: np.ndarray
        The opacity distribution in [cm**2/u.g]. If a wavelength solution has been
        provided, then it is interpolated to it
    """
    _, wl, opacity, *_ = op.readoutputfile(str(file), scat=scat)
    return interpolate_opacity(wl, opacity, wavelength_solution)


def generate_switch_string(value: Any, switch: str) -> str:
    """Generates a command string to be passed to optool

    Parameters
    ----------
    value: any
        The input value
    switch: str
        The switch to be used in the command string

    Returns
    -------
    str
        Part of the command line argument passed to optool
    """
    if value is None:
        return ""
    elif isinstance(value, Path):
        return f"-{switch} {str(value)} "
    elif isinstance(value, str) or isinstance(value, float) or isinstance(value, int):
        return f"-{switch} {value} "
    elif isinstance(value, List):
        return f"-{switch} {' '.join([str(i) for i in value])} "
    elif isinstance(value, Dict):
        return f"-{switch} {' '.join([f'{k} {v}' for k, v in value.items()])} "
    else:
        raise IOError(f"Input type '{type(value)}' for '{switch}' is not supported!")


class OptoolBackend(op.particle):
    """A wrapper for the dust distribution calculation with optool.

    Once calculated the (.dat)-files (if no other output path provided) will be stored
    internally in /data/cache/optool/ to increase the speed of future calculations

    Parameters
    ----------
    grains: Path | str | Dict[str, int]
        Path to a (.lnk)-file or string (name) of a material or dictionary containing
        materials to include in the grain. If a dictionary is passed then it needs to
        contain the name of the material as key and its mass fraction as value (e.g.,
        {"pyr-mg70": 0.87, "c": 0.13}).
        Mass fractions do not have to add up to 1, they will be renormalized.
        Up to 20 materials can be specified to build up a grain
    grain_mantels: Path | str | Dict[str, int], optional
        The same as the grains parameter, but the material will be placed into the grain
        mantle. Multiple mantle materials will be mixed using the Bruggeman rule, and
        than that mix will be added to the core using the Maxwell-Garnett rule
    porosity: int | List[int], optional
        The volume fraction of vacuum, a number smaller than 1. The default is 0.
        A single value will apply to both core and mantle, but if a list with two values
        is provided the second value will be specific for the mantle (and may be 0)
    computational_method: str, optional
    f_max: float, optional
    monomer_radius: float, optional
    scat: bool, optional
        Include the scattering matrix in the output
    nang: int, optional
        The number of evenly spaced angular grid points that cover a range between 0 and
        180 degrees. The default for nang is 180
    nsub: int, optional
        Divide the computation up into parts to produce a file for each grain size.
        Each size will be an average over a range of nsub (default 5) grains around the
        real size
    ndeg: int, optional
        Cap the first degrees of the forward scattering peak. Default 2
    cache_dir: Path, optional
        The directory to store the (.dat)-files in. Default is '/data/cache/optool/'

    Returns
    -------
    particle: op.particle
        The optool particle class containing the specified materials and attributes

    See Also
    --------
    op.particle.plot()
    op.particle.computemean(tmin=10, tmax=1500, ntemp=100)
    op.particle.scatnorm(norm='')
    op.particle.sizedist(N_of_a)
    """

    def __init__(self, grains: Path | str | Dict[str, int] = None,
                 grain_mantels: Optional[Path | str | Dict[str, int]] = None,
                 porosity: Optional[int | List[int]] = 0,
                 computational_method: Optional[str] = "dhs",
                 f_max: Optional[float] = 0.8,
                 monomer_radius: Optional[float] = 0.1,
                 scat: Optional[bool] = False,
                 nang: Optional[int] = 180,
                 nsub: Optional[int] = 5,
                 ndeg: Optional[int] = 2,
                 cache_dir: Optional[Path] = None):
        if not find_executable("optool"):
            raise RuntimeError("The 'optool' executable has not been found! "
                               "For installation instructions see:\n"
                               "https://github.com/cdominik/optool/blob/master/UserGuide.org")

        if nang > 180 or nang < 0:
            raise IOError("The number of angular grid points has to be between 0 and 180!")

        cmd = "optool "
        cmd += generate_switch_string(grains, "c")
        cmd += generate_switch_string(grain_mantels, "m")
        cmd += generate_switch_string(porosity, "p")
        cmd = cmd.strip()

        if cache_dir is None:
            cache_dir = Path(pkg_resources.resource_filename("oimodeler", "data/cache/"))

        super().__init__(cmd, str(cache_dir))

    def print_materials(self):
        ...


# particle = get_particle("pyr")
test = OptoolBackend("pyr")
breakpoint()
