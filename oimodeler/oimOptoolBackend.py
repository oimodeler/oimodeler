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
    """Interpolates the opacity to a given wavelength solution if provided,
    otherwise returns the opacity as is

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
        Either the opacity or interpolated opacity if a wavelength solution is
        provided
    """
    if wavelength_solution is not None:
        if isinstance(wavelength_solution, u.Quantity):
            wavelength_solution = wavelength_solution.value
    return interpolate.interp1d(wl, opacity)(wavelength_solution)


def load_optool_opacity(file: Path = None,
                        scat: Optional[bool] = False,
                        wavelength_solution: Optional[np.ndarray | u.Quantity] = None
                        ) -> np.ndarray:
    """Reads the absorption opacity from a (.dat)-file created with optool.

    If a wavelength solution is given, the wavelength range from the optool-file
    is interpolated to the provided wavelength solution

    Parameters
    ----------
    file: Path, optional
        Input (.dat)-file, created with optool, containing the opacity
        distribution
    scat: bool, optional
        If toggled will assume that the input (.dat)-file contains a scattering
        matrix
    wavelength_solution: np.ndarray | u.um, optional
        A wavelength solution to interpolate the opacity distribution to

    Returns
    -------
    opacity: np.ndarray
        The opacity distribution in [cm**2/u.g]. If a wavelength solution has
        been provided, then it is interpolated to it
    """
    _, wl, opacity, *_ = op.readoutputfile(str(file), scat=scat)
    return interpolate_opacity(wl, opacity, wavelength_solution)


def generate_switch_string(value: Any, switch: str,
                           bounds: Optional[List[int]] = None) -> str:
    """Generates part of the command string corresponding to the provided
    switch to be passed to optool

    Parameters
    ----------
    value: any
        The input value
    switch: str
        The switch to be used in the command string
    bounds: List[int], optional
        The bounds of the switch to be used in the command string. If value is
        out of bounds an error will be raised

    Returns
    -------
    str
        Part of the command line argument passed to optool
    """
    if value is None:
        return ""
    elif bounds is not None:
        if value < bounds[0] or value > bounds[1]:
            raise ValueError(f"Value {value} out of bounds for {bounds}")
    elif isinstance(value, (Path, str, float, int)):
        return f"-{switch} {str(value)}"
    elif isinstance(value, List):
        return f"-{switch} {' '.join([str(i) for i in value])}"
    elif isinstance(value, Dict):
        return f"-{switch} {' '.join([f'{k} {v}' for k, v in value.items()])}"
    else:
        raise IOError(f"Input type '{type(value)}' for '{switch}' is not supported!")


class oimOptoolBackend(op.particle):
    """A wrapper for the dust distribution calculation with optool that
    inherits much of the optool.particle's class functionality.

    Once calculated the (.dat)-files (if no other output path provided) will be
    stored internally in 'oimodeler/data/cache/optool/' (if no other path is
    provided) to increase the speed of future calculations

    Parameters
    ----------
    grains: Path | str | Dict[str, int]
        Path to a (.lnk)-file or string (name) of a material or dictionary
        containing materials to include in the grain. If a dictionary is passed
        then it needs to contain the name of the material as key and its mass
        fraction as value (e.g., {"pyr-mg70": 0.87, "c": 0.13}).
        Mass fractions do not have to add up to 1, they will be renormalized.
        Up to 20 materials can be specified to build up a grain
    grain_mantels: Path | str | Dict[str, int], optional
        The same as the grains parameter, but the material will be placed into
        the grain mantle. Multiple mantle materials will be mixed using the
        Bruggeman rule, and than that mix will be added to the core using the
        Maxwell-Garnett rule
    porosity: int | List[int], optional
        The volume fraction of vacuum, a number smaller than 1. The default is
        0. A single value will apply to both core and mantle, but if a list
        with two values is provided the second value will be specific for the
        mantle (and may be 0)
    dust_distribution: str, optional
        Some preset dust distributions. If toggled all other inputs will be
        ignored. Use DIANA (Woitke+2016) "diana" or DSHARP (Birnstiel+2018)
        "dsharp" compositions (or DSHARP without ice "dsharp-no-ice")
    computational_method: str, optional
        Various computational methods for the dust distribution. Possible
        choices are "dhs", "mff", "mie" and "cde". The default is "dhs" with
        f_max of 0.8
    f_max: float, optional
        The default is 0.8. A value of 0
        means to use solid spheres (Mie theory), i.e. perfectly regular grains
    monomer_radius: float, optional
        The monomer radius for the Modified Mean Field theory.
        The default is 0.1μm
    dfrac_or_fill: float, optional
        Either the fractal dimension (if > 1) or the volume filling factor
        (if < 1) for the Modified Mean Field theory. The default is 0.2
    prefactor: float, optional
        The prefactor for the Modified Mean Field theory
    scat: bool, optional
        Include the scattering matrix in the output
    gs_min: float, optional
        The minimum grain radius. Can be specified without a maximum grain
        radius nor a number of size bins. The default value is 0.05 µm
    gs_max: float, optional
        The maximum grain radius. The default value is 3000 µm
    gs_pow: float, optional
        The size distribution powerlaw. The default is 3.5
    gs_mean: float, optional
        The centroid size for a log-normal size distribution
    gs_sigma: float, optional
        The logarithmic width for a log-normal size distribution. If it is
        negative then a normal distribution with that width [µm] around gs_mean
        is created
    ngs: int, optional
        The number of size bins. The default is 15 per size decade with a fixed
        minimum of 5
    wl_min: float, optional
        The minimum wavelength. Can be specified without a maximum wavelength
        nor a the number of wavelength points. The default value is 0.05 µm
    wl_max: float, optional
        The maximum wavelength. The default value is 10000 µm
    nwl: int, optional
        The number of wavelength points for the construction of the wavelength
        grid. The default is 300
    wavelength_file: Path | str, optional
        Read the wavelength grid from a (.dat)-file. To get an example file
        'optool_lam.dat', execute this script with the option 'wgrid=True'.
        Otherwise, an (.lnk)-file could be used here as well
    nang: int, optional
        The number of evenly spaced angular grid points that cover a range
        between 0 and 180 degrees. The default for nang is 180
    nsub: int, optional
        Divide the computation up into parts to produce a file for each grain
        size. Each size will be an average over a range of nsub grains around
        the real size. The default is 5
    ndeg: int, optional
        Cap the first degrees of the forward scattering peak. The default is 2
    fits: bool, optional
        Write into a (.fits)-file instead to ASCII (.dat). With nsub, write
        files that amount to the grain size bins 'ngs'.
    radmc_label: str, optional
        If a label is provided then the file names will contain the label and
        have the extension (.inp). This makes it so the files can be used as
        input for RADMC-3D, which uses a different angular grid and scattering
        matrix normalization
    wgrid: bool, optional
        Create the additional files optool_sd.dat and optool_lam.dat with the
        grain size distribution and the wavelengths grid, respectively
    cache_dir: Path, optional
        The directory to store the (.dat)-files in so they can be read instead
        of recomputed the next time the same command is used. The cache is
        automatically cleared when CMD changes between runs. The default is
        'oimodeler/data/cache/optool/'

    Attributes
    ----------
    cmd: str
        The full command given in the particle (super.__init__()) call
    radmc: bool
        Indicates if the output is to RADMC conventions
    scat: bool
        Indicates if the scattering matrix is available
    nlam: int
        The Number of wavelength points
    lam: np.ndarray
        The wavelength grid [µm]
    nang: int
        The number of scattering angles
    scatang: np.ndarray
        The angular grid
    materials: np.ndarray
        Lists containing [location, m_{frac}, rho, material]
    np: int
        The number of particles, either 1 or (with -d) corresponds to the
        number of grain size bins 'ngs'
    fmax: np.ndarray
        The maximum volume fraction of vacuum for DHS
    pcore, pmantle: np.ndarray
        The Porosity of the core/mantle material
    amin: np.ndarray
        The minimum grain size used for each particle
    amax: np.ndarray
        The maximum grain size used for each particle
    nsub: np.ndarray
        The number of sizes averaged for each particle
    apow: np.ndarray
        The negative size distribution power law (e.g., 3.5)
    amean: np.ndarray
        The mean size for (log-) normal size distributions
    asig: np.ndarray
        The standard deviation for (log-) normal distribution
    a1: np.ndarray
        The mean grain radius [µm]
    a2: np.ndarray
        The radius of the grain with mean surface area [µm]
    a3: np.ndarray
        The radius of the grain with mean volume [µm]
    rho: np.ndarray
        The specific density of grains
    kabs: np.ndarray
        The absorption cross section [cm^2 g^-1]
    ksca: np.ndarray
        The scattering cross section [cm^2 g^-1]
    kext: np.ndarray
        The extinction cross section [cm^2 g^-1]
    gsca: np.ndarray
        The asymmetry parameter
    f11, ..., f44: float[np, nlam, nang]
        The scattering matrix elements F_11, ... ,F_44
    chop: np.ndarray
        The degrees chopped off forward scattering
    tmin: float
        The minimum temperature for mean opacities
    tmax: float
        The maximum temperature for mean opacities
    ntemp: int
        The mumber of temperatures for mean opacities
    temp: np.ndarray
        The Temperatures used to calculate the mean opacities
    kplanck: np.ndarray
        The Planck mean opacities (are only accessible after calling the
        computemean method)
    kross: np.ndarray
         Rosseland mean opacities (are only accessible after calling the
         computemean method)
    norm: str
        The current scattering matrix normalization

    Methods
    --------
    plot(minkap=1e0)
        Create interactive plots of the opacities
    select(i)
        Select just one bin from a multi-particle object.
        A multi-particle object is produced when running optool with
        a -d switch.
    sizedist(N_of_a)
        Compute opacity of a size distribution of the elements
    scatnorm(norm='')
        Check or change the normalization of the scattering matrix
    computemean(tmin=10, tmax=1500, ntemp=100)
        Compute mean opacities from the opacities

    See also
    --------
    optool.particle(cmd, cache="")
        Run optool and turn output into a python object
    optool.lnktable(file, i_lnk=[1, 2, 3], nskip=0, nlam_rho=True)
        Class to work with lnk files

    Notes
    ------
    * Computational Methods

    - The Distribution of Hollow Spheres (DHS, Min+ 2005) approach to model
      deviations from perfect spherical symmetry and low-porosity aggregates.
      Spheres with inner holes with volume fractions between 0 and f_max
      are averaged to mimic irregularities.
    - The Modified Mean Field theory (MMF, Tazaki & Tanaka 2018)
      to compute opacities of highly porous or fractal aggregates.
      The grain material, grain mantel material and the porosity determine
      the composition of monomers.
    - A standard Mie calculation for perfect spheres. This is short for DHS
      with an f_max of 0.
    - Compute CDE (continuous distribution of ellipsoids) Rayleigh limit
      opacities.

    * Available Materials

    amorph.pyroxenes  pyr-mg100/95/80/70/60/50/40
    amorph.olivines   ol-mg100/40                       (Dorschner95,Henning96)
    cryst. pyr/ol     pyr-c-mg96 ol-c-mg100/95/00     (Jäger96,Suto06,Fabian01)
    other silicates   astrosil                                       (Draine03)
    amorphous carbon  c-z    c-p                          (Zubko96,Preibisch93)
    graphite,special  c-gra  c-nano  c-org               (Dra.03,Mut.04,Hen.96)
    quartz,corundum   sio2   cor-c                         (Kitamura07,Koike95)
    iron/sulfide      fe-c   fes                                    (Henning96)
    carbides          sic                                            (Draine93)
    water ice         h2o-w  h2o-a                         (Warren08,Hudgins93)
    other ices        co2-w  nh3-m                      (Warren86,Martonchik83)
                      co-a   co2-a/c ch4-a/c ch3oh-a/c  (Palumbo06,Gerakines20)

    - The abbreviations for the optool CLI (and/or for the grains/grain_mantels
    parameters) are:

    pyr  -> pyr-mg70     c    -> c-z         iron -> fe-c      h2o  -> h2o-w
    ol   -> ol-mg50      gra  -> c-gra       qua  -> sio2      co   -> co-a
    ens  -> pyr-c-mg96   org  -> c-org       cor  -> cor-c     co2  -> co2-w
    for  -> ol-c-mg100                       tro  -> fes       nh3  -> nh3-m
    fay  -> ol-c-mg00

    * Optool Documentation
    For more information on the optool python extension see
    https://github.com/cdominik/optool/blob/b68e67e5abdf8746078d3767acbba4e28fcbc75f/UserGuide.org#sizedist

    * Optool Installation Instructions
    You can download, compile, and install optool with these simple steps,
    using the freely available GNU FORTRAN compiler 'gfortran'.

    # Clone repository
    >>> git clone https://github.com/cdominik/optool.git

    # Enter code directory and compile with multicore support
    >>> cd optool; make multi=true &&

    # Copy binaries to binary path (The 'bindir' needs to be changed and may
    vary by OS!). In the case of Linux or MacOs, if no root permissions are
    available or other errors occur try setting bindir to '/usr/local/bin/'
    >>> make install bindir=~/bin/

    In the compilation step, use 'multi=true' to add multicore support
    (recommended!), 'ifort=true' to use the Intel fortran compiler, fits=true to
    support FITS files (This requires the 'cfitsio' library to be installed on
    your system.), and 'oldio=true' if your compiler does not have the
    ISO_FORTRAN_ENV module.

    The executable is called optool. The make install step copies it and also
    optool2tex and optool-complete into bin-dir.
    For shell command line completion support, check the file optool-complete.

    Warnings
    --------
    For this backend to work the optool bin must be installed on your computer.
    See installation instructions under Notes
    """

    def __init__(self, grains: Path | str | Dict[str, int] = None,
                 grain_mantels: Optional[Path | str | Dict[str, int]] = None,
                 porosity: Optional[int | List[int]] = None,
                 dust_distribution: Optional[str] = None,
                 computational_method: Optional[str] = None,
                 f_max: Optional[float] = 0.8,
                 monomer_radius: Optional[float] = 0.1,
                 dfrac_or_fill: Optional[float] = 0.2,
                 prefactor: Optional[float] = None,
                 scat: Optional[bool] = False,
                 gs_min: Optional[float] = None,
                 gs_max: Optional[float] = None,
                 gs_pow: Optional[float] = None,
                 gs_mean: Optional[float] = None,
                 gs_sigma: Optional[float] = None,
                 ngs: Optional[int] = None,
                 grain_size_file: Optional[Path | str] = None,
                 wl_min: Optional[float] = None,
                 wl_max: Optional[float] = None,
                 nwl: Optional[int] = None,
                 wavelength_file: Optional[Path | str] = None,
                 nang: Optional[int] = None,
                 nsub: Optional[int] = None,
                 ndeg: Optional[int] = None,
                 fits: Optional[bool] = False,
                 radmc_label: Optional[str] = None,
                 wgrid: Optional[bool] = False,
                 cache_dir: Optional[Path] = None):
        """Constructs the command for the optool particle's class and then executes it"""
        if not find_executable("optool"):
            raise RuntimeError("The 'optool' executable has not been found!")

        if cache_dir is None:
            cache_dir = Path(pkg_resources.resource_filename("oimodeler",
                                                             "data/cache/optool/"))

        cmd_arguments = ["optool"]

        if dust_distribution is None:
            # NOTE: Grain composition
            cmd_arguments.append(generate_switch_string(grains, "c"))
            cmd_arguments.append(generate_switch_string(grain_mantels, "m"))
            cmd_arguments.append(generate_switch_string(porosity, "p"))

            # NOTE: Grain geometry and computational method
            if computational_method == "dhs":
                cmd_arguments.append(generate_switch_string(f_max, "dhs"))
            elif computational_method == "mmf":
                mmf_switch = generate_switch_string(monomer_radius, "mmf")
                mmf_switch += f" {dfrac_or_fill}"
                mmf_switch += f" {prefactor}" if prefactor is not None else ""
                cmd_arguments.append(mmf_switch)
            elif computational_method in ["mie", "cde"]:
                cmd_arguments.append(generate_switch_string("", computational_method))
            else:
                raise IOError("No such computational_method:"
                              f" {computational_method} is supported!")

            # NOTE: Grain size distribution and wavelength grid
            cmd_arguments.append(generate_switch_string(gs_min, "amin"))
            cmd_arguments.append(generate_switch_string(gs_max, "amax"))
            cmd_arguments.append(generate_switch_string(gs_pow, "apow"))
            cmd_arguments.append(generate_switch_string(gs_mean, "amean"))
            cmd_arguments.append(generate_switch_string(gs_sigma, "asig"))
            cmd_arguments.append(generate_switch_string(ngs, "na"))
            cmd_arguments.append(generate_switch_string(wl_min, "lmin"))
            cmd_arguments.append(generate_switch_string(wl_max, "lmax"))
            cmd_arguments.append(generate_switch_string(nwl, "nlam"))

            # NOTE: Output control
            cmd_arguments.append(generate_switch_string(nang, "s", [0, 180]))
            cmd_arguments.append(generate_switch_string(nsub, "d"))
            cmd_arguments.append(generate_switch_string(ndeg, "chop"))
            cmd_arguments.append(generate_switch_string(cache_dir, "o"))
            cmd_arguments.append(generate_switch_string("", "fits") if fits else "")
            cmd_arguments.append(generate_switch_string(radmc_label, "radmc"))
            cmd_arguments.append(generate_switch_string("", "wgrid") if wgrid else "")
        else:
            # NOTE: Uses specific, predetermined dust distribution
            cmd_arguments.append(generate_switch_string("", dust_distribution))
        cmd = " ".join([switch for switch in cmd_arguments if switch])
        print(f"[Calling] {cmd}")

        # NOTE: Call the optool executable via inheritance of the optool
        # particle class
        super().__init__(cmd, str(cache_dir))

        self.lam, self.gsca = map(lambda x: x.squeeze(), (self.lam, self.gsca))
        self.kabs, self.ksca, self.kext = map(lambda x: x.squeeze(),
                                              (self.kabs, self.ksca, self.kext))

        # TODO: Implement renaming of the files
        # TODO: Maybe use the write function for it?
