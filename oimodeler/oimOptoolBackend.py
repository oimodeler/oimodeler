import pkg_resources
import subprocess
from datetime import datetime
from distutils.spawn import find_executable
from pathlib import Path
from typing import Optional, Any, List, Dict

import numpy as np
import optool as op
import toml


def make_wavelength_file(wavelength_solution: np.ndarray):
    """Makes a (.dat)-file from a wavelength solution to be passed to optool"""
    if wavelength_solution is not None:
        wavelength_file = storage_dir / "wavelength_solution.dat"
        np.savetxt(wavelength_file, np.array(wavelength_solution))


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
        raise IOError(
            f"Input type '{type(value)}' for '{switch}' is not supported!")


# TODO: Also write number of files into (.toml)-file?
class oimOptoolBackend(op.particle):
    """A wrapper class for optool  that can either load a file or calculate
    the dust distribution.

    Inherits much of the optool.particle's class functionality.

    Once calculated the (.dat)-files (if no other output path provided) will be
    stored in a specified folder internally in 'oimodeler/data/cache/optool/'
    (if no other path is provided) to increase the speed of future calculations

    Parameters
    ----------
    grains: Path | str | Dict[str, int]
        Path to a (.lnk)-file or string (name) of a material or dictionary
        containing materials to include in the grain. If a dictionary is
        passed then it needs to contain the name of the material as key and
        its mass fraction as value (e.g., {"pyr-mg70": 0.87, "c": 0.13}).
        Mass fractions do not have to add up to 1, they will be
        renormalized.
        Up to 20 materials can be specified to build up a grain
    grain_mantels: Path | str | Dict[str, int], optional
        The same as the grains parameter, but the material will be placed
        into the grain mantle. Multiple mantle materials will be mixed
        using the Bruggeman rule, and than that mix will be added to the
        core using the Maxwell-Garnett rule
    porosity: int | List[int], optional
        The volume fraction of vacuum, a number smaller than 1. The default
        is 0. A single value will apply to both core and mantle, but if a
        list with two values is provided the second value will be specific
        for the mantle (and may be 0)
    dust_distribution: str, optional
        Some preset dust distributions. DISCLAIMER: If toggled all other
        inputs will be ignored. Use DIANA (Woitke+2016) "diana" or DSHARP
        (Birnstiel+2018) "dsharp" compositions (or DSHARP without ice
        "dsharp-no-ice")
    computational_method: str, optional
        Various computational methods for the dust distribution. Possible
        choices are "dhs", "mff", "mie" and "cde". The default is "dhs"
        with f_max of 0.8
    f_max: float, optional
        The default is 0.8. A value of 0 means to use solid spheres
        (Mie theory), i.e. perfectly regular grains
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
        negative then a normal distribution with that width [µm] around
        gs_mean is created
    ngs: int, optional
        The number of size bins. The default is 15 per size decade with a
        fixed minimum of 5
    gs_file: Path, optional
        The path to a (.dat)-file containing the grain size distribution
    wl_min: float, optional
        The minimum wavelength. Can be specified without a maximum
        wavelength nor a the number of wavelength points.
        The default value is 0.05 µm
    wl_max: float, optional
        The maximum wavelength. The default value is 10000 µm
    nwl: int, optional
        The number of wavelength points for the construction of the
        wavelength grid. The default is 300
    wavelength_file: Path | str, optional
        Read the wavelength grid from a (.dat)-file. To get an example file
        'optool_lam.dat', execute this script with the option 'wgrid=True'.
        Otherwise, an (.lnk)-file could be used here as well
    wavelength_solution: np.ndarray, optional
        A wavelength solution that will be written to a temporary file and
        used for the wavelength file parameter/switch in the optool
        executable.
        DISCLAIMER: The wavelength solution will take precendent over a
        provided wavelength file (if both is given)
    nang: int, optional
        The number of evenly spaced angular grid points that cover a range
        between 0 and 180 degrees. The default for nang is 180
    nsub: int, optional
        Divide the computation up into parts to produce a file for each
        grain size. Each size will be an average over a range of nsub
        grains around the real size. The default is 5
    ndeg: int, optional
        Cap the first degrees of the forward scattering peak.
        The default is 2
    fits: bool, optional
        Write into a (.fits)-file instead to ASCII (.dat). With nsub, write
        files that amount to the grain size bins 'ngs'.
    radmc_label: str, optional
        If a label is provided then the file names will contain the label
        and have the extension (.inp). This makes it so the files can be
        used as input for RADMC-3D, which uses a different angular grid and
        scattering matrix normalization
    wgrid: bool, optional
        Create the additional files optool_sd.dat and optool_lam.dat with
        the grain size distribution and the wavelengths grid, respectively
    cmd: str, optional
        An optool compliant command line argument that is called from within
        python. If this is provided, all other optool related arguments will
        be ignored
    opacity_file: Path | str, optional
        A (.dat)-file already containing an dust distribution in the optool
        format. If the scattering matrix is included, set the scat paramter
        to True as well
    cache_dir: Path, optional
        The directory to create a folder in to store the (.dat)-files so they
        can be read instead of recomputed the next time the same command
        is used.
        The default directory is 'oimodeler/data/cache/optool/'

    Attributes
    ----------
    cmd: str
        The full command given in the particle call
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
    cache_dir: Path
        The directory to create a folder in to store the (.dat)-files so they
        can be read instead of recomputed the next time the same command
        is used
    storage_dir: Path
        The specific directory to store the (.dat)-files for this calculation

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
        Class that runs optool and turns output into a python object as well as
        providing additional functionality
    optool.lnktable(file, i_lnk=[1, 2, 3], nskip=0, nlam_rho=True)
        Class to work with lnk files

    Notes
    ------
    * Computational Methods

    - The Distribution of Hollow Spheres (DHS, Min+ 2005) approach to model
      deviations from perfect spherical symmetry and low-porosity
      aggregates.
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

    * Optool Installation Instructions (Needed for optool calculation)
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
    For this backend to work the optool-bin must be installed on your computer
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
                 gs_file: Optional[Path | str] = None,
                 wl_min: Optional[float] = None,
                 wl_max: Optional[float] = None,
                 nwl: Optional[int] = None,
                 wavelength_file: Optional[Path | str] = None,
                 wavelength_solution: Optional[np.ndarray | List] = None,
                 nang: Optional[int] = 180,
                 nsub: Optional[int] = None,
                 ndeg: Optional[int] = None,
                 fits: Optional[bool] = False,
                 radmc_label: Optional[str] = None,
                 wgrid: Optional[bool] = False,
                 cmd: Optional[str] = None,
                 opacity_file: Optional[Path] = None,
                 cache_dir: Optional[Path] = None) -> None:
        """Instantiates a version of the class"""
        if cache_dir is None:
            self.cache_dir = Path(pkg_resources.resource_filename("oimodeler",
                                                                  "data/cache/optool/"))
        else:
            self.cache_dir = cache_dir

        self.scat = scat
        self.storage_dir = self.make_storage_dir_name()

        self.np, self.masscale = [1]*2
        self.materials, self.rho = [], []

        if wavelength_solution is not None:
            wavelength_file = make_wavelength_file(wavelength_solution)
        cmd_output = self.make_cmd(grains, grain_mantels, porosity,
                                   dust_distribution, computational_method,
                                   f_max, monomer_radius, dfrac_or_fill,
                                   prefactor, gs_min, gs_max,
                                   gs_pow, gs_mean, gs_sigma, ngs, gs_file,
                                   wl_min, wl_max, nwl, wavelength_file,
                                   wavelength_solution, nang, nsub, ndeg,
                                   fits, radmc_label, wgrid, cmd)

        self.cmd = cmd_output.split("-o")[0].strip()

        if opacity_file is not None:
            self.read_files(opacity_file)
        elif (storage_dir := self.cmd_in_cache()) is not None:
            self.storage_dir = storage_dir
            self.read_cache()
        else:
            self.run_optool(cmd_output)

    def run_optool(self, cmd_output: str) -> None:
        """Runs optool to calculate the opacities for the input parameters

        Parameters
        ----------
        cmd_output: str
            The full command containing the output path for the files
        """
        if not find_executable("optool"):
            raise RuntimeError("The 'optool' executable has not been found!")

        if not self.storage_dir.exists():
            self.storage_dir.mkdir(parents=True)

        with open(self.storage_dir / "calculation_info.toml", "w+") as toml_file:
            toml.dump({"cmd": self.cmd, "scat": self.scat}, toml_file)

        print(f"[Calling] {cmd_output}")
        subprocess.run(cmd_output, shell=True, check=True)

        self.read_cache()

    def make_cmd(self, grains: Path | str | Dict[str, int] = None,
                 grain_mantels: Optional[Path | str | Dict[str, int]] = None,
                 porosity: Optional[int | List[int]] = None,
                 dust_distribution: Optional[str] = None,
                 computational_method: Optional[str] = None,
                 f_max: Optional[float] = 0.8,
                 monomer_radius: Optional[float] = 0.1,
                 dfrac_or_fill: Optional[float] = 0.2,
                 prefactor: Optional[float] = None,
                 gs_min: Optional[float] = None,
                 gs_max: Optional[float] = None,
                 gs_pow: Optional[float] = None,
                 gs_mean: Optional[float] = None,
                 gs_sigma: Optional[float] = None,
                 ngs: Optional[int] = None,
                 gs_file: Optional[Path | str] = None,
                 wl_min: Optional[float] = None,
                 wl_max: Optional[float] = None,
                 nwl: Optional[int] = None,
                 wavelength_file: Optional[Path | str] = None,
                 wavelength_solution: Optional[np.ndarray | List] = None,
                 nang: Optional[int] = 180,
                 nsub: Optional[int] = None,
                 ndeg: Optional[int] = None,
                 fits: Optional[bool] = False,
                 radmc_label: Optional[str] = None,
                 wgrid: Optional[bool] = False,
                 cmd: Optional[str] = None) -> str:
        """Generates the optool's command line arguments from the values passed
        to the class

        Refer to the class's documentation
        """
        cmd_arguments = ["optool"]
        if dust_distribution is None:
            # NOTE: Grain composition
            cmd_arguments.append(generate_switch_string(grains, "c"))
            cmd_arguments.append(generate_switch_string(grain_mantels, "m"))
            cmd_arguments.append(generate_switch_string(porosity, "p"))

            # NOTE: Grain geometry and computational method
            if computational_method is None or computational_method == "dhs":
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
            if gs_file is not None:
                cmd_arguments.append(generate_switch_string(gs_file, "a"))
            else:
                cmd_arguments.append(generate_switch_string(gs_min, "amin"))
                cmd_arguments.append(generate_switch_string(gs_max, "amax"))
                cmd_arguments.append(generate_switch_string(gs_pow, "apow"))
                cmd_arguments.append(generate_switch_string(gs_mean, "amean"))
                cmd_arguments.append(generate_switch_string(gs_sigma, "asig"))
                cmd_arguments.append(generate_switch_string(ngs, "na"))

            if wavelength_file is not None:
                cmd_arguments.append(generate_switch_string(wavelength_file, "l"))
            else:
                cmd_arguments.append(generate_switch_string(wl_min, "lmin"))
                cmd_arguments.append(generate_switch_string(wl_max, "lmax"))
                cmd_arguments.append(generate_switch_string(nwl, "nlam"))

            # NOTE: Output control
            if self.scat:
                cmd_arguments.append(generate_switch_string(nang, "s", [0, 180]))
            cmd_arguments.append(generate_switch_string(nsub, "d"))
            cmd_arguments.append(generate_switch_string(ndeg, "chop"))
            cmd_arguments.append(generate_switch_string("", "fits") if fits else "")
            cmd_arguments.append(generate_switch_string(radmc_label, "radmc"))
            cmd_arguments.append(generate_switch_string("", "wgrid") if wgrid else "")
        else:
            # NOTE: Uses specific, predetermined dust distribution
            cmd_arguments.append(generate_switch_string("", dust_distribution))
        cmd_arguments.append(generate_switch_string(self.storage_dir, "o"))

        return " ".join([switch for switch in cmd_arguments if switch])\
            if cmd is None else cmd

    def make_storage_dir_name(self) -> Path:
        """Makes a timestamp of the calculation for the cache directory name

        Returns
        -------
        storage_dir: Path
            The directory in which the cached files are to be stored
        """
        dir_name = str(datetime.now()).replace(" ", "_").replace(":", "-")
        return self.cache_dir / dir_name

    def cmd_in_cache(self) -> Optional[Path]:
        """Checks if the command is in the cache directory and returns its
        path if found otherwise None"""
        for toml_path in self.cache_dir.rglob("*.toml"):
            with open(toml_path, "r") as toml_file:
                calculation_info = toml.load(toml_file)
                cmd = calculation_info["cmd"] if "cmd" in calculation_info else ""
                if cmd == self.cmd:
                    return toml_path.parent
        return None

    # TODO: Iterate over all cached files
    def read_cache(self) -> None:
        """Reads in all files contained in a cache directory and stores
        their values in the cache"""
        for toml_path in self.cache_dir.rglob("*.toml"):
            with open(toml_path, "r") as toml_file:
                calculation_info = toml.load(toml_file)
                self.cmd = calculation_info["cmd"]
                self.scat = calculation_info["scat"]
        if self.scat:
            filename = "*dustkapscatmat*"
        else:
            filename = "*dustkappa*"

        dust_files = list(self.storage_dir.rglob(filename))
        self.np = len(dust_files)
        self.read_files(dust_files)
        print(f"Read in cached files from '{self.storage_dir}'")

    def read_files(self, files: Path | list[Path]) -> None:
        """Reads one or more opacity files in and stores the values in the class"""
        self.header = []
        self.kabs, self.ksca,\
            self.kext, self.gsca, self.scatang = [], [], [], [], []
        self.f11, self.f12, self.f22,\
            self.f33, self.f34, self.f44 = [], [], [], [], [], []

        if isinstance(files, Path):
            files = [files]

        for file in files:
            header, *rest = op.readoutputfile(file, self.scat)
            self.header.append(header)
            if self.scat:
                lam, kabs, ksca,\
                    phase_g, scatang,\
                    f11, f12, f22, f33,\
                    f34, f44 = map(np.squeeze, rest)
                self.scatang.append(scatang)
                self.f11.append(f11)
                self.f12.append(f12)
                self.f22.append(f22)
                self.f33.append(f33)
                self.f34.append(f34)
                self.f44.append(f44)
            else:
                lam, kabs,\
                    ksca, phase_g = map(np.squeeze, rest)
            kext = kabs + ksca
            self.kabs.append(kabs)
            self.ksca.append(ksca)
            self.kext.append(kext)
            self.gsca.append(phase_g)
        self.nang = len(self.scatang) if self.scatang else 0
