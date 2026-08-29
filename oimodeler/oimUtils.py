# -*- coding: utf-8 -*-
"""Various utilities for optical interferometry"""

from __future__ import annotations

import base64
import csv
import importlib
import io
import pickle
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import toml
from astropy.coordinates import Angle
from astropy.io import fits
from astroquery.simbad import Simbad
from numpy.typing import ArrayLike
from scipy.stats import circstd

import oimodeler as oim

from .oimOptions import constants as const
from .oimOptions import oimOptions

# TODO: Should this (global variables) be moved into a configuration file?
_oimDataType = ["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI", "FLUXDATA"]
_oimDataTypeErr = [
    "VIS2ERR",
    "VISAMPERR",
    "VISPHIERR",
    "T3AMPERR",
    "T3PHIERR",
    "FLUXERR",
]
_oimDataTypeArr = ["OI_VIS2", "OI_VIS", "OI_VIS", "OI_T3", "OI_T3", "OI_FLUX"]
_oimDataAnalysisInComplex = [False, False, True, False, True, False]
_cutArr = [
    "EFF_WAVE",
    "EFF_BAND",
    "VIS2DATA",
    "VIS2ERR",
    "FLAG",
    "VISAMP",
    "VISAMPERR",
    "VISPHI",
    "VISPHIERR",
    "T3AMP",
    "T3AMPERR",
    "T3PHI",
    "T3PHIERR",
    "FLUX",
    "FLUXDATA",
    "FLUXERR",
    "FLAG",
]

OI_TARGET_KEYWORDS = [("OI_REVN", False, "Revision number")]
OI_TARGET_COLUMNS = [
    ("TARGET_ID", "I", False, "Index number. Must be >=1", None),
    ("TARGET", "16A", False, "Target name", None),
    ("RAEP0", "D", False, "RA at mean EQUINOX ", "deg"),
    ("DECEP0", "D", False, "Dec at mean EQUINOX", "deg"),
    ("EQUINOX", "E", False, "Equinox", None),
    ("RA_ERR", "D", False, "Error in RA", "deg"),
    ("DEC_ERR", "D", False, "Error in Dec", "deg"),
    ("SYSVEL", "D", False, "Systemic radial velocity", "m/s"),
    (
        "VELTYP",
        "8A",
        False,
        "Reference for radial velocity:LSR, GEOCENTR...",
        None,
    ),
    (
        "VELDEF",
        "8A",
        False,
        "Definition of radial velocity:(OPTICAL,RADIO)",
        None,
    ),
    ("PMRA", "D", False, "Proper motion in RA", "deg/yr"),
    ("PMDEC", "D", False, "Proper motion in Dec", "deg/yr"),
    ("PMRA_ERR", "D", False, "Error of proper motion in RA", "deg/yr"),
    ("PMDEC_ERR", "D", False, "Error of proper motion in Dec", "deg/yr"),
    ("PARALLAX", "E", False, "Parallax", "deg"),
    ("PARA_ERR", "E", False, "Error in parallax ", "deg"),
    ("SPECTYP", "16A", False, "Spectral type", "deg"),
    ("CATEGORY", "3A", True, "CAL or SCI", None),
]
OI_ARRAY_KEYWORDS = [
    ("OI_REVN", False, "Revision number"),
    ("ARRNAME", False, "A Array name, for cross-referencing"),
    ("FRAME", False, "A Coordinate frame"),
    ("ARRAYX", False, "Array center x coordinates (m)"),
    ("ARRAYY", False, "Array center y coordinates (m)"),
    ("ARRAYZ", False, "Array center z coordinates (m)"),
]
OI_ARRAY_COLUMNS = [
    ("TEL_NAME", "16A", False, " Telescope name", None),
    ("STA_NAME", "16A", False, "Station name", None),
    ("STA_INDEX", "I", False, "Station number. Must be >=1", None),
    ("DIAMETER", "E", False, "Element diameter", "m"),
    ("STAXYZ", "3D", False, "Station coordinates w.r.t. array center", "m"),
    ("FOV", "D", False, " Photometric field of view", "arcsec"),
    ("FOVTYPE", "6A", False, "Model for FOV: FWHM or RADIUS", None),
]
OI_WL_KEYWORDS = [
    ("OI_REVN", False, "Revision number"),
    ("INSNAME", False, "Identifies corresponding OI_WAVELENGTH table"),
]
OI_WL_COLUMNS = [
    ("EFF_WAVE", "E", False, "Effective wavelength of channel", "m"),
    ("EFF_BAND", "E", False, "Effective bandpass of channel", "m"),
]
OI_VIS2_KEYWORDS = [
    ("OI_REVN", False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME", False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME", False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True, "Identifies corresponding OI_CORR table"),
]
OI_VIS2_COLUMNS = [
    (
        "TARGET_ID",
        "I",
        False,
        "Target number as index into OI_TARGET table",
        "m",
    ),
    ("TIME", "D", False, "Zero. For backwards compatibility", None),
    ("MJD", "D", False, "Modified Julian day", "day"),
    ("INT_TIME", "D", False, "Integration time", None),
    ("VIS2DATA", "NWLD", False, "Squared Visibility", None),
    ("VIS2ERR", "NWLD", False, "Error in Squared Visibility", None),
    (
        "CORRINDX_VIS2DATA",
        "J",
        True,
        "Index into correlation matrix for 1st VIS2DATA element",
        None,
    ),
    ("UCOORD", "D", False, "U coordinate of the data", "m"),
    ("VCOORD", "D", False, "V coordinate of the data", "m"),
    (
        "STA_INDEX",
        "2I",
        False,
        "Station numbers contributing to the data",
        None,
    ),
    ("FLAG", "NWLL", False, "Flag", None),
]
OI_VIS_KEYWORDS = [
    ("OI_REVN", False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME", False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME", False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True, "Identifies corresponding OI_CORR table"),
    ("AMPTYP", True, "absolute, differential, correlated flux"),
    ("PHITYP", True, "absolute, differential"),
    (
        "AMPORDER",
        True,
        "Polynomial fit order for differential chromatic amplitudes",
    ),
    (
        "PHIORDER",
        True,
        "Polynomial fit order for differential chromatic phases",
    ),
]

# TODO: Implement RIV VISREFMAP ...
OI_VIS_COLUMNS = [
    (
        "TARGET_ID",
        "I",
        False,
        "Target number as index into OI_TARGET table",
        "m",
    ),
    ("TIME", "D", False, "Zero. For backwards compatibility", None),
    ("MJD", "D", False, "Modified Julian day", "day"),
    ("INT_TIME", "D", False, "Integration time", "s"),
    ("VISAMP", "NWLD", False, "Visibility amplitude", None),
    ("VISAMPERR", "NWLD", False, "Error in visibility amplitude", None),
    (
        "CORRINDX_VISAMP",
        "J",
        True,
        "Index into correlation matrix for 1st VISAMP element",
        None,
    ),
    ("VISPHI", "NWLD", False, "Visibility phase", "deg"),
    ("VISPHIERR", "NWLD", False, "Error in visibility Phase", "deg"),
    (
        "CORRINDX_VISPHI",
        "J",
        True,
        "Index into correlation matrix for 1st VISPHI element",
        None,
    ),
    ("UCOORD", "D", False, "U coordinate of the data", "m"),
    ("VCOORD", "D", False, "V coordinate of the data", "m"),
    (
        "STA_INDEX",
        "2I",
        False,
        "Station numbers contributing to the data",
        None,
    ),
    ("FLAG", "NWLL", False, "Flag", None),
]
OI_T3_KEYWORDS = [
    ("OI_REVN", False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME", False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME", False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True, "Identifies corresponding OI_CORR table"),
]
OI_T3_COLUMNS = [
    (
        "TARGET_ID",
        "I",
        False,
        "Target number as index into OI_TARGET table",
        "m",
    ),
    ("TIME", "D", False, "Zero. For backwards compatibility", None),
    ("MJD", "D", False, "Modified Julian day", "day"),
    ("INT_TIME", "D", False, "Integration time", "s"),
    ("T3AMP", "NWLD", False, "Triple product amplitude", None),
    ("T3AMPERR", "NWLD", False, "Error in triple product amplitude", None),
    (
        "CORRINDX_T3AMP",
        "J",
        True,
        "Index into correlation matrix for 1st T3AMP element",
        None,
    ),
    ("T3PHI", "NWLD", False, "Triple Product Phase", "deg"),
    ("T3PHIERR", "NWLD", False, "Error in Triple Product Phase", "deg"),
    (
        "CORRINDX_T3PHI",
        "J",
        True,
        "Index into correlation matrix for 1st T3PHI element",
        None,
    ),
    (
        "U1COORD",
        "D",
        False,
        "U coordinate of baseline AB of the triangle",
        "m",
    ),
    (
        "V1COORD",
        "D",
        False,
        "V coordinate of baseline AB of the triangle",
        "m",
    ),
    (
        "U2COORD",
        "D",
        False,
        "U coordinate of baseline BC of the triangle",
        "m",
    ),
    (
        "V2COORD",
        "D",
        False,
        "V coordinate of baseline BC of the triangle",
        "m",
    ),
    (
        "STA_INDEX",
        "3I",
        False,
        "Station numbers contributing to the data",
        None,
    ),
    ("FLAG", "NWLL", False, "Flag", None),
]
OI_FLUX_KEYWORDS = [
    ("OI_REVN", False, "Revision number"),
    ("DATE-OBS", False, "UTC start date of observations"),
    ("ARRNAME", False, "Identifies corresponding OI_ARRAY"),
    ("INSNAME", False, "Identifies corresponding OI_WAVELENGTH table"),
    ("CORRNAME", True, "Identifies corresponding OI_CORR table"),
    ("FOV", True, "Area on sky over which flux is integrated (arcsec)"),
    ("FOVTYPE", True, "Model for FOV: FWHM or RADIUS"),
    ("CALSTAT", True, "C: Spectrum is calibrated, U: uncalibrated"),
]
OI_FLUX_COLUMNS = [
    (
        "TARGET_ID",
        "I",
        False,
        "Target number as index into OI_TARGET table",
        "m",
    ),
    ("MJD", "D", False, "Modified Julian day", "day"),
    ("INT_TIME", "D", False, "Integration time", "s"),
    ("FLUXDATA", "NWLD", False, "Flux in units of TUNITn", "external"),
    ("FLUXERR", "NWLD", False, "Corresponding flux error", "external"),
    (
        "CORRINDX_FLUXDATA",
        "J",
        True,
        "Index into correlation matrix for 1st FLUXDATA element",
        None,
    ),
    ("STA_INDEX", "I", True, "Station number contributing to the data", None),
    ("FLAG", "NWLL", False, "Flag", None),
]


def _pickle(self, f: str | Path | io.BufferedWriter, **kwargs) -> None:
    """Save the pickled representation of the object into an already
    open file or opens a file from a string or `pathlib.Path`.
    """
    file = open(f, "wb") if isinstance(f, (str, Path)) else f
    pickle.dump(self.serialize(), file)
    if isinstance(f, (str, Path)):
        file.close()


def _unpickle(cls, f: Path | io.TextIOWrapper, **kwargs) -> object:
    """Read the pickled representation from an open file
    or reads a string or `pathlib.Path` into a file to return the reconstituted object.
    """
    file = open(f, "rb") if isinstance(f, (str, Path)) else f
    restored_object = cls.deserialize(pickle.load(file))
    if isinstance(f, (str, Path)):
        file.close()

    return restored_object


def _serialize_function(fn: Callable) -> str:
    """Serialises a function JSON safe."""
    payload = f"{fn.__module__}|{fn.__qualname__}".encode()
    return f"fn::{base64.urlsafe_b64encode(payload).decode()}"


def _deserialize_function(value: str) -> Callable:
    """Deserialises a function JSON safe."""
    payload = base64.urlsafe_b64decode(value.removeprefix("fn::")).decode()
    module_name, qualname = payload.split("|", 1)
    obj = importlib.import_module(module_name)
    for part in qualname.split("."):
        obj = getattr(obj, part)

    return obj


def load_toml(toml_file: Path) -> dict[str, Any]:
    """Loads a `toml` file into a dictionary while properly converting units
    to `astropy.units` and floats to `numpy._float`.
    """
    with open(toml_file, "r") as file:
        dictionary = toml.load(file)

    for value in dictionary.values():
        value["unit"] = u.Unit(value.get("unit", ""))
        value["mini"] = np.float16(value.get("mini", "-inf"))
        value["maxi"] = np.float16(value.get("maxi", "inf"))

    return dictionary


def attach_methods(
    functions: Callable | ArrayLike | dict[str, Callable],
) -> Callable:
    """Class decorator that attaches function(s) to a class as methods."""
    if not isinstance(functions, dict):
        if not isinstance(functions, Iterable):
            functions = [functions]

        functions = {func.__name__: func for func in functions}

    def decorator(cls):
        for name, func in functions.items():
            setattr(cls, name, func)

        return cls

    return decorator


# TODO: Think of splitting into ν/λ variants to save computation time by
# avoiding divisions outside of this function
def blackbody(
    T: float | np.ndarray, nu: float | np.ndarray
) -> float | np.ndarray:
    r"""Computes Planck's law.

    Parameters
    ----------
    T: float or numpy.ndarray
        The temperature (K).
    nu : float or numpy.ndarray
        The frequency (Hz).

    Returns
    -------
    blackbody : float or np.ndarray
        The blackbody (erg / (cm² s Hz sr)).

    Notes
    -----
    Planck's law is defined as

    .. math::

        B_\nu(\nu,T)=\frac{2h\nu^3}{c^2}\frac{1}{\exp\left(\frac{h\nu}{k_\text{B}T}\right)-1}

    This custom variant is implemented for a more efficient computation (i.e. to
    avoid the overhead of similar implementations like the astropy's
    `astropy.modeling.physical_models.BlackBody`).
    """
    return (
        2
        * const.cgs.h
        * nu**3
        / const.cgs.c**2
        / (np.exp(const.cgs.h * nu / (const.cgs.kB * T)) - 1)
    )


def spectral_index(
    data, T: float | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    r"""Computes the spectral index of a star dependent on its effective
    temperature.

    Parameters
    ----------
    data : oimData.oimData
        The data.
    T : float or numpy.ndarray
        The effective temperature of the star (K).

    Returns
    -------
    wl : numpy.ndarray
        The (unique) wavelengths of the data (m).
    spectral_index : numpy.ndarray
        The spectral index.

    Notes
    -----
    The spectral index is defined as

    .. math:: \alpha(\nu,T)=\frac{\partial\log B_\nu(T)}{\partial\log\nu}

    """
    wl = np.unique(
        np.hstack([item for sublist in data.struct_wl for item in sublist])
    )
    nu = const.c / wl
    return wl, np.gradient(np.log(blackbody(T, nu)), np.log(nu))


def pad_image(image: np.ndarray, padfact: float | None = None) -> np.ndarray:
    """Pads an image with additional zeros for Fourier transform.

    Parameters
    ----------
    image : numpy.ndarray
        The image of shape (nt, nwl, nx, ny) to be padded.
    padfact : float, optional
        The factor by which to pad (e.g. `factor=2` with original shape `4` yields `16`).
        Defaults to the `oimOptions.ft.padding`, which itself by default is `4`.

    Results
    -------
    padded_image : numpy.ndarray
        The padded image.
    """

    if padfact is None:
        padfact = oimOptions.ft.padding

    im0 = np.sum(image, axis=(0, 1))
    dimy = im0.shape[0]
    dimx = im0.shape[1]

    im0x = np.sum(im0, axis=1)
    im0y = np.sum(im0, axis=1)

    s0x = np.trim_zeros(im0x).size
    s0y = np.trim_zeros(im0y).size

    min_sizex = s0x * padfact
    min_sizey = s0y * padfact

    min_pow2x = 2 ** (min_sizex - 1).bit_length()
    min_pow2y = 2 ** (min_sizey - 1).bit_length()

    # HACK: If Image has zeros around it already then this does not work -> Rework
    if min_pow2x < dimx:
        return image

    padx = (min_pow2x - dimx) // 2
    pady = (min_pow2y - dimy) // 2

    return np.pad(
        image,
        ((0, 0), (0, 0), (padx, padx), (pady, pady)),
        "constant",
        constant_values=0,
    )


def getDataTypeIsAnalysisComplex(dataType: str) -> bool:
    """Returns if the given dataType is complex."""
    try:
        return _oimDataAnalysisInComplex[_oimDataType.index(dataType)]
    except ValueError:
        raise TypeError(f"{dataType} not a valid OIFITS2 datatype")


def getDataArrname(dataType: str) -> str:
    """Returns the dataArrname for a given datatype."""
    try:
        return _oimDataTypeArr[_oimDataType.index(dataType)]
    except ValueError:
        raise TypeError(f"{dataType} not a valid OIFITS2 datatype")


def getDataType(dataArrname: str) -> list[str]:
    """Returns the datatype for a given dataArrname."""
    return [
        datatypei
        for datatypei, arrnamei in zip(_oimDataType, _oimDataTypeArr)
        if arrnamei == dataArrname
    ]


def getDataTypeError(dataArrname: str) -> list[str]:
    """Returns the error datatype for a given dataArrname."""
    return [
        datatypei
        for datatypei, arrnamei in zip(_oimDataTypeErr, _oimDataTypeArr)
        if arrnamei == dataArrname
    ]


def getBaselineName(
    data: str | Path | fits.HDUList,
    arr: str = "OI_VIS2",
    length: bool = False,
    angle: bool = False,
    extver: int | None = None,
    squeeze: bool = True,
) -> list[str]:
    """Gets the baseline names (i.e., the telescope names connected with a dash)
    in an extension/table of a OIFITS file.

    Defaults to reading the `"OI_VIS2"` array/table.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    hduname: str, optional
        The fits array/table name. Defaults to `"OI_VIS2"`.
    length: bool, optional
        Adds baseline length to the result. Defaults to `False`.
    angle: bool, optional
        Adds baseline position angle (deg) to the result. Defaults to `False`.
    extver: int, optional
        The extension/table version. Defaults to `None`.
    squeeze: bool, optional
        If `True` and only one extension/table is found, the result is squeezed.
        Defaults to `True`.

    Returns
    -------
    names : list of str
        The array containing the baseline names (or triplet) and, optionally,
        the baseline length and orientation.
    """
    already_open = isinstance(data, fits.HDUList)
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    extnames = np.array([di.name for di in data])

    idx_arr = np.where(extnames == "OI_ARRAY")[0]
    arrnames = np.array([data[i].header["ARRNAME"] for i in idx_arr])

    if extver is not None:
        data_arrnames = [data[arr, extver].header["ARRNAME"]]
    else:
        idx = np.where(extnames == arr)[0]
        data_arrnames = [data[i].header["ARRNAME"] for i in idx]

    names = []
    for idata, data_arrname in enumerate(data_arrnames):
        iarr = idx_arr[np.where(arrnames == data_arrname)[0][0]]
        stanames = data[iarr].data["STA_NAME"]
        staindexes = data[iarr].data["STA_INDEX"]

        staidx = data[idx[idata]].data["STA_INDEX"]
        if arr == "OI_FLUX":
            staidx = staidx[:, None]
        shape = np.shape(staidx)

        namei = []
        if (length or angle) and arr not in ["OI_T3", "OI_FLUX"]:
            u = data[idx[idata]].data["UCOORD"]
            v = data[idx[idata]].data["VCOORD"]
            B, PA = np.hypot(u, v), np.rad2deg(np.arctan2(u, v))

        for i in range(shape[0]):
            namej = ""
            for j in range(shape[1]):
                namej += stanames[np.where(staindexes == staidx[i, j])[0]][0]
                if j < shape[1] - 1:
                    namej += "-"

            if arr != "OI_T3":
                if length:
                    namej += f" {B[i]:.0f}m"

                if angle:
                    namej += f" {PA[i]:.0f}$^o$"
            namei.append(namej)
        names.append(namei)

    if not already_open:
        data.close()

    return names[0] if (squeeze and len(names) == 1) else names


# TODO : Add support for multiple extensions/tables
def getConfigName(
    data: str | Path | fits.HDUList,
    arr: str = "OI_VIS2",
    extver: int | None = None,
    squeeze: bool = True,
) -> list[str]:
    """Gets the configuration names in an extension/table of a OIFITS file.

    Defaults to reading the `"OI_VIS2"` extension/table.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    hduname: str, optional
        The fits extension/table name. Defaults to `"OI_VIS2"`.
    extver: int, optional
        The extension/table version. Defaults to `None`.
    squeeze: bool, optional
        If `True` and only one extension/table is found, the result is squeezed.
        Defaults to `True`.

    Results
    -------
    names : list of str
        The array containing the configuration names.
    """
    already_open = isinstance(data, fits.HDUList)
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    extnames = np.array([di.name for di in data])

    idx_arr = np.where(extnames == "OI_ARRAY")[0]
    arrnames = np.array([data[i].header["ARRNAME"] for i in idx_arr])

    if extver is not None:
        data_arrnames = [data[arr, extver].header["ARRNAME"]]
    else:
        idx = np.where(extnames == arr)[0]
        data_arrnames = [data[i].header["ARRNAME"] for i in idx]

    names = []
    for idata, data_arrname in enumerate(data_arrnames):
        iarr = idx_arr[np.where(arrnames == data_arrname)[0][0]]

        stanames = data[iarr].data["STA_NAME"]
        staindexes = data[iarr].data["STA_INDEX"]

        staidx = np.unique(data[idx[idata]].data["STA_INDEX"].flatten())
        s = staidx.size
        namei = ""
        for i in range(s):
            namei += stanames[np.where(staindexes == staidx[i])[0]][0]
            if i < s - 1:
                namei += "-"
        names.append(namei)

    if not already_open:
        data.close()

    return names[0] if (squeeze and len(names) == 1) else names


def getBaselineLengthAndPA(
    data: str | Path | fits.HDUList,
    arr: str = "OI_VIS2",
    extver: int | None = None,
    squeeze: bool = True,
    returnUV: bool = False,
    T3Max: bool = False,
    showFlagged: bool = True,
) -> tuple[np.ndarray, ...]:
    """Return a tuple (B, PA) of the baseline lengths and orientation
    (position angles) from a fits extension/table within an opened oifits file.

    Defaults to reading the `"OI_VIS2"` extension/table.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    arr: str, optional
        The fits extension/table name. Defaults to `"OI_VIS2"`.
    extver: int, optional
        The extension/table version. Defaults to `None`.
    squeeze: bool, optional
        If True and only one extension/table is found, the result is squeezed.
        Defaults to `True`.
    returnUV : bool, optional
        If True also return the (u,v) coordinates (m). Defaults to `False`.
    T3Max : bool, optional
        If `True` and `arr="OI_T3"` then the longest baselines of the triangles
        are returned. Defaults to `False`.
    showFlagged : bool, optional
        If `True`, takes flagged (u,v) coordinates into account. Defaults to `True`.

    Returns
    -------
    B : numpy.ndarray
        The baseline lengths.
    PA : numpy.ndarray
        The baseline orientations (deg).
    ucoord : numpy.ndarray, optional
        The u coordinate (m).
    vcoord : numpy.ndarray, optional
        The v coordinate (m).
    """
    already_open = isinstance(data, fits.HDUList)
    if isinstance(data, (str, Path)):
        data = hdul = fits.open(data)

    if extver is not None:
        data = [data[arr, extver]]
    else:
        extnames = np.array([di.name for di in data])
        idx = np.where(extnames == arr)[0]
        data = [data[i] for i in idx]

    baselines, pa, ucoord, vcoord = [], [], [], []
    for datai in data:
        if arr != "OI_T3":
            u, v = datai.data["UCOORD"], datai.data["VCOORD"]
            if not showFlagged:
                flag = datai.data["FLAG"]
                if (len(flag.shape) == 2) & (u.size == flag.shape[0]):
                    flag = np.all(flag, axis=1)
                    u = u[np.logical_not(flag)]
                    v = v[np.logical_not(flag)]

            ucoord.append(u)
            vcoord.append(v)

            baselines.append(np.hypot(u, v))
            pa.append(np.rad2deg(np.arctan2(u, v)))
        else:
            u1, v1 = datai.data["U1COORD"], datai.data["V1COORD"]
            u2, v2 = datai.data["U2COORD"], datai.data["V2COORD"]
            u3, v3 = -u1 - u2, -v1 - v2
            u123 = np.array([u1, u2, u3])
            v123 = np.array([v1, v2, v3])

            ucoord.append(u123)
            vcoord.append(v123)

            b123 = np.hypot(u123, v123)
            pa123 = np.rad2deg(np.arctan2(u123, v123))

            if T3Max:
                baselines.append(np.max(b123, axis=0))
                # TODO: What's the PAS of a triangle?
                pa.append(0 * pa123)
            else:
                baselines.append(b123)
                pa.append(np.max(pa123, axis=0))

    if squeeze and len(baselines) == 1:
        baselines, pa = baselines[0], pa[0]
        ucoord, vcoord = ucoord[0], vcoord[0]

    if not already_open:
        hdul.close()

    if returnUV:
        return baselines, pa, ucoord, vcoord
    else:
        return baselines, pa


def getSpaFreq(
    data: str | Path | fits.HDUList,
    arr: str = "OI_VIS2",
    unit: str | None = None,
    extver: int | None = None,
    squeeze: bool = True,
) -> list[np.ndarray] | np.ndarray:
    """Get the spatial dimensional frequencies.

    Defaults to reading the `"OI_VIS2"` extension/table.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    arr : str, optional
        The fits extension/table name. Defaults to `"OI_VIS2"`.
    unit : str, optional
        The unit of the spatial frequency. Defaults to `None`.
    extver : int, optional
        The extension/table version. Defaults to `None`.
    squeeze : bool, optional
        If `True` and only one extension/table is found, the result is squeezed.
        Defaults to `True`.

    Returns
    -------
    spaFreq : list of numpy.ndarray or numpy.ndarray
        The Spatial frequencies.
    """
    already_open = isinstance(data, fits.HDUList)
    if isinstance(data, (str, Path)):
        data = hdul = fits.open(data)

    baselines, _ = getBaselineLengthAndPA(data, arr, extver, squeeze=False)

    if arr == "OI_T3":
        baselines = [np.max(baseline, axis=0) for baseline in baselines]

    extnames = np.array([di.name for di in data])

    if extver is not None:
        arrays = [data[arr, extver]]
        insnames = np.array([arrays.header["INSNAME"]])
    else:
        idx = np.where(extnames == arr)[0]
        insnames = [data[i].header["INSNAME"] for i in idx]
        arrays = [data[i] for i in idx]

    idx_wlarr = np.where(extnames == "OI_WAVELENGTH")[0]
    wl_insnames = np.array([data[i].header["INSNAME"] for i in idx_wlarr])

    if unit == "cycles/mas":
        mult = u.mas.to(u.rad)
    elif unit == "cycles/arcsec":
        mult = u.arcsec.to(u.rad)
    elif unit == "Mlam":
        mult = 1 / (1e6)
    else:
        mult = 1

    spaFreq = []
    for iarr, _ in enumerate(arrays):
        iwlarr = idx_wlarr[np.where(wl_insnames == insnames[iarr])[0][0]]

        lam = data[iwlarr].data["EFF_WAVE"]
        nlam = np.size(lam)
        nB = np.size(baselines[iarr])

        spaFreqi = np.ndarray([nB, nlam])
        for iB in range(nB):
            spaFreqi[iB, :] = baselines[iarr][iB] / lam * mult

        spaFreq.append(spaFreqi)

    if not already_open:
        hdul.close()

    return spaFreq[0] if (squeeze and len(spaFreq) == 1) else spaFreq


def get2DSpaFreq(
    data: str | Path | fits.HDUList,
    arr: str = "OI_VIS2",
    unit: str | None = None,
    extver: int | None = None,
    squeeze: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Get the spatial two dimensional frequencies.

    Defaults to reading the `"OI_VIS2"` extension/table.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    arr : str, optional
        The fits extension/table name. Defaults to `"OI_VIS2"`.
    unit : str, optional
        The unit of the spatial frequency. Defaults to `None`.
    extver : int, optional
        The extension/table version. Defaults to `None`.
    squeeze : bool, optional
        If `True` and only one extension/table is found, the result is squeezed.
        Defaults to `True`.

    Returns
    -------
    2DspaFreq : tuple of numpy.ndarray
        The two-dimensional spatial frequencies.
    """
    already_open = isinstance(data, fits.HDUList)
    if isinstance(data, (str, Path)):
        data = hdul = fits.open(data)

    *_, ucoord, vcoord = getBaselineLengthAndPA(
        data, arr, extver, squeeze=False, returnUV=True
    )

    if arr == "OI_T3":
        raise TypeError("get2DSpaFreq does not accept OI_T3 extension/table")

    extnames = np.array([di.name for di in data])

    if extver is not None:
        arrays = [data[arr, extver]]
        insnames = np.array([arrays.header["INSNAME"]])
    else:
        idx = np.where(extnames == arr)[0]
        insnames = [data[i].header["INSNAME"] for i in idx]
        arrays = [data[i] for i in idx]

    idx_wlarr = np.where(extnames == "OI_WAVELENGTH")[0]
    wl_insnames = np.array([data[i].header["INSNAME"] for i in idx_wlarr])

    if unit == "cycles/mas":
        mult = u.mas.to(u.rad)
    elif unit == "cycles/arcsec":
        mult = u.arcsec.to(u.rad)
    elif unit == "Mlam":
        mult = 1 / (1e6)
    else:
        mult = 1

    spaFreqU, spaFreqV = [], []
    for iarr, _ in enumerate(arrays):
        iwlarr = idx_wlarr[np.where(wl_insnames == insnames[iarr])[0][0]]

        lam = data[iwlarr].data["EFF_WAVE"]
        nlam = np.size(lam)
        nB = np.size(ucoord[iarr])

        spaFreqUi = np.ndarray([nB, nlam])
        spaFreqVi = np.ndarray([nB, nlam])
        for iB in range(nB):
            spaFreqUi[iB, :] = ucoord[iarr][iB] / lam * mult
            spaFreqVi[iB, :] = vcoord[iarr][iB] / lam * mult
        spaFreqU.append(spaFreqUi)
        spaFreqV.append(spaFreqVi)

    if squeeze and len(spaFreqU) == 1:
        spaFreqU, spaFreqV = spaFreqU[0], spaFreqV[0]

    if not already_open:
        hdul.close()

    return spaFreqU, spaFreqV


def getWlFromOifits(
    data: str | Path | fits.HDUList,
    arr: str = "OI_VIS2",
    extver: int | None = None,
    returnBand: bool = False,
) -> tuple[np.ndarray, ...]:
    """Get the wavelength.

    Defaults to reading the `"OI_VIS2"` extension/table.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    arr : str, optional
        The fits extension/table name. Defaults to `"OI_VIS2"`.
    unit : str, optional
        The unit of the spatial frequency. Defaults to `None`.
    extver : int, optional
        The extension/table version. Defaults to `None`.
    returnBand : bool, optional
        If `True` returns the bandwith. Defaults to `False`.

    Returns
    -------
    wavelength : numpy.ndarray
        The wavelength (m).
    dwl : numpy.ndarray, optional
        The bandwith.
    """
    already_open = isinstance(data, fits.HDUList)
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    try:
        if isinstance(arr, str):
            arr = data[arr, extver]
    except:
        if extver == 1:
            arr = data[arr]
        else:
            raise TypeError(f"No extver in {arr}")

    insname = arr.header["INSNAME"]
    oiwls = np.array([di for di in data if di.name == "OI_WAVELENGTH"])
    oiwls_insname = np.array([oiwli.header["INSNAME"] for oiwli in oiwls])

    iwl = np.where(oiwls_insname == insname)[0][0]
    oiwl = oiwls[iwl]

    wavelength, bandwidth = oiwl.data["EFF_WAVE"], oiwl.data["EFF_BAND"]

    if not already_open:
        data.close()

    if returnBand:
        return wavelength, bandwidth
    else:
        return wavelength


def getWlFromFitsImageCube(
    header: fits.header.Header, outputUnit: str | None = None
) -> float:
    """Returns the wavelength law from a chromatic cube image in the fits format.

    Parameters
    ----------
    header : astropy.io.fits.header
        The header of the fits cube.
    outputUnit : astropy.unit, optional
        Converts the wavelength to passed unit. Defaults to `None`.

    Returns
    -------
    wavelength : float
        The wavelength in the given unit of the fits cube or the user-specified
        if outputUnit is set.
    """
    dwl, nwl, wl0 = header["CDELT3"], header["NAXIS3"], header["CRVAL3"]
    try:
        x0 = header["CRPIX3"]
    except:
        x0 = 0

    wl = wl0 + (np.arange(nwl) - x0) * dwl
    if outputUnit:
        if "CUNIT3" in header:
            try:
                unit0 = u.Unit(header["CUNIT3"])
            except:
                unit0 = u.m
        else:
            unit0 = u.m

        wl *= unit0.to(outputUnit)
    return wl


# TODO: Check if this is needed or if copy.deepcopy does the same.
def hdulistDeepCopy(hdulist: fits.HDUList) -> fits.HDUList:
    """Deep copy of a fits HDUList."""
    res = hdulist.copy()
    res._file = hdulist._file
    for iext, exti in enumerate(res):
        res[iext] = exti.copy()
        res[iext].header = exti.header.copy()
    return res


def _createOiTab(
    extname: str,
    keywords_def: tuple[Any, ...],
    colums_def: tuple[Any, ...],
    dataTypeFromShape: str,
    **kwargs,
) -> fits.BinTableHDU:
    """Create a OIFITS table from a dictionary of data.

    Parameters
    ----------
    extname : str
        The extension/table name.
    keywords_def : tuple of typing.Any
        The keywords definition.
    colums_def : tuple of typing.Any
        The columns definition.
    dataTypeFromShape : str
        The key in keyword argument to get the data type from.

    Returns
    -------
    hdu : astropy.io.fits.BinTableHDU
    """
    keys = {}

    for el in kwargs.keys():
        keys[el.upper()] = kwargs[el]

    if "DATE_OBS" in keys:
        keys["DATE-OBS"] = keys.pop("DATE_OBS")

    if "DATEOBS" in keys:
        keys["DATE-OBS"] = keys.pop("DATEOBS")

    try:
        nb = len(keys["TARGET_ID"])
    except:
        nb = 1

    shape = np.shape(keys[dataTypeFromShape])
    dim = len(shape)

    if "TIME" not in keys:
        keys["TIME"] = np.zeros(nb)

    if "OI_REVN" not in keys:
        keys["OI_REVN"] = 2

    if dim == 2:
        nwl = shape[1]
    elif dim == 1:
        if shape[0] == nb:
            nwl = 1
        else:
            nwl = shape[0]
    cols = []
    for colname, form, optional, comment, unit in colums_def:
        if colname not in keys and not optional:
            raise TypeError(f"Missing {colname} column")

        elif colname in keys:
            if "NWL" in form:
                form = f"{nwl}{form[3:]}"

            if unit == "external":
                try:
                    unit = keys["UNIT"]
                except:
                    raise TypeError(f"missing unit for {colname}")
            cols.append(
                fits.Column(
                    name=colname,
                    format=form,
                    array=np.array(keys[colname]),
                    unit=unit,
                )
            )

    hdu = fits.BinTableHDU.from_columns(cols)
    for keyword, optional, comment in keywords_def:
        if keyword not in keys and not optional:
            raise TypeError(f"Missing {keyword} keyword")

        elif keyword in keys:
            hdu.header[keyword] = (keys[keyword], comment)

    hdu.header["EXTNAME"] = extname
    if "EXTVER" in keys:
        hdu.header["EXTVER"] = (keys["EXTVER"], f"ID number of this {extname}")
    return hdu, keys


def createOiTarget(**kwargs):
    """Create a OI_TARGET table from a dictionary of data."""
    return _createOiTab(
        "OI_TARGET",
        OI_TARGET_KEYWORDS,
        OI_TARGET_COLUMNS,
        "TARGET_ID",
        **kwargs,
    )[0]


def createOiArray(**kwargs):
    """Create a OI_ARRAY table from a dictionary of data."""
    return _createOiTab(
        "OI_ARRAY", OI_ARRAY_KEYWORDS, OI_ARRAY_COLUMNS, "STA_INDEX", **kwargs
    )[0]


def createOiWavelength(**kwargs):
    """Create a OI_WAVELENGTH table from a dictionary of data."""
    return _createOiTab(
        "OI_WAVELENGTH", OI_WL_KEYWORDS, OI_WL_COLUMNS, "EFF_WAVE", **kwargs
    )[0]


def createOiVis(**kwargs):
    """Create a OI_VIS table from a dictionary of data."""
    return _createOiTab(
        "OI_VIS", OI_VIS_KEYWORDS, OI_VIS_COLUMNS, "VISAMP", **kwargs
    )[0]


def createOiVis2(**kwargs):
    """Create a OI_VIS2 table from a dictionary of data."""
    return _createOiTab(
        "OI_VIS2", OI_VIS2_KEYWORDS, OI_VIS2_COLUMNS, "VIS2DATA", **kwargs
    )[0]


def createOiT3(**kwargs):
    """Create a OI_T3 table from a dictionary of data."""
    return _createOiTab(
        "OI_T3", OI_T3_KEYWORDS, OI_T3_COLUMNS, "T3AMP", **kwargs
    )[0]


def createOiFlux(**kwargs):
    """Create a OI_FLUX table from a dictionary of data."""
    return _createOiTab(
        "OI_FLUX", OI_FLUX_KEYWORDS, OI_FLUX_COLUMNS, "FLUXDATA", **kwargs
    )[0]


def createOiTargetFromSimbad(names: str | list[str]) -> fits.BinTableHDU:
    """Create a OI_TARGET table from a dictionary of data.

    Parameters
    ----------
    names : str or list of str
        The name of the targets.

    Results
    -------
    hdu : astropy.io.fits.BinTableHDU
    """
    customSimbad = Simbad()
    customSimbad.add_votable_fields(
        "plx", "plx_error", "propermotions", "sptype", "velocity"
    )

    if type(names) == type(""):
        names = [names]

    data = customSimbad.query_objects(names)
    ntargets = len(names)

    rad = Angle(data["RA"], unit="hourangle").deg
    dec = Angle(data["DEC"], unit="deg").deg
    ra_err = (data["COO_ERR_MAJA"].data * u.mas).to_value(unit="deg")
    dec_err = (data["COO_ERR_MINA"].data * u.mas).to_value(unit="deg")
    pmra = (data["PMRA"].data * u.mas).to_value(unit="deg")
    pmdec = (data["PMDEC"].data * u.mas).to_value(unit="deg")
    pmra_err = (data["PM_ERR_MAJA"].data * u.mas).to_value(unit="deg")
    pmdec_err = (data["PM_ERR_MINA"].data * u.mas).to_value(unit="deg")
    plx_value = (data["PLX_VALUE"].data * u.mas).to_value(unit="deg")
    plx_error = (data["PLX_ERROR"].data * u.mas).to_value(unit="deg")

    return createOiTarget(
        target_id=np.arange(1, ntargets + 1),
        target=names,
        raep0=rad,
        decep0=dec,
        equinox=np.repeat(2000, ntargets),
        ra_err=ra_err,
        dec_err=dec_err,
        sysvel=np.zeros([ntargets]),
        veltyp=np.repeat("UNKNOWN", ntargets),
        veldef=np.repeat("OPTICAL", ntargets),
        pmra=pmra,
        pmdec=pmdec,
        pmra_err=pmra_err,
        pmdec_err=pmdec_err,
        parallax=plx_value,
        para_err=plx_error,
        spectyp=data["SP_TYPE"],
    )


def cutWavelengthRange(
    data: str | Path | fits.HDUList,
    wlRange: list[float] | None = None,
    addCut: list[float] = [],
) -> fits.HDUList:
    """Cut the wavelength range of an OIFITS file.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    wlRange : list of float, optional
        The wavelength range to keep. Defaults to `None`.
    addCut : list of float, optional
        Additional columns to cut. Defaults to `[]`.

    Returns
    -------
    data : fits.HDUList
        The new oifits file with the cut wavelength range.
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    extnames = np.array([data[i].name for i in range(len(data))])
    wlRange = np.array(wlRange)

    if wlRange.ndim == 1:
        wlRange = wlRange.reshape((1, len(wlRange)))

    for i in np.where(extnames == "OI_WAVELENGTH")[0]:
        insname = data[i].header["INSNAME"]

        idx_wl_cut = []
        for wlRangei in wlRange:
            idx_wl_cut.extend(
                np.where(
                    (data[i].data["EFF_WAVE"] >= wlRangei[0])
                    & (data[i].data["EFF_WAVE"] <= wlRangei[1])
                )[0]
            )

        idx_wl_cut = np.sort(idx_wl_cut).astype("int64")
        nwl_cut = len(idx_wl_cut)
        for idata, datai in enumerate(data):
            if "INSNAME" in datai.header:
                if datai.header["INSNAME"] == insname:
                    colDefs = []
                    for col in datai.columns:
                        if np.isin(col.name, _cutArr) or np.isin(
                            col.name, addCut
                        ):
                            format0 = col.format[-1]
                            shape = datai.data[col.name].shape
                            if len(shape) == 2:
                                arr = np.take(
                                    datai.data[col.name], idx_wl_cut, axis=1
                                )
                                colDefs.append(
                                    fits.Column(
                                        name=col.name,
                                        format=f"{nwl_cut}{format0}",
                                        array=arr,
                                        unit=col.unit,
                                    )
                                )
                            else:
                                arr = np.take(datai.data[col.name], idx_wl_cut)
                                colDefs.append(
                                    fits.Column(
                                        name=col.name,
                                        format=col.format,
                                        array=arr,
                                        unit=col.unit,
                                    )
                                )

                        else:
                            colDefs.append(
                                fits.Column(
                                    name=col.name,
                                    format=col.format,
                                    array=datai.data[col.name],
                                    unit=col.unit,
                                )
                            )

                    cols = fits.ColDefs(colDefs)
                    hdu = fits.BinTableHDU.from_columns(cols)
                    hdu.header = datai.header
                    hdu.update_header()
                    data[idata] = hdu
    return data


# TODO: This current method does not allow to shift baseline of the same oifits
# individually.
def shiftWavelength(
    data: str | Path | fits.HDUList, shift: float, verbose: bool = False
) -> None:
    """Shift the wavelength of an OIFITS file.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    shift : float
        The wavelength shift to apply.
    verbose : bool, optional
        If `True` prints the tables index. Defaults to `False`.
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    extnames = np.array([data[i].name for i in range(len(data))])
    wl_idx = np.where(extnames == "OI_WAVELENGTH")[0]
    for i in wl_idx:
        if verbose:
            print(f"OI_WAVELENGTH table at index {i}")

        insname = data[i].header["INSNAME"]
        if verbose:
            print(f"INSNAME = {insname}")
        data[i].data["EFF_WAVE"] += shift


# TODO: For phases smoothing should be done in complex plane
def spectralSmoothing(
    data: str | Path | fits.HDUList,
    kernel_size: float,
    cols2Smooth: str | list[str] = "all",
    normalizeError: bool = True,
) -> None:
    """Smooth the spectral data of an OIFITS file.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    kernel_size : float
        The kernel size.
    cols2Smooth : str or list of str, optional
        The columns to smooth. Defaults to `"all"`.
    normalizeError : bool, optional
        If `True` normalize the error. Defaults to `True`.
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    tableToSmooth = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]
    kernel = np.ones(kernel_size) / kernel_size
    if not isinstance(cols2Smooth, list):
        cols2Smooth = [cols2Smooth]

    if "all" in cols2Smooth:
        cols2Smooth = [
            "VIS2DATA",
            "VIS2ERR",
            "VISAMP",
            "VISAMPERR",
            "VISPHI",
            "VISPHIERR",
            "T3AMP",
            "T3AMPERR",
            "T3PHI",
            "T3PHIERR",
            "FLUXDATA",
            "FLUXDATAERR",
        ]

        circular = [
            False,
            False,
            False,
            False,
            True,
            True,
            False,
            False,
            True,
            True,
            False,
            False,
        ]

    for i, _ in enumerate(data, start=1):
        try:
            if data[i].name in tableToSmooth:
                cols = data[i].data.columns
                for coli in cols:
                    if coli.name in cols2Smooth:
                        dims = np.shape(data[i].data[coli.name])
                        iscirc = circular[cols2Smooth.index(coli.name)]
                        if iscirc:
                            if len(dims) == 2:
                                nB = dims[0]
                                for iB in range(nB):
                                    datai = np.exp(
                                        complex(0, 1)
                                        * data[i].data[coli.name][iB, :]
                                    )
                                    datai_real = np.real(datai)
                                    datai_imag = np.imagl(datai)

                                    conv_real = np.convolve(
                                        datai_real, kernel, "same"
                                    )
                                    conv_imag = np.convolve(
                                        datai_imag, kernel, "same"
                                    )
                                    data_conv = np.complex(
                                        conv_real, conv_imag
                                    )

                                    data[i].data[coli.name][iB, :] = data_conv

                            else:
                                data[i].data[coli.name] = np.convolve(
                                    data[i].data[coli.name], kernel, "same"
                                )
                        else:
                            if len(dims) == 2:
                                nB = dims[0]
                                for iB in range(nB):
                                    data[i].data[coli.name][iB, :] = (
                                        np.convolve(
                                            data[i].data[coli.name][iB, :],
                                            kernel,
                                            "same",
                                        )
                                    )
                            else:
                                data[i].data[coli.name] = np.convolve(
                                    data[i].data[coli.name], kernel, "same"
                                )
                        if normalizeError and coli.name in _oimDataTypeErr:
                            data[i].data[coli.name] /= np.sqrt(kernel_size)
        except:
            pass


# TODO: Properly implement error propagation for circular case.
def _intpBinning(
    array: np.ndarray,
    binMasks: ArrayLike,
    binEdgeValues: ArrayLike,
    values: ArrayLike | None = None,
    nSpecChannels: float = 1.0,
    kind: str = "mean",
    **kwargs,
) -> np.ndarray:
    """Bins the given array  in the mask.

    Parameters
    ----------
    array : numpy.ndarray
        The array to be binned.
    binMasks : array_like
        Masks of the grid underlying the array that splits it into individual bins.
    binEdgeValues : array_like
        Edge points (i.e. values) of the bin windows (i.e. masks). These are included
        to make sure the edge points of the bins are always included. Without this they
        might not be the case for arbitrary values of the bin grid.
    values : array_like, optional
        If this parameters is passed the function will assume that the `array`
        provided are errors to this `values` parameter. Defaults to `None`.
    nSpecChannels : float, optional
        The number of spectral channels determined by the spectral resolution.
        Will be used to calculate the divisor within the error propagation.
        Defaults to `1.0`.

        .. math:: divisor = bin_elements / spectralChannels

    kind : bool, optional
        Specifies the kind of binning as a string. The string has to be one of
        `"mean"`, `"median"`, `"circular"`. Defaults to `"mean"`.

    Returns
    -------
    interpolation_binned_array : numpy.ndarray
        The interpolated and binned array.
    """
    bin_func = np.mean

    if kind == "median":
        bin_func = np.median
    # TODO: Properly reimplement this.
    elif kind == "circular":
        bin_func = np.mean
        # bin_func = partial(circmean, low=-180, high=180)

    res = []
    for (lower, upper), mask in zip(binEdgeValues, binMasks):
        val = np.array([lower, *array[mask], upper])
        if values is not None:
            divisor = val.size / nSpecChannels
            res.append(np.sqrt(np.sum(val**2)) / divisor)
        else:
            res.append(bin_func(val))

    return np.array(res)


# TODO: Change this to masked arrays somehow to make it even more robust?
def _interpolateBinHDU(
    hdu: fits.BinTableHDU,
    binGrid: np.ndarray,
    binMasks: np.ndarray,
    binEdgeValues: ArrayLike,
    grid: ArrayLike,
    exception: list[str] = [],
    nSpecChannels: float = 1.0,
    **kwargs,
) -> fits.BinTableHDU:
    """Bin an HDU via interpolation.

    Parameters
    ----------
    hdu : astropy.io.fits.BinTableHDU
        The HDU to re-bin.
    binGrid : numpy.ndarray
        The grid that is to be achieved/binned to.
    binMasks : numpy.ndarray
        Masks of the grid underlying the array that splits it into individual bins.
    binEdgeValues : array_like
        Edge points (i.e. values) of the bin windows (i.e. masks). These are included
        to make sure the edge points of the bins are always included. Without this they
        might not be the case for arbitrary values of the bin grid.
    grid : array_like
        The non-binned grid.
    exception : list of str
        The exceptions (i.e. table(s) that are not to be binned).
    nSpecChannels : float, optional
        The number of spectral channels determined by the spectral resolution.
        Will be used to calculate the divisor within the error propagation.
        Defaults to `1.0`.

        .. math:: divisor = bin_elements / spectralChannels

    Returns
    -------
    newhdu : astropy.io.fits.BinTableHDU
        The rebinned HDU.
    """
    if not np.all(np.diff(grid) > 0):
        indices = np.argsort(grid)
        grid = grid[indices]
    else:
        indices = grid.astype(bool)

    cols, new_cols = hdu.data.columns, []
    if 2 in [len(np.shape(hdu.data[coli.name])) for coli in cols]:
        for col in cols:
            if not np.isin(col.name, _cutArr):
                new_cols.append(col)
                continue

            kind = "circular" if "PHI" in col.name else "mean"
            newformat, shape = col.format, hdu.data[col.name].shape
            if len(shape) == 2 and (col.name not in exception):
                bini = []
                for jB in range(shape[0]):
                    if "ERR" in col.name:
                        if any(x in col.name for x in ["VIS2", "FLUX"]):
                            val_key = col.name.replace("ERR", "DATA")
                        else:
                            val_key = col.name.replace("ERR", "")

                        values = hdu.data[val_key][jB][indices]
                    else:
                        values = None

                    array = hdu.data[col.name][jB][indices]
                    binEdgeValues = np.interp(binEdgeValues, grid, array)
                    binij = _intpBinning(
                        array,
                        binMasks[:, indices],
                        binEdgeValues,
                        values,
                        nSpecChannels,
                        kind,
                        **kwargs,
                    )
                    bini.append(binij)

                bini = np.array(bini)
                newformat = f"{binGrid.shape[0]}{col.format[-1]}"
            else:
                bini = hdu.data[col.name]

            if col.name == "FLAG" and kwargs.get("resetFlags", True):
                bini = np.full(bini.shape, False)

            new_cols.append(
                fits.Column(
                    name=col.name,
                    array=bini,
                    unit=col.unit,
                    format=newformat,
                )
            )
    else:
        for col in cols:
            if not np.isin(col.name, _cutArr):
                new_cols.append(col)
                continue

            if col.name == "EFF_WAVE":
                bini = binGrid
            elif col.name == "EFF_BAND":
                diff = np.diff(binGrid)
                bini = (
                    np.full(binGrid.shape, 0)
                    if diff.size == 0
                    else np.append(diff, diff[0])
                )
            else:
                if "ERR" in col.name:
                    if any(x in col.name for x in ["VIS2", "FLUX"]):
                        val_key = col.name.replace("ERR", "DATA")
                    else:
                        val_key = col.name.replace("ERR", "")

                    values = hdu.data[val_key][indices]
                else:
                    values = None

                kind = "circular" if "PHI" in col.name else "mean"
                array = hdu.data[col.name][indices]
                binEdgeValues = np.interp(binEdgeValues, grid, array)
                bini = _intpBinning(
                    array,
                    binMasks[:, indices],
                    binEdgeValues,
                    values,
                    nSpecChannels,
                    kind,
                    **kwargs,
                )

            if col.name == "FLAG" and kwargs.get("resetFlags", True):
                bini = np.full(bini.shape, False)

            new_cols.append(
                fits.Column(
                    name=col.name,
                    array=bini,
                    unit=col.unit,
                    format=col.format,
                )
            )

    newhdu = fits.BinTableHDU.from_columns(fits.ColDefs(new_cols))
    newhdu.header = hdu.header
    newhdu.update_header()
    return newhdu


def intpBinWavelength(
    data: str | Path | fits.HDUList, binGrid: ArrayLike, **kwargs
) -> None:
    """Bin the wavelength of an OIFITS file to a specified binGrid.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fitsHDUList`.
    binGrid : numpy.ndarray
        The grid that is to be achieved/binned to.
    binWindow : array_like, optional
        The bin windows that correspond to the binGrid elements.
        If None, computes the bin windows from the distance between two
        elements in the binGrid. Defaults to None.
    resetFlags : bool, optional
        If True, resets all flags to "False" after binning. Defaults to True.
    averageError : bool, optional
        If True, forgoes the error propagation and simply averages the errors
        for each bin. Defaults to False.
    nSpecChannels : float, optional
        The number of spectral channels determined by the spectral resolution.
        Will be used to calculate the divisor within the error propagation.
        Defaults to "1.0".

        .. math:: divisor = bin_elements / spectralChannels
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    extnames = np.array([data[i].name for i in range(len(data))])
    for i in np.where(extnames == "OI_WAVELENGTH")[0]:
        insname, wl = data[i].header["INSNAME"], data[i].data["EFF_WAVE"]
        if kwargs["binWindow"] is None:
            window = np.full(binGrid.shape, np.diff(binGrid)[0])
        else:
            window = kwargs["binWindow"]

        binEdgeValues = np.array(
            [
                (bin - win / 2, bin + win / 2)
                for win, bin in zip(window, binGrid)
            ]
        )
        binMasks = np.array(
            [(wl > lower) & (wl < upper) for lower, upper in binEdgeValues]
        )

        for idata, datai in enumerate(data):
            if datai.header.get("INSNAME", None) == insname:
                data[idata] = _interpolateBinHDU(
                    datai,
                    binGrid,
                    binMasks,
                    binEdgeValues,
                    wl,
                    exception=["STA_INDEX"],
                    **kwargs,
                )


# TODO: Properly implement the circular case
def _rebin(
    array: ArrayLike,
    binSize: int,
    kind: str = "mean",
) -> np.ndarray:
    """Rebins an array.

    Parameters
    ----------
    array : numpy.ndarray
        The array to rebin.
    binSize : int, optional
        The bin size.
    kind : bool, optional
        Specifies the kind of binning as a string. The string has to be one of
        "mean", "median", "circular". Defaults to "mean".

    Returns
    -------
    rebinned_array : numpy.ndarray
        The re-binned array.
    """
    newsize = (array.shape[0] // int(binSize)) * binSize
    array = array[:newsize]

    shape = (array.shape[0] // binSize, binSize)
    array = array.reshape(shape)

    bin_func = np.mean
    if kind == "median":
        bin_func = np.median
    elif kind == "circular":
        # TODO: Should this be circmean?
        bin_func = np.mean

    return bin_func(array, axis=-1)


def _rebinHDU(
    hdu: fits.BinTableHDU,
    binSize: int,
    exception: list[str] = [],
) -> fits.BinTableHDU:
    """Rebin an HDU.

    Parameters
    ----------
    hdu : astropy.io.fits.BinTableHDU
        The HDU to rebin.
    binsize : int
        The bin size.
    exception : list of str
        The exceptions. Defaults to [].

    Returns
    -------
    newhdu : astropy.io.fits.BinTableHDU
        The rebinned HDU.
    """
    cols, newcols = hdu.data.columns, []
    if 2 in [len(np.shape(hdu.data[coli.name])) for coli in cols]:
        for coli in cols:
            newformat = coli.format
            shape = np.shape(hdu.data[coli.name])
            kind = "circular" if "PHI" in coli.name else "mean"

            if len(shape) == 2 and not (coli.name in exception):
                bini = []
                for jB in range(shape[0]):
                    binij = _rebin(hdu.data[coli.name][jB], binSize, kind=kind)
                    bini.append(binij)
                bini = np.array(bini)
                newformat = f"{shape[1]//binSize}{coli.format[-1]}"
            else:
                bini = hdu.data[coli.name]

            newcoli = fits.Column(
                name=coli.name, array=bini, unit=coli.unit, format=newformat
            )
            newcols.append(newcoli)
    else:
        for coli in cols:
            newformat = coli.format
            kind = "circular" if "PHI" in coli.name else "mean"
            bini = _rebin(hdu.data[coli.name], binSize, kind=kind)
            newcoli = fits.Column(
                name=coli.name, array=bini, unit=coli.unit, format=newformat
            )
            newcols.append(newcoli)

    newhdu = fits.BinTableHDU.from_columns(fits.ColDefs(newcols))
    newhdu.header = hdu.header
    newhdu.update_header()
    return newhdu


# TODO: For phases, binning should be done in complex plane
def binWavelength(
    data: str | Path | fits.HDUList,
    binSize: int | None = None,
    normalizeError: bool = True,
) -> None:
    """Bin the wavelength of an OIFITS file.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an oifits file or a HDUList.
    binSize : int, optional
        The bin size. Defaults to None.
    normalizeError : bool, optional
        If True normalize the error. Defaults to True.
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    to_bin = ["OI_WAVELENGTH", "OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]
    for i, _ in enumerate(data):
        if data[i].name in to_bin:
            data[i] = _rebinHDU(data[i], binSize, exception=["STA_INDEX"])
            if normalizeError:
                errname = getDataTypeError(data[i].name)
                for errnamei in errname:
                    data[i].data[errnamei] /= np.sqrt(binSize)


def oifitsFlagWithExpression(
    data,
    arr: str | list[str],
    extver0: int,
    expr: str,
    keepOldFlag: bool = False,
):
    """Flag data with an expression.

    Parameters
    ----------
    data : astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    arr : str or list of str
        The fits table/array name(s).
    extver : int
        The extension/table version.
    expr : str
        The expression to evaluate.
    keepOldFlag : bool, optional
        If True keep the old flag. Defaults to True.

    Returns
    -------
    flags : numpy.ndarray
        The flags.
    """
    if not isinstance(arr, list):
        arr = [arr]

    if arr == ["all"]:
        arr = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    for iarr in range(len(data)):
        ok = True
        if data[iarr].name in arr:
            arri = data[iarr].name
            try:
                if extver0:
                    if extver0 != data[iarr].header.get("EXTVER", 1):
                        ok = False
                else:
                    extver = data[iarr].header.get("EXTVER", 1)
            except:
                pass
        else:
            ok = False

        if ok:
            try:
                eff_wave, eff_band = getWlFromOifits(
                    data, arr=arri, extver=extver, returnBand=True
                )
                nwl = np.size(eff_wave)
                if arri != "OI_FLUX":
                    length, pa = getBaselineLengthAndPA(
                        data, arr=arri, extver=extver, T3Max=True
                    )
                    nB = np.size(length)
                else:
                    dim = data[iarr].data["FLUXDATA"].shape
                    ndim = len(dim)
                    if ndim == 2:
                        nB = dim[0]
                        length = np.ones(nB) * np.nan
                        pa = np.ones(nB) * np.nan

                EFF_WAVE = np.tile(eff_wave[None, :], (nB, 1))
                EFF_BAND = np.tile(eff_band[None, :], (nB, 1))
                LENGTH = np.tile(length[:, None], (1, nwl))

                PA = np.tile(pa[:, None], (1, nwl))
                SPAFREQ = LENGTH / EFF_WAVE

                for colname in data[arri].columns:
                    coldata = data[arri].data[colname.name]
                    s = coldata.shape

                    if len(s) == 1 and s[0] == nB:
                        coldata = np.tile(length[:, None], (1, nwl))

                    # TODO: replace globals by a dict
                    globals()[f"{colname.name}"] = coldata

                # TODO: Remove eval here as it is can be security liability
                flags = eval(expr)
                if keepOldFlag:
                    data[iarr].data["FLAG"] = np.logical_or(
                        flags, data[iarr].data["FLAG"]
                    )
                else:
                    data[iarr].data["FLAG"] = flags

            except:
                raise Warning(
                    f"oifitsFlagWithExpression: "
                    f"Couldn't resolve expression {expr} in {arri} "
                )

    return True


# TODO: Does the `baselines_to_keep` arg need to be converted to a
# list if a string is passed?
def oifitsKeepBaselines(
    data,
    arr: str | list[str],
    baselines_to_keep: str | list[str],
    extver: int | None = None,
    keepOldFlag: bool = True,
):
    """Remove all baselines except those specified by name."""
    if arr == "all" or arr == ["all"] or not arr:
        arr = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    if isinstance(arr, str) or not isinstance(arr, Iterable):
        arr = [arr]

    for arri in arr:
        try:
            baselines = getBaselineName(data, arr=arri, extver=extver)
            baselines_to_keep_ordered = []
            for Bi in baselines_to_keep:
                Bi = Bi.split("-")
                Bi.sort()
                baselines_to_keep_ordered.append("".join(Bi))

            baselines_ordered = []
            for iB, Bi in enumerate(baselines):
                Bi = Bi.split("-")
                Bi.sort()
                baselines_ordered.append("".join(Bi))

            baselines_ordered = np.array(baselines_ordered)
            baselines_to_keep_ordered = np.array(baselines_to_keep_ordered)

            idx_to_keep = []
            for Bi in baselines_to_keep_ordered:
                idx = np.where(baselines_ordered == Bi)[0]
                if len(idx != 0):
                    idx_to_keep.extend(idx)

            for iB, Bi in enumerate(baselines):
                if not (iB in idx_to_keep):
                    data[arri, extver].data["FLAG"][iB, :] = True

        except:
            pass


# TODO: Does the `baselines_to_remove` arg need to be converted to a
# list if a string is passed?
def oifitsRemoveBaselines(
    data,
    arr: str | list[str],
    baselines_to_remove: list[str],
    extver: int | None = None,
    keepOldFlag: bool = True,
):
    """Remove all baselines specified by name."""
    if arr == "all" or arr == ["all"] or not arr:
        arr = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    if isinstance(arr, str) or not isinstance(arr, Iterable):
        arr = [arr]

    for arri in arr:
        try:
            baselines = getBaselineName(data, arr=arri, extver=extver)
            baselines_to_remove_ordered = []
            for Bi in baselines_to_remove:
                Bi = Bi.split("-")
                Bi.sort()
                baselines_to_remove_ordered.append("".join(Bi))

            baselines_ordered = []
            for iB, Bi in enumerate(baselines):
                Bi = Bi.split("-")
                Bi.sort()
                baselines_ordered.append("".join(Bi))

            baselines_ordered = np.array(baselines_ordered)
            baselines_to_remove_ordered = np.array(baselines_to_remove_ordered)

            idx_to_remove = []
            for Bi in baselines_to_remove_ordered:
                idx = np.where(baselines_ordered == Bi)[0]
                if len(idx != 0):
                    idx_to_remove.extend(idx)

            for iB, Bi in enumerate(baselines):
                if iB in idx_to_remove:
                    data[arri, extver].data["FLAG"][iB, :] = True

        except:
            pass


# TODO: Does the `telescopes_to_keep` arg need to be converted to a
# list if a string is passed?
def oifitsKeepTelescopes(
    data,
    arr: str | list[str],
    telescopes_to_keep: str | list[str],
    extver: int | None = None,
    keepOldFlag: bool = True,
):
    """Remove all telescopes except those specified by name."""
    if arr == "all" or arr == ["all"] or not arr:
        arr = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    if isinstance(arr, str) or not isinstance(arr, Iterable):
        arr = [arr]

    for arri in arr:
        try:
            baselines = getBaselineName(data, arr=arri, extver=extver)

            baselines = np.array(baselines)
            telescopes_to_keep = np.array(telescopes_to_keep)

            idx_to_keep = np.where(
                [
                    set(BB.split("-")).issubset(telescopes_to_keep)
                    for BB in baselines
                ]
            )[0]

            for iB, Bi in enumerate(baselines):
                if not (iB in idx_to_keep):
                    data[arri, extver].data["FLAG"][iB, :] = True

        except:
            pass


# TODO: Does the `telescopes_to_remove` arg need to be converted to a
# list if a string is passed?
def oifitsRemoveTelescopes(
    data,
    arr: str | list[str],
    telescopes_to_remove: str | list[str],
    extver: int | None = None,
    keepOldFlag: bool = True,
):
    """Remove all telescopes specified by name."""
    if arr == "all" or arr == ["all"] or not arr:
        arr = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    if isinstance(arr, str) or not isinstance(arr, Iterable):
        arr = [arr]

    for arri in arr:
        try:
            baselines = getBaselineName(data, arr=arri, extver=extver)

            baselines = np.array(baselines)
            telescopes_to_remove = np.array(telescopes_to_remove)

            idx_to_remove = np.where(
                [
                    bool(set(BB.split("-")) & set(telescopes_to_remove))
                    for BB in baselines
                ]
            )[0]

            for iB, Bi in enumerate(baselines):
                if iB in idx_to_remove:
                    data[arri, extver].data["FLAG"][iB, :] = True

        except:
            pass


def computeDifferentialError(
    data: str | Path | fits.HDUList,
    ranges: list[list[float]] = [[0, 5]],
    excludeRange: bool = False,
    rangeType: str = "index",
    dataType: str | list[str] = "VISPHI",
    extver: list[int | None] = [None],
) -> None:
    """Compute the differential error.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    ranges : list of list of float, optional
        The ranges to compute the differential error. Defaults to `[[0, 5]]`.
    excludeRange : bool, optional
        If `True`, exclude the range. Defaults to `False`.
    rangeType : str, optional
        The range type. Defaults to `"index"`.
    dataType : str or list of str, optional
        The data type(s). Defaults to `"VISPHI"`.
    extver : list of int, optional
        The extension/table version. Defaults to `[None]`.
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    ranges = np.array(
        ranges, dtype="int64" if rangeType == "index" else "float64"
    )

    if ranges.ndim == 1:
        ranges = ranges.reshape(1, ranges.size)

    if not isinstance(dataType, list):
        dataType = [dataType]

    extnames = np.unique([getDataArrname(dti) for dti in dataType])

    for datai in data[1:]:
        if datai.name in extnames:
            if datai.header["EXTVER"] in extver or extver == [None]:

                wl = getWlFromOifits(
                    data, arr=datai.name, extver=datai.header["EXTVER"]
                )
                nwl = wl.size
                nB = np.size(datai.data["TARGET_ID"])

                idx = np.array([], dtype="int64")
                for ri in ranges:
                    if rangeType == "index":
                        idx = np.append(idx, np.arange(ri[0], ri[1]))
                    else:
                        idxi = np.where((wl >= ri[0]) & (wl <= ri[1]))[0]
                        idx = np.append(idx, idxi)

                    if excludeRange:
                        idx = np.delete(np.arange(nwl), idx)

                for dataTypei in getDataType(datai.name):
                    if dataTypei in dataType:
                        dataTypeiErr = _oimDataTypeErr[
                            _oimDataType.index(dataTypei)
                        ]
                        if getDataTypeIsAnalysisComplex(dataTypei):
                            for iB in range(nB):
                                err = circstd(
                                    datai.data[dataTypei][iB, idx],
                                    low=-180,
                                    high=180,
                                )
                                datai.data[dataTypeiErr][iB, :] = err
                        else:
                            for iB in range(nB):
                                err = np.std(datai.data[dataTypei][iB, idx])
                                datai.data[dataTypeiErr][iB, :] = err


def setMinimumError(
    data: str | Path | fits.HDUList,
    dataTypes: str | list[str],
    values: float | list[float],
    extver: int | list[int] | None = None,
) -> None:
    """Set the minimum error of a given data type to a given value.

    Parameters
    ----------
    data: str or pathlib.Path or astropy.io.fits.HDUList
        Either a path to an OIFITS file or an `astropy.io.fits.HDUList`.
    dataTypes : str or list of str
        The data types.
    values : float or list of float
        The minimum error value.
    extver : int or list of int, optional
        The extension/table version. Defaults to `None`.
    """
    if isinstance(data, (str, Path)):
        data = fits.open(data)

    if isinstance(dataTypes, str) or not isinstance(dataTypes, Iterable):
        dataTypes = [dataTypes]

    values = [values] if not isinstance(values, Iterable) else values
    extver = [extver] if not isinstance(extver, Iterable) else extver

    extnames = np.unique([getDataArrname(dti) for dti in dataTypes])
    for datai in data[1:]:
        if datai.name in extnames:
            if datai.header.get("EXTVER", 1) in extver or extver == [None]:

                for dataTypei in getDataType(datai.name):
                    if dataTypei in dataTypes:
                        dataTypeiErr = _oimDataTypeErr[
                            _oimDataType.index(dataTypei)
                        ]
                        vali = values[dataTypes.index(dataTypei)]

                        if getDataTypeIsAnalysisComplex(dataTypei):
                            # TODO: Could the astype here be changed to ``int`` or ``bool``?
                            mask = (datai.data[dataTypeiErr] < vali).astype(
                                datai.data[dataTypeiErr].dtype
                            )

                            datai.data[dataTypeiErr] = (
                                datai.data[dataTypeiErr] * (1 - mask)
                                + mask * vali
                            )
                        else:
                            vali = vali / 100
                            mask = (
                                (
                                    datai.data[dataTypeiErr]
                                    / datai.data[dataTypei]
                                )
                                < vali
                            ).astype(datai.data[dataTypeiErr].dtype)
                            datai.data[dataTypeiErr] = (
                                datai.data[dataTypeiErr] * (1 - mask)
                                + mask * vali * datai.data[dataTypei]
                            )


def _listFeatures(
    baseClass,
    featureToTextFunction,
    details: bool = False,
    save2csv: bool = False,
    header=None,
):

    list_features = []
    for obj in oim.__dict__:
        try:
            if issubclass(oim.__dict__[obj], baseClass):
                list_features.append(obj)
        except:
            pass

    table = []
    if header:
        table.append(header)

    names = []
    for cname in list_features:
        try:
            ti = featureToTextFunction(cname)
            table.append(ti)
            names.append(cname)
        except:
            print(cname)

    if details:
        # TODO: Change this to "with" generator
        if save2csv:
            f = open(save2csv, "w", newline="")
            w = csv.writer(f, delimiter="|")
            w.writerows(table)
            f.close()
        return table
    else:
        return names


def listComponents(
    details: bool = False, save2csv: bool = False, componentType: str = "all"
):

    def _componentToTextfunction(cname: str):
        c = oim.__dict__[cname]()
        p = c.params
        tab = [cname, c.name]
        txt = ""
        for pname in p:
            txt += ":abbr:`"
            txt += pname
            txt += "("
            txt += p[pname].description
            txt += ")`, "
        txt = txt[:-2]
        tab.append(txt)
        return tab

    header = ["Component Name", "Short description", "Parameters"]
    if componentType.lower() == "all":
        class0 = oim.oimComponent
    elif componentType.lower() == "fourier":
        class0 = oim.oimComponentFourier
    elif componentType.lower() == "image":
        class0 = oim.oimComponentImage
    elif componentType.lower() == "radial":
        class0 = oim.oimComponentRadialProfile

    return _listFeatures(
        class0, _componentToTextfunction, details, save2csv, header=header
    )


def listDataFilters(details: bool = False, save2csv: bool = False):
    header = ["Filter Name", "Short description", "Class Keywords"]

    def _datFilterToTextfunction(cname: str):
        filt = oim.__dict__[cname]()
        tab = [cname, filt.description]
        txt = ""
        for pname in filt.params:
            txt += pname + ", "
        tab.append(txt[:-2])
        return tab

    res = _listFeatures(
        oim.oimDataFilterComponent,
        _datFilterToTextfunction,
        details,
        save2csv,
        header=header,
    )

    return res


def listFitters(details: bool = False, save2csv: bool = False):
    header = ["Fitter Name", "Description"]

    def _fitterToTextfunction(cname: str):
        fit = oim.__dict__[cname](None, None)
        tab = [cname, getattr(fit, "description", " - ")]

        """
        txt = ""
        try:
            for pname in fit.params:
                txt += pname + ", "
            txt = txt[:-2]
        except:
            pass
        tab.append(txt)
        """
        return tab

    res = _listFeatures(
        oim.oimFitter, _fitterToTextfunction, details, save2csv, header=header
    )

    return res


def listParamInterpolators(details: bool = False, save2csv: bool = False):
    header = ["Class Name", "oimInterp macro", "Description", "parameters"]
    p = oim.oimParam()
    interp_name = list(oim._interpolators.values())
    interp_macro = list(oim._interpolators.keys())

    def _interpToTextfunction(cname):
        interp = oim.__dict__[cname](p)
        try:
            macro = interp_macro[interp_name.index(cname)]
        except:
            macro = " - "

        tab = [cname, macro, getattr(interp, "interpdescription", " - ")]
        txt = ""
        try:
            for pname in interp.interparams:
                txt += pname + ", "
            txt = txt[:-2]
        except:
            pass
        tab.append(txt)

        return tab

    res = _listFeatures(
        oim.oimParamInterpolator,
        _interpToTextfunction,
        details,
        save2csv,
        header=header,
    )

    return res


# %%
class _terminalColor:

    BACKGROUND_BLACK = "\033[40m"
    BACKGROUND_RED = "\033[41m"
    BACKGROUND_GREEN = "\033[42m"
    BACKGROUND_YELLOW = "\033[43m"  # orange on some systems
    BACKGROUND_BLUE = "\033[44m"
    BACKGROUND_MAGENTA = "\033[45m"
    BACKGROUND_CYAN = "\033[46m"
    BACKGROUND_LIGHT_GRAY = "\third-party033[47m"
    BACKGROUND_DARK_GRAY = "\033[100m"
    BACKGROUND_BRIGHT_RED = "\033[101m"
    BACKGROUND_BRIGHT_GREEN = "\033[102m"
    BACKGROUND_BRIGHT_YELLOW = "\033[103m"
    BACKGROUND_BRIGHT_BLUE = "\033[104m"
    BACKGROUND_BRIGHT_MAGENTA = "\033[105m"
    BACKGROUND_BRIGHT_CYAN = "\033[106m"
    BACKGROUND_WHITE = "\033[107m"

    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"  # orange on some systems
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    LIGHT_GRAY = "\033[37m"
    DARK_GRAY = "\033[90m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    WHITE = "\033[97m"

    def __init__(self):
        pass

    def getCode(self, text: str):
        return getattr(self, text.upper())


def colorPrint(text: str, color) -> None:
    tcol = oim.oimUtils._terminalColor()
    col_text = tcol.getCode(color)
    reset = "\033[0m"

    print(col_text + text + reset)


# %%
def oimWarning(myclass, warningName, text: str, color: str = "red") -> None:
    if myclass._firstInit and oimOptions.general.warning:
        BOLD = "\033[1m"
        colorPrint(
            BOLD + f"oimodeler {warningName} Warning ({myclass.__name__})",
            color,
        )
        colorPrint(text, color)
        myclass._firstInit = False


# %%
def oimAckWarning(myclass, text: str) -> None:
    text += (
        "\nCheck the oimodeler page for proper refrence and acknowledgment : \n"
        "https://oimodeler.readthedocs.io/en/latest/ackn.html#acknowledgment"
    )
    oimWarning(myclass, "acknowledgement", text)
