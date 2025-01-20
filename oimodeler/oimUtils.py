# -*- coding: utf-8 -*-
"""Various utilities for optical interferometry"""
from typing import Any, List, Optional, Tuple, Union

import astropy.constants as const
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.modeling import models
from scipy.stats import circmean, circstd

from .oimOptions import oimOptions
import oimodeler as oim
import csv

# TODO: Should this (global variables) be moved into a configuration file?
_oimDataType = ["VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI", "FLUXDATA"]
_oimDataTypeErr = ["VIS2ERR", "VISAMPERR",
                   "VISPHIERR", "T3AMPERR", "T3PHIERR", "FLUXERR"]
_oimDataTypeArr = ["OI_VIS2", "OI_VIS", "OI_VIS", "OI_T3", "OI_T3", "OI_FLUX"]

_oimDataAnalysisInComplex = [False, False, True, False, True, False]

_cutArr = ['EFF_WAVE', 'EFF_BAND', 'VIS2DATA', 'VIS2ERR', 'FLAG', 'VISAMP', 'VISAMPERR',
           'VISPHI', 'VISPHIERR', 'T3AMP', 'T3AMPERR', 'T3PHI', 'T3PHIERR',
           'FLUXDATA', 'FLUXERR', 'FLAG']

OI_TARGET_KEYWORDS = [
    ("OI_REVN", False, "Revision number")
    ]

OI_TARGET_COLUMNS = [
        ("TARGET_ID", "I",   False, "Index number. Must be >=1", None),
        ("TARGET",    "16A", False, "Target name", None),
        ("RAEP0",     "D",   False, "RA at mean EQUINOX ", "deg"),
        ("DECEP0",    "D",   False, "Dec at mean EQUINOX", "deg"),
        ("EQUINOX",   "E",   False, "Equinox", None),
        ("RA_ERR",    "D",   False, "Error in RA", "deg"),
        ("DEC_ERR",   "D",   False, "Error in Dec", "deg"),
        ("SYSVEL",    "D",   False, "Systemic radial velocity", "m/s"),
        ("VELTYP",   "8A",   False, "Reference for radial velocity:LSR, GEOCENTR...", None),
        ("VELDEF",   "8A",   False, "Definition of radial velocity:(OPTICAL,RADIO)", None),
        ("PMRA",     "D",    False, "Proper motion in RA", "deg/yr"),
        ("PMDEC",    "D",    False, "Proper motion in Dec", "deg/yr"),
        ("PMRA_ERR", "D",    False, "Error of proper motion in RA", "deg/yr"),
        ("PMDEC_ERR","D",    False, "Error of proper motion in Dec", "deg/yr"),
        ("PARALLAX", "E",    False, "Parallax", "deg"),
        ("PARA_ERR", "E",    False, "Error in parallax ", "deg"),
        ("SPECTYP",  "16A",  False, "Spectral type", "deg"),
        ("CATEGORY", "3A",   True,  "CAL or SCI", None)
        ]

OI_ARRAY_KEYWORDS = [
        ("OI_REVN", False, "Revision number"),
        ("ARRNAME", False, "A Array name, for cross-referencing"),
        ("FRAME",   False, "A Coordinate frame"),
        ("ARRAYX",  False, "Array center x coordinates (m)"),
        ("ARRAYY",  False, "Array center y coordinates (m)"),
        ("ARRAYZ",  False, "Array center z coordinates (m)")
        ]

OI_ARRAY_COLUMNS = [
        ("TEL_NAME", "16A", False, " Telescope name", None),
        ("STA_NAME", "16A", False, "Station name", None),
        ("STA_INDEX", "I", False, "Station number. Must be >=1", None),
        ("DIAMETER", "E", False, "Element diameter", "m"),
        ("STAXYZ", "3D", False, "Station coordinates w.r.t. array center", "m"),
        ("FOV", "D", False, " Photometric field of view", "arcsec"),
        ("FOVTYPE", "6A", False, "Model for FOV: FWHM or RADIUS", None)
        ]

OI_WL_KEYWORDS = [
        ("OI_REVN",  False, "Revision number"),
        ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table")
        ]

OI_WL_COLUMNS = [
        ("EFF_WAVE", "E", False, "Effective wavelength of channel", "m"),
        ("EFF_BAND", "E", False, "Effective bandpass of channel", "m")
        ]

OI_VIS2_KEYWORDS = [
        ("OI_REVN",  False, "Revision number"),
        ("DATE-OBS", False, "UTC start date of observations"),
        ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
        ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
        ("CORRNAME", True,  "Identifies corresponding OI_CORR table")
        ]

OI_VIS2_COLUMNS = [
        ("TARGET_ID",         "I",    False, "Target number as index into OI_TARGET table", "m"),
        ("TIME",              "D",    False, "Zero. For backwards compatibility", None),
        ("MJD",               "D",    False, "Modified Julian day", "day"),
        ("INT_TIME",          "D",    False, "Integration time", None),
        ("VIS2DATA",          "NWLD", False, "Squared Visibility", None),
        ("VIS2ERR",           "NWLD", False, "Error in Squared Visibility", None),
        ("CORRINDX_VIS2DATA", "J",    True,  "Index into correlation matrix for 1st VIS2DATA element", None),
        ("UCOORD",            "D",    False, "U coordinate of the data", "m"),
        ("VCOORD",            "D",    False, "V coordinate of the data", "m"),
        ("STA_INDEX",         "2I",   False, "Station numbers contributing to the data", None),
        ("FLAG",              "NWLL", False, "Flag", None)
        ]

OI_VIS_KEYWORDS = [
        ("OI_REVN",  False, "Revision number"),
        ("DATE-OBS", False, "UTC start date of observations"),
        ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
        ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
        ("CORRNAME", True,  "Identifies corresponding OI_CORR table"),
        ("AMPTYP",   True,  "absolute, differential, correlated flux"),
        ("PHITYP",   True,  "absolute, differential"),
        ("AMPORDER", True,  "Polynomial fit order for differential chromatic amplitudes"),
        ("PHIORDER", True,  "Polynomial fit order for differential chromatic phases")
        ]

# TODO: Implement RIV VISREFMAP ...
OI_VIS_COLUMNS = [
        ("TARGET_ID",       "I",    False, "Target number as index into OI_TARGET table", "m"),
        ("TIME",            "D",    False, "Zero. For backwards compatibility", None),
        ("MJD",             "D",    False, "Modified Julian day", "day"),
        ("INT_TIME",        "D",    False, "Integration time", "s"),
        ("VISAMP",          "NWLD", False, "Visibility amplitude", None),
        ("VISAMPERR",       "NWLD", False, "Error in visibility amplitude", None),
        ("CORRINDX_VISAMP", "J",    True,  "Index into correlation matrix for 1st VISAMP element", None),
        ("VISPHI",          "NWLD", False, "Visibility phase", "deg"),
        ("VISPHIERR",       "NWLD", False, "Error in visibility Phase", "deg"),
        ("CORRINDX_VISPHI", "J",    True,  "Index into correlation matrix for 1st VISPHI element", None),
        ("UCOORD",          "D",    False, "U coordinate of the data", "m"),
        ("VCOORD",          "D",    False, "V coordinate of the data", "m"),
        ("STA_INDEX",       "2I",   False, "Station numbers contributing to the data", None),
        ("FLAG",            "NWLL", False, "Flag", None)]

OI_T3_KEYWORDS = [
        ("OI_REVN",  False, "Revision number"),
        ("DATE-OBS", False, "UTC start date of observations"),
        ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
        ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
        ("CORRNAME", True,  "Identifies corresponding OI_CORR table")
        ]

OI_T3_COLUMNS = [
        ("TARGET_ID",      "I",    False, "Target number as index into OI_TARGET table", "m"),
        ("TIME",           "D",    False, "Zero. For backwards compatibility", None),
        ("MJD",            "D",    False, "Modified Julian day", "day"),
        ("INT_TIME",       "D",    False, "Integration time", "s"),
        ("T3AMP",          "NWLD", False, "Triple product amplitude", None),
        ("T3AMPERR",       "NWLD", False, "Error in triple product amplitude", None),
        ("CORRINDX_T3AMP", "J",    True,  "Index into correlation matrix for 1st T3AMP element", None),
        ("T3PHI",          "NWLD", False, "Triple Product Phase", "deg"),
        ("T3PHIERR",       "NWLD", False, "Error in Triple Product Phase", "deg"),
        ("CORRINDX_T3PHI", "J",    True,  "Index into correlation matrix for 1st T3PHI element", None),
        ("U1COORD",        "D",    False, "U coordinate of baseline AB of the triangle", "m"),
        ("V1COORD",        "D",    False, "V coordinate of baseline AB of the triangle", "m"),
        ("U2COORD",        "D",    False, "U coordinate of baseline BC of the triangle", "m"),
        ("V2COORD",        "D",    False, "V coordinate of baseline BC of the triangle", "m"),
        ("STA_INDEX",      "3I",   False, "Station numbers contributing to the data", None),
        ("FLAG",           "NWLL", False, "Flag", None)
        ]

OI_FLUX_KEYWORDS = [
        ("OI_REVN",  False, "Revision number"),
        ("DATE-OBS", False, "UTC start date of observations"),
        ("ARRNAME",  False, "Identifies corresponding OI_ARRAY"),
        ("INSNAME",  False, "Identifies corresponding OI_WAVELENGTH table"),
        ("CORRNAME", True,  "Identifies corresponding OI_CORR table"),
        ("FOV",      True,  "Area on sky over which flux is integrated (arcsec)"),
        ("FOVTYPE",  True,  "Model for FOV: FWHM or RADIUS"),
        ("CALSTAT",  True,  "C: Spectrum is calibrated, U: uncalibrated")
        ]

OI_FLUX_COLUMNS = [
        ("TARGET_ID",         "I",    False, "Target number as index into OI_TARGET table", "m"),
        ("MJD",               "D",    False, "Modified Julian day", "day"),
        ("INT_TIME",          "D",    False, "Integration time", "s"),
        ("FLUXDATA",          "NWLD", False, "Flux in units of TUNITn", "external"),
        ("FLUXERR",           "NWLD", False, "Corresponding flux error", "external"),
        ("CORRINDX_FLUXDATA", "J",    True,  "Index into correlation matrix for 1st FLUXDATA element", None),
        ("STA_INDEX",         "I",   True, "Station number contributing to the data", None),
        ("FLAG",              "NWLL", False, "Flag", None)
        ]


def getDataTypeIsAnalysisComplex(dataType: str) -> str:
    """Get the data type and if it is complex."""
    try:
        return _oimDataAnalysisInComplex[_oimDataType.index(dataType)]
    except Exception:
        raise TypeError(f"{dataType} not a valid OIFITS2 datatype")


def getDataArrname(dataType: str) -> str:
    """Returns the dataArrname for a given datatype."""
    try:
        return _oimDataTypeArr[_oimDataType.index(dataType)]
    except Exception:
        raise TypeError(f"{dataType} not a valid OIFITS2 datatype")


def getDataType(dataArrname: str) -> List[str]:
    """Returns the datatype for a given dataArrname."""
    return [datatypei for datatypei, arrnamei
            in zip(_oimDataType, _oimDataTypeArr)
            if arrnamei == dataArrname]


def getDataTypeError(dataArrname: str) -> List[str]:
    """Returns the error datatype for a given dataArrname."""
    return [datatypei for datatypei, arrnamei
            in zip(_oimDataTypeErr, _oimDataTypeArr)
            if arrnamei == dataArrname]


def blackbody(temperature: u.K, wavelength: u.m = None,
              frequency: u.Hz = None) -> np.ndarray:
    """Planck's law with output in cgs.

    Parameters
    ----------
    temperature : astrophy.units.K
        The temperature [K].
    wavelength : astrpy.units.m, optional
        The wavelength [cm].
    frequency : astropy.units.Hz, optional
        The frequency [Hz].

    Returns
    -------
    blackbody : np.ndarray
        The blackbody [erg/(cm² s Hz sr)].
    """
    if wavelength is not None:
        wavelength = u.Quantity(wavelength, u.m)
        frequency = (const.c.cgs/(wavelength).to(u.cm))
    else:
        frequency = u.Quantity(frequency, u.Hz)

    temperature = u.Quantity(temperature, u.K)
    return ((2*const.h.cgs*frequency**3/const.c.cgs**2)\
            * (1/(np.exp(const.h.cgs*frequency/(const.k_B.cgs*temperature))-1))).value


def compute_intensity(wavelengths: u.um, temperature: u.K,
                      pixel_size: Optional[float] = None) -> np.ndarray:
    """Calculates the blackbody_profile via Planck's law and the
    emissivity_factor for a given wavelength, temperature- and
    dust surface density profile.

    Parameters
    ----------
    wavelengths : astropy.units.um
        Wavelength value(s).
    temp_profile : astropy.units.K
        Temperature profile.
    pixSize: float, optional
        The pixel size [rad].

    Returns
    -------
    intensity : numpy.ndarray
        Intensity per pixel.
    """
    plancks_law = models.BlackBody(temperature=temperature*u.K)
    spectral_profile = []
    pixel_size *= u.rad
    for wavelength in wavelengths*u.m:
        spectral_radiance = plancks_law(wavelength).to(
                u.erg/(u.cm**2*u.Hz*u.s*u.rad**2))
        spectral_profile.append((spectral_radiance*pixel_size**2).to(u.Jy).value)
    return np.array(spectral_profile)


def compute_photometric_slope(
        data: np.ndarray, temperature: u.K) -> np.ndarray:
    """Computes the photometric slope of the data from
    the effective temperature of the star.

    Parameters
    ----------
    data : oimData.oimData
        The observed data.
    temperature : astropy.units.K
        The effective temperature of the star.

    Returns
    -------
    photometric_slope : numpy.ndarray
    """
    temperature = u.Quantity(temperature, u.K)
    struct_wl = [item for sublist in data.struct_wl for item in sublist]
    wl = (np.unique(np.hstack(struct_wl))*u.m).to(u.um)
    nu = (const.c/wl.to(u.m)).to(u.Hz)
    blackbody = models.BlackBody(temperature)(nu)
    return wl.to(u.m).value, np.gradient(np.log(blackbody.value), np.log(nu.value))


def pad_image(image: np.ndarray) -> np.ndarray:
    """Pads an image with additional zeros for Fourier transform.

    Parameters
    ----------
    image : numpy.ndarray
        The image to be padded.

    Results
    -------
    padded_image : numpy.ndarray
        The padded image.
    """
    im0 = np.sum(image, axis=(0, 1))
    dimy = im0.shape[0]
    dimx = im0.shape[1]

    im0x = np.sum(im0, axis=1)
    im0y = np.sum(im0, axis=1)

    s0x = np.trim_zeros(im0x).size
    s0y = np.trim_zeros(im0y).size

    min_sizex = s0x*oimOptions.ft.padding
    min_sizey = s0y*oimOptions.ft.padding

    min_pow2x = 2**(min_sizex - 1).bit_length()
    min_pow2y = 2**(min_sizey - 1).bit_length()

    # HACK: If Image has zeros around it already then this does not work -> Rework
    if min_pow2x < dimx:
        return image

    padx = (min_pow2x-dimx)//2
    pady = (min_pow2y-dimy)//2

    return np.pad(image, ((0, 0), (0, 0), (padx, padx), (pady, pady)),
                  'constant', constant_values=0)


def get_next_power_of_two(number: Union[int, float]) -> int:
    """Returns the next power of two for an integer or float input.

    Parameters
    ----------
    number : int or float
        An input number.

    Returns
    -------
    closest_power_of_two : int
        The, to the input, closest power of two.
    """
    return int(2**np.ceil(np.log2(number)))


def convert_angle_to_distance(radius: Union[float, np.ndarray],
                              distance: float,
                              rvalue: Optional[bool] = False) -> Union[float, np.ndarray]:
    """Converts an angle from milliarcseconds to meters.

    Parameters
    ----------
    radius : float or numpy.ndarray
        The radius of the object around the star [mas].
    distance : float
        The star's distance to the observer [pc].
    rvalue : bool, optional
        If toggled, returns the value witout u.else returns
        an astropy.u.Quantity object. The default is False.

    Returns
    -------
    radius : float or numpy.ndarray
        The radius of the object around the star [m].

    Notes
    -----
    The formula for the angular diameter small angle approximation is

    .. math:: \\delta = \\frac{d}{D}

    where d is the distance from the star and D is the distance from the star
    to the observer and ..math::`\\delta` is the angular diameter.
    """
    radius = ((radius*u.mas).to(u.arcsec).value*distance*u.pc).to(u.m)
    return radius if rvalue else radius.value


def convert_distance_to_angle(radius: Union[float, np.ndarray],
                              distance: float,
                              rvalue: Optional[bool] = False) -> Union[float, np.ndarray]:
    """Converts a distance from meters to an angle in milliarcseconds.

    Parameters
    ----------
    radius : float or numpy.ndarray
        The radius of the object around the star [au].
    distance : float
        The star's distance to the observer [pc].
    rvalue : bool, optional
        If toggled, returns the value witout u.else returns
        an astropy.u.Quantity object. The default is False.

    Returns
    -------
    radius : float or numpy.ndarray
        The radius of the object around the star [m].

    Notes
    -----
    The formula for the angular diameter small angle approximation is

    .. math:: \\delta = \\frac{d}{D}

    where d is the distance from the star and D is the distance from the star
    to the observer and ..math::`\\delta` is the angular diameter.
    """
    radius = (((radius*u.au).to(u.m)/(distance*u.pc).to(u.m))*u.rad).to(u.mas)
    return radius if rvalue else radius.value


def getBaselineName(oifits: fits.HDUList, hduname: Optional[str] = "OI_VIS2",
                    length: Optional[bool] = False, angle: Optional[bool] = False,
                    extver: Optional[int] = None, squeeze: Optional[bool] = True) -> List[str]:
    """Gets the baseline names (i.e., telescopes names
    separated by minus sign) in an extension of a oifits file.

    By default it is reading the 'OI_VIS' extension.

    Parameters
    ----------
    oifits: astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits
    hduname: str, optional
        The fits extension name. The default is "OI_VIS2".
    length: bool, optional
        Add baseline length to the returned result. The default is False.
    angle: bool, optional
        Add baseline position angle ((in deg)à=) to the returned result.
        The default is False
    extver: int, optional
        The extension version. The default is None.
    squeeze: bool, optional
        If True and only one extension is found, the result is squeezed.

    Returns
    -------
    names : list of str
        The array containing the baseline names (or triplet) and optionally
        the baseline length and orientation.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    extnames = np.array([di.name for di in data])

    idx_arr = np.where(extnames == "OI_ARRAY")[0]
    arrnames = np.array([data[i].header["ARRNAME"] for i in idx_arr])

    if extver is not None:
        data_arrnames = [data[hduname, extver].header["ARRNAME"]]
    else:
        idx = np.where(extnames == hduname)[0]
        data_arrnames = [data[i].header["ARRNAME"] for i in idx]

    names = []
    for idata, data_arrname in enumerate(data_arrnames):
        iarr = idx_arr[np.where(arrnames == data_arrname)[0][0]]
        stanames = data[iarr].data["STA_NAME"]
        staindexes = data[iarr].data["STA_INDEX"]

        staidx = data[idx[idata]].data["STA_INDEX"]
        if hduname == "OI_FLUX":
            staidx = staidx[:, None]
        shape = np.shape(staidx)

        namei = []
        if length or angle and (hduname != "OI_T3" or hduname != "OI_FLUX"):
            u = data[idx[idata]].data["UCOORD"]
            v = data[idx[idata]].data["VCOORD"]
            B, PA = np.hypot(u, v), np.rad2deg(np.arctan2(u, v))

        for i in range(shape[0]):
            namej = ""
            for j in range(shape[1]):
                namej += stanames[np.where(staindexes == staidx[i, j])[0]][0]
                if j < shape[1]-1:
                    namej += "-"

            if hduname != "OI_T3":
                if length:
                    namej += f" {B[i]:.0f}m"

                if angle:
                    namej += f" {PA[i]:.0f}$^o$"
            namei.append(namej)
        names.append(namei)

    return names[0] if (squeeze and len(names) == 1) else names


# TODO : Add support for multiple extensions
def getConfigName(oifits: fits.HDUList, hduname: Optional[str] = "OI_VIS2",
                  extver: Optional[int] = None, squeeze: Optional[bool] = True) -> List[str]:
    """Gets the configuration names in an extension of a oifits file.

    Parameters
    ----------
    oifits: astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits
    hduname: str, optional
        The fits extension name. The default is "OI_VIS2".
    extver: int, optional
        The extension version. The default is None.
    squeeze: bool, optional
        If True and only one extension is found, the result is squeezed.

    Results
    -------
    names : list of str
        The array containing the configuration names.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    extnames = np.array([di.name for di in data])

    idx_arr = np.where(extnames == "OI_ARRAY")[0]
    arrnames = np.array([data[i].header["ARRNAME"] for i in idx_arr])

    if extver is not None:
        data_arrnames = [data[hduname, extver].header["ARRNAME"]]
    else:
        idx = np.where(extnames == hduname)[0]
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
            if i < s-1:
                namei += "-"
        names.append(namei)

    return names[0] if (squeeze and len(names) == 1) else names


def getBaselineLengthAndPA(oifits: fits.HDUList, arr: Optional[str] = "OI_VIS2",
                           extver: Optional[int] = None, squeeze: Optional[bool] = True,
                           returnUV: Optional[bool] = False, T3Max: Optional[bool] = False) -> Tuple[np.ndarray]:
    """Return a tuple (B, PA) of the baseline lengths and orientation
    (position angles) from a fits extension within an opened oifits file.

    By default it is reading the OI_VIS extension.

    Parameters
    ----------
    oifits: astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr: str, optional
        The fits extension name. The default is "OI_VIS2".
    returnUV : bool, optional
        If True also return the u,v coordinates in m
        the default is False

    Returns
    -------
    B : numpy.ndarray
        The array containing the baselines length.
    PA : numpy.ndarray
        The array containing the baselines orientation (in deg).
    ucoord : numpy.ndarray
        The array containing the u coordinate (in m)(optional)
    ucoord : numpy.ndarray
        The array containing the u coordinate (in m)(optional)
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

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
            ucoord.append(u)
            vcoord.append(v)

            baselines.append(np.hypot(u, v))
            pa.append(np.rad2deg(np.arctan2(u, v)))
        else:
            u1, v1 = datai.data["U1COORD"], datai.data["V1COORD"]
            u2, v2 = datai.data["U2COORD"], datai.data["V2COORD"]
            u3, v3 = -u1-u2, -v1-v2
            u123 = np.array([u1, u2, u3])
            v123 = np.array([v1, v2, v3])

            ucoord.append(u123)
            vcoord.append(v123)

            # TODO: Check if this give the same as adding up individually
            # TODO: Shouldn't the positional angles be the other way around? -> Convention?
            b123 = np.hypot(u123, v123)
            pa123 = np.rad2deg(np.arctan2(u123, v123))

            if T3Max:
                baselines.append(np.max(b123,axis=0))
                pa.append(0*pa123) # TODO: what's the PAS of a triangle?
            else:
                baselines.append(b123)
                pa.append(np.max(pa123,axis=0))


    if squeeze and len(baselines) == 1:
        baselines, pa = baselines[0], pa[0]
        ucoord, vcoord = ucoord[0], vcoord[0]

    if returnUV:
        return baselines, pa, ucoord, vcoord
    else:
        return baselines, pa


def get2DSpaFreq(oifits: fits.HDUList, arr: Optional[str] = "OI_VIS2",
                 unit: Optional[str] = None, extver: Optional[int] = None,
                 squeeze: Optional[bool] = True) -> Tuple[np.ndarray]:
    """Get the spatial two dimensional frequencies.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr : str, optional
        The fits extension name. The default is "OI_VIS2".
    unit : str, optional
        The unit of the spatial frequency. The default is None.
    extver : int, optional
        The extension version. The default is None.
    squeeze : bool, optional
        If True and only one extension is found, the result is squeezed.

    Returns
    -------
    spaFreq : tuple of numpy.ndarray
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    _, _, ucoord, vcoord = getBaselineLengthAndPA(
            data, arr, extver, squeeze=False, returnUV=True)

    if arr == "OI_T3":
        raise TypeError("get2DSpaFreq does not accept OI_T3 extension")

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
        mult = 1/(1e6)
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
            spaFreqUi[iB, :] = ucoord[iarr][iB]/lam*mult
            spaFreqVi[iB, :] = vcoord[iarr][iB]/lam*mult
        spaFreqU.append(spaFreqUi)
        spaFreqV.append(spaFreqVi)

    if squeeze and len(spaFreqU) == 1:
        spaFreqU, spaFreqV = spaFreqU[0], spaFreqV[0]

    return spaFreqU,spaFreqV


def getSpaFreq(oifits: fits.HDUList, arr: Optional[str] = "OI_VIS2",
               unit: Optional[str] = None, extver: Optional[int] = None,
               squeeze: Optional[bool] = True) -> Tuple[np.ndarray]:
    """Get the spatial dimensional frequencies.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr : str, optional
        The fits extension name. The default is "OI_VIS2".
    unit : str, optional
        The unit of the spatial frequency. The default is None.
    extver : int, optional
        The extension version. The default is None.
    squeeze : bool, optional
        If True and only one extension is found, the result is squeezed.

    Returns
    -------
    spaFreq : tuple of numpy.ndarray
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

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
        mult = 1/(1e6)
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
            spaFreqi[iB, :] = baselines[iarr][iB]/lam*mult

        spaFreq.append(spaFreqi)

    return spaFreq[0] if (squeeze and len(spaFreq) == 1) else spaFreq


def getWlFromOifits(oifits: fits.HDUList, arr: Optional[str] = "OI_VIS2",
                    extver: Optional[int] = None,
                    returnBand: Optional[bool] = False) -> Tuple[np.ndarray]:
    """Get the wavelength

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr : str, optional
        The fits extension name. The default is "OI_VIS2".
    unit : str, optional
        The unit of the spatial frequency. The default is None.
    extver : int, optional
        The extension version. The default is None.
    returnBand : bool, optional
        If True return the bandwith. The default is False.

    Returns
    -------
    wavelength : numpy.ndarray
        The wavelength.
    dwl : numpy.ndarray
        The bandwith.
    """
    try:
        if isinstance(arr, str):
            arr = oifits[arr, extver]
    except:
        if extver == 1:
            arr = oifits[arr]
        else:
            TypeError(f"No extver in {arr}")

    insname = arr.header["INSNAME"]
    oiwls = np.array([di for di in oifits if di.name == "OI_WAVELENGTH"])
    oiwls_insname = np.array([oiwli.header["INSNAME"] for oiwli in oiwls])

    iwl = np.where(oiwls_insname == insname)[0][0]
    oiwl = oiwls[iwl]

    wavelength  = oiwl.data["EFF_WAVE"]
    if returnBand:
        return wavelength, oiwl.data["EFF_BAND"] 
    else:
        return wavelength


def hdulistDeepCopy(hdulist: fits.HDUList) -> fits.HDUList:
    """Deep copy of a fits HDUList."""
    res = hdulist.copy()
    res._file = hdulist._file
    for iext, exti in enumerate(res):
        res[iext] = exti.copy()
        res[iext].header = exti.header.copy()
    return res


def cutWavelengthRange(oifits: fits.HDUList, wlRange: Optional[List[float]] = None,
                       addCut: Optional[List[float]] = []) -> fits.HDUList:
    """Cut the wavelength range of an oifits file.
    
    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    wlRange : list of float, optional
        The wavelength range to keep. The default is None.
    addCut : list of float, optional
        Additional columns to cut. The default is [].

    Returns
    -------
    data : fits.HDUList
        The new oifits file with the cut wavelength range.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits
    extnames = np.array([data[i].name for i in range(len(data))])

    oiwl_idx = np.where(extnames == "OI_WAVELENGTH")[0]
    wlRange = np.array(wlRange)

    if wlRange.ndim == 1:
        wlRange = wlRange.reshape((1, len(wlRange)))

    for i in oiwl_idx:
        insname = data[i].header["INSNAME"]

        idx_wl_cut = []
        for wlRangei in wlRange:
            idx_wl_cut.extend(np.where((data[i].data["EFF_WAVE"] >= wlRangei[0]) &
                                       (data[i].data["EFF_WAVE"] <= wlRangei[1]))[0])

        idx_wl_cut = np.sort(idx_wl_cut).astype("int64")
        nwl_cut = len(idx_wl_cut)
        for idata, datai in enumerate(data):
            if "INSNAME" in datai.header:
                if datai.header["INSNAME"] == insname:
                    colDefs = []
                    for col in datai.columns:
                        if np.isin(col.name, _cutArr) or np.isin(col.name, addCut):
                            format0 = col.format[-1]
                            shape = datai.data[col.name].shape
                            if len(shape) == 2:
                                arr = np.take(
                                    datai.data[col.name], idx_wl_cut, axis=1)
                                colDefs.append(fits.Column(
                                    name=col.name, format=f"{nwl_cut}{format0}",
                                    array=arr, unit=col.unit))
                            else:
                                arr = np.take(datai.data[col.name], idx_wl_cut)
                                colDefs.append(fits.Column(
                                    name=col.name, format=col.format,
                                    array=arr, unit=col.unit))

                        else:
                            colDefs.append(fits.Column(
                                name=col.name, format=col.format,
                                array=datai.data[col.name], unit=col.unit))

                    cols = fits.ColDefs(colDefs)
                    hdu = fits.BinTableHDU.from_columns(cols)
                    hdu.header = datai.header
                    hdu.update_header()
                    data[idata] = hdu
    return data


def getWlFromFitsImageCube(header: fits.header.Header,
                           outputUnit: Optional[str] = None) -> float:
    """Returns the wavelength law from a chromatic cube image in the fits format.

    Parameters
    ----------
    header : astropy.io.fits.header
        The header of the fits cube.
    outputUnit : astropy.unit, optional
        If set convert the result to the proper unit. The default is None.

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

    wl = wl0+(np.arange(nwl)-x0)*dwl
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


def _createOiTab(extname: str, keywords_def: Tuple[Any], colums_def: Tuple[Any],
                 dataTypeFromShape: str, **kwargs) -> fits.BinTableHDU:
    """Create a OIFITS table from a dictionary of data.

    Parameters
    ----------
    extname : str
        The extension name.
    keywords_def : tuple of any
        The keywords definition.
    colums_def : tuple of any
        The columns definition.
    dataTypeFromShape : str
        The key in kwargs to get the data type from.

    Returns
    -------
    hdu : astropy.io.fits.BinTableHDU
    """
    keys = {}

    for el in kwargs.keys():
        keys[el.upper()] = kwargs[el]

    if "DATE_OBS"  in keys:
        keys["DATE-OBS"] = keys.pop("DATE_OBS")

    if "DATEOBS"  in keys:
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
            cols.append(fits.Column(name=colname, format=form,
                                    array=np.array(keys[colname]), unit=unit))

    hdu = fits.BinTableHDU.from_columns(cols)
    for keyword, optional, comment in keywords_def:
        if keyword not in keys and not optional:
            raise TypeError(f"Missing {keyword} keyword")

        elif keyword in keys:
            hdu.header[keyword]=(keys[keyword], comment)

    hdu.header["EXTNAME"] = extname
    if "EXTVER" in keys:
        hdu.header["EXTVER"] = (keys["EXTVER"], f"ID number of this {extname}")
    return hdu, keys


def createOiTarget(**kwargs):
    """Create a OI_TARGET table from a dictionary of data."""
    return _createOiTab("OI_TARGET", OI_TARGET_KEYWORDS,
                        OI_TARGET_COLUMNS, "TARGET_ID",**kwargs)[0]


def createOiArray(**kwargs):
    """Create a OI_ARRAY table from a dictionary of data."""
    return _createOiTab("OI_ARRAY", OI_ARRAY_KEYWORDS, 
                        OI_ARRAY_COLUMNS, "STA_INDEX", **kwargs)[0]


def createOiWavelength(**kwargs):
    """Create a OI_WAVELENGTH table from a dictionary of data."""
    return _createOiTab("OI_WAVELENGTH", OI_WL_KEYWORDS,
                        OI_WL_COLUMNS, "EFF_WAVE", **kwargs)[0]


def createOiVis(**kwargs):
    """Create a OI_VIS table from a dictionary of data."""
    return _createOiTab("OI_VIS", OI_VIS_KEYWORDS,
                        OI_VIS_COLUMNS, "VISAMP", **kwargs)[0]


def createOiVis2(**kwargs):
    """Create a OI_VIS2 table from a dictionary of data."""
    return _createOiTab("OI_VIS2", OI_VIS2_KEYWORDS,
                        OI_VIS2_COLUMNS, "VIS2DATA", **kwargs)[0]


def createOiT3(**kwargs):
    """Create a OI_T3 table from a dictionary of data."""
    return _createOiTab("OI_T3", OI_T3_KEYWORDS,
                        OI_T3_COLUMNS, "T3AMP", **kwargs)[0]


def createOiFlux(**kwargs):
    """Create a OI_FLUX table from a dictionary of data."""
    return _createOiTab("OI_FLUX", OI_FLUX_KEYWORDS,
                        OI_FLUX_COLUMNS, "FLUXDATA", **kwargs)[0]


def createOiTargetFromSimbad(names: Union[str, List[str]]) -> fits.BinTableHDU:
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
            "plx", "plx_error", "propermotions", "sptype", "velocity")

    if type(names) == type(""):
        names = [names]

    data = customSimbad.query_objects(names)
    ntargets = len(names)

    rad = Angle(data["RA"], unit="hourangle").deg
    dec = Angle(data["DEC"], unit="deg").deg
    ra_err = (data["COO_ERR_MAJA"].data*u.mas).to_value(unit="deg")
    dec_err = (data["COO_ERR_MINA"].data*u.mas).to_value(unit="deg")
    pmra = (data["PMRA"].data*u.mas).to_value(unit="deg")
    pmdec = (data["PMDEC"].data*u.mas).to_value(unit="deg")
    pmra_err = (data["PM_ERR_MAJA"].data*u.mas).to_value(unit="deg")
    pmdec_err = (data["PM_ERR_MINA"].data*u.mas).to_value(unit="deg")
    plx_value = (data["PLX_VALUE"].data*u.mas).to_value(unit="deg")
    plx_error = (data["PLX_ERROR"].data*u.mas).to_value(unit="deg")

    return createOiTarget(
            target_id=np.arange(1, ntargets+1), target=names,
            raep0=rad, decep0=dec, equinox=np.repeat(2000, ntargets),
            ra_err=ra_err, dec_err=dec_err, sysvel=np.zeros([ntargets]),
            veltyp=np.repeat("UNKNOWN", ntargets),
            veldef=np.repeat("OPTICAL", ntargets), pmra=pmra, pmdec=pmdec,
            pmra_err=pmra_err, pmdec_err=pmdec_err, parallax=plx_value,
            para_err=plx_error, spectyp=data["SP_TYPE"])


# TODO: This current method does not allow to shift baseline of the same oifits
# individually.
def shiftWavelength(oifits: fits.HDUList, shift: float,
                    verbose: Optional[bool] = False) -> None:
    """Shift the wavelength of an oifits file.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    shift : float
        The wavelength shift to apply.
    verbose : bool, optional
        If True print the tables index. The default is False.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits
    extnames = np.array([data[i].name for i in range(len(data))])

    wl_idx = np.where(extnames == "OI_WAVELENGTH")[0]
    for i in wl_idx:
        if verbose :
            print(f"OI_WAVELENGTH table at index {i}")

        insname = data[i].header["INSNAME"]
        if verbose:
            print(f"INSNAME = {insname}")
        data[i].data["EFF_WAVE"] += shift


# TODO: For phases smoothing should be done in complex plan
def spectralSmoothing(oifits: fits.HDUList, kernel_size: float,
                      cols2Smooth: Optional[Union[str, List[str]]] = "all",
                      normalizeError: Optional[bool] = True) -> None:
    """Smooth the spectral data of an oifits file.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    kernel_size : float
        The kernel size.
    cols2Smooth : str or list of str, optional
        The columns to smooth. The default is "all".
    normalizeError : bool, optional
        If True normalize the error. The default is True.
    """
    tableToSmooth = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    kernel = np.ones(kernel_size)/kernel_size

    if not isinstance(cols2Smooth, list):
        cols2Smooth = [cols2Smooth]

    if  "all" in cols2Smooth:
        cols2Smooth = ["VIS2DATA", "VIS2ERR", "VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR",
                       "T3AMP", "T3AMPERR", "T3PHI", "T3PHIERR", "FLUXDATA", "FLUXDATAERR"]
        
        circular=[False,False,False,False,True,True,False,False,True,True,False,False]

    for i, _ in enumerate(data, start=1):
        try:
            if data[i].name in tableToSmooth:
                cols = data[i].data.columns
                for coli in cols:
                    if coli.name in cols2Smooth:
                        dims = np.shape(data[i].data[coli.name])
                        iscirc=circular[cols2Smooth.index(coli.name)]
                        if iscirc:
                            if len(dims) == 2:
                                nB = dims[0]
                                for iB in range(nB):
                                    datai=np.exp(complex(0,1)*data[i].data[coli.name][iB, :])
                                    datai_real=np.real(datai)
                                    datai_imag=np.imagl(datai)
                                    
                                    conv_real=np.convolve(datai_real, kernel, "same")
                                    conv_imag=np.convolve(datai_imag, kernel, "same")  
                                    data_conv=np.complex(conv_real,conv_imag)
                                    
                                    data[i].data[coli.name][iB, :] = data_conv
                                    
                            else:
                                  data[i].data[coli.name]=np.convolve(
                                      data[i].data[coli.name], kernel, "same")
                        else:
                            if len(dims) == 2:
                                nB = dims[0]
                                for iB in range(nB):
                                    data[i].data[coli.name][iB, :] = np.convolve(
                                        data[i].data[coli.name][iB, :], kernel, "same")
                            else:
                                  data[i].data[coli.name]=np.convolve(
                                      data[i].data[coli.name], kernel, "same")
                        if normalizeError and coli.name in _oimDataTypeErr:
                            data[i].data[coli.name] /= np.sqrt(kernel_size)
        except:
            pass


def rebin_image(image: np.ndarray,
                binning_factor: Optional[int] = None,
                rdim: Optional[bool] = False) -> np.ndarray:
    """Bins a 2D-image down according.

    The down binning is according to the binning factor
    in oimOptions.ft.binning. Only accounts for
    square images.

    Parameters
    ----------
    image : numpy.ndarray
        The image to be rebinned.
    binning_factor : int, optional
        The binning factor. The default is 0
    rdim : bool
        If toggled, returns the dimension

    Returns
    -------
    rebinned_image : numpy.ndarray
        The rebinned image.
    dimension : int, optional
        The new dimension of the image.
    """
    if binning_factor is None:
        return image, image.shape[-1] if rdim else image

    new_dim = int(image.shape[-1] * 2**-binning_factor)
    binned_shape = (new_dim, int(image.shape[-1] / new_dim),
                    new_dim, int(image.shape[-1] / new_dim))
    if len(image.shape) == 4:
        shape = (image.shape[0], image.shape[1], *binned_shape)
    else:
        shape = binned_shape
    if rdim:
        return image.reshape(shape).mean(-1).mean(-2), new_dim
    return image.reshape(shape).mean(-1).mean(-2)


def _rebin(array: np.ndarray, binsize: int,
           median: Optional[bool] = True,
           circular: Optional[bool] = False) -> np.ndarray:
    """Rebin an array.

    Parameters
    ----------
    arr : numpy.ndarray
        The array to rebin.
    binsize : int
        The bin size.
    median : bool
        If True return the median.

    Returns
    -------
    res : numpy.ndarray
        The rebinned array.
    """
    
    newsize = (array.shape[0] // int(binsize)) * binsize
    array = array [:newsize]
    shape = (array.shape[0]// binsize, binsize)

    if circular:
        if median:
            res = np.median(array.reshape(shape),axis=-1)
        else:
            res = array.reshape(shape).mean(-1)
    else:
        res = np.mean(array.reshape(shape),axis=-1)

    return res


def _rebinHdu(hdu: fits.BinTableHDU, binsize: int,
              exception: Optional[List[str]] = []) -> fits.BinTableHDU:
    """Rebin an HDU.

    Parameters
    ----------
    hdu : astropy.io.fits.BinTableHDU
        The HDU to rebin.
    binsize : int
        The bin size.
    exception : list of str
        The exceptions.

    Returns
    -------
    newhdu : astropy.io.fits.BinTableHDU
        The rebinned HDU.
    """

    cols = hdu.data.columns
    newcols = []

    if 2 in [len(np.shape(hdu.data[coli.name])) for coli in cols]:
        for coli in cols:
            if "PHI" in coli.name:
                circular=True
            else:
                circular=False
            newformat = coli.format
            shape = np.shape(hdu.data[coli.name])
            if len(shape) == 2 and not(coli.name in exception):
                bini =[]
                for jB in range(shape[0]):
                    binij = _rebin(hdu.data[coli.name][jB,:],binsize,circular=circular)
                    bini.append(binij)
                bini = np.array(bini)
                newformat = f"{shape[1]//binsize}{coli.format[-1]}"
            else:
                bini = hdu.data[coli.name]

            newcoli = fits.Column(name=coli.name, array=bini,
                                  unit=coli.unit, format=newformat)
            newcols.append(newcoli)
    else:
        for coli in cols:
            if "PHI" in coli.name:
                circular=True
            else:
                circular=False
         
            
            bini = _rebin(hdu.data[coli.name],binsize,circular=circular)
            newcoli = fits.Column(name=coli.name, array=bini,
                                  unit=coli.unit, format=coli.format)
            newcols.append(newcoli)

    newhdu = fits.BinTableHDU.from_columns(fits.ColDefs(newcols))
    newhdu.header = hdu.header
    newhdu.update_header()
    return newhdu


# TODO: For phases bining should be done in complex plan
def binWavelength(oifits: fits.HDUList, binsize: int,
                  normalizeError: Optional[bool] = True) -> None:
    """Bin the wavelength of an oifits file.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    binsize : int
        The bin size.
    normalizeError : bool
        If True normalize the error
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits
    
    tobin = ["OI_WAVELENGTH", "OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]
    for i, _ in enumerate(data):
        if data[i].name in tobin:
            data[i] = _rebinHdu(data[i], binsize, exception=["STA_INDEX"])
            if normalizeError:
                errname = getDataTypeError(data[i].name)
                for errnamei in errname:
                    data[i].data[errnamei] /= np.sqrt(binsize)


def oifitsFlagWithExpression(data,arr,extver,expr,keepOldFlag = False):
    """Flag the data with an expression.

    Parameters
    ----------
    data : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    arr : str or list of str
        The fits extension name.
    extver : int
        The extension version.
    expr : str
        The expression to evaluate.
    keepOldFlag : bool
        If True keep the old flag.

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
            if extver:
                if extver != data[iarr].header["EXTVER"]:
                    ok = False
            else:
                extver = 1 #data[iarr].header["EXTVER"]
        else:
            ok = False
            
        if ok:
            try:
                eff_wave, eff_band = getWlFromOifits(data, arr=arri,
                                                     extver=extver, returnBand=True)
                nwl = np.size(eff_wave)
                length, pa = getBaselineLengthAndPA(data, arr=arri, extver=extver,T3Max=True)
                nB = np.size(length)
    
                EFF_WAVE = np.tile(eff_wave[None, :], (nB,1))
                EFF_BAND = np.tile(eff_band[None, :], (nB,1))
                LENGTH = np.tile(length[:, None], (1, nwl))
    
                PA = np.tile(pa[:, None], (1, nwl))
     
                SPAFREQ = LENGTH/EFF_WAVE
               
    
                for colname in data[arri].columns:
                    coldata = data[arri].data[colname.name]
                    s = coldata.shape
    
                    if len(s) == 1 and s[0] == nB:
                        coldata = np.tile(length[:, None], (1, nwl))
    
                    # TODO: Remove exec here as it is can be security liability
                    exec(f"{colname.name}=coldata")
                    
      
                # TODO: Remove eval here as it is can be security liability
                flags = eval(expr)
                
                if keepOldFlag:
                    data[arri].data["FLAG"] = np.logical_or(flags, data[arri].data["FLAG"])
                else:
                    data[arri].data["FLAG"] = flags
           
            except:
                pass
        

    return True


def oifitsKeepBaselines(data,arr,baselines_to_keep,extver=None,keepOldFlag = True):
    
    
    if arr == ["all"]:
        arr = ["OI_VIS", "OI_VIS2", "OI_T3", "OI_FLUX"]

    for arri in arr:
        try:
            baselines=getBaselineName(data,hduname=arr,extver=extver)
            baselines_to_keep_ordered=[]
            for Bi in baselines_to_keep:
                Bi=Bi.split("-")
                Bi.sort()
                baselines_to_keep_ordered.append(''.join(Bi))
                
            baselines_ordered=[]
            for iB,Bi in enumerate(baselines):
                Bi=Bi.split("-")
                Bi.sort()
                baselines_ordered.append(''.join(Bi))
                
            baselines_ordered=np.array(baselines_ordered)
            baselines_to_keep_ordered=np.array(baselines_to_keep_ordered)
            
            idx_to_keep=[]
            for Bi in baselines_to_keep_ordered:
                idx=np.where(baselines_ordered==Bi)[0]
                if len(idx!=0):
                    idx_to_keep.extend(idx)
            
            
            for iB,Bi in enumerate(baselines):
                if not(iB in idx_to_keep):
                    data[arr,extver].data["FLAG"][iB,:]=True
        
        except:
            pass        
    

def computeDifferentialError(
        oifits: fits.HDUList, ranges: Optional[List[int]] = [[0, 5]],
        excludeRange: Optional[bool] = False, rangeType: Optional[str] = "index",
        dataType: Optional[str] = "VISPHI", extver: Optional[List[int]]= [None]) -> None:
    """Compute the differential error.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    ranges : list of int, optional
        The ranges to compute the differential error. The default is [[0, 5]].
    excludeRange : bool, optional
        If True exclude the range. The default is False.
    rangeType : str, optional
        The range type. The default is "index".
    dataType : str, optional
        The data type. The default is "VISPHI".
    extver : list of int, optional
        The extension version. The default is [None].
    """
    if rangeType == "index":
        dtype = "int64"
    else:
        dtype = "float64"

    ranges=np.array(ranges,dtype=dtype)

    if ranges.ndim == 1:
        ranges = ranges.reshape(1,ranges.size)

    if not isinstance(dataType, list):
        dataType = [dataType]

    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    extnames = np.unique([getDataArrname(dti) for dti in dataType])

    for datai in data[1:]:
        if datai.name in extnames:
            if datai.header["EXTVER"] in extver or extver == [None]:

                wl = getWlFromOifits(oifits, arr=datai.name,
                                     extver=datai.header["EXTVER"])
                nwl = wl.size
                nB = np.size(datai.data['TARGET_ID'])


                idx = np.array([], dtype="int64")
                for ri in ranges:
                    if rangeType == "index":
                        idx = np.append(idx,np.arange(ri[0],ri[1]))
                    else:
                        idxi = np.where( (wl>= ri[0]) & (wl <= ri[1]) )[0]
                        idx = np.append(idx,idxi)


                    if excludeRange:
                        idx = np.delete(np.arange(nwl),idx)

                for dataTypei in getDataType(datai.name ):
                    if dataTypei in dataType:
                        dataTypeiErr = _oimDataTypeErr[_oimDataType.index(dataTypei)]
                        if getDataTypeIsAnalysisComplex(dataTypei):
                            for iB in range(nB):
                                err = circstd(datai.data[dataTypei][iB, idx],
                                              low=-180, high=180)
                                datai.data[dataTypeiErr][iB, :] = err
                        else:
                            for iB in range(nB):
                                err = np.std(datai.data[dataTypei][iB, idx])
                                datai.data[dataTypeiErr][iB, :] = err


def setMinimumError(oifits: fits.HDUList, dataTypes: Union[str, List[str]],
                    values: Union[float, List[float]],
                    extver: Optional[Union[int, List[int]]] = None) -> None:
    """Set the minimum error of a given data type to a given value.

    Parameters
    ----------
    oifits : astropy.io.fits.HDUList
        An oifits file structure already opened with astropy.io.fits.
    dataTypes : str or list of str
        The data types.
    values : float or list of float
        The minimum error value.
    extver : int or list of int, optional
        The extension version. The default is None.
    """
    if type(oifits) == type(""):
        data = fits.open(oifits)
    else:
        data = oifits

    if not isinstance(dataTypes, list):
        dataTypes = [dataTypes]

    if not isinstance(values, list):
        values = [values]

    if not isinstance(extver, list):
        extver = [extver]

    extnames = np.unique([getDataArrname(dti) for dti in dataTypes])
    for datai in data[1:]:
        if datai.name in extnames:
            if datai.header["EXTVER"] in extver or extver == [None]:

                for dataTypei in getDataType(datai.name):
                    if dataTypei in dataTypes:
                        dataTypeiErr = _oimDataTypeErr[_oimDataType.index(dataTypei)]
                        vali = values[dataTypes.index(dataTypei)]

                        if getDataTypeIsAnalysisComplex(dataTypei):
                            mask = (datai.data[dataTypeiErr]<vali) \
                                    .astype(datai.data[dataTypeiErr].dtype)

                            datai.data[dataTypeiErr] = datai.data[dataTypeiErr] * \
                                         (1 - mask) +  mask * vali
                        else:
                            vali = vali / 100
                            mask = ((datai.data[dataTypeiErr]/datai.data[dataTypei]) \
                                     <vali).astype(datai.data[dataTypeiErr].dtype)
                            datai.data[dataTypeiErr] = datai.data[dataTypeiErr] * \
                                (1 - mask) + mask * vali * datai.data[dataTypei]


def _listFeatures(baseClass,featureToTextFunction,details=False,
                  save2csv=None,header=None):
    
    list_features=[]
    
    for obj in oim.__dict__:
        try:
            if issubclass(oim.__dict__[obj], baseClass):
                list_features.append(obj)  
        except:  
            pass
    
    table  = []
    if header:
        table.append(header)
    names  = []
    
    for cname in list_features:
        try:
            
            ti = featureToTextFunction(cname)
            table.append(ti)
            names.append(cname)
        except:
            print(cname)
        
    if details:
        if save2csv:
            f=open(save2csv,"w")
            w=csv.writer(f,delimiter="|")
            w.writerows(table)
            f.close()
        return table
    else:
        return names


def listComponents(details=False,save2csv=None, componentType="all"):
    
    def _componentToTextfunction(cname):
        tab=[cname]
        c = oim.__dict__[cname]()
        p = c.params
        name = c.name        
        tab.append(name)
        txt=""
        for pname in p:
            txt+=":abbr:`"
            txt+=pname
            txt+="("
            txt+=p[pname].description           
            txt+=")`, "
        txt=txt[:-2]
        tab.append(txt)
        return tab
    
    header =["Component Name","Short description","Parameters"]
    if componentType.lower() == "all":
        class0 = oim.oimComponent
    elif componentType.lower() == "fourier":
        class0 = oim.oimComponentFourier
    elif componentType.lower() == "image":
        class0 = oim.oimComponentImage
    elif componentType.lower() == "radial":
        class0 = oim.oimComponentRadialProfile
    
    return _listFeatures(class0,_componentToTextfunction,details,
                         save2csv,header=header)


def listDataFilters(details=False,save2csv=None):
    header =["Filter Name","Short description","Class Keywords"]
    
    def _datFilterToTextfunction(cname):
        filt = oim.__dict__[cname]()
        tab=[cname]
        tab.append(filt.description)
        txt=""
        for pname in filt.params:
            txt+=pname+", "
        tab.append(txt[:-2])
        return tab
    
    res = _listFeatures(oim.oimDataFilterComponent,_datFilterToTextfunction,
                        details,save2csv,header=header)
    
    return res
    
def listFitters(details=False,save2csv=None):
    header =["Fitter Name","Short description","Keywords"]
    
    def _fitterToTextfunction(cname):
        #fit = oim.__dict__[cname](None,None)
        tab=[cname]
        #tab.append(fit.description)
        #txt=""
        #for pname in filt.params:
        #    txt+=pname+", "
        #tab.append(txt[:-2])
        return tab
    
    res = _listFeatures(oim.oimFitter,_fitterToTextfunction,
                        details,save2csv,header=header)
    
    return res   


def listParamInterpolators(details=False,save2csv=None):
    header =["Class Name","oimInterp macro","Description","parameters"]
    p=oim.oimParam()
    interp_name=list(oim._interpolators.values())
    interp_macro=list(oim._interpolators.keys())
    def _interpToTextfunction(cname):
        interp = oim.__dict__[cname](p)
        tab=[cname]
        try:
            macro=interp_macro[interp_name.index("oimParamMultipleGaussianTime")]
        except:
            macro="None"
        tab.append(macro)
        #txt=""
        #for pname in filt.params:
        #    txt+=pname+", "
        #tab.append(txt[:-2])
        return tab
    
    res = _listFeatures(oim.oimParamInterpolator,_interpToTextfunction,
                        details,save2csv,header=header)
    
    return res  