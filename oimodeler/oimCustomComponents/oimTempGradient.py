from typing import Dict

import astropy.units as u
import astropy.constants as const
import numpy as np
from astropy.modeling import models

from ..oimComponent import oimComponentRadialProfile
from ..oimParam import oimParam
from .oimRadialPowerLaw import oimRadialPowerLaw


def calculate_intensity(params: Dict[str, oimParam],
                        wavelengths: u.um,
                        kappa_abs: u.cm**2/u.g,
                        temp_profile: u.K,
                        sigma_profile: u.g/u.cm**2,
                        rJy: bool = False) -> np.ndarray:
    """Calculates the blackbody_profile via Planck's law and the
    emissivity_factor for a given wavelength, temperature- and
    dust surface density profile.

    Parameters
    ----------
    params : dict with keys of str and values of oimParam
        The oimParams of the oimComponent.
    wavelengths : astropy.units.um
        Wavelength value(s).
    kappa_abs : astropy.units.cm**2/astropy.units.g
        Absorption opacity
    temp_profile : astropy.units.K
        Temperature profile.
    sigma_profile : astropy.units.g/astropy.units.cm**2
        Dust surface density profile.

    Returns
    -------
    intensity : numpy.ndarray
        Intensity per pixel.
    """
    plancks_law = models.BlackBody(temperature=temp_profile*u.K)
    wavelengths = (wavelengths*u.m).to(u.um)
    spectral_radiance = plancks_law(wavelengths).to(u.erg/(u.cm**2*u.Hz*u.s*u.mas**2))
    emissivity_factor = 1-np.exp(-sigma_profile*kappa_abs)
    if rJy:
        pix = params["pixSize"].value**2*params["pixSize"].unit**2
        return ((spectral_radiance*pix).to(u.Jy)*emissivity_factor).value
    return (spectral_radiance*emissivity_factor).value


class oimTempGradient(oimComponentRadialProfile):
    """A ring defined by a radial temperature profile in r^q and a radial dust
    surface density profile in r^p.

    Parameters
    ----------
    rin : float
        Inner radius of the disk [au].
    rout : float
        Outer radius of the disk [au].
    Tin : float
        Inner radius temperature [K].
    Mdust : float
        Mass of the dusty disk [M_sun].
    q : float
        Power-law exponent for the temperature profile.
    p : float
        Power-law exponent for the dust surface density profile.
    kappa_abs : float or oimInterp
        Dust mass absorption coefficient [cm2.g-1].
    dist : float
        Distance of the star [pc].

    Attributes
    ----------
    params : dict with keys of str and values of oimParam
        Dictionary of parameters.
    _wl : array_like
        Wavelengths.
    _t : array_like
        Times.
    _r : array_like
        Radial grid.

    Methods
    -------
    _radialProfileFunction(r, wl, t)
        Calculates a radial temperature gradient profile via a dust-surface
        density- and temperature profile.
    """
    name = "Temperature Gradient"
    shortname = "TempGrad"
    elliptic = True

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0, unit=u.au,
                                      description="Inner radius of the disk")
        self.params["rout"] = oimParam(name="rout", value=0, unit=u.au,
                                       description="Outer radius of the disk")
        self.params["q"] = oimParam(name="q", value=0, unit=u.one,
                                    description="Power-law exponent for the temperature profile")
        self.params["p"] = oimParam(name="p", value=0, unit=u.one,
                                    description="Power-law exponent for the dust surface density profile")
        self.params["Mdust"] = oimParam(name="Mdust", value=0, unit=u.M_sun,
                                        description="Mass of the dusty disk")
        self.params["Tin"] = oimParam(name="Tin", value=0,
                                      unit=u.K, free=False,
                                      description="Inner radius temperature")
        self.params["kappa_abs"] = oimParam(name="kappa_abs", value=0,
                                            unit=u.cm**2/u.g, free=False,
                                            description="Dust mass absorption coefficient")
        self.params["dist"] = oimParam(name="dist", value=0,
                                       unit=u.pc, free=False,
                                       description="Distance of the star")

        self._t = np.array([0])  # constant value <=> static model
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._eval(**kwargs)

    def _radialProfileFunction(self, r: np.ndarray,
                               wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates a radial temperature gradient profile via a dust-surface
        density- and temperature profile.

        Parameters
        ----------
        r : numpy.ndarray
            Radial grid.
        wl : numpy.ndarray
            Wavelengths.
        t : numpy.ndarray
            Times.

        Results
        -------
        radial_profile : numpy.ndarray
        """
        rin, rout = map(lambda x: self.params[x](wl, t), ["rin", "rout"])
        q, p = map(lambda x: self.params[x](wl, t), ["q", "p"])
        dist, inner_temp = map(lambda x: self.params[x](wl, t), ["dist", "Tin"])
        dust_mass = self.params["Mdust"](wl, t)*const.M_sun.value*1e3
        kappa_abs = self.params["kappa_abs"](wl, t)

        rin_mas, rout_mas = map(lambda x: 1e3*x/dist, [rin, rout])
        rin_cm, rout_cm = map(lambda x: x*self.params["rin"].unit.to(u.cm), [rin, rout])

        # NOTE: Temperature radial profile
        temp_profile = inner_temp*(r / rin_mas)**q

        # NOTE: Surface density radial profile
        if p == -2:
            sigma_in = dust_mass/(2.*np.pi*np.log(rout_cm/rin_cm)*rin_cm**2)
        else:
            f = ((rout_cm/rin_cm)**(p+2)-1)/(p+2)
            sigma_in = dust_mass/(2.*np.pi*f*rin_cm**2)
        sigma_profile = sigma_in*(r / rin_mas)**p
        spectral_density = calculate_intensity(self.params, wl,
                                               kappa_abs,
                                               temp_profile,
                                               sigma_profile)
        return np.nan_to_num(np.logical_and(r > rin_mas, r < rout_mas).astype(int)*spectral_density, nan=0)

    @property
    def _r(self):
        """Gets the radius."""
        if False:
            rout = self.params["rout"](self._wl, self._t)
        else:
            rout = self.params["rout"](self._wl, 1e99)
        rmax = 4*rout
        r = np.linspace(0, 1, self.params["dim"](self._wl, self._t))*1e+3*rmax / \
            self.params["dist"](self._wl, self._t)
        return r

    @_r.setter
    def _r(self, value):
        """Sets the radius."""
        return


class oimAsymTempGradient(oimRadialPowerLaw):
    """A ring defined by a radial temperature profile in r^q and an asymmetric
    radial dust surface density profile in r^p.

    Parameters
    ----------
    rin : float
        Inner radius of the disk [mas].
    rout : float
        Outer radius of the disk [mas].
    Tin : float
        Inner radius temperature [K].
    Mdust : float
        Mass of the dusty disk [M_sun].
    a : float
        Azimuthal modulation amplitude.
    phi : float
        Azimuthal modulation angle [deg].
    q : float
        Power-law exponent for the temperature profile.
    p : float
        Power-law exponent for the dust surface density profile.
    kappa_abs : float or oimInterp
        Dust mass absorption coefficient [cm2.g-1].
    dist : float
        Distance of the star [pc].
    fov: float
        The field of view [mas].
    pa : float
        Positional angle [deg].
    elong : float
        Elongation of the disk.
    dim : float
        Dimension of the image.

    Attributes
    ----------
    params : dict with keys of str and values of oimParam
        Dictionary of parameters.
    _t : numpy.ndarray
        Array of time values.
    _wl : numpy.ndarray
        Array of wavelength values.

    Methods
    -------
    _imageFunction(xx, yy, wl, t)
        Calculates a 2D-image from a dust-surface density- and
        temperature profile.
    """
    name = "Asymmetric Temperature Gradient"
    shortname = "AsymTempGrad"
    elliptic = True
    asymmetric = True

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0, unit=u.mas,
                                      description="Inner radius of the disk")
        self.params["rout"] = oimParam(name="rout", value=0, unit=u.mas,
                                       description="Outer radius of the disk")
        self.params["q"] = oimParam(name="q", value=0, unit=u.one,
                                    description="Power-law exponent for the temperature profile")
        self.params["p"] = oimParam(name="p", value=0, unit=u.one,
                                    description="Power-law exponent for the dust surface density profile")
        self.params["Mdust"] = oimParam(name="Mdust", value=0, unit=u.M_sun,
                                        description="Mass of the dusty disk")
        self.params["Tin"] = oimParam(name="Tin", value=0,
                                      unit=u.K, free=False,
                                      description="Inner radius temperature")
        self.params["kappa_abs"] = oimParam(name="kappa_abs", value=0,
                                            unit=u.cm**2/u.g, free=False,
                                            description="Dust mass absorption coefficient")
        self.params["dist"] = oimParam(name="dist", value=0,
                                       unit=u.pc, free=False,
                                       description="Distance of the star")

        self._t = np.array([0])  # constant value <=> static model
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._eval(**kwargs)

    # FIXME: Does pixsize need to be in rad or in mas? -> Check again, probably correct as is
    def _imageFunction(self, xx: np.ndarray, yy: np.ndarray,
                       wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates a 2D-image from a dust-surface density- and
        temperature profile.

        Parameters
        ----------
        xx : numpy.ndarray
            The x-coordinate grid
        yy : numpy.ndarray
            The y-coordinate grid
        wl : numpy.ndarray
            Wavelengths.
        t : numpy.ndarray
            Times.

        Returns
        -------
        image : numpy.ndarray
        """
        r = np.sqrt(xx**2+yy**2)
        rin, rout = map(lambda x: self.params[x](wl, t), ["rin", "rout"])
        q, p = map(lambda x: self.params[x](wl, t), ["q", "p"])
        dist, inner_temp = map(lambda x: self.params[x](wl, t), ["dist", "Tin"])
        rin_cm, rout_cm = map(lambda x: x*dist*self.params["rin"].unit.to(u.arcsec)*u.au.to(u.cm), [rin, rout])
        dust_mass = self.params["Mdust"](wl, t)*const.M_sun.value*1e3
        kappa_abs = self.params["kappa_abs"](wl, t)

        # NOTE: Temperature radial profile
        temp_profile = inner_temp*(r / rin)**q

        # NOTE: Surface density radial profile
        if p == -2:
            sigma_in = dust_mass/(2.*np.pi*np.log(rout_cm/rin_cm)*rin_cm**2)
        else:
            f = ((rout_cm/rin_cm)**(p+2)-1)/(p+2)
            sigma_in = dust_mass/(2.*np.pi*f*rin_cm**2)
        sigma_profile = sigma_in*(r / rin)**p

        if self.asymmetric:
            sigma_profile = sigma_profile*(1+self._azimuthal_modulation(xx, yy, wl, t))
        spectral_density = calculate_intensity(self.params, wl,
                                               kappa_abs, temp_profile,
                                               sigma_profile, rJy=True)
        return np.nan_to_num(np.logical_and(r > rin, r < rout).astype(int)*spectral_density, nan=0)
