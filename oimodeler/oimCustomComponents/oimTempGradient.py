from typing import Optional

import astropy.units as u
import astropy.constants as const
import numpy as np
from astropy.modeling import models

from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import oimOptions
from ..oimParam import oimParam
from ..oimUtils import convert_radial_profile_to_meter
from .oimRadialPowerLaw import oimRadialPowerLaw


def calculate_intensity(wavelengths: u.um,
                        temp_profile: u.K,
                        pixSize: Optional[float] = None) -> np.ndarray:
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
    plancks_law = models.BlackBody(temperature=temp_profile*u.K)
    wavelengths = (wavelengths*u.m).to(u.um)
    spectral_radiance = plancks_law(wavelengths).to(
        u.erg/(u.cm**2*u.Hz*u.s*u.rad**2))

    if oimOptions["ModelOutput"] == "corr_flux":
        if pixSize is None:
            raise KeyError("'pixSize' needs to be directly or indirectly"
                           " (via the 'fov'-parameter) set by the user!")
        pixSize *= u.rad
        return (spectral_radiance*pixSize**2).to(u.Jy).value
    return spectral_radiance.value


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
    _r : array_like
    _wl : array_like
        Wavelengths [micron].
    _t : array_like
        Times [second].

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
            Radial grid [mas].
        wl : numpy.ndarray
            Wavelengths [micron].
        t : numpy.ndarray
            Times [second].

        Results
        -------
        radial_profile : numpy.ndarray
        """
        rin, rout = map(lambda x: self.params[x](wl, t), ["rin", "rout"])
        q, p = map(lambda x: self.params[x](wl, t), ["q", "p"])
        dist, inner_temp = map(
            lambda x: self.params[x](wl, t), ["dist", "Tin"])
        dust_mass = self.params["Mdust"](wl, t)*const.M_sun.value*1e3
        kappa_abs = self.params["kappa_abs"](wl, t)

        rin_mas, rout_mas = map(lambda x: 1e3*x/dist, [rin, rout])
        rin_cm, rout_cm = map(
            lambda x: x*self.params["rin"].unit.to(u.cm), [rin, rout])

        # NOTE: Temperature radial profile
        temp_profile = inner_temp*(r / rin_mas)**(-q)

        # NOTE: Surface density radial profile
        if p == 2:
            sigma_in = dust_mass/(2.*np.pi*np.log(rout_cm/rin_cm)*rin_cm**2)
        else:
            f = ((rout_cm/rin_cm)**(2-p)-1)/(2-p)
            sigma_in = dust_mass/(2.*np.pi*f*rin_cm**2)
        sigma_profile = sigma_in*(r / rin_mas)**(-p)
        spectral_density = (2*const.h.value*const.c.value**2/wl**5)*np.divide(
            1., np.exp((const.h.value*const.c.value)/(wl*const.k_B.value*temp_profile))-1)
        spectral_density *= 1-np.exp(-sigma_profile*kappa_abs)
        return np.nan_to_num(np.logical_and(r > rin_mas, r < rout_mas).astype(int)*spectral_density, nan=0)

    @property
    def _r(self):
        """Gets the radial profile [mas]."""
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
        """Sets the radial profile [mas]."""
        return


class oimAsymTempGradient(oimRadialPowerLaw):
    """A ring defined by a radial temperature profile in r^q
    that is multiplied by an azimuthal modulation. 
    and an asymmetric
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
    pixSize : float
        Pixel size [mas].
    _t : numpy.ndarray
        Array of time values [second].
    _wl : numpy.ndarray
        Array of wavelength values [micron].

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
    asymmetric_image = True
    asymmetric_surface_density = False
    const_temperature = False
    continuum_contribution = False

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.normalizeImage = False
        self.params["q"] = oimParam(name="q", value=0, unit=u.one,
                                    description="Power-law exponent for the temperature profile")

        if self.const_temperature:
            self.params["q"].free = False

        self.params["p"].description = "Power-law exponent for the dust surface density profile"
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
        self.params["Teff"] = oimParam(name="Teff", value=0,
                                       unit=u.K, free=False,
                                       description="The star's effective Temperature")
        self.params["lum"] = oimParam(name="lum", value=0,
                                      unit=u.Lsun, free=False,
                                      description="The star's luminosity")

        if self.continuum_contribution:
            self.params["cont_weight"] = oimParam(name="cont_weight", value=0,
                                                  unit=u.one, free=True,
                                                  description="Dust mass continuum absorption coefficient's weight")
            self.params["kappa_cont"] = oimParam(name="kappa_cont", value=0,
                                                 unit=u.cm**2/u.g, free=False,
                                                 description="Continuum dust mass absorption coefficient")

        self._t = np.array([0])  # constant value <=> static model
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._eval(**kwargs)

    def _surface_density_profile(self, xx, yy, wl, t):
        """Calculates the surface density profile.

        This can be azimuthally varied if so specified.

        Parameters
        ----------
        xx : numpy.ndarray
            The x-coordinate grid [mas].
        yy : numpy.ndarray
            The y-coordinate grid [mas].
        wl : numpy.ndarray
            Wavelengths [micron].
        t : numpy.ndarray
            Times [second].

        Returns
        -------
        surface_density_profile : np.ndarray
            The surface density profile [g/cm^2].
        """
        dist = self.params["dist"](wl, t)
        rin, rout = map(lambda x: self.params[x](wl, t), ["rin", "rout"])
        rin_cm = convert_radial_profile_to_meter(rin, dist).to(u.cm).value
        rout_cm = convert_radial_profile_to_meter(rout, dist).to(u.cm).value

        p = self.params["p"](wl, t)
        dust_mass = self.params["Mdust"](wl, t)*const.M_sun.value*1e3
        if p == 2:
            sigma_in = dust_mass/(2.*np.pi*np.log(rout_cm/rin_cm)*rin_cm**2)
        else:
            f = ((rout_cm/rin_cm)**(2-p)-1)/(2-p)
            sigma_in = dust_mass/(2.*np.pi*f*rin_cm**2)

        sigma_profile = sigma_in*(np.sqrt(xx**2+yy**2) / rin)**(-p)
        if self.asymmetric_surface_density:
            return sigma_profile*(1+self._azimuthal_modulation(xx, yy, wl, t))
        return sigma_profile

    def _temperature_profile(self, r, wl, t):
        """Calculates the temperature profile.

        Can be specified to be either as a r^q power law or an a
        constant/idealised temperature profile derived from the star's
        luminosity and the observer's distance to the star and contingent
        only on those values.

        Parameters
        ----------
        r : numpy.ndarray
            Radial grid [mas].
        wl : numpy.ndarray
            Wavelengths [micron].
        t : numpy.ndarray
            Times [second].

        Returns
        -------
        temperature_profile : numpy.ndarray
            The temperature profile [K].

        Notes
        -----
        In case of a radial power law the formula is

        .. math:: T = T_0 * (1+\\frac{r}{R_0})^\\q.

        In case of a constant grey body profile the stellar radius is
        calculated from its lumionsity via

        .. math:: R_* = \\sqrt{\\frac{L_*}{4\\pi\\sigma_sb\\T_*^4}}.

        And with this the individual grain's temperature profile is

        .. math:: T_{grain} = \\sqrt{\\frac{R_*}{2r}}\\cdot T_*.
        """
        if self.const_temperature:
            radius = convert_radial_profile_to_meter(r, self.params["dist"](wl, t))
            luminosity = (self.params["lum"](wl, t) *
                          self.params["lum"].unit).to(u.W)
            stellar_temperature = self.params["Teff"](
                wl, t)*self.params["Teff"].unit
            stellar_radius = np.sqrt(
                luminosity/(4*np.pi*const.sigma_sb*stellar_temperature**4))
            return (np.sqrt(stellar_radius/(2*radius))*stellar_temperature).value
        q, inner_temp = map(lambda x: self.params[x](wl, t), ["q", "Tin"])
        return inner_temp*(r / self.params["rin"](wl, t))**(-q)

    def _image(self, xx: np.ndarray, yy: np.ndarray,
               wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Combines the various radial profiles into an image.

        If physical output is specified, the model will produce Jansky per
        pixel else unitless intensity.

        Parameters
        ----------
        xx : numpy.ndarray
            The x-coordinate grid [mas].
        yy : numpy.ndarray
            The y-coordinate grid [mas].
        wl : numpy.ndarray
            Wavelengths [micron].
        t : numpy.ndarray
            Times [second].

        Returns
        -------
        image : numpy.ndarray
        """
        r = np.sqrt(xx**2+yy**2)
        rin, rout = map(lambda x: self.params[x](wl, t), ["rin", "rout"])

        temperature_profile = self._temperature_profile(r, wl, t)
        spectral_density = calculate_intensity(wl, temperature_profile,
                                               pixSize=self._pixSize)
        sigma_profile = self._surface_density_profile(xx, yy, wl, t)
        if self.continuum_contribution:
            optical_depth = -sigma_profile*(self.params["kappa_abs"](wl, t) +
                                            self.params["cont_weight"](wl, t) *
                                            self.params["kappa_cont"](wl, t))
        else:
            optical_depth = -sigma_profile*self.params["kappa_abs"](wl, t)
        radial_profile = np.logical_and(r > rin, r < rout).astype(int)
        return np.nan_to_num(radial_profile * spectral_density *
                             (1 - np.exp(optical_depth)), nan=0)

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
        if self.asymmetric_image:
            return self._image(xx, yy, wl, t) * \
                (1+self._azimuthal_modulation(xx, yy, wl, t))
        return self._image(xx, yy, wl, t)


class oimAsymSDTempGradient(oimAsymTempGradient):
    """A ring defined by a radial temperature profile in r^q
    and an asymmetric radial dust surface density profile in r^p.

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
    pixSize : float
        Pixel size [mas].
    _t : numpy.ndarray
        Array of time values [second].
    _wl : numpy.ndarray
        Array of wavelength values [micron].

    Methods
    -------
    _imageFunction(xx, yy, wl, t)
        Calculates a 2D-image from a dust-surface density- and
        temperature profile.
    """
    name = "Asymmetric Surface Density Temperature Gradient"
    shortname = "AsymSDTempGrad"
    asymmetric_image = False
    asymmetric_surface_density = True


class oimAsymSDGreyBody(oimAsymSDTempGradient):
    """A ring defined by a grain grey body temperature profile and
    and an asymmetric radial dust surface density profile in r^p
    with opacity curves for the dust surface density profile.

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
    pixSize : float
        Pixel size [mas].
    _t : numpy.ndarray
        Array of time values [second].
    _wl : numpy.ndarray
        Array of wavelength values [micron].

    Methods
    -------
    _imageFunction(xx, yy, wl, t)
        Calculates a 2D-image from a dust-surface density- and
        temperature profile.
    """
    name = "Asymmetric Surface Density Grey Body"
    shortname = "AsymSDGreyBody"
    const_temperature = True
    continuum_contribution = True
