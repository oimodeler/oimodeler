import astropy.units as u
import astropy.constants as const
import numpy as np
from astropy.modeling import models

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam


class oimStar(oimComponentImage):
    """A physical point source.

    Parameters
    ----------
    dist : float
        The distance to the star [pc].
    lum : float
        The star's luminosity [Lsun].
    eff_temp : float
        The star's effective temperature [K].

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
    _convert_radius_to_parallax(orbital_radius, wl, t)
        Calculates the parallax from the orbital radius.
    _calc_stellar_radius(wl, t)
        Calculates the stellar radius.
    _calc_stellar_flux(wl, t)
        Calculates the stellar flux from its distance and radius.
    _imageFunction(xx, yy, wl, t)
        Creates a 2D-image with the stellar flux in Jansky at its centre
        and only 0 for the rest.
    """
    name = "Physical Star"
    shortname = "Star"
    elliptic = False

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["dist"] = oimParam(name="dist", value=0,
                                       unit=u.pc, free=False,
                                       description="Distance of the star")
        self.params["lum"] = oimParam(name="lum", value=0,
                                      unit=u.Lsun, free=False,
                                      description="The star's luminosity")
        self.params["eff_temp"] = oimParam(name="eff_temp", value=0,
                                           unit=u.K, free=False,
                                           description="The star's effective temperature")

        self._wl = None
        self._t = np.array([0])
        self._eval(**kwargs)

    def _convert_radius_to_parallax(self, orbital_radius: u.m, wl, t) -> u.mas:
        """Calculates the parallax from the orbital radius.

        Parameters
        ----------
        orbital_radius : astropy.units.m
            The orbital radius.

        Returns
        -------
        parallax : astropy.units.mas
            The angle of the orbital radius.

        Notes
        -----
        The formula for the angular diameter $ \delta = \frac{d}{D} $ is used.
        This produces an output in radians.
        """
        distance = self.params["dist"](wl, t)*self.params["dist"].unit
        return u.rad.to(u.mas)*orbital_radius.to(u.m)/distance.to(u.m)

    def _calc_stellar_radius(self, wl, t) -> u.m:
        """Calculates the stellar radius.

        Returns
        -------
        stellar_radius : astropy.units.m
            The star's radius.
        """
        luminosity, effective_temperature = map(
            lambda x: self.params[x](wl, t), ["lum", "eff_temp"])
        return np.sqrt(luminosity*self.params["lum"].unit
                       / (4*np.pi*const.sigma_sb
                          * (effective_temperature*self.params["eff_temp"].unit)**4))

    def _calc_stellar_flux(self, wl, t) -> u.Jy:
        """Calculates the stellar flux from its distance and radius.

        Returns
        -------
        stellar_flux : astropy.units.Jy
            The star's flux.
        """
        plancks_law = models.BlackBody(
            temperature=self.params["eff_temp"](wl, t) * self.params["eff_temp"].unit)
        spectral_radiance = plancks_law(wl).to(u.erg/(u.cm**2*u.Hz*u.s*u.mas**2))
        stellar_radius = self.calc_stellar_radius(wl, t)

        # TODO: Check if that can be used in this context -> The conversion
        stellar_radius_angular = self._convert_radius_to_parallax(stellar_radius, wl, t)

        # TODO: Check this calculation for the pi in the formulas
        return (spectral_radiance*stellar_radius_angular**2).to(u.Jy)

    def _imageFunction(self, xx, yy, wl, t):
        """Creates a 2D-image with the stellar flux in Jansky at its centre
        and only 0 for the rest."""
        image = xx*0
        val = np.abs(xx)+np.abs(yy)
        idx = np.unravel_index(np.argmin(val), np.shape(val))
        image[idx] = self._calc_stellar_flux(wl, t).value
        return image