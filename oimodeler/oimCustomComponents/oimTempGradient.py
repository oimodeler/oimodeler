import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import constants as const
from ..oimOptions import oimOptions
from ..oimParam import oimParam
from ..oimUtils import blackbody, linear_to_angular


class oimTempGrad(oimComponentRadialProfile):
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
        self.params["rin"] = oimParam(
            name="rin",
            value=0,
            unit=u.au,
            description="Inner radius of the disk",
        )
        self.params["rout"] = oimParam(
            name="rout",
            value=0,
            unit=u.au,
            description="Outer radius of the disk",
        )
        self.params["q"] = oimParam(
            name="q",
            value=1,
            mini=0,
            maxi=1,
            unit=u.one,
            description="Power-law exponent for the temperature",
        )
        self.params["p"] = oimParam(
            name="p",
            value=1,
            mini=0,
            maxi=1,
            unit=u.one,
            description="Power-law exponent for the dust surface density",
        )
        self.params["dust_mass"] = oimParam(
            name="dust_mass",
            value=0,
            unit=u.M_sun,
            description="Mass of the dusty disk",
        )
        self.params["temp0"] = oimParam(
            name="temp0",
            value=0,
            unit=u.K,
            free=False,
            description="Temperature at reference radius",
        )
        self.params["kappa_abs"] = oimParam(
            name="kappa_abs",
            value=0,
            unit=u.cm**2 / u.g,
            free=False,
            description="Dust mass absorption coefficient",
        )
        self.params["dist"] = oimParam(
            name="dist",
            value=0,
            unit=u.pc,
            free=False,
            description="Distance of the star",
        )
        self.params["r0"] = oimParam(
            name="r0",
            value=1,
            unit=u.au,
            free=False,
            description="Reference radius",
        )
        self.params["f"].free = False

        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static
        self._eval(**kwargs)

    def _radialProfileFunction(
        self, r: np.ndarray, wl: np.ndarray, t: np.ndarray
    ) -> np.ndarray:
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
        # HACK: Sets the multi wavelength coordinates properly.
        # Does not account for time, improves computation time.
        # TODO: Check if this is still needed and if not remove it
        wl = np.unique(wl)
        if np.shape(wl) == np.shape(self._wl):
            dist = self.params["dist"](wl, t)
            kappa_abs = self.params["kappa_abs"](wl, t)
        else:
            dist = self.params["dist"](wl, t)
            kappa_abs = np.interp(wl,self._wl,self.params["kappa_abs"](self._wl, t))
        if len(r.shape) == 3:
            r = r[0, 0][np.newaxis, np.newaxis, :]
            wl, kappa_abs = map(
                lambda x: x[np.newaxis, :, np.newaxis], [wl, kappa_abs]
            )
        else:
            wl, kappa_abs = map(lambda x: x[:, np.newaxis], [wl, kappa_abs])
            r = r[np.newaxis, :]
        rin, rout = map(lambda x: self.params[x](wl, t), ["rin", "rout"])
        q, p = map(lambda x: self.params[x](wl, t), ["q", "p"])
        dust_mass = self.params["dust_mass"](wl, t) * self.params[
            "dust_mass"
        ].unit.to(u.g)

        r0 = linear_to_angular(self.params["r0"](wl, t), dist) * 1e3
        temp = self.params["temp0"](wl, t) * (r / r0) ** (q)
        rin_cm, rout_cm = map(
            lambda x: x * self.params["rin"].unit.to(u.cm), [rin, rout]
        )
        if p == -2:
            sigma0 = dust_mass / (
                2.0 * np.pi * np.log(rout_cm / rin_cm) * rin_cm**2
            )
        else:
            f = ((rout_cm / rin_cm) ** (2 + p) - 1) / (2 + p)
            sigma0 = dust_mass / (2.0 * np.pi * f * rin_cm**2)
        
        sigma = sigma0 * (r / r0) ** (p)
        epsilon = 1 - np.exp(-sigma * kappa_abs)
        nu=const.c/wl
        spectral_density = blackbody(temp, nu) * epsilon
        rin_mas, rout_mas = map(
            lambda x: linear_to_angular(x, dist) * 1e3, [rin, rout]
        )
        radial_profile = ((r > rin_mas) & (r < rout_mas)).astype(int)
        image = np.nan_to_num(radial_profile * spectral_density, nan=0)
        if len(r.shape) == 3:
            return image

        return image

    @property
    def _r(self):
        """Gets the radial profile (mas)."""
        rin = (
            linear_to_angular(
                self.params["rin"].value, self.params["dist"].value
            )
            * 1e3
        )
        rout = (
            linear_to_angular(
                self.params["rout"].value, self.params["dist"].value
            )
            * 1e3
        )
        if oimOptions.model.grid.type == "linear":
            return np.linspace(rin, rout, self.params["dim"].value)
        return np.logspace(
            0.0 if rin == 0 else np.log10(rin),
            np.log10(rout),
            self.params["dim"].value,
        )

    @_r.setter
    def _r(self, value):
        """Sets the radial profile [mas]."""
        return


class oimAsymTempGrad(oimTempGrad):
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

    name = "Asymmetric Temperature Gradient"
    shortname = "AsymTempGrad"
    elliptic = True
    asymmetric = True
