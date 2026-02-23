from pathlib import Path

import astropy.units as u
import numpy as np
from numpy.typing import NDArray

from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import constants as const
from ..oimOptions import oimOptions
from ..oimParam import _standardParameters, oimParam
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
    r0 : float
        Reference radius [au].
    T0 : float
        Temperature at reference radius r0 [K].
    Sigma0 : float
         Dust surface density at reference radius r0 [M_sun].
    Mdust : float
         Mass of the dusty disk [M_sun].
    q : float
        Power-law exponent for the temperature profile.
    p : float
        Power-law exponent for the dust surface density profile.

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
    compute_sigma0 = True
    elliptic = True

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(**_standardParameters["rin"])
        self.params["rout"] = oimParam(**_standardParameters["rout"])
        self.params["r0"] = oimParam(
            name="r0",
            value=1,
            unit=u.au,
            free=False,
            description="Reference radius",
        )
        self.params["temp0"] = oimParam(
            name="temp0",
            value=300,
            unit=u.K,
            description="Temperature at reference radius",
        )

        if "sigma0" in kwargs or not kwargs.get(
            "compute_sigma0", self.compute_sigma0
        ):
            self.compute_sigma0 = False
            self.params["sigma0"] = oimParam(
                name="sigma0",
                value=1e-2,
                unit=u.g / u.cm**2,
                description="Dust surface density at reference radius",
            )
        else:
            self.compute_sigma0 = True
            self.params["dust_mass"] = oimParam(
                name="dust_mass",
                value=0.2,
                unit=u.M_sun,
                description="Mass of the dusty disk",
            )

        self.params["q"] = oimParam(
            name="q",
            value=-0.5,
            unit=u.one,
            description="Power-law exponent for the temperature",
            mini=-1,
            maxi=0,
        )
        self.params["p"] = oimParam(
            name="p",
            value=-0.5,
            unit=u.one,
            description="Power-law exponent for the dust surface density",
            mini=-1,
            maxi=0,
        )
        self.params["kappa_abs"] = oimParam(
            name="kappa_abs",
            value=0,
            unit=u.cm**2 / u.g,
            description="Absorption (silicate) opacity",
            free=False,
        )

        # TODO: Potentially generalis this for many different contributions
        if "kappa_cont" in kwargs:
            self.params["kappa_cont"] = oimParam(
                name="kappa_cont",
                value=0,
                unit=u.cm**2 / u.g,
                description="Absorption continuum opacity",
                free=False,
            )
            self.params["kappa_ratio"] = oimParam(
                name="kappa_ratio",
                value=1,
                unit=u.one,
                description="Silicate to continuum ratio",
                mini=0,
                maxi=1,
            )

        self.params["dist"] = oimParam(**_standardParameters["dist"])
        self.params["f"].free = False
        self._wl, self._t = None, [0]
        self._eval(**kwargs)

    @property
    def _r(self):
        """Gets the radial profile (mas)."""
        rin, rout = self.params["rin"].value, self.params["rout"].value
        dim, dist = self.params["dim"].value
        rin = linear_to_angular(rout, dist) * 1e3
        rout = linear_to_angular(rout, dist) * 1e3
        if oimOptions.model.grid.type == "linear":
            return np.linspace(rin, rout, dim)
        return np.logspace(
            0.0 if rin == 0 else np.log10(rin), np.log10(rout), dim
        )

    @_r.setter
    def _r(self, value) -> None:
        """Sets the radial profile [mas]."""
        return None

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
            Wavelengths [m].
        t : numpy.ndarray
            Times [second].

        Results
        -------
        radial_profile : numpy.ndarray
        """
        # HACK: Sets the multi wavelength coordinates properly.
        # Does not account for time, improves computation time.
        wl = np.unique(wl)
        dist = self.params["dist"].value
        kappa_abs = self.params["kappa_abs"](wl, t)
        if "kappa_cont" in self.params:
            ratio = self.params["kappa_ratio"].value
            kappa_cont = self.params["kappa_cont"](wl, t)
            kappa_abs = (1 - ratio) * kappa_abs + ratio * kappa_cont

        if self.flat:
            elong = 1 / self.params["cosi"].value
        else:
            elong = self.params["elong"].value

        # HACK: Adding the correct dimensions to the radial grid, wavelength array and opacity table
        if len(r.shape) == 3:
            r = r[0, 0][np.newaxis, np.newaxis, :]
            wl, kappa_abs = map(
                lambda x: x[np.newaxis, :, np.newaxis], [wl, kappa_abs]
            )
        else:
            wl, kappa_abs = map(lambda x: x[:, np.newaxis], [wl, kappa_abs])
            r = r[np.newaxis, :]

        rin, rout, r0 = map(
            lambda x: self.params[x](wl, t), ["rin", "rout", "r0"]
        )
        rin_cm, rout_cm, r0_cm = map(
            lambda x: x * self.params["rin"].unit.to(u.cm), [rin, rout, r0]
        )

        q, p = map(lambda x: self.params[x](wl, t), ["q", "p"])
        r0 = linear_to_angular(r0, dist) * 1e3
        temp = self.params["temp0"](wl, t) * (r / r0) ** q

        if self.compute_sigma0:
            dust_mass = self.params["dust_mass"](wl, t) * self.params[
                "dust_mass"
            ].unit.to(u.g)
            if p == -2:
                sigma0 = dust_mass / (
                    2.0 * np.pi * np.log(rout_cm / rin_cm) * r0_cm**2
                )
            else:
                f = (rout_cm ** (2 + p) - rin_cm ** (2 + p)) / (2 + p)
                sigma0 = dust_mass / (2.0 * np.pi * f * r0_cm ** (-p))
        else:
            sigma0 = self.params["sigma0"](wl, t)

        sigma = sigma0 * (r / r0) ** p
        emissivity = 1 - np.exp(-sigma * kappa_abs * elong)
        spectral_density = (
            blackbody(temp, const.c / wl) * emissivity * (1 / elong)
        )
        rin_mas, rout_mas = map(
            lambda x: linear_to_angular(x, dist) * 1e3, [rin, rout]
        )
        radial_profile = ((r > rin_mas) & (r < rout_mas)).astype(int)
        image = np.nan_to_num(radial_profile * spectral_density, nan=0)
        if len(r.shape) == 3:
            return image

        return image
