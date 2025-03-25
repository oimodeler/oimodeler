import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentFourier
from ..oimParam import _standardParameters, oimParam


class oimGaussLorentz(oimComponentFourier):
    name = "Gauss-Lorentzian"
    shortname = "GL"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["hlr"] = oimParam(**_standardParameters["hlr"])
        self.params["flor"] = oimParam(**_standardParameters["f"])
        self.params["flor"].name = "flor"

        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        """The visibility function for the Gauss-Lorentzian model."""
        flor = self.params["flor"](wl, t)
        xx = (
            np.pi
            * self.params["hlr"](wl, t)
            * self.params["hlr"].unit.to(u.rad)
            * rho
        )
        return (1 - flor) * np.exp(-(xx**2) / np.log(2)) + flor * np.exp(
            -2 * xx / np.sqrt(3)
        )

    def _imageFunction(self, xx, yy, wl, t):
        hlr, flor = self.params["hlr"](wl, t), self.params["flor"](wl, t)
        radius = np.hypot(xx, yy)
        image_gauss = (
            np.log(2)
            / (np.pi * hlr**2)
            * np.exp(-((radius / hlr) ** 2) * np.log(2))
        )
        image_lor = (
            hlr
            / (2 * np.pi * np.sqrt(3))
            * (hlr**2 / 3 + radius**2) ** (-3 / 2)
        )
        return (1 - flor) * image_gauss + flor * image_lor
