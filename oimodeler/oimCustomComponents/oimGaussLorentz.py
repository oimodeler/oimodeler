import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentFourier
from ..oimParam import _standardParameters, oimParam


# TODO: Rename with E convnetion? For Elliptical?
class oimGaussLorentz(oimComponentFourier):
    """Gaussian-Lorentzian component defined in the Fourier space.

    Parameters
    ----------
    x : float or oimInterp
        x pos of the component (mas). Defaults to ``0``.
    y : float or oimInterp
        y pos of the component (mas). Defaults to ``0``.
    f : float or oimInterp
        Flux (ratio) of the component. Defaults to ``1``.
    hlr : float or oimInterp
    flor : float or oimInterp
    pa : float or oimInterp
        Position angle of the major axis (deg). Defaults to ``0``.
    elong : float or oimInterp
        Elongation of the major axis. Defaults to ``1``.

    Attributes
    ----------
    x : oimParam
        x pos of the component (mas).
    y : oimParam
        y pos of the component (mas).
    f : oimParam
        Flux (ratio) of the component.
    hlr : oimParam
    flor : oimParam
    pa : oimParam
        Position angle of the major axis (deg).
    elong : oimParam
        Elongation of the major axis.

    Notes
    -----
    From `2017A%26A...599A..85L <https://scixplorer.org/abs/2017A%26A...599A..85L>`_.
    """

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
