import numpy as np

from ..oimParam import _standardParameters, oimParam
from ..oimComponent import oimComponentRadialProfile


class oimRadialRing(oimComponentRadialProfile):
    """A ring defined by a radial intensity profile in r^p.

    It accounts for elongation and rotation

    Parameters
    ----------
    din : u.mas
        Inner radius of the disk
    dout : u.mas
        Outer radius of the disk
    p : u.one
        Power-law exponent for the radial profile
    pa : u.deg
        Positional angle
    elong : float
        Elongation of the disk
    dim : u.one
        Dimension of the image

    Attributes
    ----------
    params : Dict[str, oimParam]
        Dictionary of parameters
    _t : np.ndarray
        Array of time values
    _wl : np.ndarray
        Array of wavelength values
    _r : np.ndarray

    Methods
    -------
    _radialProfileFunction(xx, yy, wl, t)
        Calculates a radial power law profile
    """
    name = "Radial Ring"
    shortname = "RadRing"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self.params["p"] = oimParam(**_standardParameters["p"])

        self._t = np.array([0])  # constant value <=> static model
        self._wl = np.array([0])  # constant value <=> achromatic model
        self._eval(**kwargs)

    def _radialProfileFunction(self, r, wl, t):
        """Calculates a radial power law profile"""
        p = self.params["p"](wl, t)
        rin, rout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))
        return np.nan_to_num(np.logical_and(r > rin, r < rout).astype(int)*(r / rin)**p, nan=0)

    @property
    def _r(self):
        """Gets the radius"""
        if False:
            dout = self.params["dout"](self._wl, self._t)
        else:
            dout = self.params["dout"](1e99)
        rmax = 4*dout/2.
        return np.linspace(0, 1, self._dim)*rmax

    @_r.setter
    def _r(self, r):
        """Sets the radius"""
        pass
