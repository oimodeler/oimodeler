import numpy as np

from ..oimParam import _standardParameters, oimParam
from ..oimComponent import oimComponentRadialProfile


class oimRadialRing(oimComponentRadialProfile):
    """A ring defined by a radial intensity profile in r^p.

    It accounts for elongation and rotation.

    Parameters
    ----------
    din : float
        Inner radius of the disk [mas].
    dout : float
        Outer radius of the disk [mas].
    p : float
        Power-law exponent for the radial profile.
    pa : float
        Positional angle.
    elong : float
        Elongation of the disk.
    dim : float
        Dimension of the image.

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
    _radialProfileFunction(xx, yy, wl, t)
        Calculates a radial power law profile.
    """
    name = "Radial Ring"
    shortname = "RadRing"
    elliptic = False

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self.params["p"] = oimParam(**_standardParameters["p"])

        self._t = np.array([0])  # constant value <=> static model
        self._wl = None
        self._eval(**kwargs)

    def _radialProfileFunction(self, r: np.ndarray,
                               wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates a radial power law profile.

        Parameters
        ----------
        r : numpy.ndarray
            Radial grid.
        wl : numpy.ndarray
            Wavelengths.
        t : numpy.ndarray
            Times.

        Returns
        -------
        radial_profile : numpy.ndarray
        """
        p = self.params["p"](wl, t)
        rin, rout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))
        return np.nan_to_num(np.logical_and(r > rin, r < rout).astype(int)*(r / rin)**p, nan=0)

    @property
    def _r(self):
        """Gets the radius."""
        if False:
            dout = self.params["dout"](self._wl, self._t)
        else:
            dout = self.params["dout"](1e99)
        rmax = 4*dout/2.
        return np.linspace(0, 1, self.params["dim"](self._wl, self._t))*rmax

    @_r.setter
    def _r(self, r: np.ndarray):
        """Sets the radius."""
        pass
