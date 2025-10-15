import numpy as np

from ..oimParam import _standardParameters, oimParam
from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import oimOptions


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
    elliptic = True

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
        # HACK: Sets the multi wavelength coordinates properly. Does not account for time, improves computation time.
        wl, p = np.unique(wl), self.params["p"](wl, t)
        rin, rout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))

        if len(r.shape) == 3:
            r = r[0, 0][np.newaxis, np.newaxis, :]
            wl = wl[np.newaxis, :, np.newaxis]
        else:
            r, wl = r[np.newaxis, :], wl[np.newaxis, :]

        image = np.nan_to_num(np.logical_and(r > rin, r < rout).astype(int) * (r / rin) ** p, nan=0)
        return image * np.ones_like(wl)

    @property
    def _r(self):
        """Gets the radial profile [mas]."""
        rin, rout = map(lambda x: self.params[x].value/2, ("din", "dout"))
        if oimOptions.model.grid.type == "linear":
            return np.linspace(rin, rout, self.params["dim"].value)
        return np.logspace(0.0 if rin == 0 else np.log10(rin),
                           np.log10(rout), self.params["dim"].value)

    @_r.setter
    def _r(self, r: np.ndarray):
        """Sets the radius."""
        pass
