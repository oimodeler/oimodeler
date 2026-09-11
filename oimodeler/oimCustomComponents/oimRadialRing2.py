import numpy as np

from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import oimOptions
from ..oimParam import oimParam


class oimRadialRing2(oimComponentRadialProfile):
    """A ring defined by a radial intensity profile in r^p.

    It accounts for elongation and rotation.

    Parameters
    ----------
    din : float
        Inner diameter of the ring [mas].
    w : float
        width of the ring [mas].
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

    name = "Radial Ring2"
    shortname = "RadRing2"
    elliptic = False

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["din"] = oimParam(base="din")
        self.params["w"] = oimParam(base="w")
        self.params["p"] = oimParam(base="p")
        self._eval(**kwargs)

    @property
    def r(self):
        """Gets the radial profile [mas]."""
        rin = self.din.value / 2
        rout = rin + self.w.value
        dim, dist = self.dim.value, self.dist.value
        grid_type = oimOptions.model.grid.type

        rin, rout = rin / dist * 1e3, rout / dist * 1e3
        if grid_type == "linear":
            self._r = np.linspace(rin, rout, dim)
        else:
            if rin <= 0:
                raise ValueError("Logarithmic grid requires rin > 0.")

            self._r = np.logspace(np.log10(rin), np.log10(rout), dim)

        return self._r

    def _radialProfileFunction(
        self, r: np.ndarray, wl: np.ndarray, t: np.ndarray
    ) -> np.ndarray:
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
        rin = self.din(wl, t) / 2
        rout = rin + self.w(wl, t)
        image = ((r >= rin) & (r <= rout)) * (r / rin) ** self.p(wl, t)
        return image * np.ones_like(wl)
