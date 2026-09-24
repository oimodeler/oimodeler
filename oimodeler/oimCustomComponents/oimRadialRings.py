import numpy as np

from ..oimComponent import oimComponentRadialProfile
from ..oimOptions import oimOptions
from ..oimParam import _standardParameters, oimParam


class oimRadialPowRing(oimComponentRadialProfile):
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

    name = "Radial Pow Ring"
    shortname = "PowR"
    elliptic = True

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["din"] = oimParam(base="din")
        self.params["dout"] = oimParam(base="dout")
        self.params["p"] = oimParam(
            name="p",
            description="Power-law exponent for radial ring",
            mini=-1,
            maxi=0,
            base="amp",
        )
        self._eval(**kwargs)

    @property
    def r(self):
        """Gets the radial profile [mas]."""
        dim = self.dim.value
        rin, rout = self.din.value / 2, self.dout.value / 2
        grid_type = oimOptions.model.grid.type
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
        rin, rout = self.din.value / 2, self.dout.value / 2
        image = ((r >= rin) & (r <= rout)) * (r / rin) ** self.p(wl, t)
        return image * np.ones_like(wl)

class oimRadialPowRing2(oimComponentRadialProfile):
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

    name = "Radial Pow Ring alt"
    shortname = "PowR2"
    elliptic = False

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["din"] = oimParam(base="din")
        self.params["w"] = oimParam(base="w")
        self.params["p"] = oimParam(
            name="p",
            description="Power-law exponent for radial ring",
            mini=-1,
            maxi=0,
            base="amp",
        )
        self._eval(**kwargs)

    @property
    def r(self):
        """Gets the radial profile [mas]."""
        rin, dim = self.din.value / 2, self.dim.value
        rout = rin + self.w.value
        grid_type = oimOptions.model.grid.type
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


class oimRadialExpRing(oimComponentRadialProfile):
    name = "Radial Exponential Ring"
    shortname = "ExpR"

    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(base="d")
        self.params["fwhm"] = oimParam(base="fwhm")
        self.params["dim"] = oimParam(base="dim")
        self._eval(**kwargs)

    @property
    def r(self):
        if False:
            fwhm_max = np.max(self.fwhm(self._wl, self._t))
            r0_max = np.max(self.d(self._wl, self._t)) / 2
        else:
            fwhm_max, r0_max = self.fwhm(1e99), self.d(1e99)

        rmax = r0_max + 8 * fwhm_max
        self._r = np.linspace(0, 1, self.dim.value) * rmax
        return self._r

    def _radialProfileFunction(self, r, wl, t):
        r0, fwhm = self.d(wl, t) / 2, self.fwhm(wl, t)
        return (r > r0) * np.exp(-0.692 * np.divide(r - r0, fwhm))
