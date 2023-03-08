import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam, _standardParameters


class oimRadialPowerLaw(oimComponentImage):
    """A 2D-radial power law distribution

    It accounts for elongation and rotation
    """
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._t, self._wl = [None]*2
        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self.params["pixSize"] = oimParam(**_standardParameters["pixSize"])
        self.params["p"] = oimParam(name="p", value=0, description="Power-law exponent")
        self._eval(**kwargs)

    def _imageFunction(self, xx, yy, wl, t):
        """The function describing the component's image"""
        self._pixSize = self.params["pixSize"](
            wl, t)*self.params["pixSize"].unit.to(u.rad)
        rin, rout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))
        r, p = np.sqrt(xx**2+yy**2), self.params["p"](wl, t)
        return np.nan_to_num(np.logical_and(r > rin, r < rout).astype(int)*np.divide(r, rin)**p, nan=0)


class oimAsymmetricRadialPowerLaw(oimRadialPowerLaw):
    """An asymmetrically azimuthally modulated 2D-radial power law distribution

    It accounts for elongation and rotation as well
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["a"] = oimParam(
            name="a", value=0, description="Azimuthal modulation amplitude")
        self.params["phi"] = oimParam(
            name="phi", value=0, description="Azimuthal modulation angle")
        self._eval(**kwargs)

    def _azimuthal_modulation(self, x_arr, y_arr, wl, t):
        """Calculates the azimuthal modulation"""
        # TEST: Is it the other way around (y_arr, x_arr)?
        polar_angle = np.arctan2(x_arr, y_arr)
        return (1 + self.params["a"](wl, t)*np.cos(polar_angle-self.params["phi"](wl, t)))

    def _imageFunction(self, xx, yy, wl, t):
        img = super()._imageFunction(xx, yy, wl, t)
        return self._azimuthal_modulation(xx, yy, wl, t)*img
