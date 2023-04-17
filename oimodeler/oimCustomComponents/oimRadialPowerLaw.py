import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam, _standardParameters


class oimRadialPowerLaw(oimComponentImage):
    """A 2D-radial power law distribution.

    It accounts for elongation and rotation

    Parameters
    ----------
    din : u.mas
        Inner radius of the disk
    dout : u.mas
        Outer radius of the disk
    pixSize: u.mas/px
        Pixel size
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
    _pixSize: u.mas/px
        Pixel size
    _t : np.ndarray
        Array of time values
    _wl : np.ndarray
        Array of wavelength values

    Methods
    -------
    _imageFunction(xx, yy, wl, t)
        Calculates a radial power law
    """
    name = "Radial Power Law"
    shortname = "RadPowerLaw"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(**_standardParameters["din"])
        self.params["rout"] = oimParam(**_standardParameters["dout"])
        self.params["fov"] = oimParam(**_standardParameters["fov"])
        self.params["pixSize"] = oimParam(**_standardParameters["pixSize"])
        self.params["p"] = oimParam(**_standardParameters["p"])
        self._t = np.array([0])  # constant value <=> static model
        self._wl = np.array([0])  # constant value <=> achromatic model
        self._eval(**kwargs)

    @property
    def _pixSize(self):
        self.params["pixSize"].value = pix = self.params["fov"].value/self.params["dim"].value
        return pix*self.params["fov"].unit.to(u.rad)

    @_pixSize.setter
    def _pixSize(self, value: u.mas):
        pass

    def _imageFunction(self, xx, yy, wl, t):
        """Calculates a radial power law"""
        din, dout = map(lambda x: self.params[x](wl, t)/2, ("din", "dout"))
        r, q = np.sqrt(xx**2+yy**2), self.params["p"](wl, t)
        return np.nan_to_num(np.logical_and(r > din, r < dout).astype(int)*(r / din)**q, nan=0)


class oimAsymRadialPowerLaw(oimRadialPowerLaw):
    """An asymmetrically azimuthally modulated 2D-radial power law distribution.

    It accounts for elongation and rotation, as well

    Parameters
    ----------
    din : u.mas
        Inner radius of the disk
    dout : u.mas
        Outer radius of the disk
    pixSize: u.mas/px
        Pixel size
    p : u.one
        Power-law exponent for the radial profile
    a : u.one
        Azimuthal modulation amplitude
    phi : float
        Azimuthal modulation angle
    pa : u.deg
        Positional angle
    elong : u.one
        Elongation of the disk
    dim : u.one
        Dimension of the image [px]

    Attributes
    ----------
    params : Dict[str, oimParam]
        Dictionary of parameters
    _t : np.array
        Array of time values
    _wl : np.ndarray
        Array of wavelength values
    _r : np.ndarray

    Methods
    -------
    _azimuthal_modulation(xx, yy, wl, t)
        Calculates the azimuthal modulation
    _imageFunction(xx, yy, wl, t)
        Calculates a radial power law with an azimuthal asymmetry
    """
    name = "Asymmetric Radial Power Law"
    shortname = "AsymRadPowerLaw"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["a"] = oimParam(name="a", value=0, unit=u.one,
                                    description="Azimuthal modulation amplitude")
        self.params["phi"] = oimParam(name="phi", value=0, unit=u.deg,
                                      description="Azimuthal modulation angle")
        self._eval(**kwargs)

    def _azimuthal_modulation(self, xx, yy, wl, t):
        """Calculates the azimuthal modulation"""
        # TEST: Is it the other way around (y_arr, x_arr)?
        polar_angle = np.arctan2(xx, yy)
        return self.params["a"](wl, t)*np.cos(polar_angle-self.params["phi"](wl, t))

    def _imageFunction(self, xx, yy, wl, t):
        """Calculates a radial power law with an azimuthal asymmetry"""
        img = super()._imageFunction(xx, yy, wl, t)
        return self._azimuthal_modulation(xx, yy, wl, t)*(1 + img)
