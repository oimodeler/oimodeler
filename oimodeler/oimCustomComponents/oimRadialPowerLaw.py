import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentImage
from ..oimParam import oimParam, _standardParameters
from ..oimUtils import get_next_power_of_two


class oimRadialPowerLaw(oimComponentImage):
    """A 2D-radial power law distribution.

    It accounts for elongation and rotation.

    Parameters
    ----------
    rin : float
        Inner radius of the disk [mas].
    rout : float
        Outer radius of the disk [mas].
    pixSize: float
        Pixel size [mas/px].
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
    _pixSize
    params : dict with keys of str and values of oimParam
        Dictionary of parameters.
    _wl : array_like
        Wavelengths.
    _t : array_like
        Times.

    Methods
    -------
    _azimuthal_modulation(xx, yy, wl, t)
        Calculates the azimuthal modulation.
    _imageFunction(xx, yy, wl, t)
        Calculates the image from a radial power law.
    """
    name = "Radial Power Law"
    shortname = "RadPowerLaw"
    elliptic = True
    asymmetric = False

    def __init__(self, **kwargs):
        """The class's constructor."""
        super().__init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0, unit=u.mas,
                                      description="Inner radius of the disk")
        self.params["rout"] = oimParam(name="rout", value=0, unit=u.mas,
                                       description="Outer radius of the disk")
        self.params["p"] = oimParam(**_standardParameters["p"])
        self.params["fov"] = oimParam(**_standardParameters["fov"])
        self.params["pixSize"] = oimParam(**_standardParameters["pixSize"])

        if self.asymmetric:
            self.params["a"] = oimParam(name="a", value=0, unit=u.one,
                                        description="Azimuthal modulation amplitude")
            self.params["phi"] = oimParam(name="phi", value=0, unit=u.deg,
                                          description="Azimuthal modulation angle")

        self._t = np.array([0])  # constant value <=> static model
        self._wl = np.array([0])  # constant value <=> achromatic model
        self._eval(**kwargs)

    @property
    def _pixSize(self):
        """Returns the pixel size [mas/px].

        Additionally this property also calculates the pixel size from
        the field of view [mas] and the dimension [px], if the pixel size
        [mas/px] is not provided. Otherwise, it will automatically adjust the
        dimension [px] depending on the field of view [mas] and the pixel
        size [mas/px].
        """
        if self.params["fov"].value == 0:
            raise ValueError("The field of view must be specified as a keyword argument.")
        if self.params["pixSize"].value == 0:
            self.params["pixSize"].value = pix = self.params["fov"].value/self.params["dim"].value
            return pix*self.params["fov"].unit.to(u.rad)
        dim = get_next_power_of_two(self.params["fov"].value/\
                                    self.params["pixSize"].value)
        if dim != self.params["dim"].value:
            self.params["dim"].value = dim
        return self.params["pixSize"].value*self.params["pixSize"].unit.to(u.rad)

    @_pixSize.setter
    def _pixSize(self, value: u.mas):
        """The pixel size [mas/px]."""
        pass

    def _azimuthal_modulation(self, xx: np.ndarray, yy: np.ndarray,
                              wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates the azimuthal modulation.

        Parameters
        ----------
        xx : numpy.ndarray
            The x-coordinate grid
        yy : numpy.ndarray
            The y-coordinate grid
        wl : numpy.ndarray
            Wavelengths.
        t : numpy.ndarray
            Times.

        Returns
        -------
        azimuthal_modulation : numpy.ndarray
        """
        phi = self.params["phi"](wl, t) * self.params["phi"].unit.to(u.rad)
        return self.params["a"](wl, t)*np.cos(np.arctan2(yy, xx)-phi)

    def _image(self, xx: np.ndarray, yy: np.ndarray,
               wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        rin, rout = map(lambda x: self.params[x](wl, t), ("rin", "rout"))
        r, p = np.sqrt(xx**2+yy**2), self.params["p"](wl, t)
        return np.nan_to_num(np.logical_and(r > rin, r < rout).astype(int)*(r / rin)**p, nan=0)

    def _imageFunction(self, xx: np.ndarray, yy: np.ndarray,
                       wl: np.ndarray, t: np.ndarray) -> np.ndarray:
        """Calculates a 2D-image from a radial power law.

        Parameters
        ----------
        xx : numpy.ndarray
            The x-coordinate grid
        yy : numpy.ndarray
            The y-coordinate grid
        wl : numpy.ndarray
            Wavelengths.
        t : numpy.ndarray
            Times.

        Returns
        -------
        image : numpy.ndarray
        """
        if self.asymmetric:
            return self._image(xx, yy, wl, t) * \
                (1-self._azimuthal_modulation(xx, yy, wl, t))
        return self._image(xx, yy, wl, t)


class oimAsymRadialPowerLaw(oimRadialPowerLaw):
    """An asymmetrically azimuthally modulated 2D-radial power law distribution.

    It accounts for elongation and rotation, as well

    Parameters
    ----------
    rin : float
        Inner radius of the disk [mas].
    rout : float
        Outer radius of the disk [mas].
    pixSize: float
        Pixel size [mas/px].
    p : float
        Power-law exponent for the radial profile.
    a : float
        Azimuthal modulation amplitude.
    phi : float
        Azimuthal modulation angle [deg].
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
    _wl : numpy.ndarray
        Wavelengths.
    _t : numpy.array
        Times.
    """
    name = "Asymmetric Radial Power Law"
    shortname = "AsymRadPowerLaw"
    asymmetric = True
