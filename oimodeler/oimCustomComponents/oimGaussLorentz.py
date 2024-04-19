import astropy.units as u
import numpy as np

from ..oimComponent import oimComponentFourier
from ..oimParam import _standardParameters, oimParam


class oimGaussLorentz(oimComponentFourier):
    name = "Gauss-Lorentzian"
    shortname = "GL"
    elliptic = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self.params["flor"] = oimParam(**_standardParameters["f"])
        self.params["flor"].name = "flor"
        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        flor = self.params["flor"](wl, t)
        xx = np.pi*self.params["fwhm"](wl, t)*self.params["fwhm"].unit.to(u.rad)*rho
        return (1-flor)*np.exp(-xx**2/(4*np.log(2))) + flor*np.exp(-xx/np.sqrt(3))

    def _imageFunction(self, xx, yy, wl, t):
        fwhm, radius = self.params["fwhm"](wl, t), np.hypot(xx, yy)
        flor = self.params["flor"](wl, t)
        image_gauss = 4*np.log(2)/(np.pi*fwhm**2)*np.exp(-(radius/fwhm)**2*4*np.log(2))
        image_lor = fwhm/(4*np.pi*np.sqrt(3))*(fwhm**2/12+radius**2)**(-1.5)
        return (1-flor)*image_gauss + flor*image_lor
