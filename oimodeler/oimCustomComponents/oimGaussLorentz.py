import astropy.units as u
import numpy as np


from ..oimBasicFourierComponents import oimEGauss
from ..oimParam import _standardParameters, oimParam


class oimGaussLorentz(oimEGauss):
    name = "Gauss-Lorentzian component"
    shorname = "gaussLor"

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["flor"] = oimParam(**_standardParameters["f"])
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        vis_gauss = super()._visFunction(xp, yp, rho, wl, t)
        xx = self.params["fwhm"](wl, t)*self.params["fwhm"].unit.to(u.rad)*rho/2
        vis_lor = np.exp(-2*np.pi*xx/np.sqrt(3))
        return (1-self.params["flor"](wl, t))*vis_gauss + self.params["flor"](wl, t)*vis_lor

    # def _imageFunction(self, xx, yy, wl, t):
    #     image_gauss = self.gauss._imageFunction(xx, yy, wl, t)
    #     image_lor = self.lor._imageFunction(xx, yy, wl, t)
    #     return (1-self.params["flor"](wl, t))*image_gauss + self.params["flor"](wl, t)*image_lor
