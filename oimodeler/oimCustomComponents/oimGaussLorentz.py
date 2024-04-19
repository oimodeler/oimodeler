import astropy.units as u
import numpy as np

from ..oimBasicFourierComponents import oimEGauss
from ..oimParam import _standardParameters, oimParam


class oimGaussLorentz(oimEGauss):
    name = "Gauss-Lorentzian component"
    shorname = "gaussLor"
    elliptic = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["flor"] = oimParam(**_standardParameters["f"])
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        flor = self.params["flor"](wl, t)
        fwhm = self.params["fwhm"](wl, t)*self.params["fwhm"].unit.to(u.rad)
        vis_gauss = np.exp(-(np.pi*fwhm*rho)**2/(4*np.log(2)))
        vis_lor = np.exp(-np.pi*fwhm*rho/np.sqrt(3))
        return (1-flor)*vis_gauss + flor*vis_lor

    def _imageFunction(self, xx, yy, wl, t):
        fwhm, radius = self.params["fwhm"](wl, t), np.hypot(xx, yy)
        flor = self.params["flor"](wl, t)
        image_gauss = 4*np.log(2)/(np.pi*fwhm**2)*np.exp(-4*(radius/fwhm)**2*np.log(2))
        image_lor = fwhm/(4*np.pi*np.sqrt(3))*(fwhm**2/12+radius**2)**(-1.5)
        return (1-flor)*image_gauss + flor*image_lor
