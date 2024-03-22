import astropy.units as u
import numpy as np
from scipy.special import j0, j1

from .oimGaussLorentz import oimGaussLorentz
from ..oimOptions import oimOptions
from ..oimParam import _standardParameters, oimParam


class oimStarHaloGaussLorentz(oimGaussLorentz):
    name = "Star and Halo component with Gauss-Lorentzian disk"
    shorname = "SHGL"
    elliptic = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["fs"] = oimParam(**_standardParameters["f"])
        self.params["fc"] = oimParam(**_standardParameters["f"])
        self.params["kc"] = oimParam(**_standardParameters["p"])
        self.params["kc"].description = "Photometric slope of the disk component"
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])

        self.params["ks"] = oimParam(**_standardParameters["p"])
        self.params["ks"].description = "Photometric slope of the star component"
        self.params["ks"].free = False
        self.params["wl0"] = oimParam(**_standardParameters["wl"])
        self.params["wl0"].description = "Central wavelength"
        self.params["wl0"].free = False

        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh, kc, ks = 1-(fs+fc), self.params["kc"](wl, t), self.params["ks"](wl, t)
        wavelength_ratio = self.params["wl0"](wl, t)/wl

        vis_star = fs*wavelength_ratio**ks
        divisor = (fh+fs)*wavelength_ratio**ks + fc*wavelength_ratio**kc
        vis_disk = super()._visFunction(xp, yp, rho, wl, t)
        vis_comp = fc*vis_disk*wavelength_ratio**kc
        return (vis_star+vis_comp)/divisor

    def _imageFunction(self, xx, yy, wl, t):
        ...


class oimStarHaloRing(oimGaussLorentz):
    name = "Star and Halo component with Gauss-Lorentzian ring convolved disk"
    shorname = "SHGLR"
    elliptic = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0, unit=u.au,
                                      description="Inner radius of the disk")
        self.params["a"] = oimParam(name="a", value=0, unit=u.one,
                                    description="Azimuthal modulation amplitude")
        self.params["phi"] = oimParam(name="phi", value=0, unit=u.deg,
                                      description="Azimuthal modulation angle")
        self.params["w"] = oimParam(name="w", value=0, unit=u.mas,
                                    description="Azimuthal modulation angle")

        self.params["fs"] = oimParam(**_standardParameters["f"])
        self.params["fc"] = oimParam(**_standardParameters["f"])
        self.params["kc"] = oimParam(**_standardParameters["p"])
        self.params["kc"].description = "Photometric slope of the disk component"

        self.params["ks"] = oimParam(**_standardParameters["p"])
        self.params["ks"].description = "Photometric slope of the star component"
        self.params["ks"].free = False
        self.params["wl0"] = oimParam(**_standardParameters["wl"])
        self.params["wl0"].description = "Central wavelength"
        self.params["wl0"].free = False

        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh, kc, ks = 1-(fs+fc), self.params["kc"](wl, t), self.params["ks"](wl, t)
        wavelength_ratio = self.params["wl0"](wl, t)/wl

        rin, width = self.params["rin"](wl, t), self.params["w"](wl, t)
        a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
        baseline_angle = np.arctan2(yp, xp)

        if oimOptions.model.grid.type == "linear":
            radius = np.linspace(rin, rin+width, self.params["dim"](wl, t))
        else:
            radius = np.logspace(0.0 if rin == 0 else np.log10(rin),
                                 np.log10(rin+width), self.params["dim"](wl, t))

        vis_ring = np.trapz(j0(2*np.pi*radius.to(u.rad)*rho), radius)/width
        factor_ring = -1j*a*np.cos(baseline_angle-phi*phi.unit.to(u.rad))
        modulation_ring = (factor_ring*np.trapz(j1(2*np.pi*radius.to(u.rad)*rho), radius))/width

        vis_gausslor = super()._visFunction(xp, yp, rho, wl, t)
        vis_disk = vis_gausslor * (vis_ring + modulation_ring)

        vis_star = fs*wavelength_ratio**ks
        divisor = (fh+fs)*wavelength_ratio**ks + fc*wavelength_ratio**kc
        vis_comp = fc*vis_disk*wavelength_ratio**kc
        return (vis_star+vis_comp)/divisor

    def _imageFunction(self, xx, yy, wl, t):
        ...
