import astropy.units as u
import numpy as np
from scipy.special import j0, j1
from scipy.signal import fftconvolve

from .oimGaussLorentz import oimGaussLorentz
from ..oimOptions import oimOptions
from ..oimParam import _standardParameters, oimParam


class oimStarHaloGaussLorentz(oimGaussLorentz):
    name = "Star and Halo component with Gauss-Lorentzian disk"
    shortname = "SHGL"

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

    # TODO: Maybe include the wavelength dependent ratios (with the exponents here)
    def _imageFunction(self, xx, yy, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh = 1-(fs+fc)
        val = np.abs(xx)+np.abs(yy)
        idx = np.unravel_index(np.argmin(val), np.shape(val))
        image_disk = super()._imageFunction(xx, yy, wl, t)
        image_disk[idx] += fs+fh
        return image_disk


class oimStarHaloIRing(oimGaussLorentz):
    name = "Star and Halo component with Gauss-Lorentzian ring convolved disk"
    shortname = "SHGLR"
    thin = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0, unit=u.mas,
                                      description="Inner radius of the disk")
        self.params["a"] = oimParam(name="a", value=0, unit=u.one, mini=0, maxi=1,
                                    description="Azimuthal modulation amplitude")
        self.params["phi"] = oimParam(name="phi", value=0, unit=u.deg, mini=-180, maxi=180, 
                                      description="Azimuthal modulation angle")
        self.params["w"] = oimParam(name="w", value=0, unit=u.mas,
                                    description="Azimuthal modulation angle")
        if self.thin:
            self.params["w"].free = False

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

        factor_ring = -1j*a*np.cos(baseline_angle-phi*self.params["phi"].unit.to(u.rad))
        if self.thin:
            xx = rin*self.params["rin"].unit.to(u.rad)*rho
            vis_ring = j0(2*np.pi*xx) + factor_ring*j1(2*np.pi*xx) 
        else:
            if oimOptions.model.grid.type == "linear":
                radius = np.linspace(rin, rin+width, self.params["dim"](wl, t))
            else:
                radius = np.logspace(0.0 if rin == 0 else np.log10(rin),
                                     np.log10(rin+width), self.params["dim"](wl, t))

            xx = radius*self.params["rin"].unit.to(u.rad)*rho
            vis_ring = np.trapz(j0(2*np.pi*xx) + factor_ring*j1(2*np.pi*xx), radius)/width

        vis_gausslor = super()._visFunction(xp, yp, rho, wl, t)
        vis_star = fs*wavelength_ratio**ks
        divisor = (fh+fs)*wavelength_ratio**ks + fc*wavelength_ratio**kc
        vis_comp = fc*vis_gausslor*vis_ring*wavelength_ratio**kc
        return (vis_star+vis_comp)/divisor

    # TODO: Maybe include the wavelength dependent ratios (with the exponents here)
    # FIXME:Right now not quite correctly implemented
    def _imageFunction(self, xx, yy, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh, rin = 1-(fs+fc), self.params["rin"](wl, t)
        a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
        phi *= self.params["phi"].unit.to(u.rad)
        c, s, polar_angle = a*np.cos(phi), a*np.sin(phi), np.arctan2(yy, xx)

        radius, val = np.hypot(xx, yy), np.abs(xx)+np.abs(yy)
        idx = np.unravel_index(np.argmin(val), np.shape(val))
        image_gausslor = super()._imageFunction(xx, yy, wl, t)

        dx = np.max([np.abs(1.*(xx[0, 0, 0, 1]-xx[0, 0, 0, 0])),
                     np.abs(1.*(yy[0, 0, 1, 0]-yy[0, 0, 0, 0]))])

        radial_profile = (radius <= rin) & (radius >= rin+dx)
        image_ring = 1/(2*np.pi)*radius*radial_profile
        image_ring *= 1 + c*np.cos(polar_angle) + s*np.sin(polar_angle)
        image_disk = fftconvolve(image_gausslor, image_ring, mode="same")
        image_disk[idx] += fs+fh
        return image_disk


class oimStarHaloRing(oimGaussLorentz):
    name = "Star and Halo component with Gauss-Lorentzian ring convolved disk"
    shortname = "SHGLR"
    thin = False
