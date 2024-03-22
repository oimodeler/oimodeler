import astropy.units as u
import numpy as np
from scipy.special import j0, j1

from ..oimComponent import oimComponentFourier
from ..oimOptions import oimOptions
from ..oimParam import oimParam


class oimAEIRing(oimComponentFourier):
    name = "Asymmetrical Elliptical Infinitesimal Ring"
    shorname = "AEIR"
    elliptic = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["rin"] = oimParam(name="rin", value=0, unit=u.au,
                                      description="Inner radius of the disk")
        self.params["a"] = oimParam(name="a", value=0, unit=u.one,
                                    description="Azimuthal modulation amplitude")
        self.params["phi"] = oimParam(name="phi", value=0, unit=u.deg,
                                      description="Azimuthal modulation angle")
        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
        rin, psi = self.params["rin"](wl, t), np.arctan2(yp, xp)
        factor = -1j*a*np.cos(psi-phi*phi.unit.to(u.rad))
        modulation = factor*j1(2*np.pi*rho*rin*rin.unit.to(u.rad))
        return j0(2*np.pi*rin*rin.unit.to(u.rad)*rho) + modulation

    def _imageFunction(self, xx, yy, wl, t):
        a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
        c, s, polar_angle = a*np.cos(phi), a*np.sin(phi), np.arctan2(yy, xx)

        rin, radius = self.params["rin"](wl, t), np.hypot(xx, yy)
        dx = np.max([np.abs(1.*(xx[0, 0, 0, 1]-xx[0, 0, 0, 0])),
                     np.abs(1.*(yy[0, 0, 1, 0]-yy[0, 0, 0, 0]))])

        radial_profile = (radius <= rin) & (radius >= rin+dx)
        image = 1/(2*np.pi)*radius*radial_profile
        modulation = c*np.cos(polar_angle) + s*np.sin(polar_angle)
        return image*(1 + modulation)


class oimAERing(oimAEIRing):
    name = "Asymmetrical Elliptical Ring"
    shorname = "AER"
    elliptic = True

    def __init__(self, **kwargs):
        super(). __init__(**kwargs)
        self.params["w"] = oimParam(name="w", value=0, unit=u.mas,
                                    description="Azimuthal modulation angle")
        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        rin, width = self.params["rin"](wl, t), self.params["w"](wl, t)
        a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
        baseline_angle = np.arctan2(yp, xp)

        if oimOptions.model.grid.type == "linear":
            radius = np.linspace(rin, rin+width, self.params["dim"](wl, t))
        else:
            radius = np.logspace(0.0 if rin == 0 else np.log10(rin),
                                 np.log10(rin+width), self.params["dim"](wl, t))

        vis = np.trapz(j0(2*np.pi*radius.to(u.rad)*rho), radius)/width
        factor = -1j*a*np.cos(baseline_angle-phi*phi.unit.to(u.rad))
        modulation = (factor*np.trapz(j1(2*np.pi*radius.to(u.rad)*rho), radius))/width
        return vis + modulation

    def _imageFunction(self, xx, yy, wl, t):
        a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
        c, s, polar_angle = a*np.cos(phi), a*np.sin(phi), np.arctan2(yy, xx)

        rin, radius = self.params["rin"](wl, t), np.hypot(xx, yy)
        radial_profile = (radius <= rin/2) & (radius >= rin-self.params["w"](wl, t))
        image = 1/(2*np.pi)*radius*radial_profile
        modulation = c*np.cos(polar_angle) + s*np.sin(polar_angle)
        return image*(1 + modulation)
