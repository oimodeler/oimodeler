import astropy.units as u
import numpy as np
from scipy.special import j0, j1

from ..oimBasicFourierComponents import oimIRing, oimRing
from ..oimOptions import oimOptions
from ..oimParam import _standardParameters, oimParam


class oimAEIRing(oimIRing):
    name = "Asymmetrical Elliptical Infinitesimal Ring"
    shorname = "AEIR"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["skw"] = oimParam(**_standardParameters["skw"])
        self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])
        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        vis = super()._visFunction(xp, yp, rho, wl, t)
        skw, skwPa = self.params["skw"](wl, t), self.params["skwPa"](wl, t)
        skwPa *= self.params["skwPa"].unit.to(u.rad)
        d = self.params["d"](wl, t) * self.params["d"].unit.to(u.rad)
        return vis + -1j * skw * np.cos(np.arctan2(xp, yp) - skwPa) * j1(
            np.pi * d * rho
        )

    def _imageFunction(self, xx, yy, wl, t):
        img = super()._imageFunction(xx, yy, wl, t)
        skw, skwPa = self.params["skw"](wl, t), self.params["skwPa"](wl, t)
        skwPa *= self.params["skwPa"].unit.to(u.rad)
        c, s, polar_angle = (
            skw * np.sin(skwPa),
            skw * np.cos(skwPa),
            np.arctan2(xx, yy),
        )
        return img * (1 + c * np.cos(polar_angle) + s * np.sin(polar_angle))


# TODO: Adapt this to the oimRing
# class oimAERing(oimRing):
#     name = "Asymmetrical Elliptical Ring"
#     shorname = "AER"
#
#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
#         self.params["skw"] = oimParam(**_standardParameters["skw"])
#         self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])
#         self._t = np.array([0])
#         self._wl = None
#         self._eval(**kwargs)
#
#     def _visFunction(self, xp, yp, rho, wl, t):
#         rin, width = self.params["d"](wl, t) / 2, self.params["w"](wl, t)
#         skw, skwPa = self.params["skw"](wl, t), self.params["skwPa"](wl, t)
#         skwPa *= self.params["skwPa"].unit.to(u.rad)
#         baseline_angle = np.arctan2(xp, yp)
#
#         if oimOptions.model.grid.type == "linear":
#             radius = np.linspace(rin, rin + width, self.params["dim"](wl, t))
#         else:
#             radius = np.logspace(
#                 0.0 if rin == 0 else np.log10(rin),
#                 np.log10(rin + width),
#                 self.params["dim"](wl, t),
#             )
#
#         xx = 2 * np.pi * radius * self.params["rin"].unit.to(u.rad) * rho
#         return (
#             np.trapezoid(
#                 j0(xx) + -1j * skw * np.cos(baseline_angle - skwPa) * j1(xx),
#                 radius,
#             )
#             / width
#         )
#
#     def _imageFunction(self, xx, yy, wl, t):
#         a, phi = self.params["a"](wl, t), self.params["phi"](wl, t)
#         c, s, polar_angle = (
#             a * np.cos(phi),
#             a * np.sin(phi),
#             np.arctan2(yy, xx),
#         )
#
#         rin, radius = self.params["d"](wl, t) / 2, np.hypot(xx, yy)
#         radial_profile = (radius >= rin) & (
#             radius <= rin - self.params["w"](wl, t)
#         )
#         image = 1 / (2 * np.pi) * radius * radial_profile
#         return image * (1 + c * np.cos(polar_angle) + s * np.sin(polar_angle))
