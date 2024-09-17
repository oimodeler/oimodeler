import astropy.units as u
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import j0, j1

from ..oimComponent import oimComponentFourier
from ..oimParam import _standardParameters, oimParam


class oimStarHaloGaussLorentz(oimComponentFourier):
    name = "Star and Halo component with Gauss-Lorentzian disk"
    shortname = "SHGL"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._ar, self._ak = None, None
        self.ar, self.ak = None, None
        self.params["la"] = oimParam(name="la", value=0, unit=u.one, mini=-1, maxi=1.5,
                                     description="Logarithm of the half-light/-flux radius")
        self.params["flor"] = oimParam(name="flor", value=0, unit=u.one, mini=0, maxi=1,
                                       free=True, description="The Gauss to Lorentz ratio")
        self.params["fs"] = oimParam(name="fs", value=0, unit=u.one, mini=0, maxi=1,
                                     free=True, description="The star's flux ratio")
        self.params["fc"] = oimParam(name="fc", value=0, unit=u.one, mini=0, maxi=1,
                                     free=True, description="The component's flux ratio")
        self.params["kc"] = oimParam(name="kc", value=0, unit=u.one, mini=-10, maxi=10,
                                     free=True, description="Photometric slope of the disk component")
        self.params["ks"] = oimParam(name="ks", value=0, unit=u.one,
                                     free=False, description="Photometric slope of the star component")
        self.params["wl0"] = oimParam(name="wl0", value=0, unit=u.one, free=False,
                                      description="Central wavelength of the instrument")

        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    def _vis_gauss_lorentz(self, xp, yp, rho, wl, t):
        flor = self.params["flor"](wl, t)
        hlr = 10 ** self.params["la"](wl, t) if self.ak is None else self.ak
        xx = np.pi * hlr * u.mas.to(u.rad) * rho
        vis_gauss = np.exp(-(xx**2) / np.log(2))
        vis_lor = np.exp(-2 * xx / np.sqrt(3))
        return (1 - flor) * vis_gauss  + flor * vis_lor

    def _visFunction(self, xp, yp, rho, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh = 1 - (fs + fc)
        kc, ks = self.params["kc"](wl, t), self.params["ks"](wl, t)
        wavelength_ratio = self.params["wl0"](wl, t) / wl

        vis_star = fs * wavelength_ratio**ks
        vis_comp = fc * self._vis_gauss_lorentz(xp, yp, rho, wl, t) * wavelength_ratio**kc
        divisor = (fs + fh) * wavelength_ratio**ks + fc * wavelength_ratio**kc
        return (vis_star + vis_comp) / divisor

    def _image_gauss_lorentz(self, xx, yy, wl, t):
        flor = self.params["flor"](wl, t)
        hlr = 10 ** self.params["la"](wl, t) if self.ak is None else self.ak
        radius = np.hypot(xx, yy)
        image_gauss = np.log(2) / (np.pi * hlr**2) * np.exp(-(radius / hlr) ** 2 * np.log(2))
        image_lor = hlr / (2 * np.pi * np.sqrt(3)) * (hlr**2 / 3 + radius**2) ** (-3 / 2)
        return (1 - flor) * image_gauss + flor * image_lor

    def _imageFunction(self, xx, yy, wl, t):
        fh = 1 - (self.params["fs"](wl, t) + self.params["fc"](wl, t))
        val = np.abs(xx) + np.abs(yy)
        idx = np.unravel_index(np.argmin(val), np.shape(val))
        image_disk = self._image_gauss_lorentz(xx, yy, wl, t)
        image_disk[idx] += self.params["fs"](wl, t)
        return image_disk + fh


class oimStarHaloIRing(oimStarHaloGaussLorentz):
    name = "Star and Halo component with Gauss-Lorentzian ring convolved disk"
    shortname = "SHGLR"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["lkr"] = oimParam(name="lkr", value=0, unit=u.one, mini=-1, maxi=1,
                                      description="Logarithm of the kernel/ring radius")
        self.params["skw"] = oimParam(**_standardParameters["skw"])
        self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])

        self._t = np.array([0])
        self._wl = None
        self._eval(**kwargs)

    @property
    def ar(self):
        la, lkr = self.params["la"].value, self.params["lkr"].value
        if self._ar is None:
            self._ar = np.sqrt(10 ** (2 * la) / (1 + 10 ** (2 * lkr)))
        return self._ar

    @ar.setter
    def ar(self, value):
        self._ar = value

    @property
    def ak(self):
        la, lkr = self.params["la"].value, self.params["lkr"].value
        if self._ak is None:
            self._ak = np.sqrt(10 ** (2 * la) / (1 + 10 ** (-2 * lkr)))
        return self._ak

    @ak.setter
    def ak(self, value):
        self._ak = value

    def _visFunction(self, xp, yp, rho, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh = 1 - (fs + fc)
        kc, ks = self.params["kc"](wl, t), self.params["ks"](wl, t)
        wavelength_ratio = self.params["wl0"](wl, t) / wl

        skw, skwPa = self.params["skw"](wl, t), self.params["skwPa"](wl, t)
        skwPa = (skwPa+90) * self.params["skwPa"].unit.to(u.rad)
        baseline_angle = np.arctan2(yp, xp)

        xx = 2 * np.pi * self.ar * u.mas.to(u.rad) * rho
        vis_ring = j0(xx) + -1j * skw * np.cos(baseline_angle - skwPa) * j1(xx)
        vis_gauss_lor = self._vis_gauss_lorentz(xp, yp, rho, wl, t)
        vis_star = fs * wavelength_ratio**ks
        vis_comp = fc * vis_gauss_lor * vis_ring * wavelength_ratio**kc
        divisor = (fs + fh) * wavelength_ratio**ks + fc * wavelength_ratio**kc
        return (vis_star + vis_comp) / divisor

    def _imageFunction(self, xx, yy, wl, t):
        fs, fc = self.params["fs"](wl, t), self.params["fc"](wl, t)
        fh = 1 - (fs + fc)
        skw, skwPa = self.params["skw"](wl, t), self.params["skwPa"](wl, t)
        skwPa = (skwPa+90) * self.params["skwPa"].unit.to(u.rad)
        c, s = skw * np.cos(skwPa), skw * np.sin(skwPa)
        polar_angle = np.arctan2(yy, xx)

        radius, val = np.hypot(xx, yy), np.abs(xx) + np.abs(yy)
        idx = np.unravel_index(np.argmin(val), np.shape(val))
        
        dx = np.max(
            [
                np.abs(1.0 * (xx[0, 0, 0, 1] - xx[0, 0, 0, 0])),
                np.abs(1.0 * (yy[0, 0, 1, 0] - yy[0, 0, 0, 0])),
            ]
        )

        radial_profile = (radius >= self.ar) & (radius <= (self.ar + dx))
        image_ring = 1 / (2 * np.pi) * radial_profile
        image_ring *= 1 + c * np.cos(polar_angle) + s * np.sin(polar_angle)
        image = fftconvolve(
            self._image_gauss_lorentz(xx, yy, wl, t), image_ring, mode="same")
        image[idx] += fs
        return image + fh
