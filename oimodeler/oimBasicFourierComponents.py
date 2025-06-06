# p -*- coding: utf-8 -*-
"""
Basic model-components defined in the Fourier plan
"""
import operator
from functools import reduce

import astropy.units as u
import numpy as np
from scipy.signal import convolve2d
from scipy.special import gamma, j0, j1, jn, jv

from .oimComponent import oimComponentFourier
from .oimParam import _standardParameters, oimParam


class oimPt(oimComponentFourier):
    """Point Source component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    """

    name = "Point source"
    shortname = "Pt"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        return 1

    def _imageFunction(self, xx, yy, wl, t):
        if len(xx.shape) != 1:
            image = xx * 0
            val = np.abs(xx) + np.abs(yy)
            nwl = xx.shape[1]
            nt = xx.shape[0]
            # TODO rewrite without loop
            for it in range(nt):
                for iwl in range(nwl):
                    val = np.abs(xx[it, iwl, :, :]) + np.abs(yy[it, iwl, :, :])
                    idx = np.unravel_index(np.argmin(val), np.shape(val))
                    image[it, iwl, idx[0], idx[1]] = 1
            return image
        else:
            return (xx == 0) & (yy == 0)


class oimBackground(oimComponentFourier):
    """Background component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    """

    name = "Background"
    shortname = "Bckg"

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
        self._eval(**kwargs)

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        vc = rho * 0
        idx = np.where(rho == 0)[0]
        if np.size(idx) != 0:
            vc[idx] = 1
        return vc

    def _imageFunction(self, xx, yy, wl, t):
        return xx * 0 + 1


class oimUD(oimComponentFourier):
    """Uniform Disk component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the disk (in mas). The default is 0.
    """

    name = "Uniform Disk"
    shortname = "UD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self._eval(**kwargs)

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )
        return np.nan_to_num(np.divide(2 * j1(xx), xx), nan=1)

    def _imageFunction(self, xx, yy, wl, t):
        return (
            (xx**2 + yy**2) <= (self.params["d"](wl, t) / 2) ** 2
        ).astype(float)


class oimEllipse(oimUD):
    """Uniform Ellipse component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        major-axis diameter of the ellipse (in mas). The default is 0.
    pa: u.deg | oimInterp
        position angle of the major axis of the ellipse (in deg).
        The default is 0.
    elong : u.dimensionless_unscaled | oimInterp
        elongation of the ellipse. The default is 1.
    """

    name = "Uniform Ellipse"
    shortname = "eUD"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)


class oimGauss(oimComponentFourier):
    """Gaussian Disk component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    fwhm: u.mas | oimInterp
        FWHM of the Gaussian (in mas). The default is 0.
    """

    name = "Gaussian Disk"
    shortname = "GD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            self.params["fwhm"](wl, t)
            * self.params["fwhm"].unit.to(u.rad)
            * rho
        )
        return np.exp(-((np.pi * xx) ** 2) / (4 * np.log(2)))

    def _imageFunction(self, xx, yy, wl, t):
        r2 = xx**2 + yy**2
        return (
            4
            * np.log(2)
            / self.params["fwhm"](wl, t) ** 2
            / np.pi
            * np.exp(-4 * np.log(2) * r2 / self.params["fwhm"](wl, t) ** 2)
        )


class oimEGauss(oimGauss):
    """Elliptical Gaussian component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    fwhm: u.mas | oimInterp
        FWHM of the Gaussian (in mas). The default is 0.
    pa: u.deg | oimInterp
        position angle of the major axis of the Gaussian (in deg).
        The default is 0.
    elong: u.dimensionless_unscaled | oimInterp
        elongation of the Gaussian. The default is 1.
    """

    name = "Gaussian Ellipse"
    shortname = "EG"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)


class oimIRing(oimComponentFourier):
    """Infinitesimal Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    """

    name = "Infinitesimal Ring"
    shortname = "IR"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )
        return j0(xx)

    def _imageFunction(self, xx, yy, wl, t, minPixSize=None):
        r2 = xx**2 + yy**2
        dx = np.max(
            [
                np.abs(1.0 * (xx[0, 0, 0, 1] - xx[0, 0, 0, 0])),
                np.abs(1.0 * (yy[0, 0, 1, 0] - yy[0, 0, 0, 0])),
            ]
        )
        return (
            (r2 <= (self.params["d"](wl, t) / 2 + dx) ** 2)
            & (r2 >= (self.params["d"](wl, t) / 2) ** 2)
        ).astype(float)


class oimEIRing(oimIRing):
    """Infinitesimal Elliptical Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    pa: u.deg | oimInterp
        position angle of the major axis of the ring (in deg).
        The default is 0.
    elong: u.dimensionless_unscaled | oimInterp
        elongation of the ring. The default is 1.
    """

    name = "Ellitical Infinitesimal Ring"
    shortname = "EIR"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)


class oimRing(oimComponentFourier):
    """Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    din: u.mas | oimInterp
        inner diameter of the ring (in mas). The default is 0.
    dout: u.mas | oimInterp
        outer diameter of the ring (in mas). The default is 0.
    """

    name = "Ring"
    shortname = "R"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xxin = (
            np.pi
            * self.params["din"](wl, t)
            * self.params["din"].unit.to(u.rad)
            * rho
        )
        xxout = (
            np.pi
            * self.params["dout"](wl, t)
            * self.params["dout"].unit.to(u.rad)
            * rho
        )

        fin = (self.params["din"](wl, t)) ** 2
        fout = (self.params["dout"](wl, t)) ** 2

        return np.nan_to_num(
            2
            * (j1(xxout) / xxout * fout - j1(xxin) / xxin * fin)
            / (fout - fin),
            nan=1,
        )

    def _imageFunction(self, xx, yy, wl, t):

        r2 = xx**2 + yy**2
        return (
            (r2 <= (self.params["dout"](wl, t) / 2) ** 2)
            & (r2 >= (self.params["din"](wl, t) / 2) ** 2)
        ).astype(float)


class oimRing2(oimComponentFourier):
    """Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    width: u.mas | oimInterp
        width of the ring (in mas). The default is 0.
    """

    name = "IRing convolved with UD"
    shortname = "R2"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["w"] = oimParam(**_standardParameters["w"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        d = self.params["d"](wl, t) * self.params["d"].unit.to(u.rad)
        w = self.params["w"](wl, t) * self.params["w"].unit.to(u.rad)

        xx = np.pi * (d) * rho
        dxx = np.pi * w * rho

        return j0(xx) * np.nan_to_num(np.divide(2 * j1(dxx), dxx), nan=1)

    def _imageFunction(self, xx, yy, wl, t):

        r2 = xx**2 + yy**2
        return (
            (
                (
                    r2
                    <= (
                        self.params["d"](wl, t) / 2
                        + self.params["w"](wl, t) / 2
                    )
                    ** 2
                )
                & (
                    r2
                    >= (
                        self.params["d"](wl, t) / 2
                        - self.params["w"](wl, t) / 2
                    )
                    ** 2
                )
            )
        ).astype(float)


class oimERing(oimRing):
    """Elliptical Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    width: u.mas | oimInterp
        width of the ring (in mas). The default is 0.
    pa: u.deg | oimInterp
        position angle of the major axis of the ring (in deg).
        The default is 0.
    elong: u.dimensionless_unscaled | oimInterp
        elongation of the ring. The default is 1.
    """

    name = "Elliptical Ring"
    shortname = "ER"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)


class oimERing2(oimRing2):
    """Elliptical Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    width: u.mas | oimInterp
        width of the ring (in mas). The default is 0.
    pa: u.deg | oimInterp
        position angle of the major axis of the ring (in deg).
        The default is 0.
    elong: u.dimensionless_unscaled | oimInterp
        elongation of the ring. The default is 1.
    """

    name = "Elliptical Ring2"
    shortname = "ER2"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)


class oimESKIRing(oimComponentFourier):
    """Skewed Elliptical Infinitesimal Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    skw: u.dimensionless_unscaled | oimInterp
        skew of the ring. The default is 0.
    skwPa: u.deg | oimInterp
        elongation of the ring. The default is 1.
    """

    name = "Skewed Elliptical Infinitesimal Ring"
    shortname = "SKEIR"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["skw"] = oimParam(**_standardParameters["skw"])
        self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])
        self._eval(**kwargs)

    # TODO change definition of skwPA
    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        phi = (
            self.params["skwPa"](wl, t) - self.params["pa"](wl, t)
        ) * self.params["skwPa"].unit.to(u.rad) + np.arctan2(yp, xp)

        return np.nan_to_num(
            j0(xx) - 1j * np.sin(phi) * j1(xx) * self.params["skw"](wl, t),
            nan=1,
        )

    def _imageFunction(self, xx, yy, wl, t):
        r2 = xx**2 + yy**2
        # dr=np.sqrt(np.abs(np.roll(r2,(-1,-1),(0,1))-np.roll(r2,(1,1),(0,1))))
        phi = np.arctan2(yy, xx) + self.params["skwPa"](wl, t) * self.params[
            "skwPa"
        ].unit.to(u.rad)
        dx = np.abs(
            1
            * (
                xx[0, 0, 0, 1]
                - xx[0, 0, 0, 0]
                + xx[0, 0, 1, 0]
                - xx[0, 0, 0, 0]
            )
            * self.params["elong"](wl, t)
        )
        # dx=np.abs(1*(xx[0,1]-xx[0,0]+xx[1,0]-xx[0,0]))*3
        F = 1 + self.params["skw"](wl, t) * np.cos(phi)
        return (
            (r2 <= (self.params["d"](wl, t) / 2 + dx / 2) ** 2)
            & (r2 >= (self.params["d"](wl, t) / 2 - dx / 2) ** 2)
        ).astype(float) * F


class oimESKGRing(oimComponentFourier):
    """Skewed Elliptical Gaussian Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    fwhm: u.mas | oimInterp
        outer diameter of the ring (in mas). The default is 0.
    skw: u.dimensionless_unscaled | oimInterp
        skew of the ring. The default is 0.
    skwPa: u.deg | oimInterp
        elongation of the ring. The default is 1.
    """

    name = "Skewed Elliptical Ring"
    shortname = "SKER"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self.params["skw"] = oimParam(**_standardParameters["skw"])
        self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        xx2 = (
            self.params["fwhm"](wl, t)
            * self.params["fwhm"].unit.to(u.rad)
            * rho
        )

        phi = (
            self.params["skwPa"](wl, t) - self.params["pa"](wl, t)
        ) * self.params["skwPa"].unit.to(u.rad) + np.arctan2(yp, xp)

        res = (
            j0(xx) - 1j * np.sin(phi) * j1(xx) * self.params["skw"](wl, t)
        ) * np.nan_to_num(
            np.exp(-((np.pi * xx2) ** 2) / (4 * np.log(2))), nan=1
        )

        return res


class oimESKRing(oimComponentFourier):
    """Skewed Elliptical Ring component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    din: u.mas | oimInterp
        inner diameter of the ring (in mas). The default is 0.
    dout: u.mas | oimInterp
        outer diameter of the ring (in mas). The default is 0.
    skw: u.dimensionless_unscaled | oimInterp
        skew of the ring. The default is 0.
    skwPa: u.deg | oimInterp
        elongation of the ring. The default is 1.
    """

    name = "Skewed Elliptical Ring"
    shortname = "SKER"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["din"] = oimParam(**_standardParameters["din"])
        self.params["dout"] = oimParam(**_standardParameters["dout"])
        self.params["skw"] = oimParam(**_standardParameters["skw"])
        self.params["skwPa"] = oimParam(**_standardParameters["skwPa"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xxin = (
            np.pi
            * self.params["din"](wl, t)
            * self.params["din"].unit.to(u.rad)
            * rho
        )
        xxout = (
            np.pi
            * self.params["dout"](wl, t)
            * self.params["dout"].unit.to(u.rad)
            * rho
        )

        xx = (xxin + xxout) / 2
        xx2 = (xxout - xxin) / 2

        phi = (
            self.params["skwPa"](wl, t) - self.params["pa"](wl, t)
        ) * self.params["skwPa"].unit.to(u.rad) + np.arctan2(yp, xp)

        res = (
            j0(xx) - 1j * np.sin(phi) * j1(xx) * self.params["skw"](wl, t)
        ) * np.nan_to_num(np.divide(2 * j1(xx2), xx2), nan=1)

        return res

    # def _imageFunction(self, xx, yy, wl, t):

    # r2 = (xx**2+yy**2)

    # phi = (self.params["skwPa"](wl, t)-self.params["pa"](wl, t)) * \
    #     self.params["skwPa"].unit.to(u.rad) + np.arctan2(yy, xx)

    # F = (1+self.params["skw"](wl, t)*np.sin(phi)) / \
    #     (1+self.params["skw"](wl, t))

    # return ((r2 <= (self.params["dout"](wl, t)/2)**2) &
    #         (r2 >= (self.params["din"](wl, t)/2)**2)).astype(float)*F


# TODO
class oimLorentz(oimComponentFourier):
    """Pseudo-Lorentzian component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    fwhm: u.mas | oimInterp
        FWHM of the Lorentzian (in mas). The default is 0.
    """

    # NOTE: From Lazareff 2017 A&A 599, 85
    # TODO : Small difference between images using direct formula or inverse of vis function
    name = "Pseudo Lorentzian"
    shortname = "LZ"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):

        # 1.13 factor was computed to transform the Half-light Radius originally
        # in Lazareff paper into FWHM

        xx = (
            self.params["fwhm"](wl, t)
            * self.params["fwhm"].unit.to(u.rad)
            * rho
            / 1.13
        )
        return np.exp(-2 * np.pi * xx / 3**0.5)

    def _imageFunction(self, xx, yy, wl, t):
        r2 = xx**2 + yy**2
        a = (
            self.params["fwhm"](wl, t)
            * self.params["fwhm"].unit.to(u.mas)
            / 1.13
        )
        return a / (2 * np.pi * 3**0.5) * (a**2 / 3 + r2) ** (-1.5)


class oimELorentz(oimLorentz):
    """Elliptical-Lorentzian component defined in the fourier space

    Parameters
    ----------
    x : u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y : u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f : u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    fwhm : u.mas | oimInterp
        FWHM of the Lorentzian (in mas). The default is 0.
    pa : u.deg | oimInterp
        position angle of the major axis of the Lorentzian (in deg).
        The default is 0.
    elong : u.dimensionless_unscaled | oimInterp
        elongation of the Lorentzian. The default is 1.
    """

    name = "Elliptical Pseudo Lorentzian"
    shortname = "ELZ"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._eval(**kwargs)


class oimLinearLDD(oimComponentFourier):
    """Linear Limb Darkened Disk component defined in the fourier space

    Parameters
    ----------
    x: u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y: u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f: u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    a: u.dimensionless_unscaled | oimInterp
        linear limb darkening coefficient

    I(mu)/I(1) = 1  - a(1-mu)
    """

    name = "Linear Limb Darkened Disk "
    shortname = "LLDD"

    # NOTE: From Domiciano de Souza 2003 (phd thesis) and 2021
    # https://www.aanda.org/articles/aa/pdf/2021/10/aa40478-21.pdf

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a"] = oimParam(
            name="a",
            value=0,
            description="Linear LDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):

        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        a = self.params["a"](wl, t)

        c1 = 2 * np.divide(j1(xx), xx)
        c2 = 1.5 * (np.pi * 2) ** 0.5 * np.divide(jv(1.5, xx), xx**1.5)
        return np.nan_to_num((1 - a) * c1 + a * c2, nan=1)


class oimQuadLDD(oimComponentFourier):
    """Quadratic Limb Darkened Disk component defined in the fourier space

    Parameters
    ----------
    x : u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y : u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f : u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    a1: u.dimensionless_unscaled | oimInterp
        first quadratic limb darkening coefficient
    a2: u.dimensionless_unscaled | oimInterp
        second quadratic limb darkening coefficient

    I(mu)/I(1) = 1  - a1(1-mu) - a2(1 - mu)**2
    """

    # NOTE: From Domiciano de Souza 2003 (phd thesis) and 2021
    # https://www.aanda.org/articles/aa/pdf/2021/10/aa40478-21.pdf
    name = "Quadratic Limb Darkened Disk "
    shortname = "QLDD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a1"] = oimParam(
            name="a1",
            value=0,
            description="1st QLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )
        self.params["a2"] = oimParam(
            name="a2",
            value=0,
            description="2nd QLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        a1 = self.params["a1"](wl, t)
        a2 = self.params["a2"](wl, t)

        c1 = np.divide(j1(xx), xx)
        c2 = (np.pi / 2) ** 0.5 * np.divide(jv(1.5, xx), xx**1.5)
        c3 = 2 * np.divide(jv(2.0, xx), xx**2.0)
        s = (6 - 2 * a1 - a2) / 12
        return np.nan_to_num(
            ((1 - a1 - a2) * c1 + (a1 + 2 * a2) * c2 - a2 * c3) / s, nan=1
        )


class oimPowerLawLDD(oimComponentFourier):
    """Power Law Limb Darkened Disk component defined in the fourier space
    Parameters
    ----------
    x : u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y : u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f : u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    a: u.dimensionless_unscaled | oimInterp
        power law limb darkening exponent

    I(mu)/I(1) = mu**a
    """

    # NOTE: From Domiciano de Souza 2003 (phd thesis) and 2021
    # https://www.aanda.org/articles/aa/pdf/2021/10/aa40478-21.pdf

    name = "Power Law Limb Darkened Disk "
    shortname = "PLLDD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a"] = oimParam(
            name="a",
            value=0,
            description="Power Law LDD coeff",
            unit=u.one,
            mini=0,
            maxi=3,
        )

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):

        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        a = self.params["a"](wl, t)

        nu = a / 2 + 1

        return np.nan_to_num(
            nu * gamma(nu) * 2**nu * jn(nu, xx) / xx**nu, nan=1
        )


class oimSqrtLDD(oimComponentFourier):
    """Square-root Limb Darkened Disk component defined in the fourier space

    Parameters
    ----------
    x : u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y : u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f : u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    a1: u.dimensionless_unscaled | oimInterp
        first square-root limb darkening coefficient
    a2: u.dimensionless_unscaled | oimInterp
        second square-root darkening coefficient

    I(mu)/I(1) = 1  - a1 (1-mu) - a2 (1 - sqrt(mu))
    """

    # NOTE: From Domiciano de Souza 2003 (phd thesis) and 2021
    # https://www.aanda.org/articles/aa/pdf/2021/10/aa40478-21.pdf
    name = "square-root Limb Darkened Disk "
    shortname = "SLDD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a1"] = oimParam(
            name="a1",
            value=0,
            description="1st SLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )
        self.params["a2"] = oimParam(
            name="a2",
            value=0,
            description="2nd SLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        a1 = self.params["a1"](wl, t)
        a2 = self.params["a2"](wl, t)

        c1 = np.divide(j1(xx), xx)
        c2 = gamma(5 / 2) * (2 ** (3 / 2)) * np.divide(jv(1.5, xx), xx**1.5)
        c3 = (gamma(9 / 4)) * (2**1.25) * np.divide(jv(1.25, xx), xx**1.25)
        s = (15 - 5 * a1 - 3 * a2) / 30
        return np.nan_to_num(
            ((1 - a1 - a2) * c1 + (2 * a1 / 6) * c2 + (4 * a2 / 10) * c3) / s,
            nan=1,
        )


class oim4CLDD(oimComponentFourier):
    """4 Coefficients Limb Darkened Disk component defined in the fourier space

    Parameters
    ----------
    x : u.mas | oimInterp
        x pos of the component (in mas). The default is 0.
    y : u.mas | oimInterp
        y pos of the component (in mas). The default is 0.
    f : u.dimensionless_unscaled | oimInterp
        flux of the component. The default is 1.
    d: u.mas | oimInterp
        diameter of the ring (in mas). The default is 0.
    a1: u.dimensionless_unscaled | oimInterp
        first 4 Coefficients limb darkening coefficient
    a2: u.dimensionless_unscaled | oimInterp
        second 4 Coefficients limb darkening coefficient
    a3: u.dimensionless_unscaled | oimInterp
        third 4 Coefficients limb darkening coefficient
    a4: u.dimensionless_unscaled | oimInterp
        forth 4 Coefficients limb darkening coefficient

    I(mu)/I(1) = 1 - a1(1-mu**0.5) - a2(1 - mu)
                   - a3(1 - mu**1.5) - a4(1 - mu**2)
    """

    # NOTE: From Domiciano de Souza 2003 (phd thesis) and 2021
    # https://www.aanda.org/articles/aa/pdf/2021/10/aa40478-21.pdf
    name = "4 Coefficients Limb Darkened Disk "
    shortname = "4CLDD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a1"] = oimParam(
            name="a1",
            value=0,
            description="1st 4CLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )
        self.params["a2"] = oimParam(
            name="a2",
            value=0,
            description="2nd 4CLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )
        self.params["a3"] = oimParam(
            name="a3",
            value=0,
            description="3rd 4CLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )
        self.params["a4"] = oimParam(
            name="a4",
            value=0,
            description="4th 4CLDD coeff",
            unit=u.one,
            mini=-1,
            maxi=1,
        )

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = (
            np.pi
            * self.params["d"](wl, t)
            * self.params["d"].unit.to(u.rad)
            * rho
        )

        a1 = self.params["a1"](wl, t)
        a2 = self.params["a2"](wl, t)
        a3 = self.params["a3"](wl, t)
        a4 = self.params["a4"](wl, t)

        c0 = np.divide(j1(xx), xx)
        c1 = (
            2
            / 5
            * (gamma(9 / 4))
            * (2**1.25)
            * np.divide(jv(1.25, xx), xx**1.25)
        )
        c2 = (np.pi / 2) ** 0.5 * np.divide(jv(1.5, xx), xx**1.5)
        c3 = (
            2
            / 7
            * (gamma(11 / 4))
            * (2**1.75)
            * np.divide(jv(1.75, xx), xx**1.75)
        )
        c4 = 2 * np.divide(jv(2.0, xx), xx**2.0)
        s = (210 - 42 * a1 - 70 * a2 - 90 * a3 - 105 * a4) / 420
        return np.nan_to_num(
            (
                (1 - a1 - a2 - a3 - a4) * c0
                + a1 * c1
                + a2 * c2
                + a3 * c3
                + a4 * c4
            )
            / s,
            nan=1,
        )


class oimConvolutor(oimComponentFourier):
    """Convolves two components.

    Parameters
    ----------
    component1 : oimComponentFourier
        first fourier component of the convolution
    component2 : oimComponentFourier
        first fourier component of the convolution

    Attributes
    ----------
    components : list of oimComponentFourier
        The components that are to be convolved.
    params : dict
        Dictionary with the convolutors parameters.

    Warnings
    -----
    This component overloads the methods "getComplexCoherentFlux" and "getImage"
    and does not use "_visFunction" and "_imageFunction", beware during usage.
    """

    name = "Convolution Component"
    shortname = "Conv"

    def __init__(
        self,
        component1: oimComponentFourier,
        component2: oimComponentFourier,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.components = [component1, component2]
        for index, component in enumerate(self.components, start=1):
            for key in component.params:
                self.params[f"c{index}_{key}"] = component.params[key]

        self._eval(**kwargs)

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        raise ValueError(
            f"vis function not implemented for {self.shortname}."
            " This component overloads the 'getComplexCoherentFlux' method."
        )

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        vcs = []
        for index, component in enumerate(self.components, start=1):
            fxp, fyp = ucoord.copy(), vcoord.copy()
            if component.elliptic:
                pa_rad = (self.params[f"c{index}_pa"](wl, t)) * self.params[
                    f"c{index}_pa"
                ].unit.to(u.rad)
                co, si = np.cos(pa_rad), np.sin(pa_rad)
                fxpt = (fxp * co - fyp * si) / self.params[f"c{index}_elong"](
                    wl, t
                )
                fypt = fxp * si + fyp * co
            else:
                fxpt, fypt = fxp, fyp

            vcs.append(
                self.params[f"c{index}_f"](wl, t)
                * component._visFunction(
                    fxpt, fypt, np.hypot(fxpt, fypt), wl, t
                )
            )

        return (
            self.params["f"](wl, t)
            * self._ftTranslateFactor(ucoord, vcoord, wl, t)
            * reduce(operator.mul, vcs)
        )

    def _imageFunction(self, xx, yy, wl, t):
        raise ValueError(
            f"image function not implemented for {self.shortname}."
            " This component overloads the 'getImage' method."
        )

    def getImage(self, dim, pixSize, wl=None, t=None):
        t, wl = np.array(t).flatten(), np.array(wl).flatten()
        nt, nwl = t.size, wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim, endpoint=False)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr * pixSize * dim).flatten()
        y_arr = (vy_arr * pixSize * dim).flatten()
        t_arr, wl_arr = t_arr.flatten(), wl_arr.flatten()
        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)

        images = []
        for index, component in enumerate(self.components, start=1):
            xp, yp = x_arr.copy(), y_arr.copy()
            if component.elliptic:
                pa_rad = (
                    self.params[f"c{index}_pa"](wl_arr, t_arr)
                ) * self.params[f"c{index}_pa"].unit.to(units.rad)

                co, si = np.cos(pa_rad), np.sin(pa_rad)
                xpt = (xp * co - yp * si) * self.params[f"c{index}_elong"](
                    wl_arr, t_arr
                )
                ypt = xp * si + yp * co
            else:
                xpt, ypt = xp, yp

            img = component._imageFunction(
                xpt.reshape(dims),
                ypt.reshape(dims),
                wl_arr.reshape(dims),
                t_arr.reshape(dims),
            )

            tot = np.sum(img, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        img[it, iwl] = (
                            img[it, iwl]
                            / tot[it, iwl]
                            * self.params[f"c{index}_f"](wli, ti)
                        )
            images.append(img)

        img = np.zeros_like(images[0].shape)
        for iwl in range(nwl):
            for it in range(nt):
                img[it, iwl] = convolve2d(
                    *[img[it, iwl] for img in images],
                    mode="same",
                    boundary="fill",
                    fillvalue=0,
                )

        return img
