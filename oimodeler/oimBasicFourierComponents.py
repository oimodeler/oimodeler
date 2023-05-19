# -*- coding: utf-8 -*-
"""
basic model-components defined in the Fourier plan
"""
import astropy.units as u
import numpy as np
from scipy.special import j0, j1, jv
from scipy.signal import convolve2d
from .oimComponent import oimComponentFourier
from .oimParam import oimParam, _standardParameters


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
        image = xx*0
        val = np.abs(xx)+np.abs(yy)
        idx = np.unravel_index(np.argmin(val), np.shape(val))
        image[idx] = 1
        return image


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
        vc = rho*0
        idx = np.where(rho == 0)[0]
        if np.size(idx) != 0:
            vc[idx] = 1
        return vc

    def _imageFunction(self, xx, yy, wl, t):
        return xx*0+1


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
        xx = np.pi*self.params["d"](wl, t)*self.params["d"].unit.to(u.rad)*rho
        return np.nan_to_num(np.divide(2*j1(xx), xx), nan=1)

    def _imageFunction(self, xx, yy, wl, t):
        return ((xx**2+yy**2) <= (self.params["d"](wl, t)/2)**2).astype(float)


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
        return np.exp(-1*(np.pi*self.params["fwhm"](wl, t) *
                          self.params["fwhm"].unit.to(u.rad)*rho)**2/(4*np.log(2)))

    def _imageFunction(self, xx, yy, wl, t):
        r2 = (xx**2+yy**2)
        return np.sqrt(4*np.log(2*self.params["fwhm"](wl, t))/np.pi) * \
            np.exp(-4*np.log(2)*r2/self.params["fwhm"](wl, t)**2)


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
        xx = np.pi*self.params["d"](wl, t)*self.params["d"].unit.to(u.rad)*rho
        return j0(xx)

    def _imageFunction(self, xx, yy, wl, t, minPixSize=None):
        r2 = (xx**2+yy**2)
        dx = np.max([np.abs(1.*(xx[0, 0, 0, 1]-xx[0, 0, 0, 0])),
                     np.abs(1.*(yy[0, 0, 1, 0]-yy[0, 0, 0, 0]))])
        return ((r2 <= (self.params["d"](wl, t)/2+dx)**2) &
                (r2 >= (self.params["d"](wl, t)/2)**2)).astype(float)


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
        xxin = np.pi*self.params["din"](wl, t) * \
            self.params["din"].unit.to(u.rad)*rho
        xxout = np.pi*self.params["dout"](wl, t) * \
            self.params["dout"].unit.to(u.rad)*rho

        fin = (self.params["din"](wl, t))**2
        fout = (self.params["dout"](wl, t))**2

        return np.nan_to_num(2*(j1(xxout)/xxout*fout-j1(xxin)/xxin*fin)/(fout-fin), nan=1)

    def _imageFunction(self, xx, yy, wl, t):

        r2 = (xx**2+yy**2)
        return ((r2 <= (self.params["dout"](wl, t)/2)**2) &
                (r2 >= (self.params["din"](wl, t)/2)**2)).astype(float)


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
    name = "Ring2"
    shortname = "R2"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["w"] = oimParam(**_standardParameters["d"])
        self.params["w"].name = "w"
        self.params["w"].description = "width of the ring"
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        d = self.params["d"](wl, t) * self.params["d"].unit.to(u.rad)
        w = self.params["w"](wl, t) * self.params["w"].unit.to(u.rad)

        xx = np.pi*(d+w/2)*rho
        dxx = np.pi*w*rho

        return j0(xx)*np.nan_to_num(np.divide(2*j1(dxx), dxx), nan=1)

    def _imageFunction(self, xx, yy, wl, t):

        r2 = (xx**2+yy**2)
        return (((r2 <= (self.params["d"](wl, t)/2+self.params["w"](wl, t)/2)**2) &
                (r2 >= (self.params["d"](wl, t)/2-self.params["w"](wl, t)/2)**2))).astype(float)


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

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = np.pi*self.params["d"](wl, t)*self.params["d"].unit.to(u.rad)*rho

        phi = (self.params["skwPa"](wl, t)-self.params["pa"](wl, t)) * \
            self.params["skwPa"].unit.to(u.rad) + np.arctan2(yp, xp)

        return np.nan_to_num(j0(xx)-1j*np.sin(phi)*j1(xx)*self.params["skw"](wl, t), nan=1)

    def _imageFunction(self, xx, yy, wl, t):
        r2 = (xx**2+yy**2)
        # dr=np.sqrt(np.abs(np.roll(r2,(-1,-1),(0,1))-np.roll(r2,(1,1),(0,1))))
        phi = np.arctan2(yy, xx) + self.params["skwPa"](wl, t) * \
            self.params["skwPa"].unit.to(u.rad)
        dx = np.abs(1*(xx[0, 0, 0, 1]-xx[0, 0, 0, 0]
                    + xx[0, 0, 1, 0]-xx[0, 0, 0, 0])*self.params["elong"](wl, t))
        # dx=np.abs(1*(xx[0,1]-xx[0,0]+xx[1,0]-xx[0,0]))*3
        F = 1+self.params["skw"](wl, t)*np.cos(phi)
        return ((r2 <= (self.params["d"](wl, t)/2+dx/2)**2) &
                (r2 >= (self.params["d"](wl, t)/2-dx/2)**2)).astype(float)*F


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
        xxin = np.pi*self.params["din"](wl, t) * \
            self.params["din"].unit.to(u.rad)*rho
        xxout = np.pi*self.params["dout"](wl, t) * \
            self.params["dout"].unit.to(u.rad)*rho

        xx = (xxin+xxout)/2
        xx2 = (xxout-xxin)/2

        phi = (self.params["skwPa"](wl, t)-self.params["pa"](wl, t)) * \
            self.params["skwPa"].unit.to(u.rad) + np.arctan2(yp, xp)

        res = (j0(xx)-1j*np.sin(phi)*j1(xx)*self.params["skw"](wl, t))  \
            * np.nan_to_num(np.divide(2*j1(xx2), xx2), nan=1)

        return res

    def _imageFunction(self, xx, yy, wl, t):
        r2 = (xx**2+yy**2)

        phi = (self.params["skwPa"](wl, t)-self.params["pa"](wl, t)) * \
            self.params["skwPa"].unit.to(u.rad) + np.arctan2(yy, xx)

        F = (1+self.params["skw"](wl, t)*np.sin(phi)) / \
            (1+self.params["skw"](wl, t))

        return ((r2 <= (self.params["dout"](wl, t)/2)**2) &
                (r2 >= (self.params["din"](wl, t)/2)**2)).astype(float)*F


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
    name = "Pseudo-Lorentzian"
    shortname = "LO"
    elliptic = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["fwhm"] = oimParam(**_standardParameters["fwhm"])
        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):

        xx = np.pi*self.params["fwhm"](wl, t) * \
            self.params["fwhm"].unit.to(u.rad)*rho
        return np.exp(-2*np.pi*xx/3**0.5)

    def _imageFunction(self, xx, yy, wl, t):
        r2 = (xx**2+yy**2)
        a = np.pi*self.params["fwhm"](wl, t) * \
            self.params["fwhm"].unit.to(u.mas)
        return a/(2*np.pi*3**0.5)*(a**2/3+r2)**(-1.5)


class oimELorent(oimLorentz):
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
    name = "Elliptical Pseudo-Lorentzian"
    shortname = "ELO"
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
    """
    name = "Linear Limb Darkened Disk "
    shortname = "LLDD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a"] = oimParam(name="a", value=0, description="Linear LDD coeff",
                                    unit=u.one, mini=-1, maxi=1)

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):

        xx = np.pi*self.params["d"](wl, t) * \
            self.params["d"].unit.to(u.rad)*rho

        a = self.params["a"](wl, t)

        c1 = 2*np.divide(j1(xx), xx)
        c2 = 1.5*(np.pi*2)**0.5*np.divide(jv(1.5, xx), xx**1.5)
        return np.nan_to_num((1-a)*c1+a*c2, nan=1)


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
    """
    # NOTE: From Domiciano de Souza 2003 (phd thesis)
    name = "Quadratic Limb Darkened Disk "
    shortname = "QLDD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params["d"] = oimParam(**_standardParameters["d"])
        self.params["a1"] = oimParam(name="a1", value=0, description="1st QLDD coeff",
                                     unit=u.one, mini=-1, maxi=1)
        self.params["a2"] = oimParam(name="a2", value=0, description="2nd QLDD coeff",
                                     unit=u.one, mini=-1, maxi=1)

        self._eval(**kwargs)

    def _visFunction(self, xp, yp, rho, wl, t):
        xx = np.pi*self.params["d"](wl, t) * \
            self.params["d"].unit.to(u.rad)*rho

        a1 = self.params["a1"](wl, t)
        a2 = self.params["a2"](wl, t)

        c1 = np.divide(j1(xx), xx)
        c2 = (np.pi/2)**0.5*np.divide(jv(1.5, xx), xx**1.5)
        c3 = 2*np.divide(jv(2., xx), xx**2.)
        s = (6-2*a1-a2)/12
        return np.nan_to_num(((1-a1-a2)*c1+(a1+2*a2)*c2-a2*c3)/s, nan=1)


class oimConvolutor(oimComponentFourier):
    """

    Parameters
    ----------
    component1: oimComponentFourier
        first fourier component of the convolution
    component2: oimComponentFourier
        first fourier component of the convolution
    """
    def __init__(self, component1, component2, **kwargs):
        super().__init__(**kwargs)

        self.component1 = component1
        self.component2 = component2
        self.name = "Convolution Component"
        self.shortname = "Conv"

        for key in component2.params:
            self.params["c2_"+key] = component2.params[key]

        for key in component1.params:
            self.params["c1_"+key] = component1.params[key]

        self._eval(**kwargs)

    def _visFunction(self,xp,yp,rho,wl,t):     
    
        return  self.component1.getComplexCoherentFlux(xp,yp,wl,t)* \
                self.component2._visFunction(xp,yp,rho,wl,t)
     
    def _imageFunction(self,xx,yy,wl,t):
            
        im1=self.component1._imageFunction(xx,yy,wl,t)
        im2=self.component2._imageFunction(xx,yy,wl,t)
        nt,nwl,nx,ny=im1.shape

        #TODO : no loop
        imConv=np.ndarray([nt,nwl,nx,ny])
        for iwl in range(nwl):
            for it in range(nt):
                imConv[it,iwl,:,:] = convolve2d(im1[it,iwl,:,:],im2[it,iwl,:,:],
                                     mode='same', boundary='fill', fillvalue=0)
                
        
        return imConv

