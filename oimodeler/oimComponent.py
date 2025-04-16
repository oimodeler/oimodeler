# -*- coding: utf-8 -*-
"""Components defined in Fourier or image planes"""
from typing import Any

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy import units as units
from scipy import interpolate, integrate
from scipy.special import j0
from pathlib import Path

from . import __dict__ as oimDict
from .oimOptions import oimOptions
from .oimParam import (
    oimInterp,
    oimParam,
    oimParamInterpolator,
    _standardParameters,
    oimParamNorm,
    oimParamLinker,
)
from .oimUtils import (
    getWlFromFitsImageCube,
    pad_image,
    rebin_image,
)


# TODO: Move somewhere else
def getFourierComponents():
    """A function to get the list of all available components deriving from the
    oimComponentFourier class

    Returns
    -------
    res : list
        list of all available components deriving from the oimComponentFourier class.
    """
    fnames = dir()
    res = []
    for f in fnames:
        try:
            if issubclass(oimDict[f], oimComponentFourier):
                res.append(f)
        except:
            pass
    return res


# TODO: Should elliptical parameters be moved to here?
class oimComponent:
    """The OImComponent class is the parent abstract class for all types of
    components that can be added to a OImModel.

    It has a similar interface than the oimModel and allow to compute images
    (or image cubes fore wavelength dependent models or time dependent models)
    and complex coherentFluxes for a vector of u,v,wl, and t coordinates

    Variables
    ----------
    name:
        The name of the component
    shortname:
        Short name for the component
    description:
        Detailed description of the component
    params:
        The dictionary of the component parameters
    """

    name = "Generic component"
    shortname = "Gen comp"
    description = "This is the class from which all components derived"

    def __init__(self, **kwargs):
        """Create and initiliaze a new instance of the oimComponent class.

        All components have at least three parameters the position
        x and y and their flux f
        """
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static

        self.params = {}
        self.params["x"] = oimParam(**_standardParameters["x"])
        self.params["y"] = oimParam(**_standardParameters["y"])
        self.params["f"] = oimParam(**_standardParameters["f"])
        # self.params["dim"] = oimParam(**_standardParameters["dim"])
        self._eval(**kwargs)

    def _paramstr(self):
        txt = []
        for paramname, param in self.params.items():
            if isinstance(param, oimParam):
                if isinstance(param, oimParamInterpolator):
                    # TODO: Have a string for each oimParamInterpolator
                    txt.append(f"{param.name}={param.__class__.__name__}")
                else:
                    txt.append(f"{param.name}={param.value:.2f}")
            elif isinstance(param, oimParamNorm) or isinstance(
                param, oimParamLinker
            ):
                txt.append(f"{paramname}={param.__class__.__name__}")

        return " ".join(txt)

    def __str__(self):
        return self.name + ": " + self._paramstr()

    def __repr__(self):
        return (
            f"{self.__class__.__name__} at "
            f"{str(hex(id(self)))}: {self._paramstr()}"
        )

    @property
    def _wl(self) -> np.ndarray:
        """Gets the wavelengths."""
        return self.__wl

    @_wl.setter
    def _wl(self, value: Any) -> np.ndarray:
        """Sets the wavelengths."""
        if value is None:
            self.__wl = None
            return

        if isinstance(value, (np.ndarray, tuple, list)):
            value = value
        elif isinstance(value, u.Quantity):
            if not isinstance(value.value, (np.ndarray, tuple, list)):
                value = [value]
        else:
            value = [value]
        self.__wl = np.array(value)

    @property
    def _t(self) -> np.ndarray:
        """Gets the times."""
        return self.__t

    @_t.setter
    def _t(self, value: Any) -> np.ndarray:
        """Sets the times."""
        if value is None:
            self.__t = None
            return

        if isinstance(value, (np.ndarray, tuple, list)):
            value = value
        elif isinstance(value, u.Quantity):
            if not isinstance(value.value, (np.ndarray, tuple, list)):
                value = [value]
        else:
            value = [value]
        self.__t = np.array(value)

    def _eval(self, **kwargs):
        for key, value in kwargs.items():
            if key in self.params.keys():
                if isinstance(value, oimInterp):
                    if not isinstance(self.params[key], oimParamInterpolator):
                        self.params[key] = value.type(
                            self.params[key], **value.kwargs
                        )
                else:
                    self.params[key].value = value

    def getComplexCoherentFlux(self, u, v, wl=None, t=None):
        """Compute and return the complex coherent flux for an array of u,v
        (and optionally wavelength and time ) coordinates

        Parameters
        ----------
        u : list or numpy array
            spatial coordinate u (in cycles/rad)
        v : list or numpy array
            spatial coordinate vu (in cycles/rad) .
        wl : list or numpy array, optional
            wavelength(s) in meter. The default is None.
        t :  list or numpy array, optional
            time in s (mjd). The default is None.

        Returns
        -------
        A numpy array of  the same size as u & v
            The complex coherent flux.
        """
        return np.array(u) * 0

    def getImage(self, dim, pixSize, wl=None, t=None):
        """Compute and return an image or and image cube (if wavelength and time
        are given).

        The returned image as the x,y dimension dim in pixel with
        an angular pixel size pixSize in rad. Image is returned as a numpy
        array unless the keyword fits is set to True. In that case the image is
        returned as an astropy.io.fits hdu

        Parameters
        ----------
        dim : integer
            image x & y dimension in pixels
        pixSize : float
            pixel angular size in rad
        wl : integer, list or numpy array, optional
             wavelength(s) in meter. The default is None
        t :  integer, list or numpy array, optional
            time in s (mjd). The default is None
        fits : bool, optional
            if True returns result as a fits hdu. The default is False

        Returns
        -------
        a numpy 2D array (or 3 or 4D array if wl, t or both are given) or an
        astropy.io.fits hdu. image hdu if fits=True.
            The image of the component with given size in pixels and rad
        """
        return np.zeros((dim, dim))

    def _ftTranslateFactor(self, ucoord, vcoord, wl, t):
        x = self.params["x"](wl, t) * self.params["x"].unit.to(units.rad)
        y = self.params["y"](wl, t) * self.params["y"].unit.to(units.rad)
        return np.exp(-2 * 1j * np.pi * (ucoord * x + vcoord * y))

    def _directTranslate(self, x, y, wl, t):
        x = x - self.params["x"](wl, t)
        y = y - self.params["y"](wl, t)
        return x, y
    
    def getNonRegularImage(self, xx, yy, wl=None, t=None):
        """ Compute and return a non-regular image function at the xx, yy and 
        optional wl and t coordinates)"""
    
        return 0*xx


class oimComponentFourier(oimComponent):
    """Class for all component analytically defined in the Fourier plan.
    Inherit from the oimComponent.

    Implements translation in direct and Fourier space, getImage from the
    Fourier definition of the object, ellipticity (i.e. flatening)
    Children classes should only implement the _visFunction and _imageFunction
    functions.
    """

    elliptic = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # NOTE: Add ellipticity if either elong or pa is specified in kwargs
        if ("elong" in kwargs) or ("pa" in kwargs) or self.elliptic == True:
            self.params["elong"] = oimParam(**_standardParameters["elong"])
            self.params["pa"] = oimParam(**_standardParameters["pa"])
            self.elliptic = True

        self._eval(**kwargs)

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        if self.elliptic:
            pa_rad = (self.params["pa"](wl, t)) * self.params["pa"].unit.to(
                units.rad
            )
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = (ucoord * co - vcoord * si) / self.params["elong"](wl, t)
            fyp = ucoord * si + vcoord * co
        else:
            fxp, fyp = ucoord, vcoord

        vc = self._visFunction(fxp, fyp, np.hypot(fxp, fyp), wl, t)
        return (
            vc
            * self._ftTranslateFactor(ucoord, vcoord, wl, t)
            * self.params["f"](wl, t)
        )

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        return ucoord * 0

    def getImage(self, dim, pixSize, wl=None, t=None):
        t = np.array(t).flatten()
        nt = t.size
        wl = np.array(wl).flatten()
        nwl = wl.size

        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim, endpoint=False)

        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr * pixSize * dim).flatten()
        y_arr = (vy_arr * pixSize * dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)

        if self.elliptic:

            pa_rad = (self.params["pa"](wl_arr, t_arr)) * self.params[
                "pa"
            ].unit.to(units.rad)

            xp = x_arr * np.cos(pa_rad) - y_arr * np.sin(pa_rad)
            yp = x_arr * np.sin(pa_rad) + y_arr * np.cos(pa_rad)

            x_arr = xp * self.params["elong"](wl_arr, t_arr)
            y_arr = yp

        image = self._imageFunction(
            x_arr.reshape(dims),
            y_arr.reshape(dims),
            wl_arr.reshape(dims),
            t_arr.reshape(dims),
        )

        tot = np.sum(image, axis=(2, 3))

        # TODO: No loop for normalization
        for it, ti in enumerate(t):
            for iwl, wli in enumerate(wl):
                if tot[it, iwl] != 0:
                    image[it, iwl, :, :] = (
                        image[it, iwl, :, :]
                        / tot[it, iwl]
                        * self.params["f"](wli, ti)
                    )
        return image
    
    
    def getNonRegularImage(self, xx, yy, wl=None, t=None):
        
        xx, yy = self._directTranslate(xx, yy, wl, t)

        if self.elliptic:

            pa_rad = (self.params["pa"](wl, t)) * \
                self.params["pa"].unit.to(units.rad)

            xp = xx*np.cos(pa_rad)-yy*np.sin(pa_rad)
            yp = xx*np.sin(pa_rad)+yy*np.cos(pa_rad)

            xx = xp*self.params["elong"](wl, t)
            yy = yp

        image = self._imageFunction(xx, yy, wl, t)

        """
        tot = np.sum(image, axis=(2, 3))

        # TODO: No loop for normalization
        
        for it, ti in enumerate(t):
            for iwl, wli in enumerate(wl):
                if tot[it, iwl] != 0:
                    image[it, iwl, :, :] = image[it, iwl, :, :]  \
                        / tot[it, iwl]*self.params["f"](wli, ti)
        """
        return image
        

    def _imageFunction(self, xx, yy, wl, t):
        raise ValueError(
            f"image function not implemented for {self.shortname}."
            "Use the fromFT=True option to get a model image"
            " from the inverse Fourier Transform"
        )


class oimComponentImage(oimComponent):
    """Base class for components define in 2D : x,y (regular grid) in the image plan"""

    elliptic = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._pixSize = 0  # NOTE: In rad
        self._allowExternalRotation = True
        self.normalizeImage = True
        self.params["pa"] = oimParam(**_standardParameters["pa"])
        self.params["dim"] = oimParam(**_standardParameters["dim"])

        # NOTE: Add ellipticity
        if self.elliptic:
            self.params["elong"] = oimParam(**_standardParameters["elong"])

        if "FTBackend" in kwargs:
            self.FTBackend = kwargs["FTBackend"]()
        else:
            self.FTBackend = oimOptions.ft.backend.active()

        self.FTBackendData = None
        self._eval(**kwargs)

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        if wl is None:
            wl = ucoord * 0
        if t is None:
            t = ucoord * 0

        im = self.getInternalImage(wl, t)

        if oimOptions.ft.binning is not None:
            im = rebin_image(im, oimOptions.ft.binning)
 
        im = pad_image(im)
        
        if self._pixSize!=0 :
            pix = self._pixSize
        else:
            pix = self.getPixelSize()

        tr = self._ftTranslateFactor(
            ucoord, vcoord, wl, t
        )  # ♣*self.params["f"](wl, t)

        if self._allowExternalRotation == True:
            pa_rad = (self.params["pa"](wl, t)) * self.params["pa"].unit.to(
                units.rad
            )
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = ucoord * co - vcoord * si
            fyp = ucoord * si + vcoord * co
            vcoord = fyp

            if self.elliptic:
                ucoord = fxp / self.params["elong"](wl, t)
            else:
                ucoord = fxp

        if self._wl is None:
            wl0 = np.sort(np.unique(wl))
        else:
            wl0 = self._wl

        if self._t is None:
            t0 = np.sort(np.unique(t))
        else:
            t0 = self._t

        if (
            self.FTBackend.check(
                self.FTBackendData, im, pix, wl0, t0, ucoord, vcoord, wl, t
            )
            == False
        ):

            self.FTBackendData = self.FTBackend.prepare(
                im, pix, wl0, t0, ucoord, vcoord, wl, t
            )

        vc = self.FTBackend.compute(
            self.FTBackendData, im, pix, wl0, t0, ucoord, vcoord, wl, t
        )

        return vc * tr * self.params["f"](wl, t)

    def getImage(self, dim, pixSize, wl=None, t=None):
        if wl is None:
            wl = 0
        if t is None:
            t = 0

        t = np.array(t).flatten()
        nt = t.size
        wl = np.array(wl).flatten()
        nwl = wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim, endpoint=False)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr * pixSize * dim).flatten()
        y_arr = (vy_arr * pixSize * dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)

        if self._allowExternalRotation == True:
            pa_rad = (self.params["pa"](wl_arr, t_arr)) * self.params[
                "pa"
            ].unit.to(units.rad)

            xp = x_arr * np.cos(pa_rad) - y_arr * np.sin(pa_rad)
            yp = x_arr * np.sin(pa_rad) + y_arr * np.cos(pa_rad)
            y_arr = yp

            if self.elliptic == True:
                x_arr = xp * self.params["elong"](wl_arr, t_arr)
            else:
                x_arr = xp

        im0 = self._internalImage()

        if im0 is None:
            im = self._imageFunction(x_arr, y_arr, wl_arr, t_arr)
        else:
            im0 = np.swapaxes(im0, -2, -1)
            grid = self._getInternalGrid()
            coord = np.transpose(np.array([t_arr, wl_arr, x_arr, y_arr]))

            im = interpolate.interpn(
                grid, im0, coord, bounds_error=False, fill_value=0)
            f0 = np.sum(im0)
            f = np.sum(im)
            im = im / f * f0

        im = im.reshape(dims)

        if self.normalizeImage == True:
            # TODO: No loop for normalization
            tot = np.sum(im, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        im[it, iwl, :, :] = (
                            im[it, iwl, :, :]
                            / tot[it, iwl]
                            * self.params["f"](wli, ti)
                        )
        return im

    def getInternalImage(self, wl, t):
        res = self._internalImage()

        if res is None:
            t_arr, wl_arr, x_arr, y_arr = self._getInternalGrid(
                simple=False, wl=wl, t=t
            )
            res = self._imageFunction(x_arr, y_arr, wl_arr, t_arr)

        if self.normalizeImage == True:
            for it in range(res.shape[0]):
                for iwl in range(res.shape[1]):
                    res[it, iwl, :, :] = res[it, iwl, :, :] / np.sum(
                        res[it, iwl, :, :]
                    )

        return res

    def _internalImage(self):
        return

    def _imageFunction(self, xx, yy, wl, t):
        image = xx * 0 + 1
        return image

    def _getInternalGrid(self, simple=True, flatten=False, wl=None, t=None):
        if self._wl is None:
            wl0 = np.sort(np.unique(wl))
        else:
            wl0 = self._wl

        if self._t is None:
            t0 = np.sort(np.unique(t))
        else:
            t0 = self._t

        dim = self.params["dim"](wl, t)

        
        #pix = self._pixSize * units.rad.to(units.mas)
        
        if self._pixSize!=0 :
            pix = self._pixSize* units.rad.to(units.mas)
        else:
            pix = self.getPixelSize()* units.rad.to(units.mas)
        
        
        v = np.linspace(-0.5, 0.5, dim)
        xy = v * pix * dim

        if simple:
            return t0, wl0, xy, xy

        else:
            t = np.array(t0).flatten()
            nt = t.size
            wl = np.array(wl0).flatten()
            nwl = wl.size
            xx, yy = np.meshgrid(xy, xy)
            x_arr = np.tile(xx[None, None, :, :], (nt, nwl, 1, 1))
            y_arr = np.tile(yy[None, None, :, :], (nt, nwl, 1, 1))
            wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
            t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

            if flatten == True:
                return (
                    t_arr.flatten(),
                    wl_arr.flatten(),
                    x_arr.flatten(),
                    y_arr.flatten(),
                )
            else:
                return t_arr, wl_arr, x_arr, y_arr

    def getPixelSize(self):
        raise ValueError("setPixelSize Method not implemented"
                         " while self._pixSize = None")

class oimComponentRadialProfile(oimComponent):
    """Base class for components define by their radial profile"""

    elliptic = False
    asymmetric = False
    modulation_order = 1

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._r = None
        self._wl = None  # None value <=> All wavelengths (from Data)
        self._t = [0]  # This component is static
        self.normalizeImage = True
        self.precision = None  # Precision for the Hankel transform

        # CHECK: Is this not redundant as oimComponent is already ellpitical?
        # NOTE: Add ellipticity
        if self.elliptic:
            self.params["pa"] = oimParam(**_standardParameters["pa"])
            self.params["elong"] = oimParam(**_standardParameters["elong"])

        # NOTE: Add asymmetry
        if self.asymmetric:
            for i in range(1, self.modulation_order + 1):
                self.params[f"skw{i}"] = oimParam(**_standardParameters["skw"])
                self.params[f"skwPa{i}"] = oimParam(
                    **_standardParameters["skwPa"]
                )

        self.params["dim"] = oimParam(**_standardParameters["dim"])

        self._eval(**kwargs)

    def _getInternalGrid(self, simple=True, flatten=False, wl=None, t=None):

        wl0 = np.sort(np.unique(wl)) if self._wl is None else self._wl
        t0 = np.sort(np.unique(t)) if self._t is None else self._t

        if self._r is None:
            pix = self._pixSize * units.rad.to(units.mas)
            r = (
                np.linspace(
                    0, self.params["dim"].value - 1, self.params["dim"].value
                )
                * pix
            )
        else:
            r = self._r

        if simple:
            return r, wl, t
        else:
            nt = np.array(t0).flatten().size
            nwl = np.array(wl0).flatten().size
            nr = r.flatten().size

            r_arr = np.tile(r[None, None, :], (nt, nwl, 1))
            wl_arr = np.tile(wl[None, :, None], (nt, 1, nr))
            t_arr = np.tile(t[:, None, None], (1, nwl, nr))

            if flatten:
                return t_arr.flatten(), wl_arr.flatten(), r_arr.flatten()
            else:
                return t_arr, wl_arr, r_arr

    def _internalRadialProfile(self):
        return None

    def _radialProfileFunction(self, r=None, wl=None, t=None):
        return 0

    def getInternalRadialProfile(self, wl, t):
        res = self._internalRadialProfile()
        if res is None:
            t_arr, wl_arr, r_arr = self._getInternalGrid(
                simple=False, wl=wl, t=t
            )
            res = self._radialProfileFunction(r_arr, wl_arr, t_arr)
        return res

    # FIXME: Not working!
    # @staticmethod
    # def fht(Ir, r, wlin, tin, sfreq, wl, t):
    #     nfreq, nwl = len(sfreq), len(wlin)
    #     r, Ir = (
    #         r[np.newaxis, np.newaxis, :, np.newaxis],
    #         Ir[:, :, :, np.newaxis],
    #     )
    #     sfreq = sfreq.reshape(nfreq // nwl, nwl).T[
    #         np.newaxis, :, np.newaxis, :
    #     ]
    #     num_hankel = (
    #         2.0
    #         * np.pi
    #         * integrate.trapezoid(
    #             r * Ir * j0(2 * np.pi * r * sfreq), r, axis=2
    #         )
    #     )
    #     norm = integrate.trapezoid(2 * np.pi * r * Ir, r, axis=2)
    #     return (num_hankel / norm).T.reshape(nfreq).astype(complex)

    # TODO: Convert this to non-statimethod for the asymmetric case
    @staticmethod
    def hankel(Ir, r, wlin, tin, sfreq, wl, t, precision=None):
        if precision is None:
            sfreq0 = np.unique(sfreq)
        else:
            sfreq0 = np.linspace(0, np.max(sfreq), num=precision)
        r1D = r[np.newaxis, np.newaxis, :]
        r2D = r[np.newaxis, np.newaxis, :, np.newaxis]
        Ir2D = Ir[:, :, :, np.newaxis]
        sf2D = sfreq0[np.newaxis, np.newaxis, np.newaxis, :]

        res0 = integrate.trapezoid(
            2.0 * np.pi * r2D * Ir2D * j0(2.0 * np.pi * r2D * sf2D),
            r2D,
            axis=2,
        )

        flux = 2.0 * np.pi * integrate.trapezoid(r1D * Ir, r1D, axis=2)
        flux_r = flux[:, :, np.newaxis]
        res0 /= flux_r

        grid = (tin, wlin, sfreq0)
        coord = np.transpose([t, wl, sfreq])

        real = interpolate.interpn(
            grid, np.real(res0), coord, bounds_error=False, fill_value=None
        )
        imag = interpolate.interpn(
            grid, np.imag(res0), coord, bounds_error=False, fill_value=None
        )
        return real + imag * 1j, flux

    def getImage(self, dim, pixSize, wl=None, t=None):
        wl, t = 0 if wl is None else wl, 0 if t is None else t
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
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)
        if self.elliptic:
            pa_rad = (self.params["pa"](wl_arr, t_arr)) * self.params[
                "pa"
            ].unit.to(units.rad)

            xp = x_arr * np.cos(pa_rad) - y_arr * np.sin(pa_rad)
            yp = x_arr * np.sin(pa_rad) + y_arr * np.cos(pa_rad)
            y_arr = yp
            x_arr = xp * self.params["elong"](wl_arr, t_arr)

        r_arr = np.hypot(x_arr, y_arr)
        im = self._radialProfileFunction(r_arr, wl_arr, t_arr)
        im = im.reshape(dims)

        if self.normalizeImage:
            # TODO: No loop for normalization
            tot = np.sum(im, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        im[it, iwl, :, :] = (
                            im[it, iwl, :, :]
                            / tot[it, iwl]
                            * self.params["f"](wli, ti)
                        )

        return im

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):

        wl = ucoord * 0 if wl is None else wl
        t = ucoord * 0 if t is None else t

        if self.elliptic:
            pa_rad = (self.params["pa"](wl, t)) * self.params["pa"].unit.to(
                units.rad
            )
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = ucoord * co - vcoord * si
            fyp = ucoord * si + vcoord * co
            vcoord, ucoord = fyp, fxp / self.params["elong"](wl, t)

        spf = np.hypot(ucoord, vcoord)
        psi = np.arctan2(vcoord, ucoord) if self.asymmetric else None
        wl0 = np.sort(np.unique(wl)) if self._wl is None else self._wl
        t0 = np.sort(np.unique(t)) if self._t is None else self._t

        Ir = self.getInternalRadialProfile(wl0, t0)
        vc, ftot = self.hankel(
            Ir,
            self._r * units.mas.to(units.rad),
            wl0,
            t0,
            spf,
            wl,
            t,
            precision=self.precision,
        )
        nwl0 = np.size(wl0)
        ftot = ftot.reshape(nwl0)
        ftot_Jy_interp = np.interp(wl, wl0, ftot * 1e23)
        if self.shortname == "TempGrad":
            return (
                vc
                * self._ftTranslateFactor(ucoord, vcoord, wl, t)
                * ftot_Jy_interp
            )
        else:
            return (
                vc
                * self._ftTranslateFactor(ucoord, vcoord, wl, t)
                * self.params["f"](wl, t)
            )


class oimComponentFitsImage(oimComponentImage):
    """Component load load images or chromatic-cubes from fits files"""

    elliptic = False
    name = "Fits Image Component"
    shortname = "Fits_Comp"

    def __init__(self, fitsImage=None, useinternalPA=False, **kwargs):
        super().__init__(**kwargs)
        if fitsImage:
            self.loadImage(fitsImage, useinternalPA=useinternalPA)
        self.params["pa"] = oimParam(**_standardParameters["pa"])
        self.params["scale"] = oimParam(**_standardParameters["scale"])

        self._eval(**kwargs)

    def loadImage(self, fitsImage, useinternalPA=False):
        if isinstance(fitsImage, str) or isinstance(fitsImage, Path):
            try:
                im = fits.open(fitsImage)[0]
            except:
                raise TypeError("Not a valid fits file")
        elif isinstance(fitsImage, fits.hdu.hdulist.HDUList):
            im = fitsImage[0]
        elif isinstance(fitsImage, fits.hdu.image.PrimaryHDU):
            im = fitsImage

        self._header = im.header

        dims = self._header["NAXIS"]
        if dims < 2:
            raise TypeError(
                "oimComponentFitsImage require 2D images or "
                "3D chromatic-image-cubes"
            )

        dimx = self._header["NAXIS1"]
        dimy = self._header["NAXIS2"]
        if dimx != dimy:
            raise TypeError("Current version only works with square images")
        self._dim = dimx

        pixX = self._header["CDELT1"]
        pixY = self._header["CDELT2"]
        if pixX != pixY:
            raise TypeError(
                "Current version only works with the same pixels"
                " scale in x and y dimension"
            )
        if "CUNIT1" in self._header:
            try:
                unit0 = units.Unit(self._header["CUNIT1"])
            except:
                unit0 = units.rad

        else:
            unit0 = units.rad
        self._pixSize0 = pixX * unit0.to(units.rad)

        if "CROTA1" in self._header:
            pa0 = self._header["CROTA1"]
        elif "CROTA2" in self._header:
            pa0 = self._header["CROTA2"]

        # TODO: Check units
        if useinternalPA:
            self.params["pa"].value = pa0

        if dims == 3:
            self._wl = getWlFromFitsImageCube(self._header, units.m)
            # NOTE: Adding the time dimension (nt,nwl,ny,nx)
            self._image = im.data[None, :, :, :]
        else:
            # NOTE: Adding the wl and time dimensions (nt,nwl,ny,nx)
            self._image = im.data[None, None, :, :]
            self._wl = np.array([0])

    def _internalImage(self):
        self.params["dim"].value = self._dim
        self._pixSize = self._pixSize0 * self.params["scale"].value
        return self._image
