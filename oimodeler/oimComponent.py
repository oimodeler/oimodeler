# -*- coding: utf-8 -*-
"""Components defined in Fourier or image planes"""
from typing import Any

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy import units as units
from scipy import interpolate, integrate
from scipy.special import j0

from . import __dict__ as oimDict
from .oimOptions import oimOptions, standard_parameters
from .oimParam import oimInterp, oimParam, oimParamInterpolator
from .oimUtils import getWlFromFitsImageCube, pad_image, rebin_image


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


class oimComponent(object):
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
    """
    name = "Generic component"
    shortname = "Gen comp"
    description = "This is the class from which all components derived"

    def __init__(self, **kwargs):
        """Create and initiliaze a new instance of the oimComponent class.

        All components have at least three parameters the position
        x and y and their flux f
        """
        self._wl = None    # None value <=> All wavelengths (from Data) 
        self._t = [0]      # This component is static

        self.x = oimParam(**standard_parameters.x)
        self.y = oimParam(**standard_parameters.y)
        self.f = oimParam(**standard_parameters.f)
        self.dim = oimParam(**standard_parameters.dim)
        self._eval(**kwargs)

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
            if not isinstance(
                    value.value, (np.ndarray, tuple, list)):
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
            if not isinstance(
                    value.value, (np.ndarray, tuple, list)):
                value = [value]
        else:
            value = [value]
        self.__t = np.array(value)
        

    # TODO: Should the interpolator not take the attributes of the overhead param?
    # See how to do this…
    def _eval(self, **kwargs):
        """Evaluates the parameters from the keyword arguments and
        if the parameter exists for the model it will set its value
        accordingly.
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                attribute = getattr(self, key)
                if isinstance(value, oimInterp):
                    if not isinstance(attribute, oimParamInterpolator):
                        setattr(self, key, value.type(attribute, **value.kwargs))
                else:
                    attribute.value = value

    # def _paramstr(self):
    #     txt = ""
    #     for _, param in self.params.items():
    #         if isinstance(param, oimParam):
    #             if 'value' in param.__dict__:
    #                 txt += " {0}={1:.2f}".format(param.name, param.value)
    #             else:
    #                 # TODO: Have a string for each oimParamInterpolator
    #                 txt += " {0}={1}".format(param.name,
    #                                          param.__class__.__name__)
    #     return txt

    # def __str__(self):
    #     """String representation of the class."""
    #     return self.name + self._paramstr()

    # def __repr__(self):
    #     """Console representation of the class."""
    #     return self.__class__.__name__ + " at " + str(hex(id(self)))\
    #             + " : " + self._paramstr()

    def _ftTranslateFactor(self, ucoord, vcoord, wl, t):
        x = self.x(wl, t)*self.x.unit.to(units.rad)
        y = self.y(wl, t)*self.y.unit.to(units.rad)
        return np.exp(-2*1j*np.pi*(ucoord*x+vcoord*y))

    def _directTranslate(self, x, y, wl, t):
        """Translate the x,y coordinates to the direct space
        (i.e. the image space) using the x and y parameters of the component."""
        return x-self.x(wl, t), y-self.y(wl, t)

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
        return np.array(u)*0

    def getImage(self, dim, pixSize, wl=None, t=None):
        """ Compute and return an image or and image cube (if wavelength and time
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
        return np.zeros(dim, dim)


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
        if ("elong" in kwargs) or ("pa" in kwargs) or self.elliptic:
            self.elong = oimParam(**standard_parameters.elong)
            self.pa = oimParam(**standard_parameters.pa)
            self.elliptic = True

        self._eval(**kwargs)

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        if self.elliptic:
            pa_rad = (self.pa(wl, t)) * self.pa.unit.to(units.rad)
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = (ucoord*co-vcoord*si)/self.elong(wl, t)
            fyp = ucoord*si+vcoord*co
            rho = np.hypot(fxp, fyp)
        else:
            fxp, fyp = ucoord, vcoord
            rho = np.hypot(fxp, fyp)

        vc = self._visFunction(fxp, fyp, rho, wl, t)

        return vc*self._ftTranslateFactor(ucoord, vcoord, wl, t) * self.f(wl, t)

    def _visFunction(self, ucoord, vcoord, rho, wl, t):
        return ucoord*0

    def getImage(self, dim, pixSize, wl=None, t=None):
        t = np.array(t).flatten()
        wl = np.array(wl).flatten()
        nt, nwl = t.size, wl.size

        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim)

        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr*pixSize*dim).flatten()
        y_arr = (vy_arr*pixSize*dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)

        if self.elliptic:
            pa_rad = (self.pa(wl_arr, t_arr)) * self.pa.unit.to(units.rad)

            xp = x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp = x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)

            x_arr = xp*self.elong(wl_arr, t_arr)
            y_arr = yp

        image = self._imageFunction(x_arr.reshape(dims), y_arr.reshape(dims),
                                    wl_arr.reshape(dims), t_arr.reshape(dims))

        tot = np.sum(image, axis=(2, 3))

        # TODO: No loop for normalization
        for it, ti in enumerate(t):
            for iwl, wli in enumerate(wl):
                if tot[it, iwl] != 0:
                    image[it, iwl, :, :] = image[it, iwl, :, :]  \
                        / tot[it, iwl]*self.f(wli, ti)
        return image

    def _imageFunction(self, xx, yy, wl, t):
        raise ValueError(f"image function not implemented for {self.shortname}."
                         "Use the fromFFT=True option to get a model image"
                         " from the inverse Fourier Transform")


class oimComponentImage(oimComponent):
    """Base class for components define in 2D : x,y (regular grid) in the image plan"""
    elliptic = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._pixSize = 0   # NOTE: In rad
        self._allowExternalRotation = True
        self.normalizeImage = True
        self.dim = oimParam(**standard_parameters.dim)
        self.pa = oimParam(**standard_parameters.pa)

        # NOTE: Add ellipticity
        if self.elliptic:
            self.elong = oimParam(**standard_parameters.elong)

        if 'FTBackend' in kwargs:
            self.FTBackend = kwargs['FTBackend']()
        else:
            self.FTBackend = oimOptions.ft.backend.active()

        self.FTBackendData = None
        self._eval(**kwargs)

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        if wl is None:
            wl = ucoord*0
        if t is None:
            t = ucoord*0

        im = self.getInternalImage(wl, t)

        if oimOptions.ft.binning is not None:
            im = rebin_image(im, oimOptions.ft.binning)

        im = pad_image(im)
        pix = self._pixSize

        tr = self._ftTranslateFactor(
            ucoord, vcoord, wl, t)#♣*self.params["f"](wl, t)

        if self._allowExternalRotation:
            pa_rad = (self.pa(wl, t)) * \
                self.pa.unit.to(units.rad)
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = ucoord*co-vcoord*si
            fyp = ucoord*si+vcoord*co
            vcoord = fyp

            if self.elliptic == True:
                ucoord = fxp/self.elong(wl, t)
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

        if self.FTBackend.check(self.FTBackendData, im, pix, wl0,
                                t0, ucoord, vcoord, wl, t) == False:

            self.FTBackendData = self.FTBackend.prepare(
                    im, pix, wl0, t0, ucoord, vcoord, wl, t)

        vc = self.FTBackend.compute(self.FTBackendData, im, pix, wl0,
                                    t0, ucoord, vcoord, wl, t)

        return vc*tr*self.f(wl, t)

    def getImage(self, dim, pixSize, wl=None, t=None):
        if wl is None:
            wl = 0
        if t is None:
            t = 0

        t = np.array(t).flatten()
        wl = np.array(wl).flatten()
        nt, nwl = t.size, wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr*pixSize*dim).flatten()
        y_arr = (vy_arr*pixSize*dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)

        if self._allowExternalRotation:
            pa_rad = self.pa(wl_arr, t_arr) * self.pa.unit.to(units.rad)

            xp = x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp = x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)
            y_arr = yp

            if self.elliptic:
                x_arr = xp*self.elong(wl_arr, t_arr)
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
                grid, im0, coord, bounds_error=False, fill_value=None)
            f0 = np.sum(im0)
            f = np.sum(im)
            im = im/f*f0

        im = im.reshape(dims)

        if self.normalizeImage:
            # TODO: No loop for normalization
            tot = np.sum(im, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        im[it, iwl, :, :] = im[it, iwl, :, :]  \
                            / tot[it, iwl]*self.f(wli, ti)
        return im

    def getInternalImage(self, wl, t):
        res = self._internalImage()

        if res is None:
            t_arr, wl_arr, x_arr, y_arr = self._getInternalGrid(
                simple=False, wl=wl, t=t)
            res = self._imageFunction(x_arr, y_arr, wl_arr, t_arr)

        if self.normalizeImage:
            for it in range(res.shape[0]):
                for iwl in range(res.shape[1]):
                    res[it, iwl, :, :] = res[it, iwl, :, :] / \
                        np.sum(res[it, iwl, :, :])

        return res

    def _internalImage(self):
        return

    def _imageFunction(self, xx, yy, wl, t):
        image = xx*0+1
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

        dim = self.dim(wl, t)

        pix = self._pixSize*units.rad.to(units.mas)
        v = np.linspace(-0.5, 0.5, dim)
        xy = v*pix*dim

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

            if flatten:
                return t_arr.flatten(), wl_arr.flatten(), \
                    x_arr.flatten(), y_arr.flatten()
            else:
                return t_arr, wl_arr, x_arr, y_arr


class oimComponentRadialProfile(oimComponent):
    """Base class for components define by their radial profile"""
    elliptic = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._r = None
        self.normalizeImage = True
        self.hankel = self.fht

        self._t = None
        self._wl = None
        self._r = None

        # FIXME: Is this not redundant as oimComponent is already ellpitical?
        # NOTE: Add ellipticity
        if self.elliptic:
            self.pa = oimParam(**standard_parameters.pa)
            self.elong = oimParam(**standard_parameters.elong)

        self._eval(**kwargs)

    def _getInternalGrid(self, simple=True, flatten=False, wl=None, t=None):
        if self._wl is None:
            wl0 = np.sort(np.unique(wl))
        else:
            wl0 = self._wl

        if self._t is None:
            t0 = np.sort(np.unique(t))
        else:
            t0 = self._t

        if self._r is None:
            pix = self._pixSize*units.rad.to(units.mas)
            r = np.linspace(0, self.dim.value-1, self.dim.value)*pix
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
            t_arr = np.tile(t[:, None, None],  (1, nwl, nr))

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
                simple=False, wl=wl, t=t)
            res = self._radialProfileFunction(r_arr, wl_arr, t_arr)
        for it in range(res.shape[0]):
            for iwl in range(res.shape[1]):
                res[it, iwl, :] = res[it, iwl, :]/np.sum(res[it, iwl, :])

        return res

    @staticmethod
    def fht(Ir, r, wlin, tin, sfreq, wl, t):
        pad = oimOptions.ft.padding
        nr = r.size
        ntin = tin.size
        nwlin = wlin.size

        fov = r[-1]-r[0]
        dsfreq0 = 1/(fov*pad)
        sfreq0 = np.linspace(0, pad*nr-1, pad*nr)*dsfreq0

        r2D = np.tile(r[None, None, :, None], (ntin, nwlin, 1, pad*nr))
        Ir2D = np.tile(Ir[:, :, :, None], (1, 1, 1, pad*nr))
        sf2D = np.tile(sfreq0[None, None, None, :], (ntin, nwlin, nr, 1))

        res0 = integrate.trapezoid(
            2.*np.pi*r2D*Ir2D*j0(2.*np.pi*r2D*sf2D), r2D, axis=2)

        for it in range(ntin):
            for iwl in range(nwlin):
                res0[it,iwl,:] /= integrate.trapezoid(2.*np.pi*r*Ir[it,iwl,:], r)
                # res0[it,iwl,:]/=np.sum(r*Ir[it,iwl,:]*dr)

        grid = (tin, wlin, sfreq0)
        coord = np.transpose([t, wl, sfreq])

        real = interpolate.interpn(grid, np.real(res0), coord,
                                   bounds_error=False, fill_value=None)
        imag = interpolate.interpn(grid, np.imag(res0), coord,
                                   bounds_error=False, fill_value=None)
        return real+imag*1j

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

        v = np.linspace(-0.5, 0.5, dim)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        x_arr = (vx_arr*pixSize*dim).flatten()
        y_arr = (vy_arr*pixSize*dim).flatten()
        wl_arr = wl_arr.flatten()
        t_arr = t_arr.flatten()

        x_arr, y_arr = self._directTranslate(x_arr, y_arr, wl_arr, t_arr)

        if self.elliptic:
            pa_rad = self.pa(wl_arr, t_arr) * self.pa.unit.to(units.rad)

            xp = x_arr*np.cos(pa_rad)-y_arr*np.sin(pa_rad)
            yp = x_arr*np.sin(pa_rad)+y_arr*np.cos(pa_rad)
            y_arr = yp
            x_arr = xp*self.elong(wl_arr, t_arr)

        r_arr = np.sqrt(x_arr**2+y_arr**2)
        im = self._radialProfileFunction(r_arr, wl_arr, t_arr)
        im = im.reshape(dims)

        if self.normalizeImage == True:
            # TODO: No loop for normalization
            tot = np.sum(im, axis=(2, 3))
            for it, ti in enumerate(t):
                for iwl, wli in enumerate(wl):
                    if tot[it, iwl] != 0:
                        im[it, iwl, :, :] = im[it, iwl, :, :]  \
                            / tot[it, iwl]*self.f(wli, ti)

        return im

    def getComplexCoherentFlux(self, ucoord, vcoord, wl=None, t=None):
        if wl is None:
            wl = ucoord*0
        if t is None:
            t = ucoord*0

        if self.elliptic:
            pa_rad = self.pa(wl, t) * self.pa.unit.to(units.rad)
            co, si = np.cos(pa_rad), np.sin(pa_rad)
            fxp = ucoord*co-vcoord*si
            fyp = ucoord*si+vcoord*co
            vcoord = fyp
            ucoord = fxp/self.elong(wl, t)

        spf = np.hypot(ucoord, vcoord)

        if self._wl is None:
            wl0 = np.sort(np.unique(wl))
        else:
            wl0 = self._wl

        if self._t is None:
            t0 = np.sort(np.unique(t))
        else:
            t0 = self._t

        Ir = self.getInternalRadialProfile(wl0, t0)

        vc = self.hankel(Ir, self._r*units.mas.to(units.rad),
                         wl0, t0, spf, wl, t)
        return vc*self._ftTranslateFactor(ucoord, vcoord, wl, t)*self.f(wl, t)


class oimComponentFitsImage(oimComponentImage):
    """Component load load images or chromatic-cubes from fits files"""
    elliptic = False
    name = "Fits Image Component"
    shortname = "Fits_Comp"

    def __init__(self, fitsImage, useinternalPA=False, **kwargs):
        super().__init__(**kwargs)
        self.loadImage(fitsImage, useinternalPA=useinternalPA)
        self.pa = oimParam(**standard_parameters.pa)
        self.scale = oimParam(**standard_parameters.scale)

        self._eval(**kwargs)

    def loadImage(self, fitsImage, useinternalPA=False):
        if isinstance(fitsImage, str):
            try:
                im = fits.open(fitsImage)[0]
            except:
                raise TypeError("Not a valid fits file")
        elif isinstance(fitsImage, fits.hdu.hdulist.HDUList):
            im = fitsImage[0]
        elif isinstance(fitsImage, fits.hdu.image.PrimaryHDU):
            im = fitsImage

        self._header = im.header

        dims = self._header['NAXIS']
        if dims < 2:
            raise TypeError("oimComponentFitsImage require 2D images or "
                            "3D chromatic-image-cubes")

        dimx = self._header['NAXIS1']
        dimy = self._header['NAXIS2']
        if dimx != dimy:
            raise TypeError("Current version only works with square images")
        self._dim = dimx

        pixX = self._header["CDELT1"]
        pixY = self._header["CDELT2"]
        if pixX != pixY:
            raise TypeError("Current version only works with the same pixels"
                            " scale in x and y dimension")
        if "CUNIT1" in self._header:
            try:
                unit0 = units.Unit(self._header["CUNIT1"])
            except:
                unit0 = units.rad
                    
        else:
            unit0 = units.rad
        self._pixSize0 = pixX*unit0.to(units.rad)

        if "CROTA1" in self._header:
            pa0 = self._header["CROTA1"]
        elif "CROTA2" in self._header:
            pa0 = self._header["CROTA2"]

        # TODO: Check units
        if useinternalPA:
            self.pa.value = pa0

        if dims == 3:
            self._wl = getWlFromFitsImageCube(self._header, units.m)
            # NOTE: Adding the time dimension (nt,nwl,ny,nx)
            self._image = im.data[None, :, :, :]
        else:
            # NOTE: Adding the wl and time dimensions (nt,nwl,ny,nx)
            self._image = im.data[None, None, :, :]

    def _internalImage(self):
        self.dim.value = self._dim
        self._pixSize = self._pixSize0*self.scale.value
        return self._image
