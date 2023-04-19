# -*- coding: utf-8 -*-
"""Creation of models"""
from typing import Union, Optional, Tuple, Dict, List

import astropy.units as u
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.io.fits import PrimaryHDU
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike

from .oimComponent import oimComponent
from .oimParam import oimParam, oimParamLinker, oimParamInterpolator


###############################################################################
class oimModel:
    """The oimModel class hold a model made of one or more components (derived
    from the oimComponent class).

    It allows to compute images (or image cubes for wavelength or time
    dependent models) and complex coherent fluxes for a vector of u,v,wl,
    and t coordinates.

    Parameters
    ----------
    *components : list of oimComponent
       The components of the model.

    Attributes
    ----------
    components : list of oimComponent
       The components of the model.
    """

    def __init__(self, *components: List[oimComponent]) -> None:
        """Constructor of the class"""
        if len(components) == 1 and type(components[0]) == list:
            self.components = components[0]
        else:
            self.components = components

    def getComplexCoherentFlux(self, ucoord: ArrayLike, vcoord: ArrayLike,
                               wl: Optional[ArrayLike] = None,
                               t: Optional[ArrayLike] = None) -> np.ndarray:
        """Compute and return the complex coherent flux for an array of u,v
        (and optionally wavelength and time) coordinates.

        Parameters
        ----------
        ucoord : array_like
            Spatial coordinate u (in cycles/rad).
        vcoord : array_like
            Spatial coordinate vu (in cycles/rad).
        wl : array_like, optional
            Wavelength(s) in meter. The default is None.
        t :  array_like, optional
            Time in s (mjd). The default is None.

        Returns
        -------
        numpy.ndarray
            The complex coherent flux. The same size as u & v
        """
        res = complex(0, 0)
        for c in self.components:
            res += c.getComplexCoherentFlux(ucoord, vcoord, wl, t)
        return res

    def getParameters(self, free: Optional[bool] = False) -> Dict[str, oimParam]:
        """Get the Model paramters (or free parameters)

        Parameters
        ----------
        free : bool, optional
            If True retrieve the free parameters of the models only.
            The default is False.

        Returns
        -------
        params : Dict of oimParam
            Dictionary of the model's parameters (or free parameters).
        """
        params = {}
        for i, c in enumerate(self.components):
            for name, param in c.params.items():
                if param not in params.values():
                    if isinstance(param, oimParamInterpolator):
                        for iparam, parami in enumerate(param.params):
                            if parami not in params.values():
                                if (parami.free or not free):
                                    params["c{0}_{1}_{2}_interp{3}".format(i+1, c.shortname.replace(" ", "_"), name, iparam+1)] = parami
                    elif isinstance(param, oimParamLinker):
                        pass
                    else:
                        if (param.free or not free):
                            params["c{0}_{1}_{2}".format(i+1, c.shortname.replace(" ", "_"), name)] = param
        return params

    def getFreeParameters(self) -> Dict[str, oimParam]:
        """Get the Model free paramters

        Returns
        -------
        Dict of oimParam
            A Dictionary of the model's free parameters.
        """
        return self.getParameters(free=True)

    def getImage(self, dim: int, pixSize: float,
                 wl: Optional[Union[float, ArrayLike]] = None,
                 t: Optional[Union[float, ArrayLike]] = None,
                 toFits: Optional[bool] = False,
                 fromFT: Optional[bool] = False,
                 squeeze: Optional[bool] = True,
                 normalize: Optional[bool] = False) -> Union[np.ndarray, PrimaryHDU]:
        """Compute and return an image or and image cube (if wavelength and time
        are given).

        The returned image as the x,y dimension dim in pixel with
        an angular pixel size pixSize in rad. Image is returned as a numpy
        array unless the keyword fits is set to True. In that case the image is
        returned as an astropy.io.fits hdu.

        Parameters
        ----------
        dim : int
            Image x & y dimension in pixels.
        pixSize : float
            Pixel angular size in mas.
        wl : int or array_like, optional
            Wavelength(s) in meter. The default is None.
        t :  int or array_like, optional
            Time in s (mjd). The default is None.
        toFits : bool, optional
            If True returns result as a fits hdu. The default is False.
        fromFT : bool, optional
            If True compute the image using FT formula when available.
            The default is False.
        squeeze : bool, optional
            If False returns a (nt,nwl,dim,dim) array even if nt and/or nwl equal 1.
            The default is True.
        normalize: bool, optional
            If True normalizes the image.

        Returns
        -------
        numpy.ndarray or astropy.io.fits.hdu
             A numpy 2D array (or 3D/4D array if wl, t or both are given) or an
             astropy.io.fits hdu.imagehdu if fits=True.
             The image of the component with given size in pixels and mas or rad
        """
        # TODO : maybe we should change all None to zero as default values
        if wl is None:
            wl = 0
        if t is None:
            t = 0

        t, wl = map(lambda x: np.array(x).flatten(), [t, wl])
        nt, nwl = t.size, wl.size
        dims = (nt, nwl, dim, dim)

        if fromFT:
            v = np.linspace(-0.5, 0.5, dim)
            vx, vy = np.meshgrid(v, v)

            vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
            vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
            wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
            t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

            spfx_arr = (vx_arr/pixSize/u.mas.to(u.rad)).flatten()
            spfy_arr = (vy_arr/pixSize/u.mas.to(u.rad)).flatten()
            wl_arr, t_arr = map(lambda x: x.flatten(), [wl_arr, t_arr])

            ft = self.getComplexCoherentFlux(spfx_arr, spfy_arr, wl_arr, t_arr).reshape(dims)
            image = np.abs(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(
                ft, axes=[-2, -1]), axes=[-2, -1]), axes=[-2, -1]))

        else:
            image = np.zeros(dims)
            for c in self.components:
                image += c.getImage(dim, pixSize, wl, t)

        if normalize:
            for it in range(nt):
                for iwl in range(nwl):
                    image[it, iwl, :, :] /= np.max(image[it, iwl, :, :])

        # Always squeeze dim which are equal to one if exported to fits format
        if squeeze or toFits:
            image = np.squeeze(image)

        if toFits:
            hdu = fits.PrimaryHDU(image)
            hdu.header['CDELT1'] = pixSize*u.mas.to(u.rad)
            hdu.header['CDELT2'] = pixSize*u.mas.to(u.rad)
            hdu.header['CRVAL1'] = 0
            hdu.header['CRVAL2'] = 0
            hdu.header['CRPIX1'] = dim/2
            hdu.header['CRPIX2'] = dim/2
            hdu.header['CUNIT1'] = "rad"
            hdu.header['CUNIT2'] = "rad"
            hdu.header['CROTA1'] = 0
            hdu.header['CROTA2'] = 0

            naxis = 3
            if nwl != 1:
                dwl = (np.roll(wl, -1)-wl)[:-1]

                if np.all(np.abs(dwl-dwl[0]) < 1e-12):
                    dwl = dwl[0]
                    hdu.header[f'CDELT{naxis}'] = dwl
                    hdu.header[f'CRPIX{naxis}'] = 1
                    hdu.header[f'CRVAL{naxis}'] = wl[0]
                    hdu.header[f'CUNIT{naxis}'] = "m"
                    naxis += 1

                else:
                    raise TypeError("Wavelength vector is not regular. Fit image"
                                    " with irregular grid not yet implemented")

            if nt != 1:
                dt = (np.roll(t, -1)-t)[:-1]

                if np.all(np.abs(dt-dt[0]) < 1e-12):
                    dt = dt[0]
                    hdu.header[f'CDELT{naxis}'] = dt
                    hdu.header[f'CRPIX{naxis}'] = 1
                    hdu.header[f'CRVAL{naxis}'] = t[0]
                    hdu.header[f'CUNIT{naxis}'] = "day"

                else:
                    raise TypeError("Time vector is not regular. Fit image"
                                    " with irregular grid not yet implemented")
            return hdu
        return image

    def saveImage(self, filename: str, dim: int, pixSize: float,
                  wl: Optional[Union[int, ArrayLike]] = None,
                  t: Optional[Union[int, ArrayLike]] = None,
                  fromFT: Optional[bool] = False,
                  normalize: Optional[bool] = False) -> Union[np.ndarray, PrimaryHDU]:
        """Save the model image

        Parameters
        ----------
        filename: str
            The name the file is to be saved as.
        dim : int
            Image x & y dimension in pixels.
        pixSize : float
            Pixel angular size in mas.
        wl : int or array_like, optional
            Wavelength(s) in meter. The default is None.
        t :  int or array_like, optional
            Time in s (mjd). The default is None.
        fromFT : bool, optional
            If True compute the image using FT formula when available.
            The default is False.
        normalize: bool, optional
            If True normalizes the image.

        Returns
        -------
        numpy.ndarray or astropy.io.fits.hdu
             A numpy 2D array (or 3D/4D array if wl, t or both are given).
             The image of the component with given size in pixels and mas or rad
        """
        im = self.getImage(dim, pixSize, wl=wl, t=t, toFits=True,
                           fromFT=fromFT, normalize=normalize)

        im.writeto(filename, overwrite=True)
        return im

    def showModel(self, dim: int, pixSize: float,
                  wl: Optional[Union[int, ArrayLike]] = None,
                  t: Optional[Union[int, ArrayLike]] = None,
                  fromFT: Optional[bool] = False,
                  axe: Optional[Axes] = None,
                  normPow: Optional[float] = 0.5,
                  figsize: Optional[Tuple[float]] = (3.5, 2.5),
                  savefig: Optional[str] = None,
                  colorbar: Optional[bool] = True,
                  legend: Optional[bool] = False,
                  swapAxes: Optional[bool] = True,
                  kwargs_legend: Dict = {},
                  normalize: Optional[bool] = False,
                  **kwargs: Dict) -> Tuple[Figure, Axes, np.ndarray]:
        """Show the mode Image or image-Cube

        Parameters
        ----------
        dim : int
            Image x & y dimension in pixels.
        pixSize : float
            Pixel angular size in mas.
        wl : int or array_like, optional
            Wavelength(s) in meter. The default is None.
        t :  int or array_like, optional
            Time in s (mjd). The default is None.
        fromFT : bool, optional
            If True compute the image using FT formula when available.
            The default is False.
        axe : matplotlib.axes.Axes, optional
            If provided the image will be shown in this axe. If not a new figure
            will be created. The default is None.
        normPow : float, optional
            Exponent for the Image colorscale powerLaw normalisation.
            The default is 0.5.
        figsize : tuple of float, optional
            The Figure size in inches. The default is (8., 6.).
        savefig : str, optional
            Name of the files for saving the figure If None the figure is not saved.
            The default is None.
        colorbar: bool, optional
            Add a colobar to the Axe. The default is True.
        legend : bool, optional
            If True displays a legend. Default is False.
        swapAxes : bool, optional
            If True swaps the axes of the wavelength and time.
            Default is True.
        kwargs_legend: dict, optional
        normalize : bool, optional
            If True normalizes the image.
        **kwargs : dict
            Arguments to be passed to the plt.imshow function.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The Figure created if needed.
        axe : matplotlib.axes.Axes
            The Axes instances, created if needed.
        im  : numpy.array
            The image(s).
        """
        im = self.getImage(dim, pixSize, wl, t, fromFT=fromFT,
                           squeeze=False, normalize=normalize)

        t, wl = map(lambda x: np.array(x).flatten(), [t, wl])

        if swapAxes:
            t, wl = wl, t

        nt, nwl = t.size, wl.size

        if axe is None:
            fig, axe = plt.subplots(nwl, nt, figsize=(
                figsize[0]*nt, figsize[1]*nwl), sharex=True, sharey=True, subplot_kw=dict(projection='oimAxes'))
        else:
            try:
                fig = axe.get_figure()
            except Exception:
                fig = axe.flatten()[0].get_figure()

        axe = np.array(axe).flatten().reshape((nwl, nt))

        if 'norm' not in kwargs:
            kwargs['norm'] = colors.PowerNorm(gamma=normPow)

        for iwl, wli in enumerate(wl):
            for it, ti in enumerate(t):
                if not swapAxes:
                    cb = axe[iwl, it].imshow(im[it, iwl, :, :],
                                             extent=[-dim/2*pixSize, dim/2*pixSize,
                                                     -dim/2*pixSize, dim/2*pixSize],
                                             origin='lower', **kwargs)
                else:
                    cb = axe[iwl, it].imshow(im[iwl, it, :, :],
                                             extent=[-dim/2*pixSize, dim/2*pixSize,
                                                     -dim/2*pixSize, dim/2*pixSize],
                                             origin='lower', **kwargs)

                axe[iwl, it].set_xlim(dim/2*pixSize, -dim/2*pixSize)

                if iwl == nwl-1:
                    axe[iwl, it].set_xlabel("$\\alpha$(mas)")
                if it == 0:
                    axe[iwl, it].set_ylabel("$\\delta$(mas)")

                if legend:
                    txt = ""
                    if not swapAxes:

                        if wl[0] is not None:
                            txt += "wl={:.4f}$\mu$m\n".format(wli*1e6)
                        if t[0] is not None:
                            txt += "Time={}".format(ti)
                        if 'color' not in kwargs_legend:
                            kwargs_legend['color'] = "w"
                    else:
                        if t[0] is not None:
                            txt += "wl={:.4f}$\mu$m\n".format(ti*1e6)
                        if wl[0] is not None:
                            txt += f"Time={wli}"
                        if 'color' not in kwargs_legend:
                            kwargs_legend['color'] = "w"
                    axe[iwl, it].text(0, 0.95*dim/2*pixSize, txt,
                                      va='top', ha='center', **kwargs_legend)

        if colorbar:
            fig.colorbar(cb, ax=axe, label="Normalized Intensity")

        if savefig is not None:
            plt.savefig(savefig)

        return fig, axe, im

    def showFourier(self, dim: int, pixSize: float,
                    wl: Optional[Union[int, ArrayLike]] = None,
                    t: Optional[Union[int, ArrayLike]] = None,
                    axe: Optional[Axes] = None,
                    normPow: Optional[float] = 0.5,
                    figsize: Optional[Tuple[float]] = (3.5, 2.5),
                    savefig: Optional[str] = None,
                    colorbar: Optional[bool] = True,
                    legend: Optional[bool] = False,
                    swapAxes: Optional[bool] = True,
                    display_mode: Optional[str] = "amp",
                    kwargs_legend: Optional[Dict] = {},
                    normalize: Optional[bool] = False,
                    **kwargs: Dict):
        """Show the amplitude and phase of the Fourier space

        Parameters
        ----------
        dim : int
            Image x & y dimension in pixels.
        pixSize : float
            Pixel angular size in mas.
        wl : int or array_like, optional
            Wavelength(s) in meter. The default is None.
        t :  int or array_like, optional
            Time in s (mjd). The default is None.
        axe : matplotlib.axes.Axes, optional
            If provided the image will be shown in this axe. If not a new figure
            will be created. The default is None.
        normPow : float, optional
            Exponent for the Image colorscale powerLaw normalisation.
            The default is 0.5.
        figsize : tuple of float, optional
            The Figure size in inches. The default is (8., 6.).
        savefig : str, optional
            Name of the files for saving the figure If None the figure is not saved.
            The default is None.
        colorbar : bool, optional
            Add a colobar to the Axe. The default is True.
        legend : bool, optional
            If True displays a legend. Default is False.
        swapAxes : bool, optional
            If True swaps the axes of the wavelength and time.
            Default is True.
        display_mode : str, optional
            Displays either the amplitude "amp" or the phase "phase".
            Default is "amp".
        kwargs_legend: dict, optional
        normalize : bool, optional
            If True normalizes the image.
        **kwargs : dict
            Arguments to be passed to the plt.imshow function.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The Figure created if needed
        axe : matplotlib.axes.Axes
            The Axes instances, created if needed.
        im  : numpy.ndarray
            The image(s).
        """
        t, wl = map(lambda x: np.array(x).flatten(), [t, wl])

        if swapAxes:
            t, wl = wl, t

        nt, nwl = t.size, wl.size
        dims = (nt, nwl, dim, dim)

        v = np.linspace(-0.5, 0.5, dim)
        vx, vy = np.meshgrid(v, v)

        vx_arr = np.tile(vx[None, None, :, :], (nt, nwl, 1, 1))
        vy_arr = np.tile(vy[None, None, :, :], (nt, nwl, 1, 1))
        wl_arr = np.tile(wl[None, :, None, None], (nt, 1, dim, dim))
        t_arr = np.tile(t[:, None, None, None], (1, nwl, dim, dim))

        spfx_arr, spfy_arr = map(lambda x: (x/pixSize/u.mas.to(u.rad)).flatten(), [vx_arr, vy_arr])
        wl_arr, t_arr = map(lambda x: x.flatten(), [wl_arr, t_arr])
        spfx_extent = spfx_arr.max()

        ft = self.getComplexCoherentFlux(spfx_arr, spfy_arr, wl_arr, t_arr).reshape(dims)
        if display_mode == "amp":
            im = np.abs(ft)
            im /= im.max()
        elif display_mode == "phase":
            im = np.angle(ft, deg=True)
        else:
            raise NameError("Only 'amp' and 'phase' are valid choice for the display_mode!")

        if normalize:
            for it in range(nt):
                for iwl in range(nwl):
                    im[it, iwl, :, :] /= np.max(im[it, iwl, :, :])

        if axe is None:
            fig, axe = plt.subplots(nwl, nt,
                                    figsize=(figsize[0]*nt, figsize[1]*nwl),
                                    sharex=True, sharey=True,
                                    subplot_kw={"projection": 'oimAxes'})
        else:
            try:
                fig = axe.get_figure()
            except Exception:
                fig = axe.flatten()[0].get_figure()

        axe = np.array(axe).flatten().reshape((nwl, nt))

        if 'norm' not in kwargs:
            kwargs['norm'] = colors.PowerNorm(gamma=normPow)

        for iwl, wli in enumerate(wl):
            for it, ti in enumerate(t):
                if not swapAxes:
                    cb = axe[iwl, it].imshow(im[it, iwl, :, :],
                                             extent=[-spfx_extent, spfx_extent,
                                                     -spfx_extent, spfx_extent],
                                             origin='lower', **kwargs)
                else:
                    cb = axe[iwl, it].imshow(im[iwl, it, :, :],
                                             extent=[-spfx_extent, spfx_extent,
                                                     -spfx_extent, spfx_extent],
                                             origin='lower', **kwargs)

                axe[iwl, it].set_xlim(-spfx_extent, spfx_extent)

                if iwl == nwl-1:
                    axe[iwl, it].set_xlabel("sp. freq. (cycles/rad)")
                if it == 0:
                    axe[iwl, it].set_ylabel("sp. freq. (cycles/rad)")

                if legend:
                    txt = ""
                    if not swapAxes:
                        if wl[0] is not None:
                            txt += "wl={:.4f}$\mu$m\n".format(wli*1e6)
                        if t[0] is not None:
                            txt += "Time={}".format(ti)
                        if 'color' not in kwargs_legend:
                            kwargs_legend['color'] = "w"
                    else:
                        if t[0] is not None:
                            txt += "wl={:.4f}$\mu$m\n".format(ti*1e6)
                        if wl[0] is not None:
                            txt += f"Time={wli}"
                        if 'color' not in kwargs_legend:
                            kwargs_legend['color'] = "w"
                    axe[iwl, it].text(0, 0.95*dim/2*pixSize, txt,
                                      va='top', ha='center', **kwargs_legend)
        if colorbar:
            fig.colorbar(cb, ax=axe, label="Normalized Intensity")
           
        if savefig is not None:
            plt.savefig(savefig)
        return fig, axe, im
