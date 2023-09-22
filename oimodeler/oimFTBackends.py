# -*- coding: utf-8 -*-
"""
The oimBackends.py module contains the backends for Fast-Fourier Transform
computation used in all model components derived from the oimImageComponent
semi-abstract component. This include most of the 2D-image-defined component.
The varoius FFT backends are static classes with three methods:

-   ``check`` : is called each time the component :func:`getComplexCoherentFlux
    <oimodeler.oimComponent.oimComponent.getComplexCoherentFlux>` method
    is called. It checks if the backend is ready to perform the FFT
-   ``pepare`` : is called only if the `check` method returned False.It
    prepares the backends and the return various data relative to the
    preparation
-   ``compute`` : is finally called when the backend is ready to compute
    the FFT and return the complexCoherentFlux.
"""
from typing import Tuple

import numpy as np
import pyfftw
from numpy.typing import ArrayLike
from scipy import interpolate


try:
    # NOTE: Check if `FFTW` backend is properly installed.
    test_array = pyfftw.empty_aligned(100, dtype="complex128")
    fft_object = pyfftw.FFTW(test_array, test_array, axes=(0,))
    transformed = fft_object(test_array)
    FFTW_CHECK = True
except Exception:
    FFTW_CHECK = False
    def _err() -> None:
        """Import error function that gets called when the FFTWBackend
        dependendy pyfftw has not been properly imported."""
        raise ImportError("This `FFTWBackend`-class has not been defined for "
                          "`oimodeler`. This means that the `FFTW` library is "
                          "either missing or has not been properly "
                          "installed.")


class numpyFFTBackend:
    """
    Default FFT backend using the numpy np.fft.fft2 function.

    Empty ``check`` and ``prepare`` methods. interpolation of the FFT done 
    with scipy.interpolate.interpn in the ``compute`` method.
    """
    @staticmethod
    def check(backendPreparation: bool, im: np.ndarray,
              pix: float, wlin: np.ndarray, tin: np.ndarray,
              ucoord: ArrayLike, vcoord: ArrayLike,
              wl: ArrayLike, t: ArrayLike) -> bool:
        """Checks if the backend is ready to compute the FFT.

        In the case of the simple numpy FFT backend no preparation are needed.
        Always return True

        Parameters
        ----------
        backendPreparation : bool
            The FFTBackendPreparation structure. Always True for the
            numpy backend.
        im : numpy.ndarray
            4D image (t,wl,x,y) to be FFTed.
        pix : float
            pixel size of the image in rad.
        wlin : numpy.ndarray
            the input wavelength vector of the image.
        tin : numpy.ndarray
            the input time vector of the image.
        ucoord : array_like
            the u coordinate of the baselines.
        vcoord : array_like
            the v coordinate of the baselines.
        wl : array_like
            the wl coordinate of the baselines.
        t : array_like
            the t coordinate of the baselines.

        Returns
        -------
        bool
            In the case of the numpy FFT backend it is always equal to True.
        """
        return True

    @staticmethod
    def prepare(im: np.ndarray, pix: float,
                wlin: np.ndarray, tin: np.ndarray,
                ucoord: ArrayLike, vcoord: ArrayLike,
                wl: ArrayLike, t: ArrayLike) -> bool:
        """Prepares the backend to compute the FFT if not ready.

        In the case of the simple numpy FFT backend no preparation are needed.
        Always return True.

        Parameters
        ----------
        im : numpy.ndarray
            4D image (t,wl,x,y) to be FFTed.
        pix : float
            pixel size of the image in rad.
        wlin : numpy.ndarray
            the input wavelength vector of the image.
        tin : numpy.ndarray
            the input time vector of the image.
        ucoord : array_like
            the u coordinate of the baselines.
        vcoord : array_like
            the v coordinate of the baselines.
        wl : array_like
            the wl coordinate of the baselines.
        t : array_like
            the t coordinate of the baselines.

        Returns
        -------
        bool
            The FFTBackendPreparation structure. Always True for the numpy
            backend.
        """
        return True

    @staticmethod
    def compute(backendPreparation: bool, im: np.ndarray,
                pix: float, wlin: np.ndarray, tin: np.ndarray,
                ucoord: ArrayLike, vcoord: ArrayLike,
                wl: ArrayLike, t: ArrayLike) -> np.ndarray:
        """Computes the FFT and interpolate the results at the
        required coordinates.

        It computes the FFT of the 4D image (t,wl,x,y) using the
        `numpy.fft.fft2` function and interpolate the results at the proper
        spatial, spectral and temporal coordinates using the
        `scipy.interpolate.interpn` function.

        Parameters
        ----------
        backendPreparation : bool
            The FFTBackendPreparation structure. Always True for the
            numpy backend.
        im : numpy.ndarray
            4D image (t,wl,x,y) to be FFTed.
        pix : float
            pixel size of the image in rad.
        wlin : numpy.ndarray
            the input wavelength vector of the image.
        tin : numpy.ndarry
            the input time vector of the image.
        ucoord : array_like
            the u coordinate of the baselines.
        vcoord : array_like
            the v coordinate of the baselines.
        wl : array_like
            the wl coordinate of the baselines.
        t : array_like
            the t coordinate of the baselines.

        Returns
        -------
        numpy.ndarray (complex)
           The computed and interpolated complex FFT of the image at the the
           proper spatial, spectral and temporal coordinates.
        """
        dim = im.shape[3]
        fft2D = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(
            im, axes=[-2, -1]), axes=[-2, -1]), axes=[-2, -1])

        freqVectYX = np.fft.fftshift(np.fft.fftfreq(dim, pix))
        grid = (tin, wlin, freqVectYX, freqVectYX)

        coord = np.transpose([t, wl, vcoord, ucoord])
        real = interpolate.interpn(grid, np.real(
            fft2D), coord, bounds_error=False, fill_value=None)
        imag = interpolate.interpn(grid, np.imag(
            fft2D), coord, bounds_error=False, fill_value=None)
        return real+imag*1j


class FFTWBackend:
    """FFT backend based on the python implementation of FFTW library.

    The ``prepare`` method create three FFTW objects: The IN and
    OUT arrays, at the format `pyfftw.empty_aligned`, and the
    `pyfftw.FFTW` object for the transformation.

    The ``check`` method only check the size of the arrays in x and y (dim)
    and also in wavelength and time

    After the FFT computation in the ``compute`` method, interpolation to
    the proper coordinates are done with the
    `scipy.interpolate.interpn` method.

    The backendPreparation contains the following six elements: fft_in,
    (pyfftw.empty_aligned), fft_out (pyfftw.empty_aligned), fft_object
    (pyfftw.FFTW), dim (int), nwl (int), nt (int).
    """
    @staticmethod
    def check(backendPreparation: bool, im: np.ndarray,
              pix: float, wlin: np.ndarray, tin: np.ndarray,
              ucoord: ArrayLike, vcoord: ArrayLike,
              wl: ArrayLike, t: ArrayLike) -> bool:
        """Checks if the backend is ready to compute the FFT.

        Parameters
        ----------
        backendPreparation : bool
            The FFTBackendPreparation structure. Always True for the
            numpy backend.
        im : numpy.ndarray
            4D image (t,wl,x,y) to be FFTed.
        pix : float
            pixel size of the image in rad.
        wlin : numpy.ndarray
            the input wavelength vector of the image.
        tin : numpy.ndarray
            the input time vector of the image.
        ucoord : array_like
            the u coordinate of the baselines.
        vcoord : array_like
            the v coordinate of the baselines.
        wl : array_like
            the wl coordinate of the baselines.
        t : array_like
            the t coordinate of the baselines.

        Returns
        -------
        bool
            In the case of the numpy FFT backend it is always equal to True.
        """
        if not FFTW_CHECK:
            _err()
            return
        try:
            _, _, _, dim0, nwl0, nt0 = backendPreparation
        except Exception:
            return False
        nwl1, nt1, dim1 = wlin.size, tin.size, im.shape[3]
        return (dim0, nwl0, nt0) == (dim1, nwl1, nt1)

    @staticmethod
    def prepare(im: np.ndarray, pix: float,
                wlin: np.ndarray, tin: np.ndarray,
                ucoord: ArrayLike, vcoord: ArrayLike,
                wl: ArrayLike, t: ArrayLike) -> Tuple:
        """Prepares the backend  to compute the FFT if not ready.

        The preparation is done by initializing two `pyfftw.empty_aligned`
        arrays and one `pyfftw.FFTW` on these arrays.

        Parameters
        ----------
        im : numpy.ndarray
            4D image (t,wl,x,y) to be FFTed.
        pix : float
            pixel size of the image in rad.
        wlin : numpy.ndarray
            the input wavelength vector of the image.
        tin : numpy.ndarray
            the input time vector of the image.
        ucoord : array_like
            the u coordinate of the baselines.
        vcoord : array_like
            the v coordinate of the baselines.
        wl : array_like
            the wl coordinate of the baselines.
        t : array_like
            the t coordinate of the baselines.

        Returns
        -------
        tuple
            A tuple of six elements containing the dimension of the image,
            the number of wavelength and time, and  three FFTW objects:
            IN and OUT arrays at `pyfftw.empty_aligned` and its transform
            as `pyfftw.FFTW`.
        """
        if not FFTW_CHECK:
            _err()
            return
        nwl, nt, dim = wlin.size, tin.size, im.shape[3]
        fft_in = pyfftw.empty_aligned(
            (nt, nwl, dim, dim), dtype='complex128')
        fft_out = pyfftw.empty_aligned(
            (nt, nwl, dim, dim), dtype='complex128')
        fft_object = pyfftw.FFTW(fft_in, fft_out, axes=(2, 3))
        return fft_in, fft_out, fft_object, dim, nwl, nt

    @staticmethod
    def compute(backendPreparation: Tuple, im: np.ndarray,
                pix: float, wlin: np.ndarray, tin: np.ndarray,
                ucoord: ArrayLike, vcoord: ArrayLike,
                wl: ArrayLike, t: ArrayLike) -> np.ndarray:
        """Computes the FFT and interpolate the results at the
        required coordinates.

        It computes the FFT of the 4D image (t,wl,x,y) using FFTW
        function and interpolate the results at the proper spatial,
        spectral and temporal coordinates using the
        `scipy.interpolate.interpn` function.

        Parameters
        ----------


        Returns
        -------
        numpy.ndarray (complex)
           The computed and interpolated complex FFT of the image at the
           proper spatial, spectral and temporal coordinates.

        """
        if not FFTW_CHECK:
            _err()
            return

        fft_in, fft_out, fft_object, dim, _, _ = backendPreparation
        fft_in[:] = np.fft.fftshift(im, axes=[-2, -1])
        fft_object()

        fft2D = np.fft.ifftshift(fft_out, axes=[-2, -1])

        freqVectYX = np.fft.fftshift(np.fft.fftfreq(dim, pix))
        grid = (tin, wlin, freqVectYX, freqVectYX)
        coord = np.transpose([t, wl, vcoord, ucoord])

        real = interpolate.interpn(grid, np.real(
            fft2D), coord, bounds_error=False, fill_value=None)
        imag = interpolate.interpn(grid, np.imag(
            fft2D), coord, bounds_error=False, fill_value=None)
        return real+imag*1j
