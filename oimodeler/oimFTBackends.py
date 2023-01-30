# -*- coding: utf-8 -*-
"""
Backends for Fourier Transform computation
"""
import numpy as np
from scipy import interpolate

import oimodeler as oim


class numpyFFTBackend():
    def check(backendPreparation, im, pix, wlin, tin, ucoord, vcoord, wl, t):
        return True

    def prepare(im, pix, wlin, tin, ucoord, vcoord, wl, t):
        return True

    def compute(backendPreparation, im, pix, wlin, tin, ucoord, vcoord, wl, t):

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
        vc = real+imag*1j

        return vc


oim.oimOptions["FTBackend"] = numpyFFTBackend
oim.oimOptions["AvailableFTBackends"] = [numpyFFTBackend]

try:
    import pyfftw

    class FFTWBackend():
        def check(backendPreparation, im, pix, wlin, tin, ucoord, vcoord, wl, t):
            try:
                fft_in, fft_out, fft_object, dim0, nwl0, nt0 = backendPreparation
            except:
                return False

            nwl1 = wlin.size
            nt1 = tin.size
            dim1 = im.shape[3]

            if (dim0, nwl0, nt0) == (dim1, nwl1, nt1):
                return True
            else:
                return False

        def prepare(im, pix, wlin, tin, ucoord, vcoord, wl, t):

            nwl = wlin.size
            nt = tin.size
            dim = im.shape[3]

            fft_in = pyfftw.empty_aligned(
                (nt, nwl, dim, dim), dtype='complex128')
            fft_out = pyfftw.empty_aligned(
                (nt, nwl, dim, dim), dtype='complex128')
            fft_object = pyfftw.FFTW(fft_in, fft_out, axes=(2, 3))

            return fft_in, fft_out, fft_object, dim, nwl, nt

        def compute(backendPreparation, im, pix, wlin, tin, ucoord, vcoord, wl, t):

            fft_in, fft_out, fft_object, dim, nwl, nt = backendPreparation

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
            vc = real+imag*1j

            return vc

    oim.oimOptions["AvailableFTBackends"].append(FFTWBackend)

except:
    pass
