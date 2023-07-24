# -*- coding: utf-8 -*-
"""Set global options of the oimodeler software."""
import pyfftw

from .oimFTBackends import numpyFFTBackend, FFTWBackend

# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.
oimOptions = {}
oimOptions["FTpaddingFactor"] = 8
oimOptions["FTbinningFactor"] = None
oimOptions["FTBackend"] = numpyFFTBackend
oimOptions["AvailableFTBackends"] = [numpyFFTBackend, FFTWBackend]

try:
    # NOTE: Check if `FFTW` backend works.
    test_array = pyfftw.empty_aligned(100, dtype="complex128")
    fft_object = pyfftw.FFTW(test_array, test_array, axes=(0,))
    transformed = fft_object(test_array)
    oimOptions["FTfftw"] = True
except Exception:
    oimOptions["FTfftw"] = False
