# -*- coding: utf-8 -*-
"""Set global options of the oimodeler software."""
from .oimFTBackends import numpyFFTBackend

# NOTE: The dictionary oimOption contains all the customizable option of oimodeler
oimOptions = {}
oimOptions["ModelType"] = "non-physical"
oimOptions["FTpaddingFactor"] = 8
oimOptions["FTBackend"] = numpyFFTBackend

# TODO: Should this be a dictionary?
oimOptions["AvailableFTBackends"] = [numpyFFTBackend]

try:
    # NOTE: Only append the `FFTWBackend` if `fftw` is installed and imported
    import pyfftw
    from .oimFTBackends import FFTWBackend
    oimOptions["AvailableFTBackends"].append(FFTWBackend)
except ImportError:
    pass
