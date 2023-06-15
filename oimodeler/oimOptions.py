"""Set global options of the oimodeler software."""
from .oimFTBackends import numpyFFTBackend

# NOTE: The dictionary oimOption contains all the customizable option of oimodeler
oimOptions = {}

# NOTE: Fourier transform settings
oimOptions["AvailableFTBackends"] = [numpyFFTBackend]
oimOptions["FTpaddingFactor"] = 8
oimOptions["FTbinningFactor"] = None
oimOptions["FTBackend"] = numpyFFTBackend
oimOptions["ModelOutput"] = "vis"

try:
    # NOTE: Only append the `FFTWBackend` if `fftw` is installed and imported
    import pyfftw
    from .oimFTBackends import FFTWBackend
    oimOptions["AvailableFTBackends"].append(FFTWBackend)
except ImportError:
    pass

# TODO: Maybe implement this as self.params of Fitter?
# NOTE: Fitter settings
oimOptions["FittingNCores"] = 6
