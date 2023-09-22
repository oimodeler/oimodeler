"""Set global options of the oimodeler software."""
from .oimFTBackends import numpyFFTBackend, FFTWBackend

# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.
oimOptions = {}

# NOTE: Fourier transform settings
oimOptions["AvailableFTBackends"] = [numpyFFTBackend]
oimOptions["FTPaddingFactor"] = None
oimOptions["FTBinningFactor"] = None
oimOptions["FTBackend"] = numpyFFTBackend
oimOptions["GridType"] = "linear"

try:
    # NOTE: Only append the `FFTWBackend` if `fftw` is installed and imported
    import pyfftw
    from .oimFTBackends import FFTWBackend
    oimOptions["AvailableFTBackends"].append(FFTWBackend)
except ImportError:
    pass
