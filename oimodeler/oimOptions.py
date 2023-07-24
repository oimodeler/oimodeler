# -*- coding: utf-8 -*-
"""Set global options of the oimodeler software."""
from .oimFTBackends import numpyFFTBackend, FFTWBackend

# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.
oimOptions = {}
oimOptions["FTpaddingFactor"] = 8
oimOptions["FTbinningFactor"] = None
oimOptions["FTBackend"] = numpyFFTBackend
oimOptions["AvailableFTBackends"] = [numpyFFTBackend, FFTWBackend]
