"""Set global options of the oimodeler software."""
# NOTE: The dictionary oimOption contains all the customizable option
# of `oimodeler`.
oimOptions = {}

# NOTE: Fourier transform settings
oimOptions["AvailableFTBackends"] = []
oimOptions["FTPaddingFactor"] = None
oimOptions["FTBinningFactor"] = None
oimOptions["FTBackend"] = None
oimOptions["GridType"] = "linear"
oimOptions["FFTW_Initialized"] = False
